#include "singlesetup.hpp"
#include "utils.hpp"
#include "border.hpp"
#include <fstream>
#include <filesystem>


static std::vector<std::vector<int>> split_at_most_N_into_M_slots(int N, int M) {
    if (M == 0)
        return {};
    if (N == 0)
        return {std::vector<int>(M, 0)};
    if (M == 1) {
        std::vector<std::vector<int>> result;
        for (int i = 0; i <= N; i++) {
            result.push_back({i});
        }
        return result;
    }

    std::vector<std::vector<int>> result;
    for (int i = 0; i <= N; i++) {
        auto subresult = split_at_most_N_into_M_slots(N - i, M - 1);
        for (auto& sub: subresult) {
            sub.push_back(i);
            result.push_back(sub);
        }
    }
    std::sort(result.begin(), result.end(), 
        [](const auto& x, const auto& y) {
            int sum_x = 0, sum_y = 0;
            for (auto& px: x) sum_x += px;
            for (auto& py: y) sum_y += py;
            return sum_x < sum_y;
        }
    );
    return result;
}


integral_system::integral_system(integralfamily&                       family, 
                                 const std::vector<std::vector<int>>&  integrals, 
                                 const std::vector<std::vector<int>>&  preferred, 
                                 const YAML::Node&                     config, 
                                 const char*                           workdir)
    : pfamily(&family), psymbols(family.psymbols), pconfig(&config), reducer(config, family), workdirstr(workdir) {
    
    for (auto& idx: integrals)
        targets.push_back(integral(family, idx));
    for (auto& idx: preferred)
        preferred_masters.push_back(integral(family, idx));
    
    std::filesystem::create_directory(std::filesystem::path(workdirstr));
    tempdirstr = std::filesystem::path(workdirstr).append("KIRATEMP");
    std::filesystem::create_directory(std::filesystem::path(tempdirstr));

    eta_mode = {Prescription, Mass, Propagator};
    try {
        if (has_non_null_key(config, "Run")) {
            auto run_options = config["Run"];
            if (has_non_null_key(run_options, "AMFMode")) {
                auto mode_strings = run_options["AMFMode"].as<std::vector<std::string>>();
                eta_mode = {};
                for (auto& mode: mode_strings) {
                    if (mode == "Prescription") {
                        eta_mode.push_back(Prescription);
                    } else if (mode == "Mass") {
                        eta_mode.push_back(Mass);
                    } else if (mode == "Propagator") {
                        eta_mode.push_back(Propagator);
                    } else if (mode == "Branch") {
                        eta_mode.push_back(Branch);
                    } else if (mode == "Loop") {
                        eta_mode.push_back(Loop);
                    } else if (mode == "All") {
                        eta_mode.push_back(All);
                    } else {
                        std::cerr << "SingleSetup: encountered unrecognized eta insertion rules\n";
                        exit(1);
                    }
                }
            }
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "SingleSetup: failed while reading eta insertion rules\n";
        exit(1);
    }
    std::vector<std::string> set_modes;
    for (auto& m: eta_mode) {
        switch (m) {
            case Prescription:
                set_modes.push_back("Prescription");
                break;
            case Mass:
                set_modes.push_back("Mass");
                break;
            case Propagator:
                set_modes.push_back("Propagator");
                break;
            case Branch:
                set_modes.push_back("Branch");
                break;
            case Loop:
                set_modes.push_back("Loop");
                break;
            case All:
                set_modes.push_back("All");
                break;
        }
    }
    std::cerr << "SingleSetup: set eta insertion rules [" << set_modes << "]\n";
}


void integral_system::reduce_targets() {
    std::cerr << "SingleSetup: reducing " << targets.size() << " target integrals\n";
    auto reduce_result = reducer.reduce(targets, preferred_masters, tempdirstr.c_str());
    std::ofstream out(std::filesystem::path(workdirstr).append("REDUCE"));
    out << "Targets: "      << targets                             << "\n";
    out << "Masters: "      << reduce_result.first                 << "\n";
    out << "Coefficients: " << matrix_to_lst(reduce_result.second) << "\n";
    out.close();
    masters = reduce_result.first;
}


void integral_system::build_diffeq() {
    std::cerr << "SingleSetup: building differential equation for " << masters.size() << " master integrals\n";
    auto diffeq_result = reducer.diffeq(masters, tempdirstr.c_str());
    std::ofstream out(std::filesystem::path(workdirstr).append("DIFFEQ"));
    out << "Masters: " << diffeq_result.first                 << "\n";
    out << "Diffeq: "  << matrix_to_lst(diffeq_result.second) << "\n";
    out.close();
    masters = diffeq_result.first;
}


void integral_system::determine_boundaries() {
    master_sector = top_sector(masters);
    auto all_regions = boundary_region::all_boundary_regions(*pfamily, master_sector);
    for (auto& region: all_regions) {
        std::cerr << "SingleSetup: considering region " << region.to_string() << "\n";
        if (!region.is_trivial_region(master_sector))
            regions.push_back(region);
        else
            std::cerr << "SingleSetup: region " << region.to_string() << " is a trivial region\n";
    }
    std::ofstream out(std::filesystem::path(workdirstr).append("BOUNDARIES"));
    for (auto& region: regions) {
        out << region.to_string() << "\n";
    }
    out.close();
}


void integral_system::determine_border() {
    std::ifstream in(std::filesystem::path(workdirstr).append("DIFFEQ"));
    std::string diffeq_str;
    for (int i = 0; i < 4; i++)
        in >> diffeq_str;
    in.close();

    auto solver = border_solve::from_string(diffeq_str, *pconfig);
    auto result = solver.get_boundary_form();
    
    // match result from differential equation and from boundary
    GiNaC::ex eps_value = GiNaC::ex(1) / 10000;
    if (has_non_null_key(*pconfig, "ExpansionOption")) {
        auto expansion_opt = (*pconfig)["ExpansionOption"];
        if (has_non_null_key(expansion_opt, "TryEpsilon")) {
            eps_value = GiNaC::parser()(expansion_opt["TryEpsilon"].as<std::string>());
        }
    }
    d0 = 4;
    if (has_non_null_key(*pconfig, "Run")) {
        auto run_opt = (*pconfig)["Run"];
        if (has_non_null_key(run_opt, "D0"))
            d0 = GiNaC::parser()(run_opt["D0"].as<std::string>());
    }

    int n_masters = masters.size();
    for (auto& region: regions) {
        borders.push_back({});
        for (int i = 0; i < n_masters; i++) {
            auto numeric_exp = region.exponent(masters[i])
                .subs((*psymbols)["d"] == (d0 - 2 * eps_value), GiNaC::subs_options::algebraic);
            bool exist = false;
            int  border_one = 0;
            for (auto& presult: result) {
                auto diff = GiNaC::ex_to<GiNaC::numeric>(-presult.first - numeric_exp);
                if (GiNaC::is_integer(diff)) {
                    exist = true;
                    border_one = presult.second[i] - diff.to_int();
                    break;
                }
            }
            if (!exist) {
                std::cerr << "SingleSetup: region " << region.to_string() << " is forbidden by differential equation\n";
                exit(1);
            }
            borders.back().push_back(border_one);
        }
    }

    std::ofstream out(std::filesystem::path(workdirstr).append("BORDERS"));
    for (auto& border: borders) {
        out << "[" << border << "]\n";
    }
    out.close();
}


void integral_system::determine_direction() {
    master_sector  = top_sector(masters);
    auto prop_pres = pfamily->propagator_prescription();
    int  n_props   = master_sector.size();
    std::vector<int> interested_props;
    for (int i = 0; i < n_props; i++) {
        if (master_sector[i] && (pfamily->propagators[i].has((*psymbols)["eta"]))) {
            interested_props.push_back(i);
        }
    }
    int n_pos = 0, n_neg = 0, n_zero = 0, n_failed = 0;
    for (auto& prop: interested_props) {
        if (prop_pres[prop] == 1) {
            n_pos++;
        } else if (prop_pres[prop] == -1) {
            n_neg++;
        } else if (prop_pres[prop] == 0) {
            n_zero++;
        } else if (prop_pres[prop] == integralfamily::PRES_FAILED) {
            n_failed++;
        }
    }
    if ((n_failed != 0) || (n_pos != 0 && n_neg != 0)) {
        std::cerr << "SingleSetup: cannot define a self-consistent direction for continuation\n";
        exit(1);
    }

    std::ofstream out(std::filesystem::path(workdirstr).append("DIRECTION"));
    if (n_neg != 0) {
        out << "Im\n";
    } else {
        out << "NegIm\n";
    }
    out.close();
}


std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>> 
integral_system::setup_subfamilies() {
    auto subsysdir = std::filesystem::path(workdirstr).append("SUBSYSTEMS");
    std::filesystem::create_directory(subsysdir);
    auto boundaryconddir = std::filesystem::path(workdirstr).append("BOUNDARY_CONDITION");
    std::filesystem::create_directory(boundaryconddir);

    int nregions = regions.size(), familycnt = 0, startingcnt = 0, nmasters = masters.size();
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>> subsystem_ptrs;
    for (int r = 0; r < nregions; r++) {
        startingcnt = familycnt;
        auto subintgr = regions[r].subintegrals(master_sector, masters, borders[r]);
        for (auto& subf: subintgr.first) {
            auto cursubsysdir = std::filesystem::path(subsysdir).append(std::to_string(familycnt));
            std::filesystem::create_directory(cursubsysdir);
            std::vector<std::vector<int>> indices;
            for (auto& idx: subf.second)
                indices.push_back(idx.indices);
            integral_system subsystem(*subf.first, indices, {}, *pconfig, cursubsysdir.c_str());
            subsystem.reduce_targets();
            subsystem_ptrs.push_back({subf.first, subsystem});
            familycnt++;
        }

        std::ofstream out(std::filesystem::path(boundaryconddir).append(std::to_string(r)));
        for (int i = 0; i < nmasters; i++) {
            out << masters[i] << " = " << regions[r].factor(masters[i])
                    .subs((*psymbols)["d"] == (d0 - 2 * (*psymbols)["eps"]), GiNaC::subs_options::algebraic)
                << " * (\n";
            int order = subintgr.second[i].size();
            if (order == 0)
                out << "    0\n";
            else
                for (int o = 0; o < order; o++)
                    for (auto& p: subintgr.second[i][o])
                        out << "    + eta^(-" << o << ")*"
                            << "subintegral[" << (p.first.first + startingcnt) 
                            << ", "           << p.first.second << "]*("
                            << p.second       << ")\n";
            out << ");\n\n";
        }
        out.close();
    }
    return subsystem_ptrs;
}


std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
integral_system::setup_eta() {
    master_sector = top_sector(masters);
    auto newfamily = pfamily->insert_eta(master_sector, eta_mode);
    if (newfamily == nullptr) {
        std::cerr << "SingleSetup: the current system is an ending system\n";
        return {};
    }

    auto eta_dir = std::filesystem::path(workdirstr).append("ETA_SYSTEM");
    std::filesystem::create_directory(eta_dir);

    std::vector<std::vector<int>> newtargets;
    for (auto& mi: masters)
        newtargets.push_back(mi.indices);
    
    integral_system etasystem(*newfamily, newtargets, {}, *pconfig, eta_dir.c_str());
    etasystem.reduce_targets();
    return {{newfamily, etasystem}};
}


std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
integral_system::setup_singlemass_vacuum() {
    master_sector = top_sector(masters);
    // do nothing if current integral family is not a single-mass
    // vacuum type ending family
    if (!pfamily->is_singlemass_vacuum_ending(master_sector))
        return {};

    auto vacuum_subsysdir = std::filesystem::path(workdirstr).append("VACUUM_SUBSYSTEMS");
    std::filesystem::create_directory(vacuum_subsysdir);

    // factorize the integral family into connected components
    int n_master_prop = master_sector.size();
    bool failed = false;
    for (auto& mi: masters) {
        // we would like each master integral to be fully factorized into
        // integrals in its connected components
        for (int i = 0; i < n_master_prop; i++) {
            if (!master_sector[i] && mi.indices[i] != 0) {
                failed = true;
                break;
            }
        }
        if (failed)
            break;
    }
    if (failed) {
        std::cerr << "SingleSetup: trying to regenerate master integrals\n";
        // suggest a new list of preferred masters
        int sector_slots = 0;
        for (int i = 0; i < n_master_prop; i++)
            if (master_sector[i])
                sector_slots++;

        int max_pos_sum = 0;
        std::vector<std::vector<int>> target_indices;
        for (auto& ti: targets) {
            int ti_pos_sum = ti.dot() + sector_slots;
            if (ti_pos_sum > max_pos_sum)
                max_pos_sum = ti_pos_sum;
            target_indices.push_back(ti.indices);
        }
        for (auto& mi: masters) {
            int mi_pos_sum = mi.dot() + sector_slots;
            if (mi_pos_sum > max_pos_sum)
                max_pos_sum = mi_pos_sum;
        }

        auto new_preferred = split_at_most_N_into_M_slots(max_pos_sum, sector_slots);
        std::vector<std::vector<int>> suggested_preferred;
        for (auto& prfrrd: new_preferred) {
            std::vector<int> pad_prfrrd;
            int idx = 0;
            for (int i = 0; i < n_master_prop; i++) {
                if (master_sector[i])
                    pad_prfrrd.push_back(prfrrd[idx++]);
                else
                    pad_prfrrd.push_back(0);
            }
            suggested_preferred.push_back(pad_prfrrd);
        }

        integral_system newsystem(*pfamily, target_indices, suggested_preferred, *pconfig, workdirstr.c_str());
        newsystem.reduce_targets();
        // try again
        for (auto& mi: newsystem.masters) {
            for (int i = 0; i < n_master_prop; i++) {
                if (!master_sector[i] && mi.indices[i] != 0) {
                    std::cerr << "SingleSetup: master integral " << mi << " cannot be fully factorized into integrals in its components\n";
                    std::cerr << "Please consider another set of master integrals for this family\n";
                    exit(1);
                }
            }
        }
        masters = newsystem.masters;
    }
    
    pfamily->generate_components(master_sector);
    auto mass_list = pfamily->masses();
    auto loop_list = pfamily->loop_momenta();
    GiNaC::ex component_mass;
    int subsyscnt = 0;
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>> subsystem_ptrs;
    for (auto& comp: pfamily->all_components) {
        auto comp_props = comp.propagator_indices;
        // for every component, choose a new set of loop momenta
        // such that the propagator with non-zero mass has loop
        // momenta "l1"
        int nonzero_p = -1;
        for (auto& p: comp_props) {
            if (!(bool)(mass_list[p].expand() == 0)) {
                nonzero_p = p;
                component_mass = mass_list[p].expand();
                break;
            }
        }
        if (nonzero_p == -1) {
            std::cerr << "SingleSetup: detected a zero component\n";
            exit(1);
        }

        auto nonzero_loop = loop_list[nonzero_p];
        int  raw_loop_num = pfamily->loops.nops();
        GiNaC::matrix expand_nonzero_loop(raw_loop_num, 1);
        int norm_one_idx = -1, pos = 0;
        for (int i = 0; i < raw_loop_num; i++) {
            expand_nonzero_loop(i, 0) = nonzero_loop.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->loops[i]));
            if ((bool)(expand_nonzero_loop(i, 0).expand() == 1)) {
                norm_one_idx = i;
                pos = 1;
            } else if ((bool)(expand_nonzero_loop(i, 0).expand() == -1)) {
                norm_one_idx = i;
                pos = -1;
            }
        }

        GiNaC::symbol newl1("newl1");
        GiNaC::ex     newrule;
        if (pos == 1) {
            GiNaC::ex newl1_expand = newl1;
            for (int i = 0; i < raw_loop_num; i++)
                if (i != norm_one_idx && !(bool)(expand_nonzero_loop(i, 0).expand() == 0))
                    newl1_expand -= (expand_nonzero_loop(i, 0).expand() * pfamily->loops[i]);
            newrule = (pfamily->loops[norm_one_idx] == newl1_expand);
        } else if (pos == -1) {
            GiNaC::ex newl1_expand = -newl1;
            for (int i = 0; i < raw_loop_num; i++)
                if (i != norm_one_idx && !(bool)(expand_nonzero_loop(i, 0).expand() == 0))
                    newl1_expand += (expand_nonzero_loop(i, 0).expand() * pfamily->loops[i]);
            newrule = (pfamily->loops[norm_one_idx] == newl1_expand);
        }

        GiNaC::lst new_propagators;
        int new_p = 0, new_nonzero_p = 0;
        for (auto& p: comp_props) {
            new_propagators.append(pfamily->propagators[p].subs(newrule, GiNaC::subs_options::algebraic));
            if (p == nonzero_p)
                new_nonzero_p = new_p;
            new_p++;
        }

        GiNaC::lst new_loops;
        new_loops.append(newl1);
        for (auto& l: pfamily->loops) {
            for (auto& prop: new_propagators) {
                if (prop.has(l)) {
                    new_loops.append(l);
                    break;
                }
            }
        }

        integralfamily newf;
        newf.name               = pfamily->name;
        newf.psymbols           = psymbols;
        newf.loops              = new_loops;
        newf.indeplegs          = {};
        newf.propagators        = new_propagators;
        newf.cut                = {};
        newf.prescription       = {};
        newf.conservation_rules = {};
        newf.numeric_rules      = pfamily->numeric_rules;
        newf.sps_numeric_rules  = {};

        auto complete_info      = newf.sps_in_complete_propagator();
        newf.propagators        = complete_info.prop;

        // target integrals to reduce
        std::vector<std::vector<int>> subtargets;
        for (auto& mi: masters) {
            std::vector<int> atarget;
            for (auto& p: comp_props)
                atarget.push_back(mi.indices[p]);
            while (atarget.size() != newf.propagators.nops())
                atarget.push_back(0);
            
            bool exist = false;
            for (auto& existing: subtargets) {
                if (existing == atarget) {
                    exist = true;
                    break;
                }
            }
            if (!exist)
                subtargets.push_back(atarget);
        }
        auto cursubsysdir = std::filesystem::path(vacuum_subsysdir).append(std::to_string(subsyscnt));
        std::filesystem::create_directory(cursubsysdir);
        integral_system subsys(newf, subtargets, {}, *pconfig, cursubsysdir.c_str());
        subsys.reduce_targets();

        // further reduce the component to p-integral family
        GiNaC::lst pint_loops;
        for (auto& l: new_loops) {
            if (!(bool)((l - newl1).expand() == 0)) {
                pint_loops.append(l);
            }
        }

        GiNaC::lst pint_propagators;
        int remove_pos = 0;
        pos = 0;
        for (auto& prop: newf.propagators) {
            if (!(bool)((prop - new_propagators[new_nonzero_p]).expand() == 0)) {
                pint_propagators.append(prop);
            } else {
                remove_pos = pos;
            }
            pos++;
        }

        int       L = newf.loops.nops();
        GiNaC::ex D = (*psymbols)["d"];
        GiNaC::lst                    prefactor_list;
        std::vector<std::vector<int>> pint_targets;
        for (auto& mi: subsys.masters) {
            int mi_len = mi.indices.size();
            int nu = 0, nu1 = mi.indices[remove_pos];
            std::vector<int> newtarget;
            for (int i = 0; i < mi_len; i++) {
                nu += mi.indices[i];
                if (i != remove_pos)
                    newtarget.push_back(mi.indices[i]);
            }

            prefactor_list.append(GiNaC::pow(component_mass, L * D / 2 - nu)
                                * GiNaC::tgamma(nu - L * D / 2)
                                * GiNaC::tgamma(L * D / 2 - nu + nu1)
                                / GiNaC::pow(-1, nu1)
                                / GiNaC::tgamma(nu1)
                                / GiNaC::tgamma(D / 2));

            bool exist = false;
            for (auto& existing: pint_targets) {
                if (existing == newtarget) {
                    exist = true;
                    break;
                }
            }
            if (!exist)
                pint_targets.push_back(newtarget);
        }
        auto pintsysdir = std::filesystem::path(vacuum_subsysdir).append(std::to_string(subsyscnt) + "_reduce");
        std::filesystem::create_directory(pintsysdir);
        std::ofstream prefactor_out(std::filesystem::path(pintsysdir).append("PREFACTOR"));
        prefactor_out << prefactor_list << "\n";
        prefactor_out.close();

        if (L > 1) {
            std::shared_ptr<integralfamily> pintf = std::make_shared<integralfamily>();
            pintf->name               = pfamily->name;
            pintf->psymbols           = psymbols;
            pintf->loops              = pint_loops;
            pintf->indeplegs          = {newl1};
            pintf->propagators        = pint_propagators;
            pintf->cut                = {};
            pintf->prescription       = {};
            pintf->conservation_rules = {};
            pintf->numeric_rules      = pfamily->numeric_rules;
            pintf->sps_numeric_rules  = {newl1 * newl1 == -1};

            integral_system pint_sys(*pintf, pint_targets, {}, *pconfig, pintsysdir.c_str());
            pint_sys.reduce_targets();
            subsystem_ptrs.push_back({pintf, pint_sys});
        } else {
            std::shared_ptr<integralfamily> pintf = std::make_shared<integralfamily>();
            integral_system pint_sys(*pintf, {}, {}, *pconfig, pintsysdir.c_str());
            subsystem_ptrs.push_back({nullptr, pint_sys});
        }

        subsyscnt++;
    }

    return subsystem_ptrs;
}


std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
integral_system::setup_purely_phasespace() {
    master_sector = top_sector(masters);
    // do nothing if current integral family is not a purely
    // phase-space type ending family
    if (!pfamily->is_purely_phasespace_ending(master_sector))
        return {};

    int n_cut_loops = 0;
    for (auto& comp: pfamily->all_components) {
        if (comp.is_purely_phasespace()) {
            n_cut_loops = comp.num_loops();
            break;
        }
    }
    
    auto phase_subsysdir = std::filesystem::path(workdirstr).append("PHASE_SUBSYSTEMS");
    std::filesystem::create_directory(phase_subsysdir);

    std::shared_ptr<integralfamily> newf = std::make_shared<integralfamily>();
    newf->name               = pfamily->name;
    newf->psymbols           = psymbols;
    newf->loops              = pfamily->loops;
    newf->indeplegs          = pfamily->indeplegs;
    newf->propagators        = pfamily->propagators;
    newf->cut                = {};
    newf->prescription       = {};
    newf->conservation_rules = pfamily->conservation_rules;
    newf->numeric_rules      = pfamily->numeric_rules;
    newf->sps_numeric_rules  = pfamily->sps_numeric_rules;

    std::vector<std::vector<int>> newtargets;
    for (auto& mi: masters) {
        newtargets.push_back(mi.indices);
    }

    int       L = n_cut_loops;
    GiNaC::ex D = (*psymbols)["d"];
    GiNaC::ex prefactor = 2 * GiNaC::pow(GiNaC::Pi, L * D / 2) 
                            / GiNaC::pow(2 * GiNaC::Pi, L * D) 
                            * GiNaC::pow(-1, L + 1);

    std::ofstream prefactor_out(std::filesystem::path(phase_subsysdir).append("PREFACTOR"));
    prefactor_out << prefactor << "\n";
    prefactor_out << "after taking imaginary part\n";
    prefactor_out.close();

    integral_system newsys(*newf, newtargets, {}, *pconfig, phase_subsysdir.c_str());
    newsys.reduce_targets();
    return {{newf, newsys}};
}

