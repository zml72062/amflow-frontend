#include "singlesetup.hpp"
#include "utils.hpp"
#include "border.hpp"
#include <fstream>
#include <filesystem>


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
                if (GiNaC::is_nonneg_integer(diff)) {
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


void integral_system::setup_subfamilies() {
    auto subsysdir = std::filesystem::path(workdirstr).append("SUBSYSTEMS");
    std::filesystem::create_directory(subsysdir);
    auto boundaryconddir = std::filesystem::path(workdirstr).append("BOUNDARY_CONDITION");
    std::filesystem::create_directory(boundaryconddir);

    int nregions = regions.size(), familycnt = 0, startingcnt = 0, nmasters = masters.size();
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
}

