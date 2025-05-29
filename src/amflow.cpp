#include "amflow.hpp"
#include "reduction.hpp"
#include "boundary.hpp"
#include "singlesetup.hpp"
#include "utils.hpp"
#include "solve.hpp"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>


#define EVALUATE_IF_ZERO_ELSE(system_list, asystem, thing) do {         \
    if ((asystem).masters.size() == 0) {                                \
        std::ofstream out(std::filesystem::path((asystem).workdirstr)   \
                                        .append("RESULTS"));            \
        for (auto& target: (asystem).targets) {                         \
            out << target << " = 0\n";                                  \
        }                                                               \
        out.close();                                                    \
        (system_list).push_back({system_type::ZERO,                     \
                                 reduction_type::NO_REDUCE,             \
                                 (asystem), {}});                       \
    } else {                                                            \
        thing                                                           \
    }                                                                   \
} while (0)


#define DETERMINE_AND_PUSH(system_list, asystem, pfamily) do {          \
    auto master_sector = top_sector((asystem).masters);                 \
    if ((pfamily)->is_singlemass_vacuum_ending(master_sector)) {        \
        (system_list).push_back({system_type::VACUUM,                   \
                                 reduction_type::VACUUM_REDUCE,         \
                                 (asystem), {}});                       \
    } else if ((pfamily)->is_purely_phasespace_ending(master_sector)) { \
        (system_list).push_back({system_type::CUT,                      \
                                 reduction_type::CUT_REDUCE,            \
                                 (asystem), {}});                       \
    } else {                                                            \
        (system_list).push_back({system_type::NORMAL,                   \
                                 reduction_type::INSERT_ETA,            \
                                 (asystem), {}});                       \
    }                                                                   \
} while (0)


#define MASTER_TO_TARGETS(systemid, master_values) do {                 \
    std::string cursys_ibp_string = std::filesystem::path(base_workdir) \
                                      .append(std::to_string(systemid)) \
                                      .append("REDUCE");                \
    std::ifstream cursys_ibp_file(cursys_ibp_string);                   \
    std::string target_str, master_str, coeff_str;                      \
    cursys_ibp_file >> target_str;                                      \
    cursys_ibp_file >> target_str;                                      \
    cursys_ibp_file >> master_str;                                      \
    cursys_ibp_file >> master_str;                                      \
    cursys_ibp_file >> coeff_str;                                       \
    cursys_ibp_file >> coeff_str;                                       \
    cursys_ibp_file.close();                                            \
    std::vector<std::vector<int>> target_list =                         \
        string_to_integral_list(target_str);                            \
    GiNaC::matrix coeffmatrix  = read_matrix(coeff_str, globals);       \
    int nmaster_values = master_values.size();                          \
    GiNaC::matrix mastermatrix(nmaster_values, 1);                      \
    for (int q = 0; q < nmaster_values; q++) {                          \
        mastermatrix(q, 0) = master_values[q];                          \
    }                                                                   \
    auto targetmatrix = coeffmatrix.mul(mastermatrix);                  \
    int ntarget_values = target_list.size();                            \
    std::string cursys_result_string =                                  \
        std::filesystem::path(base_workdir)                             \
            .append(std::to_string(systemid)).append("RESULTS");        \
    std::ofstream cursys_result_file(cursys_result_string);             \
    for (int q = 0; q < ntarget_values; q++) {                          \
        cursys_result_file << family.name << "["                        \
                           << target_list[q] << "] = "                  \
                           << GiNaC::ex_to<GiNaC::numeric>(             \
                            targetmatrix(q, 0)                          \
                            .subs(globals["eta"] == 0,                  \
                                  GiNaC::subs_options::algebraic)       \
                            .subs(eps_subs_rules,                       \
                                  GiNaC::subs_options::algebraic)       \
                           ).evalf() << "\n";                           \
    }                                                                   \
    cursys_result_file.close();                                         \
} while (0)


static std::vector<int> string_to_integral(const std::string& input) {
    std::size_t start = input.find('[') + 1,
                end   = input.find(']');
    std::string id = input.substr(start, end - start);
    std::vector<int> nums;
    const char* ptr = id.c_str();
    char* endptr;
    while (true) {
        nums.push_back(std::strtol(ptr, &endptr, 10));
        if (*endptr == '\0')
            return nums;
        ptr = endptr + 1;
    }
}


static std::vector<std::vector<int>> string_to_integral_list(const std::string& input) {
    std::size_t startptr = 0;
    std::vector<std::vector<int>> results;

    while (true) {
        if (input.find('[', startptr) == std::string::npos)
            break;
        std::size_t start = input.find('[', startptr) + 1,
                    end   = input.find(']', startptr);
        std::string id = input.substr(start, end - start);
        std::vector<int> nums;
        const char* ptr = id.c_str();
        char* endptr;
        while (true) {
            nums.push_back(std::strtol(ptr, &endptr, 10));
            if (*endptr == '\0') {
                results.push_back(nums);
                break;
            }
            ptr = endptr + 1;
        }
        startptr = end + 1;
    }

    return results;
}


static GiNaC::lst read_list(const std::string& str, const GiNaC::symtab& tab) {
    GiNaC::parser parser(tab);
    return GiNaC::ex_to<GiNaC::lst>(parser(str));
}


static GiNaC::matrix read_matrix(const std::string& str, const GiNaC::symtab& tab) {
    GiNaC::parser parser(tab);
    auto lst_out = parser(str);
    int rows = lst_out.nops();
    if (rows == 0)
        return GiNaC::matrix(0, 0);
    int cols = lst_out[0].nops();
    GiNaC::matrix final_out(rows, cols);
    int i = 0, j;
    for (auto r = lst_out.begin(); r != lst_out.end(); ++r, i++) {
        j = 0;
        for (auto c = r->begin(); c != r->end(); ++c, j++)
            final_out(i, j) = *c;
    }
    return final_out;
}


amflow::amflow(const char* config_path,
               const std::vector<std::vector<int>>& initial_targets,
               const std::vector<std::vector<int>>& initial_preferred,
               const char* work_dir) {
    START_TIME(initialize);
    startup();

    config = YAML::LoadFile(config_path);
    checkdep();

    base_workdir = work_dir;
    std::filesystem::create_directory(std::filesystem::path(base_workdir));

    globals["d"]   = GiNaC::symbol("d");        // space-time dimension
    globals["eps"] = GiNaC::symbol("eps");      // d = d0 - 2 * eps
    globals["eta"] = GiNaC::symbol("eta");      // auxiliary mass

    d0 = 4;
    if (has_non_null_key(config, "Run")) {
        auto run_opt = config["Run"];
        if (has_non_null_key(run_opt, "D0"))
            d0 = GiNaC::parser()(run_opt["D0"].as<std::string>());
    }

    family  = integralfamily::from_yaml(config, &globals);
    nextnum = 0;
    enroute.push(nextnum);
    END_TIME(initialize);
    PRINT_TIME(initialize);
    
    START_TIME(initial_reduction);
    if (initial_targets.size() == 0) {
        std::cerr << "AMFlow: no input target integrals\n";
        exit(1);
    }
    std::string     initial_sysdir = std::filesystem::path(base_workdir).append(std::to_string(nextnum));
    integral_system initial_system(family, initial_targets, initial_preferred, config, initial_sysdir.c_str());
    initial_system.reduce_targets();
    
    if (initial_system.masters.size() == 0) {
        std::cerr << "AMFlow: the input system is trivial\n";
        std::ofstream out(std::filesystem::path(initial_sysdir).append("RESULTS"));
        for (auto& target: initial_system.targets) {
            out << target << " = 0\n";
        }
        out.close();
        systems.push_back({system_type::ZERO,
                           reduction_type::NO_REDUCE,
                           initial_system, {}});
    } else {
        DETERMINE_AND_PUSH(systems, initial_system, &family);
    }
    nextnum++;
    report_forward_status();
    END_TIME(initial_reduction);
    PRINT_TIME(initial_reduction);
}


void amflow::forward() {
    std::cerr << "AMFlow: starting forward pass\n";
    START_TIME(forward_pass);
    
    while (!enroute.empty()) {
        int cursysnum = enroute.front();
        enroute.pop();
        std::cerr << "AMFlow: reducing system " << cursysnum << "\n";

        auto transtype = systems[cursysnum].transtype;
        if (transtype == reduction_type::NO_REDUCE)
            continue;
        
        switch (transtype) {
            case reduction_type::INSERT_ETA: {
                auto subsystem_ptrs = systems[cursysnum].system.setup_eta();
                if (subsystem_ptrs.size() == 0) {
                    std::cerr << "AMFlow: unexpected ending system, please try setting a different \"AMFMode\", especially try including \"Propagator\" or \"All\"\n";
                    exit(1);
                }

                std::string old_dir = std::filesystem::path(base_workdir).append(std::to_string(cursysnum)).append("ETA_SYSTEM"),
                            new_dir = std::filesystem::path(base_workdir).append(std::to_string(nextnum));
                std::filesystem::create_directory(std::filesystem::path(new_dir));
                std::filesystem::copy(std::filesystem::path(old_dir).append("REDUCE"),
                                      std::filesystem::path(new_dir).append("REDUCE"));

                subfamilies.push_back(subsystem_ptrs[0].first);
                auto subsystem = subsystem_ptrs[0].second;
                subsystem.workdirstr = new_dir;
                subsystem.tempdirstr = std::filesystem::path(new_dir).append("KIRATEMP");
                systems.push_back({system_type::ETA, 
                                   reduction_type::BOUNDARY_REDUCE,
                                   subsystem, {}});
                systems[cursysnum].subsystems.push_back(nextnum);
                enroute.push(nextnum);
                nextnum++;
                break;
            }
            case reduction_type::BOUNDARY_REDUCE: {
                auto cursys = systems[cursysnum].system;
                cursys.build_diffeq();
                cursys.determine_boundaries();
                cursys.determine_border();
                cursys.determine_direction();
                auto subsystem_ptrs = cursys.setup_subfamilies();
                
                int nsubsystems = subsystem_ptrs.size();
                for (int i = 0; i < nsubsystems; i++) {
                    std::string old_dir = std::filesystem::path(base_workdir).append(std::to_string(cursysnum))
                                                                             .append("SUBSYSTEMS")
                                                                             .append(std::to_string(i)),
                                new_dir = std::filesystem::path(base_workdir).append(std::to_string(nextnum));
                    std::filesystem::create_directory(std::filesystem::path(new_dir));
                    std::filesystem::copy(std::filesystem::path(old_dir).append("REDUCE"),
                                          std::filesystem::path(new_dir).append("REDUCE"));
                    
                    subfamilies.push_back(subsystem_ptrs[i].first);
                    auto subsystem = subsystem_ptrs[i].second;
                    subsystem.workdirstr = new_dir;
                    subsystem.tempdirstr = std::filesystem::path(new_dir).append("KIRATEMP");
                    EVALUATE_IF_ZERO_ELSE(systems, subsystem, { 
                        DETERMINE_AND_PUSH(systems, subsystem, subsystem_ptrs[i].first); 
                    });
                    systems[cursysnum].subsystems.push_back(nextnum);
                    enroute.push(nextnum);
                    nextnum++;
                }
                break;
            }
            case reduction_type::VACUUM_REDUCE: {
                auto subsystem_ptrs = systems[cursysnum].system.setup_singlemass_vacuum();

                int nsubsystems = subsystem_ptrs.size();
                for (int i = 0; i < nsubsystems; i++) {
                    std::string old_dir = std::filesystem::path(base_workdir).append(std::to_string(cursysnum))
                                                                             .append("VACUUM_SUBSYSTEMS")
                                                                             .append(std::to_string(i) + "_reduce"),
                                new_dir = std::filesystem::path(base_workdir).append(std::to_string(nextnum));
                    std::filesystem::create_directory(std::filesystem::path(new_dir));
                    if (std::filesystem::exists(std::filesystem::path(old_dir).append("REDUCE"))) {
                        std::filesystem::copy(std::filesystem::path(old_dir).append("REDUCE"),
                                              std::filesystem::path(new_dir).append("REDUCE"));
                    }

                    auto subsystem = subsystem_ptrs[i].second;
                    subsystem.workdirstr = new_dir;
                    subsystem.tempdirstr = std::filesystem::path(new_dir).append("KIRATEMP");
                    if (subsystem_ptrs[i].first != nullptr) {
                        subfamilies.push_back(subsystem_ptrs[i].first);
                        EVALUATE_IF_ZERO_ELSE(systems, subsystem, { 
                            DETERMINE_AND_PUSH(systems, subsystem, subsystem_ptrs[i].first); 
                        });
                    } else {
                        subsystem.pfamily = nullptr;
                        systems.push_back({system_type::TRIVIAL,
                                           reduction_type::NO_REDUCE,
                                           subsystem, {}});
                    }

                    systems[cursysnum].subsystems.push_back(nextnum);
                    enroute.push(nextnum);
                    nextnum++;
                }
                break;
            }
            case reduction_type::CUT_REDUCE: {
                auto subsystem_ptrs = systems[cursysnum].system.setup_purely_phasespace();
                if (subsystem_ptrs.size() != 0) {
                    std::string old_dir = std::filesystem::path(base_workdir).append(std::to_string(cursysnum)).append("PHASE_SUBSYSTEMS"),
                                new_dir = std::filesystem::path(base_workdir).append(std::to_string(nextnum));
                    std::filesystem::create_directory(std::filesystem::path(new_dir));
                    std::filesystem::copy(std::filesystem::path(old_dir).append("REDUCE"),
                                          std::filesystem::path(new_dir).append("REDUCE"));

                    subfamilies.push_back(subsystem_ptrs[0].first);
                    auto subsystem = subsystem_ptrs[0].second;
                    subsystem.workdirstr = new_dir;
                    subsystem.tempdirstr = std::filesystem::path(new_dir).append("KIRATEMP");
                    EVALUATE_IF_ZERO_ELSE(systems, subsystem, { 
                        DETERMINE_AND_PUSH(systems, subsystem, subsystem_ptrs[0].first); 
                    });
                    systems[cursysnum].subsystems.push_back(nextnum);
                    enroute.push(nextnum);
                    nextnum++;
                }
                break;
            }
            default:
                break;
        }
        report_forward_status();
    }

    END_TIME(forward_pass);
    PRINT_TIME(forward_pass);
}


void amflow::backward(const GiNaC::numeric& eps_value) {
    START_TIME(backward_pass);
    int digits = 100;
    if (has_non_null_key(config, "ExpansionOption")) {
        auto expansion_opt = config["ExpansionOption"];
        if (has_non_null_key(expansion_opt, "WorkingPre")) {
            digits = expansion_opt["WorkingPre"].as<int>();
        }
    }
    GiNaC::Digits = digits;
    int nsystems = systems.size();
    GiNaC::lst eps_subs_rules;
    eps_subs_rules.append(globals["eps"] == eps_value);
    eps_subs_rules.append(globals["d"] == d0 - 2 * eps_value);

    for (int vid = 0; vid < nsystems; vid++) {
        int id = nsystems - vid - 1;
        std::cerr << "AMFlow: solving system " << id << "\n";

        // from subsystem target integrals to 
        // parent system master integrals
        auto cursystype = systems[id].systype;
        switch (cursystype) {
            case system_type::ZERO:
            case system_type::TRIVIAL:
            case system_type::ETA:
                break;
            case system_type::VACUUM: {
                std::vector<GiNaC::ex> master_values;
                for (auto& mi: systems[id].system.masters) {
                    // decompose "mi" as a product of target
                    // integrals in subcomponents
                    std::vector<int>              subsysnums;
                    std::vector<std::vector<int>> subtargets;
                    int subsysidx = 0;
                    for (auto& comp: systems[id].system.pfamily->all_components) {
                        int subsysnum = systems[id].subsystems[subsysidx++];
                        subsysnums.push_back(subsysnum);

                        std::vector<int> asubtarget;
                        for (auto& p: comp.propagator_indices)
                            asubtarget.push_back(mi.indices[p]);
                        int subsystem_nprops;
                        if (systems[subsysnum].systype == system_type::TRIVIAL
                         || systems[subsysnum].systype == system_type::ZERO)
                            subsystem_nprops = 1;
                        else
                            subsystem_nprops = systems[subsysnum].system.pfamily->propagators.nops() + 1;
                        while ((int)asubtarget.size() < subsystem_nprops)
                            asubtarget.push_back(0);
                        subtargets.push_back(asubtarget);
                    }

                    // write every subtarget as a linear
                    // combination of submasters
                    int nsubsystems = subsysnums.size();
                    // subtarget[i] = sum(submasters[i] * subcoeffs[i])
                    std::vector<std::vector<std::vector<int>>> submasters;
                    std::vector<std::vector<GiNaC::ex>>        subcoeffs;
                    for (int i = 0; i < nsubsystems; i++) {
                        auto subsystype = systems[subsysnums[i]].systype;
                        if (subsystype == system_type::ZERO) {
                            submasters.push_back({});
                            subcoeffs.push_back({});
                            continue;
                        }
                        std::string reduce_filename = 
                            std::filesystem::path(base_workdir)
                                .append(std::to_string(id))
                                .append("VACUUM_SUBSYSTEMS")
                                .append(std::to_string(i))
                                .append("REDUCE");
                        std::ifstream reduce_file(reduce_filename);
                        std::string targetstr, masterstr, coeffstr;
                        reduce_file >> targetstr;
                        reduce_file >> targetstr;
                        reduce_file >> masterstr;
                        reduce_file >> masterstr;
                        reduce_file >> coeffstr;
                        reduce_file >> coeffstr;
                        reduce_file.close();

                        std::vector<std::vector<int>> 
                            targetlist = string_to_integral_list(targetstr),
                            masterlist = string_to_integral_list(masterstr);
                        GiNaC::matrix coeffmatrix = read_matrix(coeffstr, globals);
                        submasters.push_back(masterlist);
                        int len_targetlist = targetlist.size();
                        int target_id = -1;
                        for (int q = 0; q < len_targetlist; q++) {
                            if (targetlist[q] == subtargets[i]) {
                                target_id = q;
                                break;
                            }
                        }
                        std::vector<GiNaC::ex> sub_coefflist;
                        int nc = coeffmatrix.cols();
                        for (int q = 0; q < nc; q++) {
                            sub_coefflist.push_back(coeffmatrix(target_id, q));
                        }
                        subcoeffs.push_back(sub_coefflist);
                    }

                    // write submasters in terms of p-integrals
                    // submasters[i][j] = subsubtargets[i][j] * prefactors[i][j]
                    std::vector<std::vector<std::vector<int>>> subsubtargets;
                    std::vector<std::vector<GiNaC::ex>>        prefactors;
                    for (int i = 0; i < nsubsystems; i++) {
                        auto subsystype = systems[subsysnums[i]].systype;
                        if (subsystype == system_type::ZERO) {
                            subsubtargets.push_back({});
                            prefactors.push_back({});
                            continue;
                        }
                        std::vector<std::vector<int>> subsubtargetlist;
                        std::vector<GiNaC::ex>        subprefactors;
                        std::string nonzero_pos_string = 
                            std::filesystem::path(base_workdir)
                                .append(std::to_string(id))
                                .append("VACUUM_SUBSYSTEMS")
                                .append(std::to_string(i))
                                .append("MASSIVE_PROPAGATOR");
                        std::ifstream nonzero_pos_file(nonzero_pos_string);
                        int nonzero_pos = -1;
                        nonzero_pos_file >> nonzero_pos;
                        nonzero_pos_file.close();
                        std::string prefactor_path_string = 
                            std::filesystem::path(base_workdir)
                                .append(std::to_string(id))
                                .append("VACUUM_SUBSYSTEMS")
                                .append(std::to_string(i) + "_reduce")
                                .append("PREFACTOR");
                        std::ifstream prefactor_file(prefactor_path_string);
                        std::string prefactor_string;
                        prefactor_file >> prefactor_string;
                        prefactor_file.close();
                        GiNaC::lst prefactor_list = read_list(prefactor_string, globals);
                        int nsubmasters = submasters[i].size();
                        for (int p = 0; p < nsubmasters; p++) {
                            auto submi = submasters[i][p];
                            std::vector<int> newtarget;
                            int submi_len = submi.size();
                            for (int q = 0; q < submi_len; q++) {
                                if (q != nonzero_pos)
                                    newtarget.push_back(submi[q]);
                            }
                            subsubtargetlist.push_back(newtarget);
                            subprefactors.push_back(prefactor_list[p]);
                        }
                        subsubtargets.push_back(subsubtargetlist);
                        prefactors.push_back(subprefactors);
                    }

                    // read results from subsystems
                    std::vector<std::vector<GiNaC::ex>> submaster_values;
                    for (int i = 0; i < nsubsystems; i++) {
                        auto subsystype = systems[subsysnums[i]].systype;
                        if (subsystype == system_type::ZERO) {
                            submaster_values.push_back({});
                            continue;
                        }
                        std::vector<GiNaC::ex> submaster_subvalues;
                        int n_submasters = submasters[i].size();
                        for (int j = 0; j < n_submasters; j++) {
                            auto subresult    = get_subsystem_result(subsysnums[i], subsubtargets[i][j]);
                            auto subprefactor = prefactors[i][j].subs(eps_subs_rules, GiNaC::subs_options::algebraic);
                            submaster_subvalues.push_back(
                                GiNaC::ex_to<GiNaC::numeric>(subprefactor * subresult).evalf()
                            );
                        }
                        submaster_values.push_back(submaster_subvalues);
                    }

                    std::vector<GiNaC::ex> subtarget_values;
                    for (int i = 0; i < nsubsystems; i++) {
                        auto subsystype = systems[subsysnums[i]].systype;
                        if (subsystype == system_type::ZERO) {
                            subtarget_values.push_back(0);
                            continue;
                        }
                        GiNaC::ex value = 0;
                        int n_submasters = submaster_values[i].size();
                        for (int j = 0; j < n_submasters; j++) {
                            value += (submaster_values[i][j] * subcoeffs[i][j]);
                        }
                        subtarget_values.push_back(
                            GiNaC::ex_to<GiNaC::numeric>(
                                value.subs(eps_subs_rules, 
                                           GiNaC::subs_options::algebraic)
                            ).evalf()
                        );
                    }

                    GiNaC::ex mvalue = 1;
                    for (auto& sub: subtarget_values)
                        mvalue *= sub;

                    master_values.push_back(GiNaC::ex_to<GiNaC::numeric>(mvalue).evalf());
                }
                MASTER_TO_TARGETS(id, master_values);
                break;
            }
            case system_type::CUT: {
                std::string prefactor_path_string = std::filesystem::path(base_workdir)
                                .append(std::to_string(id))
                                .append("PHASE_SUBSYSTEMS")
                                .append("PREFACTOR");
                std::ifstream prefactor_file(prefactor_path_string);
                std::string prefactor_string;
                prefactor_file >> prefactor_string;
                prefactor_file.close();
                
                GiNaC::ex prefactor = GiNaC::parser(globals)(prefactor_string)
                    .subs(eps_subs_rules, GiNaC::subs_options::algebraic);
                std::vector<GiNaC::ex> master_values;
                for (auto& mi: systems[id].system.masters) {
                    master_values.push_back(
                        GiNaC::ex_to<GiNaC::numeric>(
                            get_subsystem_result(systems[id].subsystems[0], mi.indices)
                                .imag_part() * prefactor
                        ).evalf()
                    );
                }
                MASTER_TO_TARGETS(id, master_values);
                break;
            }
            case system_type::NORMAL: {
                // collect boundary condition for every region
                int eta_id = systems[id].subsystems[0];
                // the most singular exponent of a region
                std::vector<GiNaC::ex> region_exponents;
                // region_expansion[r][i][o] means:
                // in region "r", the "o"-th order of master integral "i"
                std::vector<std::vector<GiNaC::lst>> region_expansion;
                int r = 0;
                std::string boundary_dir = std::filesystem::path(base_workdir)
                                            .append(std::to_string(eta_id))
                                            .append("BOUNDARY_CONDITION");
                // use a symbol table to store subintegral values
                int n_boundary_families = systems[eta_id].subsystems.size();
                GiNaC::symtab subintegral_values;
                for (int i = 0; i < n_boundary_families; i++) {
                    int subsysnum   = systems[eta_id].subsystems[i];
                    int nsubtargets = systems[subsysnum].system.targets.size();
                    for (int j = 0; j < nsubtargets; j++) {
                        subintegral_values["subfamily" + std::to_string(i)
                                         + "integral"  + std::to_string(j)]
                            = get_subsystem_result(subsysnum, 
                                systems[subsysnum].system.targets[j].indices);
                    }
                }
                // copy existing key-value pairs to the symbol table
                for (auto& p: globals) {
                    subintegral_values[p.first] = p.second;
                }
                // read spurious orders
                std::vector<GiNaC::lst> spurious_orders;
                std::string spuriouspath = std::filesystem::path(base_workdir)
                                            .append(std::to_string(eta_id))
                                            .append("SPURIOUS_ORDERS");
                std::ifstream spuriousfile(spuriouspath);
                std::string spuriousarray;
                while (true) {
                    spuriousfile >> spuriousarray;
                    if (spuriousfile.eof())
                        break;
                    spurious_orders.push_back(read_list(spuriousarray, {}));
                }
                spuriousfile.close();
                while (true) {
                    if (!std::filesystem::exists(std::filesystem::path(boundary_dir).append(std::to_string(r))))
                        break;
                    
                    std::string curbounddir = std::filesystem::path(boundary_dir).append(std::to_string(r) + "_machine");
                    
                    // read exponent list
                    GiNaC::lst exponentlist;
                    std::string exponentpath = std::filesystem::path(curbounddir).append("EXPONENT");
                    std::string exponentstr;
                    std::ifstream exponentfile(exponentpath);
                    exponentfile >> exponentstr;
                    exponentfile.close();
                    exponentlist = read_list(exponentstr, globals);
                    
                    // deal with spurious orders
                    int n_exponents = exponentlist.nops();
                    for (int j = 0; j < n_exponents; j++) {
                        exponentlist[j] -= spurious_orders[r][j];
                    }

                    // find most singular order
                    GiNaC::ex most_singular = exponentlist[0];
                    for (int j = 1; j < n_exponents; j++) {
                        int highj = GiNaC::ex_to<GiNaC::numeric>(
                            (exponentlist[j] - most_singular).expand()
                        ).to_int();
                        if (highj > 0)
                            most_singular = exponentlist[j];
                    }
                    region_exponents.push_back(most_singular);

                    // find "gratis" zero orders
                    std::vector<int> gratis_orders;
                    for (int j = 0; j < n_exponents; j++)
                        gratis_orders.push_back(
                            GiNaC::ex_to<GiNaC::numeric>(
                                (most_singular - exponentlist[j]).expand()
                            ).to_int()
                        );
                    
                    // read expansion terms
                    std::string expansionpath = std::filesystem::path(curbounddir).append("EXPANSION");
                    std::string expansionstr;
                    std::ifstream expansionfile(expansionpath);
                    expansionfile >> expansionstr;
                    expansionfile.close();
                    GiNaC::lst expansionlist = read_list(expansionstr, subintegral_values);

                    // read jacobian
                    std::string jacobianpath = std::filesystem::path(curbounddir).append("JACOBIAN");
                    std::string jacobianstr;
                    std::ifstream jacobianfile(jacobianpath);
                    jacobianfile >> jacobianstr;
                    jacobianfile.close();
                    auto jacobian = GiNaC::parser(globals)(jacobianstr);

                    // build final expansion
                    std::vector<GiNaC::lst> final_expansion;
                    for (int j = 0; j < n_exponents; j++) {
                        GiNaC::lst final_expj;
                        for (int q = 0; q < gratis_orders[j]; q++)
                            final_expj.append(0);

                        int given_terms = expansionlist[j].nops();
                        int leading_useless_terms = 
                            GiNaC::ex_to<GiNaC::numeric>(
                                spurious_orders[r][j]
                            ).to_int();
                        for (int q = leading_useless_terms; q < given_terms; q++)
                            final_expj.append(GiNaC::ex_to<GiNaC::numeric>(
                                (expansionlist[j][q] * jacobian)
                                    .subs(eps_subs_rules, GiNaC::subs_options::algebraic)
                                ).evalf()
                            );
                        
                        final_expansion.push_back(final_expj);
                    }
                    region_expansion.push_back(final_expansion);
                    r++;
                }
                // gather information from different regions
                std::vector<GiNaC::ex>        classified_exponents;
                std::vector<std::vector<int>> classified_indices;
                std::vector<std::vector<int>> classified_gratis;
                int nregions = region_exponents.size();
                std::vector<int>              classified(nregions, 0);
                for (int i = 0; i < nregions; i++) {
                    if (classified[i])
                        continue;
                    GiNaC::ex representative_exponent = region_exponents[i];
                    classified[i] = 1;
                    classified_indices.push_back({i});
                    for (int j = i + 1; j < nregions; j++) {
                        if (classified[j])
                            continue;
                        auto diff = GiNaC::ex_to<GiNaC::numeric>(
                            (region_exponents[j] - representative_exponent)
                                .subs(eps_subs_rules, GiNaC::subs_options::algebraic)
                                .expand()
                        );
                        if (GiNaC::is_integer(diff)) {
                            classified[j] = 1;
                            classified_indices.back().push_back(j);
                            if (diff.to_int() > 0)
                                representative_exponent = region_exponents[j];
                        }
                    }
                    classified_exponents.push_back(representative_exponent);
                    classified_gratis.push_back({});
                    for (auto& ind: classified_indices.back()) {
                        classified_gratis.back().push_back(
                            GiNaC::ex_to<GiNaC::numeric>(
                                (representative_exponent - region_exponents[ind])
                                    .subs(eps_subs_rules, GiNaC::subs_options::algebraic)
                                    .expand()
                        ).to_int());
                    }
                }
                int nclasses = classified_exponents.size();
                for (int i = 0; i < nclasses; i++) {
                    int ngratis = classified_gratis[i].size();
                    for (int j = 0; j < ngratis; j++) {
                        int grat = classified_gratis[i][j];
                        int ind  = classified_indices[i][j];
                        if (grat > 0) {
                            int nmas = region_expansion[ind].size();
                            for (int k = 0; k < nmas; k++) {
                                auto raw_list = region_expansion[ind][k];
                                GiNaC::lst newlist;
                                for (int q = 0; q < grat; q++)
                                    newlist.append(0);
                                for (auto& raw: raw_list)
                                    newlist.append(raw);
                                region_expansion[ind][k] = newlist;
                            }
                        }
                    }
                }

                // systematically build boundary conditions
                std::vector<boundary_condition> conditions;
                for (int i = 0; i < nclasses; i++) {
                    boundary_condition bc;
                    bc.exponent = classified_exponents[i].subs(globals["d"] == d0 - 2 * globals["eps"],
                                                               GiNaC::subs_options::algebraic);
                    int ngratis = classified_gratis[i].size();
                    int nmas = region_expansion[classified_indices[i][0]].size();
                    bc.expansion = std::vector<GiNaC::lst>(nmas, GiNaC::lst());
                    for (int j = 0; j < ngratis; j++) {
                        int ind = classified_indices[i][j];
                        for (int k = 0; k < nmas; k++) {
                            auto onepart = region_expansion[ind][k];
                            int  curlen = bc.expansion[k].nops(),
                                 newlen = onepart.nops();
                            if (curlen >= newlen) {
                                for (int q = 0; q < newlen; q++) {
                                    bc.expansion[k][q] += onepart[q];
                                }
                            } else {
                                for (int q = 0; q < curlen; q++) {
                                    bc.expansion[k][q] += onepart[q];
                                }
                                for (int q = curlen; q < newlen; q++) {
                                    bc.expansion[k].append(onepart[q]);
                                }
                            }
                        }
                    }
                    conditions.push_back(bc);
                }

                // read differential equation
                std::string de_path = std::filesystem::path(base_workdir)
                                        .append(std::to_string(eta_id))
                                        .append("DIFFEQ");
                std::ifstream de_file(de_path);
                std::string de_string;
                de_file >> de_string; // read
                de_file >> de_string; // four
                de_file >> de_string; // times
                de_file >> de_string;
                de_file.close();

                // read direction
                std::string dir_path = std::filesystem::path(base_workdir)
                                        .append(std::to_string(eta_id))
                                        .append("DIRECTION");
                std::ifstream dir_file(dir_path);
                std::string dir_string;
                dir_file >> dir_string;
                dir_file.close();
                
                // solve master integrals at zero
                auto solver = run_solve::from_string(de_string, eps_value, dir_string, config);
                auto masters = solver.link(conditions);
                int  nmasters = masters.rows();
                std::vector<GiNaC::ex> master_values;
                for (int i = 0; i < nmasters; i++)
                    master_values.push_back(masters(i, 0));
                MASTER_TO_TARGETS(eta_id, master_values);

                // solve target integrals of original system
                std::vector<GiNaC::ex> original_master_values;
                for (auto& mi: systems[id].system.masters) {
                    original_master_values.push_back(
                        get_subsystem_result(eta_id, mi.indices)
                    );
                }
                MASTER_TO_TARGETS(id, original_master_values);
                break;
            }
        }
    }
    END_TIME(backward_pass);
    PRINT_TIME(backward_pass);
}


void amflow::report_forward_status() {
    std::ofstream out(std::filesystem::path(base_workdir).append("STATUS"));
    out << "  SYSTEM #  |  SYSTEM TYPE  |  TRANSITION TYPE  |          SUBSYSTEMS          \n"
        << "------------|---------------|-------------------|------------------------------\n";
    int nsystems = systems.size();
    for (int i = 0; i < nsystems; i++) {
        out << "  " << std::left << std::setw(10) << i << "|";
        switch (systems[i].systype) {
            case system_type::ZERO:
                out << "     ZERO      |";
                break;
            case system_type::TRIVIAL:
                out << "    TRIVIAL    |";
                break;
            case system_type::NORMAL:
                out << "    NORMAL     |";
                break;
            case system_type::ETA:
                out << "      ETA      |";
                break;
            case system_type::VACUUM:
                out << "    VACUUM     |";
                break;
            case system_type::CUT:
                out << "      CUT      |";
                break;
        }
        switch (systems[i].transtype) {
            case reduction_type::NO_REDUCE:
                out << "     NO_REDUCE     |";
                break;
            case reduction_type::INSERT_ETA:
                out << "    INSERT_ETA     |";
                break;
            case reduction_type::BOUNDARY_REDUCE:
                out << "  BOUNDARY_REDUCE  |";
                break;
            case reduction_type::VACUUM_REDUCE:
                out << "   VACUUM_REDUCE   |";
                break;
            case reduction_type::CUT_REDUCE:
                out << "    CUT_REDUCE     |";
                break;
        }
        out << "  " << systems[i].subsystems << "\n";
    }
    out << "------------|---------------|-------------------|------------------------------\n";
    out.close();
}


GiNaC::ex amflow::get_subsystem_result(int subsysnum, const std::vector<int>& indices) {
    if (systems[subsysnum].systype == system_type::ZERO)
        return 0;
    if (systems[subsysnum].systype == system_type::TRIVIAL)
        return 1;
    if (subsystem_results.find(subsysnum) != subsystem_results.end())
        return subsystem_results[subsysnum][indices];

    std::string result_path = std::filesystem::path(base_workdir)
                                            .append(std::to_string(subsysnum))
                                            .append("RESULTS");
    if (!std::filesystem::exists(result_path)) {
        std::cerr << "AMFlow: subsystem " << subsysnum << " not solved yet\n";
        exit(1);
    }

    subsystem_results[subsysnum] = std::map<std::vector<int>, GiNaC::ex>();

    std::ifstream result_file(result_path);
    std::string lhs, equals, rhs;
    while (true) {
        result_file >> lhs;
        result_file >> equals;
        result_file >> rhs;
        if (result_file.eof())
            break;
        
        auto aind    = string_to_integral(lhs);
        auto aresult = GiNaC::parser()(rhs);
        subsystem_results[subsysnum][aind] = aresult;
    }
    result_file.close();

    return subsystem_results[subsysnum][indices];
}


void amflow::startup() {
    std::cerr << "--------------------------------------------------\n";
    std::cerr << "AMFlow: package loaded\n";
    std::cerr << "Author: "       << author  << "\n";
    std::cerr << "Version: "      << version << "\n";
    std::cerr << "Release date: " << release << "\n";
    std::cerr << "--------------------------------------------------\n";
}


void amflow::usage() {
    std::cerr << "Usage: amflow <config-file>\n";
    std::cerr << "Check examples to see config file format.\n";
    exit(0);
}


void amflow::checkdep() {
    try {
        kira   = config["Install"]["KiraExecutable"].as<std::string>();
        fermat = config["Install"]["FermatExecutable"].as<std::string>();
    } catch (YAML::BadConversion& ex) {
        std::cerr << "AMFlow: failed to find Kira or Fermat executable\n";
        exit(1);
    }

    if (!std::filesystem::exists(kira) || !std::filesystem::exists(fermat)) {
        std::cerr << "AMFlow: failed to find Kira or Fermat executable\n";
        exit(1);
    }

    std::cerr << "Kira executable path:   " << kira   << "\n";
    std::cerr << "Fermat executable path: " << fermat << "\n";
}

