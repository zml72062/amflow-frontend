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
    START_TIME(forward_pass);
    
    while (!enroute.empty()) {
        int cursysnum = enroute.front();
        enroute.pop();

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

