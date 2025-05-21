#include "kira.hpp"
#include "utils.hpp"
#include <iostream>
#include <numeric>
#include <bitset>
#include <sstream>
#include <fstream>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <filesystem>


class pretty_lst {
public:
    pretty_lst(const GiNaC::lst& lst): plst(&lst) { }

    friend std::ostream& operator<<(std::ostream& out, const pretty_lst& lst) {
        int nops = lst.plst->nops();
        for (int i = 0; i < nops; i++)
            out << (i ? "," : "") << (*lst.plst)[i];
        return out;
    }
private:
    const GiNaC::lst* plst;
};



kira_agent::kira_agent(const YAML::Node& _node, const integralfamily& _family)
    : helper(_family, *_family.psymbols) {
    kira    = _node["Install"]["KiraExecutable"].as<std::string>();
    fermat  = _node["Install"]["FermatExecutable"].as<std::string>();
    pfamily = &_family;

    integral_order = 5;
    reduction_mode = Kira;
    d0             = 4;
    nthreads       = 1;

    try {
        if (has_non_null_key(_node, "ReduceOption")) {
            auto reduce_options = _node["ReduceOption"];
            if (has_non_null_key(reduce_options, "IntegralOrder")) {
                integral_order = reduce_options["IntegralOrder"].as<int>();
                if (integral_order < 1 || integral_order > 8) {
                    std::cerr << "KiraAgent: " << integral_order
                              << " is not a valid integral ordering\n";
                    exit(1);
                }
            }
            if (has_non_null_key(reduce_options, "ReductionMode")) {
                auto mode_str = reduce_options["ReductionMode"].as<std::string>();
                if (mode_str == "Masters")
                    reduction_mode = Masters;
                else if (mode_str == "Kira")
                    reduction_mode = Kira;
                else if (mode_str == "FireFly")
                    reduction_mode = FireFly;
                else if (mode_str == "Mixed")
                    reduction_mode = Mixed;
                else if (mode_str == "NoFactorScan")
                    reduction_mode = NoFactorScan;
                else {
                    std::cerr << "KiraAgent: " << mode_str
                              << " is not a valid reduction mode\n";
                    exit(1);
                }
            }
        }
        if (has_non_null_key(_node, "Run")) {
            auto run_options = _node["Run"];
            if (has_non_null_key(run_options, "D0"))
                d0 = GiNaC::parser()(run_options["D0"].as<std::string>());
            if (has_non_null_key(run_options, "NThread")) {
                nthreads = run_options["NThread"].as<int>();
                if (nthreads <= 0) {
                    std::cerr << "KiraAgent: invalid thread number\n";
                    exit(1);
                }
            }
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "KiraAgent: failed while reading reduction options\n";
        exit(1);
    }

    std::cerr << "KiraAgent: Set IntegralOrder = " << integral_order << "\n";
    switch (reduction_mode) {
        case Masters:
            std::cerr << "KiraAgent: Set ReductionMode = Masters\n";
            break;
        case Kira:
            std::cerr << "KiraAgent: Set ReductionMode = Kira\n";
            break;
        case FireFly:
            std::cerr << "KiraAgent: Set ReductionMode = FireFly\n";
            break;
        case Mixed:
            std::cerr << "KiraAgent: Set ReductionMode = Mixed\n";
            break;
        case NoFactorScan:
            std::cerr << "KiraAgent: Set ReductionMode = NoFactorScan\n";
            break;
    }
}



void kira_agent::write_integralfamilies_yaml(std::ostream& stream, unsigned long top_sector) {
    stream << "integralfamilies:\n"
           << "  - name: \""             << name()              << "\"\n"
           << "    loop_momenta: ["      << pretty_lst(loops()) << "]\n"
           << "    top_level_sectors: [" << top_sector          << "]\n"
           << "    propagators:\n";
    for (auto& prop: propagators())
        stream << "      - [\"" << prop << "\", 0]\n";

    if (cut().size() > 0 && std::accumulate(cut().begin(), cut().end(), 0) > 0) {
        stream << "    cut_propagators: ";
        int nprops = cut().size(), ncut = 0;
        for (int i = 0; i < nprops; i++)
            if (cut()[i])
                stream << ",["[ncut++ == 0] << i;
        stream << "]\n";
    }
}


void kira_agent::write_kinematics_yaml(std::ostream& stream) {
    GiNaC::symbol auxleg(name() + "AuxLeg");
    GiNaC::ex     auxconservation = 0;
    std::for_each(indeplegs().begin(), indeplegs().end(), 
                  [&auxconservation](const auto& ex){ auxconservation -= ex; });

    stream << "kinematics:\n"
           << "  incoming_momenta: ["      << auxleg << "," << pretty_lst(indeplegs()) << "]\n"
           << "  outgoing_momenta: []\n"
           << "  momentum_conservation: [" << auxleg << "," << auxconservation         << "]\n"
           << "  kinematic_invariants:\n";

    for (auto& prop: propagators()) {
        if (prop.diff(GiNaC::ex_to<GiNaC::symbol>(psymbols()->at("eta"))) != 0) {
            // there is eta dependence in propagators
            stream << "    - [eta, 2]\n";
            break;
        }
    }

    stream << "  scalarproduct_rules:\n";

    int num_indeplegs = indeplegs().nops(), n = 0;
    for (int i = 0; i < num_indeplegs; i++)
        for (int j = i; j < num_indeplegs; j++)
            stream << "    - [[" << indeplegs()[i] << "," << indeplegs()[j] << "],"
                   << sps_numeric_rules()[n++].op(1) << "]\n";
}


void kira_agent::write_jobs_yaml(std::ostream& stream, unsigned long top_sector, 
                                 unsigned rank, unsigned dot) {
    unsigned maxr = std::bitset<64>(top_sector).count() + dot;

    stream << "jobs:\n"
           << "  - reduce_sectors:\n"
           << "      reduce:\n"
           << "        - {topologies: [" << name()     << "], "
           <<            "sectors: ["    << top_sector << "], "
           <<            "r: "           << maxr       << ", "
           <<            "s: "           << rank       << "}\n";
    
    if (reduction_mode == Masters) {
        stream << "      select_integrals:\n"
               << "        select_mandatory_recursively:\n"
               << "          - {topologies: [" << name()     << "], "
               <<              "sectors: ["    << top_sector << "], "
               <<              "r: "           << maxr       << ", "
               <<              "s: "           << rank       << ", "
               <<              "d: "           << dot        << "}\n";
    } else {
        stream << "      select_integrals:\n"
               << "        select_mandatory_list:\n"
               << "          - [" << name() << ", target]\n";
    }

    stream << "      preferred_masters: preferred\n"
           << "      integral_ordering: " << integral_order << "\n";

    switch (reduction_mode) {
        case Masters:
            stream << "      run_initiate: masters\n";
            break;
        case Kira:
            stream << "      run_initiate: true\n"
                   << "      run_triangular: true\n"
                   << "      run_back_substitution: true\n";
            break;
        case FireFly:
            stream << "      run_initiate: true\n"
                   << "      run_firefly: true\n";
            break;
        case Mixed:
            stream << "      run_initiate: true\n"
                   << "      run_triangular: true\n"
                   << "      run_firefly: back\n";
            break;
        case NoFactorScan:
            stream << "      run_initiate: true\n"
                   << "      run_triangular: true\n"
                   << "      run_firefly: back\n"
                   << "      factor_scan: false\n";
            break;
    }

    if (reduction_mode != Masters)
        stream << "  - kira2math:\n"
               << "      target:\n"
               << "        - [" << name() << ", target]\n";
}


void kira_agent::write_preferred(std::ostream& stream, const integral_list& preferred) {
    for (auto& intgr: preferred)
        stream << intgr.to_string() << "\n\n";
}


void kira_agent::write_target(std::ostream& stream, const integral_list& target) {
    for (auto& intgr: target)
        stream << intgr.to_string() << "\n\n";
}


int kira_agent::run_kira() {
    pid_t pid = fork();

    if (pid == 0) { // run Kira in child process
        char kira_exec[1024] = {0};
        strcpy(kira_exec, kira.c_str());

        char parallel_opt[128] = {0};
        strcpy(parallel_opt, ("-p" + std::to_string(nthreads)).c_str());

        char kira_jobs[128] = {0};
        strcpy(kira_jobs, "jobs.yaml");

        char fermat_exec[1024] = {0};
        strcpy(fermat_exec, ("FERMATPATH=" + fermat).c_str());

        char* argv[] = {kira_exec, parallel_opt, kira_jobs, 0};
        char* envp[] = {fermat_exec, 0};

        std::cerr << fermat_exec  << " "
                  << kira_exec    << " "
                  << parallel_opt << " "
                  << kira_jobs    << std::endl;
        
        // suppress Kira output
        if (freopen("/dev/null", "w", stdout) == nullptr) {
            std::cerr << "KiraAgent: error occurred when trying to redirect output\n";
            exit(1);
        }

        if (execve(kira_exec, argv, envp) < 0) {
            std::cerr << "KiraAgent: error occurred when trying to run Kira in a subprocess\n";
            exit(1);
        }
    } else { // wait for child process
        if (waitpid(pid, 0, 0) < 0) {
            std::cerr << "KiraAgent: error occurred when waiting for subprocess\n";
            return -1;
        }
    }
    return 0;
}


std::vector<unsigned long> kira_agent::trivial_sectors(const char* workdir, unsigned long top_sector) {
    pid_t pid = fork();

    if (pid == 0) { // work in child process
        std::filesystem::create_directory(std::filesystem::path(workdir));
        if (chdir(workdir) < 0) {
            std::cerr << "KiraAgent: error occurred when trying to chdir() in a subprocess\n";
            exit(1);
        }

        // remove previous results
        std::filesystem::remove_all(std::filesystem::path(workdir).append("tmp"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("firefly_saves"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("sectormappings"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("results"));

        
        // write "integralfamilies.yaml"
        std::filesystem::create_directory(std::filesystem::path(workdir).append("config"));
        std::ofstream outf(std::filesystem::path(workdir).append("config").append("integralfamilies.yaml"));
        write_integralfamilies_yaml(outf, top_sector);
        outf.close();

        // write "kinematics.yaml"
        std::ofstream outk(std::filesystem::path(workdir).append("config").append("kinematics.yaml"));
        write_kinematics_yaml(outk);
        outk.close();

        // write "jobs.yaml"
        std::ofstream outj(std::filesystem::path(workdir).append("jobs.yaml"));
        outj << "jobs:\n"
             << "  - reduce_sectors:\n"
             << "      reduce:\n"
             << "        - {topologies: [" << name()     << "], "
             <<            "sectors: ["    << top_sector << "], "
             <<            "r: "           << 100        << ", "
             <<            "s: "           << 100        << "}\n"
             << "      run_symmetries: true\n";
        outj.close();

        // run Kira
        run_kira();
        exit(0);
    } else {
        if (waitpid(pid, 0, 0) < 0) {
            std::cerr << "KiraAgent: error occurred when waiting for subprocess\n";
            return {};
        }

        std::string sectors_path = std::filesystem::path(workdir)
                                        .append("sectormappings")
                                        .append(name())
                                        .append("trivialsector");
        std::ifstream intrivial(sectors_path);
        std::vector<unsigned long> result;
        unsigned long sect;
        while (intrivial >> sect) {
            result.push_back(sect);
            intrivial.get();
        }
        intrivial.close();
        return result;
    }
}


std::pair<integral_list, GiNaC::matrix>
kira_agent::reduce(const char* workdir, const integral_list& target, 
                   const integral_list& preferred, unsigned long top_sector, 
                   unsigned rank, unsigned dot, bool masters_only) {
    pid_t pid = fork();

    if (pid == 0) { // work in child process
        std::filesystem::create_directory(std::filesystem::path(workdir));
        if (chdir(workdir) < 0) {
            std::cerr << "KiraAgent: error occurred when trying to chdir() in a subprocess\n";
            exit(1);
        }

        // remove previous results
        std::filesystem::remove_all(std::filesystem::path(workdir).append("tmp"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("firefly_saves"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("sectormappings"));
        std::filesystem::remove_all(std::filesystem::path(workdir).append("results"));
        
        // write "integralfamilies.yaml"
        std::filesystem::create_directory(std::filesystem::path(workdir).append("config"));
        std::ofstream outf(std::filesystem::path(workdir).append("config").append("integralfamilies.yaml"));
        write_integralfamilies_yaml(outf, top_sector);
        outf.close();

        // write "kinematics.yaml"
        std::ofstream outk(std::filesystem::path(workdir).append("config").append("kinematics.yaml"));
        write_kinematics_yaml(outk);
        outk.close();

        // write "jobs.yaml"
        if (masters_only)
            reduction_mode = Masters;
        std::ofstream outj(std::filesystem::path(workdir).append("jobs.yaml"));
        write_jobs_yaml(outj, top_sector, rank, dot);
        outj.close();

        // write "preferred"
        std::ofstream outp(std::filesystem::path(workdir).append("preferred"));
        write_preferred(outp, preferred);
        outp.close();

        // write "target"
        if (!masters_only) {
            std::ofstream outt(std::filesystem::path(workdir).append("target"));
            write_target(outt, target);
            outt.close();
        }

        // run Kira
        run_kira();
        exit(0);
    } else {
        if (waitpid(pid, 0, 0) < 0) {
            std::cerr << "KiraAgent: error occurred when waiting for subprocess\n";
            return std::make_pair(integral_list(), GiNaC::matrix(0, 0));
        }

        std::string masters_path = std::filesystem::path(workdir)
                                                .append("results")
                                                .append(name())
                                                .append("masters");
        if (masters_only) {
            helper.load_masters(masters_path.c_str());
            return std::make_pair(helper.get_master_integral_list(), GiNaC::matrix(0, 0));
        } else {
            // try reading masters again
            ibphelper newhelper(*pfamily, *psymbols());
            newhelper.load_masters(masters_path.c_str());
            for (auto& mi: newhelper.master_table) {
                if (helper.master_table.find(mi.first) == helper.master_table.end()) {
                    std::cerr << "KiraAgent: inconsistent masters from Kira\n";
                    exit(1);
                }
            }

            auto master_list = helper.get_master_integral_list();
            int ntarget = target.size(), nmaster = master_list.size();

            if (nmaster == 0)
                return std::make_pair(integral_list(), GiNaC::matrix(0, 0));

            if (ntarget == 0)
                return std::make_pair(master_list, GiNaC::matrix(0, 0));

            std::string ibp_path = std::filesystem::path(workdir)
                                                .append("results")
                                                .append(name())
                                                .append("kira_target.m");
            helper.load_ibps(ibp_path.c_str());

            GiNaC::matrix coeff(ntarget, nmaster);
            for (int i = 0; i < ntarget; i++) {
                for (int j = 0; j < nmaster; j++) {
                    coeff(i, j) = helper.ibp_table.at(target[i].to_string())
                        .diff(GiNaC::ex_to<GiNaC::symbol>(
                            helper.master_table.at(master_list[j].to_string())
                        ))
                        .subs(psymbols()->at("d") == d0 - 2 * psymbols()->at("eps"), 
                              GiNaC::subs_options::algebraic);
                }
            }
            return std::make_pair(master_list, coeff);
        }
    }
}

