#include "amflow.hpp"
#include "reduction.hpp"
#include "boundary.hpp"
#include "singlesetup.hpp"
#include "utils.hpp"
#include "solve.hpp"
#include <fstream>
#include <sstream>
#include <filesystem>


amflow::amflow(const char* config_path) {
    START_TIME(initialize);
    startup();

    config = YAML::LoadFile(config_path);
    checkdep();

    globals["d"]   = GiNaC::symbol("d");        // space-time dimension
    globals["eps"] = GiNaC::symbol("eps");      // d = d0 - 2 * eps
    globals["eta"] = GiNaC::symbol("eta");      // auxiliary mass

    family = integralfamily::from_yaml(config, &globals);
    std::cout << family.to_string() << "\n";
    std::cout << family.propagator_prescription() << "\n";
    // family.generate_components({1,1,1,1,1,1,1});
    // auto newfamily = family.insert_eta({1,1,1,1,1,1,1}, {Prescription});
    // if (newfamily == nullptr)
    //     std::cout << "null\n";
    // else
    //     std::cout << newfamily->to_string() << "\n";

    integral_system sys(family, {{1,1,1,1,1,1,1},{1,2,1,1,1,1,1}},{},config, "/root/amflow_cpp_test3");
    sys.reduce_targets();
    sys.build_diffeq();
    sys.determine_boundaries();
    sys.determine_border();
    sys.determine_direction();
    sys.setup_subfamilies();

    // std::ifstream in(std::filesystem::path("/root/amflow_cpp_test2").append("DIFFEQ"));
    // std::string diffeq_str;
    // for (int i = 0; i < 4; i++)
    //     in >> diffeq_str;
    // in.close();

    // auto solver = run_solve::from_string(diffeq_str, GiNaC::numeric(1)/100, "NegIm", config);
    // std::vector<boundary_condition> v;
    // boundary_condition a1;
    // a1.exponent = 1-globals["eps"];
    // a1.expansion = {{}, {1}, {}, {}, {}, {}};
    // v.push_back(a1);
    // std::cout << solver.link(v) << std::endl;
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

