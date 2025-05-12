#include "amflow.hpp"
#include "reduction.hpp"
#include "utils.hpp"
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

    
    std::cout << family.to_string();

    integral_reducer reducer(config, family);
    // kira_agent agent(config, family);
    auto out = reducer.diffeq(integral_list{integral(family, {2,2,1,1,0,0,0,0,0}),
                                   integral(family, {2,1,1,1,0,0,0,0,0}),
                                   integral(family, {1,1,1,1,0,0,0,0,0}),
                                   integral(family, {1,1,1,0,0,0,0,0,0})}, "/root/amflow_cpp_test");
    // auto out = reducer.reduce(integral_list{integral(family, {1,2,2,2}),
    //                                integral(family, {2,2,2,2})},
                    //  {},"/root/amflow_cpp_test");
    std::cout << out.first << "\n";
    std::cout << out.second << "\n";
// out = agent.reduce("/root/amflow_cpp_test", integral_list{integral(family, {1,2,2,2,0,0,0,0,0}),
//                                    integral(family, {2,2,2,2,0,0,0,-1,0})},
//                      integral_list{integral(family, {2,2,2,0,0,0,0,0,0}),
//                                    integral(family, {1,1,1,1,0,0,0,0,0})},
//                 15, 3, 4, false);
//     std::cout << out.first << "\n";
//     std::cout << out.second << "\n";
    

// BlackBoxReduce[{j[tt,1,2,2,2,1,2,1,0,0], j[tt,2,2,2,2,-1,2,-1,-2,-1]}, 
//                {j[tt,1,1,1,1,1,1,1,0,0], j[tt,1,0,0,1,0,0,1,0,0], j[tt,1,1,1,0,0,0,0,-1,0]}, "/root/amflow_test"]
//     Print[BlackBoxReduce[{j[box1,1,2,2,2], j[box1,2,2,2,2]}, 
//  {j[box1,2,2,2,0], j[box1,1,1,1,1]}, "/root/amflow_test"]] 
    
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

