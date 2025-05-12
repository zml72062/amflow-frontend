#ifndef AMFLOW_HPP
#define AMFLOW_HPP


#include "family.hpp"


class amflow {
public:
    amflow(const char* config_path);

    static void usage();
private:
    // message
    const std::string author  = "Xiao Liu and Yan-Qing Ma";
    const std::string version = "1.1";
    const std::string release = "5-Jun-2022";
    void startup();

    YAML::Node config;

    // dependency
    std::string kira;
    std::string fermat;
    void checkdep();
    
    // global symbol table
    GiNaC::symtab globals;

    // integral family
    integralfamily family;
};


#endif // AMFLOW_HPP
