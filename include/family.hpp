#ifndef FAMILY_HPP
#define FAMILY_HPP


/**
 * family.hpp - Define a Feynman integral family.
 */


#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>


struct integralfamily {
    /**
     * Create an integral family from the configure information stored in a
     * YAML node.
     * 
     * @param _node configure YAML node
     * @param _psymbols pointer to an external symbol table
     */
    static integralfamily from_yaml(const YAML::Node& _node, GiNaC::symtab* _psymbols);

    std::string to_string() const;      // for debug

    std::string         name;
    GiNaC::symtab*      psymbols;

    GiNaC::lst          loops;
    GiNaC::lst          indeplegs;
    GiNaC::lst          propagators;
    std::vector<int>    cut;
    std::vector<int>    prescription;

    GiNaC::lst      conservation_rules;     // momentum conservation
    GiNaC::lst      numeric_rules;          // numeric kinematics
    GiNaC::lst      sps_numeric_rules;      // numeric scalar products
};


#endif // FAMILY_HPP
