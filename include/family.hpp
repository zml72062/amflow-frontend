#ifndef FAMILY_HPP
#define FAMILY_HPP


/**
 * family.hpp - Define a Feynman integral family.
 */


#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>


struct sps_struct {
    GiNaC::lst       prop;
    GiNaC::lst       sps;
    GiNaC::matrix    coeff;
    GiNaC::matrix    bias;
    std::vector<int> cut;
};


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

    /**
     * Decide whether the current list of propagators is complete.
     */
    bool       is_complete_propagator() const;

    /**
     * First, extend the current list of propagators into a complete
     * basis. Then represent all [S = L(L+1)/2 + NL] scalar products 
     * into linear combinations of the basis.
     * 
     * @returns An `sps_struct`, consisting of
     *      - `prop`:  the list of complete propagators
     *      - `sps`:   the list of all scalar products
     *      - `coeff`: coefficient matrix (S by S)
     *      - `bias`:  bias matrix (S by 1)
     *      - `cut`:   whether each of the complete propagators is
     *                 cut, if the original propagators have cuts
     * `coeff` and `bias` satisfy `coeff * prop + bias = sps`.
     */
    sps_struct sps_in_complete_propagator() const;

    /**
     * Determine the momentum that flows on each propagator.
     */
    GiNaC::lst momenta() const;

    /**
     * Determine the loop momentum that flows on each propagator.
     */
    GiNaC::lst loop_momenta() const;
};


#endif // FAMILY_HPP
