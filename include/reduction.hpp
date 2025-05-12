#ifndef REDUCTION_HPP
#define REDUCTION_HPP


/**
 * reduction.hpp - High-level interface for integral reduction.
 */


#include "kira.hpp"


class integral_reducer {

public:
    /**
     * Instantiate an integral reducer based on the information stored
     * in a YAML node and a given integral family.
     * 
     * @param _node configure YAML node
     * @param _family an integral family
     */
    integral_reducer(const YAML::Node& _node, const integralfamily& _family);

    /**
     * Reduce a list of target integrals to master integrals.
     * 
     * @param target a list of target integrals
     * @param preferred a list of preferred masters
     * @param workdir Kira working directory
     * 
     * @returns a pair of objects:
     *    (1) an integral list consisting of all master integrals
     *    (2) a matrix whose (i, j) element is the expansion coefficient
     *        of the i-th target integral into the j-th master integral
     */
    std::pair<integral_list, GiNaC::matrix>
    reduce(const integral_list& target, const integral_list& preferred, const char* workdir);

    /**
     * Derive differential equations of master integrals with respect to
     * auxiliary mass "eta".
     * 
     * @param preferred a list of preferred masters
     * @param workdir Kira working directory
     * 
     * @returns a pair of objects:
     *    (1) an integral list consisting of all master integrals
     *    (2) coefficient matrix of differential equation
     */
    std::pair<integral_list, GiNaC::matrix>
    diffeq(const integral_list& preferred, const char* workdir);
private:
    kira_agent            agent;
    const integralfamily* pfamily;
    int                   mrank;
    int                   mdot;
};


#endif // REDUCTION_HPP
