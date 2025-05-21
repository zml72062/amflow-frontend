#ifndef BORDER_HPP
#define BORDER_HPP


#include <yaml-cpp/yaml.h>
#include "symdiffeq.h"


class border_solve: public symdiffeq {
public:
    /**
     * Initialize a differential equation [d/dx Y(x, eps) = 
     * -x^2 A(1/x, eps) Y(x, eps)] given the input A(x, eps).
     * 
     * @param _coeff the input coefficient matrix A(x, eps)
     * @param _var the differential variable x
     * @param _eps the auxiliary variable eta
     * @param _config configure YAML file
     */
    border_solve(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps, const YAML::Node& _config);

    /**
     * Initialize a differential equation [d/dx Y(x, eps) = 
     * -x^2 A(1/x, eps) Y(x, eps)] from the string representation of
     * A(x, eps).
     * 
     * @param _str the string representation of A(x, eps) in 
     * Mathematica form
     * @param _config configure YAML file
     */
    static border_solve from_string(const std::string& _str, const YAML::Node& _config);

    /**
     * Get the form of boundary condition.
     * 
     * @returns a vector, each element being a pair
     *      (a) an exponent, representing a specific analytic behavior
     *      (b) a list of integers, representing how many terms should
     *          be determined for the analytic factor on this exponent,
     *          for every master integral
     */
    std::vector<std::pair<GiNaC::ex, std::vector<int>>> get_boundary_form();

private:
    int           try_order;
    GiNaC::ex     eps_value;
    GiNaC::matrix reciprocal_coeff;
};



#endif // BORDER_HPP
