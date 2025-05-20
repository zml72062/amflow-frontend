#ifndef APART_HPP
#define APART_HPP


/**
 * apart.hpp - Partial fraction expansion.
 */


#include <ginac/ginac.h>


struct parfrac {
    GiNaC::lst        denoms;
    std::vector<int>  npowers;
    GiNaC::ex         coeff;
};


/**
 * Compute polynomial division `numer / denom` with respect to variable
 * `var`. Returns a pair (quotient, remainder).
 */
GiNaC::lst generalized_div(const GiNaC::ex& numer, const GiNaC::ex& denom, const GiNaC::symbol& var);


/**
 * Expand `frac` into partial fractions, with respect to variable `var`.
 * 
 * @param frac a rational function
 * @param var variable to expand over
 */
std::vector<parfrac> parfrac_expansion(const GiNaC::ex& frac, const GiNaC::symbol& var);


/**
 * Expand `frac` into partial fractions, with respect to variables `vars`.
 * 
 * @param frac a rational function
 * @param vars a list of variables to expand over
 */
std::vector<parfrac> parfrac_expansion(const GiNaC::ex& frac, const GiNaC::lst& vars);


#endif // APART_HPP
