#ifndef SOLVE_HPP
#define SOLVE_HPP


#include "border.hpp"


struct boundary_condition {
    GiNaC::ex               exponent;
    std::vector<GiNaC::lst> expansion;
};


class public_diffeq: public symdiffeq {
public:
    public_diffeq(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps)
        : symdiffeq(_coeff, _var, _eps) { }
    friend class run_solve;
};


class finite_solver {
public:
    finite_solver(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order);
    GiNaC::matrix solution();
private:
    GiNaC::ex     var;
    std::vector<GiNaC::ex>     denom;
    std::vector<GiNaC::matrix> numer;
    std::vector<GiNaC::matrix> sol;
};


class run_solve {
public:
    /**
     * Given the coefficient matrix A(eta, eps), run variable eta from
     * infinity to zero at a specific eps value.
     * 
     * @param _coeff coefficient matrix A(eta, eps)
     * @param _var differentiation variable eta
     * @param _eps dimensional parameter eps
     * @param _epsvalue specific value for eps
     * @param _rundir running direction (must be either "Im" or "NegIm")
     * @param _config configure YAML file
     */
    run_solve(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps, const GiNaC::numeric& _epsvalue, const std::string& _rundir, const YAML::Node& _config);

    /**
     * Given the string representation of coefficient matrix A(eta, eps),
     * initialize a `run_solve` object.
     * 
     * @param _str string representation of A(eta, eps)
     * @param _epsvalue specific value for eps
     * @param _rundir running direction (must be either "Im" or "NegIm")
     * @param _config configure YAML file
     */
    static run_solve from_string(const std::string& _str, const GiNaC::numeric& _epsvalue, const std::string& _rundir, const YAML::Node& _config);

    /**
     * Give a list of evaluation points.
     */
    GiNaC::lst get_eval_points();

    /**
     * Given the boundary condition at infinity, find the numeric solution
     * at point `first_pnt`.
     * 
     * @param bc boundary condition of the following form: a vector of 
     * `boundary_condition` structs, each containing a leading exponent
     * over eta, and an expansion in negative eta powers that should be
     * multiplied on the leading power
     * @param first_pnt point to evaluate
     */
    GiNaC::matrix  infinity_solve(const std::vector<boundary_condition>& bc, const GiNaC::numeric& first_pnt);

    /**
     * Given the value at `first_pnt`, find the numeric solution at point
     * `next_pnt`.
     * 
     * @param first_pnt point of initial condition
     * @param next_pnt next point to evaluate
     * @param init_val value at `first_pnt`
     */
    GiNaC::matrix  finite_solve(const GiNaC::numeric& first_pnt, const GiNaC::numeric& next_pnt, const GiNaC::matrix& init_val);

    /**
     * Given the value at `last_pnt`, find the numeric solution at origin.
     * 
     * @param last_pnt point of initial condition
     * @param boundary_val value at `last_pnt`
     */
    GiNaC::matrix  zero_solve(const GiNaC::numeric& last_pnt, const GiNaC::matrix& boundary_val);

    /**
     * Given the boundary condition at infinity, solve result at origin.
     */
    GiNaC::matrix  link(const std::vector<boundary_condition>& bc);
private:
    GiNaC::matrix  coeff_zero;
    GiNaC::matrix  coeff_inf;
    GiNaC::symbol  raw_eps;
    GiNaC::numeric eps_value;

    int            expansion_order;
    int            working_precision;
    int            run_length;
    enum {
        IM,
        NEGIM
    }              run_direction;
    const int      run_candidates = 10;

    public_diffeq  eq_zero;
    border_solve   eq_inf;

    GiNaC::matrix  general_solution_zero;
    GiNaC::matrix  analytic_behavior_zero;
};


#endif // SOLVE_HPP
