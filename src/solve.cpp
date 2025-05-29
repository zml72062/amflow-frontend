#include "solve.hpp"
#include "apart.hpp"
#include "utils.hpp"
#include "wrapper.h"
#include "ratsolver.h"
#include <sstream>


#define SUBS(m, r)  (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).subs((r), GiNaC::subs_options::algebraic)))
#define NORMAL(m)   (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).normal()))
#define EXPAND(m)   (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).expand()))


static GiNaC::ex gnq2q(const GiNaC::ex& gnq) {
    GiNaC::parser parser;
    if (GiNaC::ex_to<GiNaC::numeric>(gnq).is_real()) {
        std::ostringstream ostr;
        ostr << gnq.normal();
        std::size_t pkt, E;
        int power = 0;
        std::string ostrstr;
        ostrstr = ostr.str();
        GiNaC::ex fq;

        DETERMINE:
        if (((pkt = ostrstr.find('.')) == std::string::npos) &
            ((E = ostrstr.find('E')) == std::string::npos))
            fq = parser(ostrstr.c_str());
        else if (E == std::string::npos)
            fq = parser(ostrstr.substr(0, pkt) + ostrstr.substr(pkt + 1) 
                 + "/1" + std::string(ostrstr.size() - pkt - 1, '0'));
        else {
            power = std::stoi(ostrstr.substr(E + 1));
            ostrstr = ostrstr.substr(0, E);
            goto DETERMINE;
        }

        if (power != 0) {
            fq *= GiNaC::pow(10, power);
        }

        return fq;
    } else {
        auto real = gnq.real_part(), imag = gnq.imag_part();
        return gnq2q(real) + GiNaC::I * gnq2q(imag);
    }
}


run_solve::run_solve(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps, const GiNaC::numeric& _epsvalue, const std::string& _rundir, const YAML::Node& _config)
    : coeff_zero(_coeff), coeff_inf(NORMAL(SUBS(_coeff, _var == 1 / _var).mul_scalar(-GiNaC::pow(_var, -2)))),
      raw_eps(_eps), eps_value(_epsvalue),
      eq_zero(SUBS(coeff_zero, _eps == _epsvalue), _var, _eps),
      eq_inf (SUBS(coeff_zero, _eps == _epsvalue), _var, _eps, _config) {

    expansion_order   = 100;
    working_precision = 100;
    run_length        = 1000;
    try {
        if (has_non_null_key(_config, "ExpansionOption")) {
            auto expansion_opt = _config["ExpansionOption"];
            if (has_non_null_key(expansion_opt, "SolvingOrder")) {
                expansion_order = expansion_opt["SolvingOrder"].as<int>();
            }
            if (has_non_null_key(expansion_opt, "WorkingPre")) {
                working_precision = expansion_opt["WorkingPre"].as<int>();
            }
            if (has_non_null_key(expansion_opt, "RunLength")) {
                run_length = expansion_opt["RunLength"].as<int>();
            }
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "RunSolve: failed while reading expansion options\n";
        exit(1);
    }

    GiNaC::Digits = working_precision;
    std::cerr << "RunSolve: solving for eps = "       << _epsvalue         << "\n";
    std::cerr << "RunSolve: set expansion order = "   << expansion_order   << "\n";
    std::cerr << "RunSolve: set working precision = " << working_precision << "\n";

    if (_rundir == "Im") {
        run_direction = IM;
    } else if (_rundir == "NegIm") {
        run_direction = NEGIM;
    } else {
        std::cerr << "RunSolve: invalid running direction\n";
        exit(1);
    }
}


run_solve run_solve::from_string(const std::string& _str, const GiNaC::numeric& _epsvalue, const std::string& _rundir, const YAML::Node& _config) {
    GiNaC::symbol eta("eta"), eps("eps");
    GiNaC::symtab tab;
    tab["eta"] = eta;
    tab["eps"] = eps;
    GiNaC::parser parser(tab);
    auto result = parser(_str);
    int rows = result.nops();
    GiNaC::matrix coeff(rows, rows);
    int i = 0, j;
    for (auto r = result.begin(); r != result.end(); ++r, i++) {
        if ((int)r->nops() != rows)
            throw std::runtime_error("Input not a square matrix!");
        j = 0;
        for (auto c = r->begin(); c != r->end(); ++c, j++) {
            coeff(i, j) = *c;
        }
    }
    return run_solve(coeff, eta, eps, _epsvalue, _rundir, _config);
}


GiNaC::lst run_solve::get_eval_points() {
    GiNaC::Digits = working_precision;
    std::cerr << "RunSolve: solving around zero\n";
    auto result = eq_zero.solve(0, expansion_order);
    auto J = eq_zero.reg_struct.J;
    if (J.rows() == 1 && J.cols() == 1 && J(0, 0) == 0) {
        analytic_behavior_zero = GiNaC::matrix(coeff_zero.rows(), coeff_zero.rows());
    } else {
        analytic_behavior_zero = J;
    }
    general_solution_zero = result.first;
    
    auto eq_singularities = result.second[0].expand();
    flint::ca_ctx  ctx;
    flint::ca_poly poly(ctx);
    poly.import_from_ginac(eq_singularities, eq_zero.x);
    GiNaC::lst  singularities;
    
    GiNaC::ex closest = GiNaC::pow(10, 20), farthest = GiNaC::pow(10, -20);
    int length = poly.ptr()->length;
    if (length == 0) {
        std::cerr << "RunSolve: invalid singularity list\n";
        exit(1);
    }
    if (length > 1) {
        int degree = length - 1;
        ca_vec_t roots;
        std::vector<ulong> mults(degree);
        ca_vec_init(roots, degree, ctx.ctx);
        if (ca_poly_roots(roots, mults.data(), poly.ptr(), ctx.ctx) == 0) {
            std::cerr << "RunSolve: fail to find singularities\n";
            ca_vec_clear(roots, ctx.ctx);
            exit(1);
        }
        int nroots = ca_vec_length(roots, ctx.ctx);
        GiNaC::Digits = 10;
        for (int i = 0; i < nroots; i++) {
            auto sing = gnq2q(fc2g(ca_vec_entry(roots, i), 10, ctx.ctx));
            if ((bool)(sing == 0))
                continue;
            if (GiNaC::abs(sing) < GiNaC::abs(closest))
                closest = GiNaC::abs(sing);
            if (GiNaC::abs(sing) > GiNaC::abs(farthest))
                farthest = GiNaC::abs(sing);
            singularities.append(sing);
        }
        ca_vec_clear(roots, ctx.ctx);
    }

    std::cerr << "RunSolve: list of singularities is " << singularities << "\n";

    GiNaC::ex dir_real = 100000000, dir_imag;
    int failed = 0;
    switch (run_direction) {
        case IM:
            dir_imag = 1;
            for (auto& sing: singularities) {
                if ((bool)(sing.real_part() == 0) && (bool)(sing.imag_part() > 0)) {
                    failed = 1;
                    break;
                }
            }
            if (!failed) {
                dir_real = 0;
            } else {
                for (int test = 1; test <= run_candidates; test++) {
                    int failed_pos = 0, failed_neg = 0;
                    for (auto& sing: singularities) {
                        if ((bool)(sing.imag_part() > 0) && 
                            (bool)((sing.real_part() / sing.imag_part()) 
                                == (GiNaC::numeric(test) / run_candidates))) {
                            failed_pos = 1;
                            break;
                        }
                    }
                    if (!failed_pos) {
                        dir_real = GiNaC::numeric(test) / run_candidates;
                        break;
                    }
                    for (auto& sing: singularities) {
                        if ((bool)(sing.imag_part() > 0) && 
                            (bool)((sing.real_part() / sing.imag_part()) 
                                == -(GiNaC::numeric(test) / run_candidates))) {
                            failed_neg = 1;
                            break;
                        }
                    }
                    if (!failed_neg) {
                        dir_real = -GiNaC::numeric(test) / run_candidates;
                        break;
                    }
                }
            }

            if ((bool)(dir_real == 100000000)) {
                std::cerr << "RunSolve: failed to find a running direction\n";
                exit(1);
            }
            break;
        case NEGIM:
            dir_imag = -1;
            for (auto& sing: singularities) {
                if ((bool)(sing.real_part() == 0) && (bool)(sing.imag_part() < 0)) {
                    failed = 1;
                    break;
                }
            }
            if (!failed) {
                dir_real = 0;
            } else {
                for (int test = 1; test <= run_candidates; test++) {
                    int failed_pos = 0, failed_neg = 0;
                    for (auto& sing: singularities) {
                        if ((bool)(sing.imag_part() < 0) && 
                            (bool)((sing.real_part() / sing.imag_part()) 
                                == (GiNaC::numeric(test) / run_candidates))) {
                            failed_pos = 1;
                            break;
                        }
                    }
                    if (!failed_pos) {
                        dir_real = -GiNaC::numeric(test) / run_candidates;
                        break;
                    }
                    for (auto& sing: singularities) {
                        if ((bool)(sing.imag_part() < 0) && 
                            (bool)((sing.real_part() / sing.imag_part()) 
                                == -(GiNaC::numeric(test) / run_candidates))) {
                            failed_neg = 1;
                            break;
                        }
                    }
                    if (!failed_neg) {
                        dir_real = GiNaC::numeric(test) / run_candidates;
                        break;
                    }
                }
            }

            if ((bool)(dir_real == 100000000)) {
                std::cerr << "RunSolve: failed to find a running direction\n";
                exit(1);
            }
            break;
    }

    auto dir_c = dir_real + GiNaC::I * dir_imag;
    std::cerr << "RunSolve: choose a running direction " << dir_c << "\n";

    GiNaC::Digits = 20;
    GiNaC::ex last_pnt  = gnq2q(closest  / 2 / GiNaC::abs(dir_c)) * dir_c;
    GiNaC::ex first_pnt = gnq2q(farthest * 2 / GiNaC::abs(dir_c)) * dir_c;
    
    if ((bool)(GiNaC::abs(last_pnt) >= GiNaC::abs(first_pnt))) {
        if ((bool)(closest  == GiNaC::pow(10, 20)) &&
            (bool)(farthest == GiNaC::pow(10, -20)))
            return {1};
        else
            return {first_pnt};
    }

    GiNaC::lst points;
    points.append(first_pnt);
    GiNaC::ex  current_nearest = first_pnt;
    int npoints = 1;
    while ((bool)(GiNaC::abs(current_nearest) > GiNaC::abs(last_pnt))) {
        npoints++;
        GiNaC::ex min_dist = GiNaC::abs(current_nearest);
        for (auto& sing: singularities) {
            auto dist = GiNaC::abs(sing - current_nearest);
            if (dist < min_dist)
                min_dist = dist;
        }
        
        if ((bool)(GiNaC::abs(current_nearest) > min_dist / 2)) {
            auto next_pnt = gnq2q((GiNaC::abs(current_nearest) - min_dist / 2) / GiNaC::abs(dir_c)) * dir_c;
            points.append(next_pnt);
            current_nearest = next_pnt;
        } else {
            auto next_pnt = last_pnt;
            points.append(next_pnt);
            current_nearest = 0;
        }

        if (npoints >= run_length) {
            std::cerr << "RunSolve: run exceeded designated number of steps\n";
            exit(1);
        }
    }
    
    return points;
}


GiNaC::matrix run_solve::infinity_solve(const std::vector<boundary_condition>& bc, const GiNaC::numeric& first_pnt) {
    GiNaC::Digits = working_precision;
    std::cerr << "RunSolve: solving around infinity\n";
    eq_inf.clear_all_reduction();
    eq_inf.moser_reduction(0);

    if (!(bool)(eq_inf.get_moser_struct(0).moser_rank() == 1))
        throw std::runtime_error("not implemented: infinity is an irregular singularity");

    GiNaC::matrix A0(eq_inf.N(), eq_inf.N()), A_ana(eq_inf.N(), eq_inf.N());
    for (int i = 0; i < eq_inf.N(); i++)
        for (int j = 0; j < eq_inf.N(); j++)
            A0(i, j) = eq_inf.coeff(i, j).series(eq_inf.x, 1).coeff(eq_inf.x, -1).normal();

    eq_inf.reg_struct.initialize(A0);
    eq_inf.regular_reduction(0);
    auto J = eq_inf.reg_struct.J;

    // find all exponents (we do not care logarithms)
    GiNaC::lst exponent_list;
    std::vector<std::vector<int>> columns;
    for (int j = 0; j < eq_inf.N(); j++) {
        auto exponent = J(j, j);
        bool exist = false;
        int current_l = exponent_list.nops();
        for (int q = 0; q < current_l; q++) {
            if ((bool)((exponent_list[q] - exponent).expand() == 0)) {
                exist = true;
                columns[q].push_back(j);
                break;
            }
        }
        if (!exist) {
            exponent_list.append(exponent);
            columns.push_back({j});
        }
    }
    std::cerr << "RunSolve: allowed exponents are " << exponent_list << "\n";

    // now we can remove spurious singular part of coefficients
    for (int i = 0; i < eq_inf.N(); i++) {
        for (int j = 0; j < eq_inf.N(); j++) {
            std::cerr << "RunSolve: setting up coefficient matrix (" << i << ", " << j << ")          \r";
            A_ana(i, j) = (eq_inf.coeff(i, j) - J(i, j) / eq_inf.x).normal();
        }
    }
    std::cerr << "\nRunSolve: done setting up coefficient matrix\n";

    ratsolver_regular solver(A_ana, eq_inf.x, expansion_order, eq_inf.reg_struct);
    GiNaC::symbol l("l");
    auto out = SUBS(EXPAND(NORMAL(eq_inf.transform.mul(solver.solution()))), GiNaC::log(eq_inf.x) == l);
    
    int num_behaviors = columns.size();
    GiNaC::matrix integral_constants(eq_inf.N(), 1);
    for (int i = 0; i < num_behaviors; i++) {
        auto exponent = exponent_list[i];
        auto span     = columns[i];
        int  crank    = 0, 
             lspan    = span.size(),
             order    = 0;
        
        boundary_condition thisbc;
        int  powerdiff = 0;
        bool found = false;
        for (auto& current_bc: bc) {
            std::ostringstream outf;
            outf << current_bc.exponent;
            GiNaC::symtab tab;
            tab["eps"] = raw_eps;
            GiNaC::parser parser(tab);
            GiNaC::ex boundary_determined_exponent 
                = parser(outf.str()).subs(raw_eps == eps_value).expand();
            auto diff = 
                GiNaC::ex_to<GiNaC::numeric>(-boundary_determined_exponent - exponent);
            if (GiNaC::is_nonneg_integer(diff)) {
                thisbc = current_bc;
                powerdiff = diff.to_int();
                found = true;
                break;
            }
        }
        if (!found) { // this analytic behavior is not relevant
            for (auto& c: span)
                integral_constants(c, 0) = 0;
            continue;
        }

        GiNaC::lst order_eqs;
        GiNaC::lst variables;
        for (int p = 0; p < lspan; p++)
            variables.append(GiNaC::symbol("temp" + std::to_string(p)));
        GiNaC::matrix coeff(0, lspan);
        while (crank < lspan) {
            for (int j = 0; j < eq_inf.N(); j++) {
                GiNaC::matrix next_row(1, lspan);
                for (int q = 0; q < lspan; q++)
                    next_row(0, q) = (out(j, span[q]) / GiNaC::pow(eq_inf.x, exponent))
                        .expand().collect(l).coeff(l, 0)
                                 .collect(eq_inf.x).coeff(eq_inf.x, order);
                
                auto new_coeff = append_row(coeff, next_row);
                if ((int)new_coeff.rank() > crank) {
                    GiNaC::ex lhs = 0, rhs;
                    int order_index = order - powerdiff;
                    if (order_index < 0 || order_index >= (int)thisbc.expansion[j].nops())
                        rhs = 0;
                    else
                        rhs = thisbc.expansion[j][order_index];
                    for (int q = 0; q < lspan; q++)
                        lhs += (variables[q] * next_row(0, q));
                    order_eqs.append(lhs == rhs);
                    coeff = new_coeff;
                    crank++;
                }
                if (crank == lspan)
                    break;
            }
            if (crank == lspan)
                break;
                
            order++;
            if (order >= expansion_order) {
                std::cerr << "RunSolve: inconsistent boundary condition\n";
                exit(1);
            }
        }
        
        auto csol = GiNaC::lsolve(order_eqs, variables);
        for (int q = 0; q < lspan; q++)
            integral_constants(span[q], 0) = variables[q].subs(csol, GiNaC::subs_options::algebraic);
    }

    auto result = SUBS(out, eq_inf.x == 1 / first_pnt).mul(integral_constants);
    for (int i = 0; i < eq_inf.N(); i++)
        result(i, 0) = GiNaC::ex_to<GiNaC::numeric>(result(i, 0)).evalf();
    
    return result;
}


finite_solver::finite_solver(const GiNaC::matrix& A, const GiNaC::ex& x, const GiNaC::matrix& Y0, int order)
    : var(x), sol(order + 1) {
    int N = A.rows();
    if ((int)A.cols() != N)
        throw std::runtime_error("A is not a square matrix!");
    if ((int)Y0.rows() != N)
        throw std::runtime_error("Y0 and A have inconsistent shapes!");
    
    GiNaC::ex common_denom = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto current_denom = A(i, j).normal().denom();
            common_denom = generalized_lcm(common_denom, current_denom, GiNaC::ex_to<GiNaC::symbol>(x));
        }
    }

    GiNaC::ex expanded_common_denom = common_denom.expand().collect(var);
    if (expanded_common_denom.ldegree(var) > 0)
        throw std::runtime_error("x=0 is not an analytic point of A");

    int d_degree = expanded_common_denom.degree(var);
    denom = std::vector<GiNaC::ex>(d_degree + 1);
    for (int i = 0; i <= d_degree; i++)
        denom[i] = GiNaC::ex_to<GiNaC::numeric>(expanded_common_denom.coeff(var, i)).evalf();

    GiNaC::matrix numerator_polys(N, N);
    int n_degree = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            auto numer_denom_A = A(i, j).numer_denom();
            numerator_polys(i, j) = (
                generalized_div(common_denom, numer_denom_A[1], 
                    GiNaC::ex_to<GiNaC::symbol>(var))[0] 
                * numer_denom_A[0]
            ).expand().collect(var);
            int current_n_degree = numerator_polys(i, j).degree(var);
            if (current_n_degree > n_degree)
                n_degree = current_n_degree;
        }
    }

    numer = std::vector<GiNaC::matrix>(n_degree + 1, GiNaC::matrix(N, N));
    for (int d = 0; d <= n_degree; d++)
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                numer[d](i, j) = GiNaC::ex_to<GiNaC::numeric>(numerator_polys(i, j).coeff(var, d)).evalf();

    for (int i = 0; i <= order; i++)
        sol[i] = GiNaC::matrix(N, Y0.cols());
    
    int C = Y0.cols();
    for (int i = 0; i < N; i++)
        for (int j = 0; j < C; j++)
            sol[0](i, j) = GiNaC::ex_to<GiNaC::numeric>(Y0(i, j)).evalf();
}


GiNaC::matrix finite_solver::solution() {
    int N = sol[0].rows(), C = sol[0].cols();
    int order = sol.size() - 1;
    int L_D = denom.size() - 1, L_N = numer.size() - 1;

    for (int k = 0; k < order; k++) {
        std::cerr << "DiffEqSolver: solving order " << k << "\r";
        auto n_sum = GiNaC::matrix(N, C);
        int n_terms = std::min(k, L_N);
        for (int j = 0; j <= n_terms; j++)
            n_sum = NORMAL(n_sum.add(NORMAL(numer[j].mul(sol[k - j]))));
        
        auto d_sum = GiNaC::matrix(N, C);
        int d_terms = std::min(k, L_D);
        for (int j = 1; j <= d_terms; j++)
            d_sum = NORMAL(d_sum.add(NORMAL(sol[k - j + 1].mul_scalar(denom[j] * (k - j + 1)))));
        
        sol[k + 1] = NORMAL(n_sum.sub(d_sum).mul_scalar(GiNaC::ex(1) / denom[0] / (k + 1)));
    }
    std::cerr << "\nDiffEqSolver: solving done\n";

    GiNaC::matrix M(N, C);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < C; j++) {
            M(i, j) = 0;
            for (int k = 0; k <= order; k++) {
                M(i, j) += (sol[k](i, j) * GiNaC::pow(var, k));
            }
        }
    }

    return M;
}


GiNaC::matrix run_solve::finite_solve(const GiNaC::numeric& first_pnt, const GiNaC::numeric& next_pnt, const GiNaC::matrix& init_val) {
    GiNaC::Digits = working_precision;
    std::cerr << "RunSolve: solving around " << first_pnt.evalf() << "\n";

    eq_zero.clear_all_reduction();
    GiNaC::symbol t("t");
    finite_solver fnsolver(SUBS(eq_zero.coeff, eq_zero.x == first_pnt + t), t, init_val, expansion_order);
    auto result = SUBS(fnsolver.solution(), t == eq_zero.x - first_pnt.evalf());
    auto matresult = SUBS(result, eq_zero.x == next_pnt.evalf());
    for (int i = 0; i < eq_zero.N(); i++)
        matresult(i, 0) = GiNaC::ex_to<GiNaC::numeric>(matresult(i, 0)).evalf();
  
    return matresult;
}


GiNaC::matrix run_solve::zero_solve(const GiNaC::numeric& last_pnt, const GiNaC::matrix& boundary_val) {
    GiNaC::Digits = working_precision;
    std::cerr << "RunSolve: solving around zero\n";
    auto coeff = SUBS(general_solution_zero, eq_zero.x == last_pnt);
    auto integral_constants = coeff.inverse().mul(boundary_val);

    
    // select only analytic part
    int N_cols = analytic_behavior_zero.rows();
    GiNaC::matrix final_result(N_cols, 1), temp(N_cols, 1);
    for (int i = 0; i < N_cols; i++) {
        if (GiNaC::is_nonneg_integer(GiNaC::ex_to<GiNaC::numeric>(-analytic_behavior_zero(i, i)))
         && (i == 0 || (bool)(analytic_behavior_zero(i - 1, i) == 0))) {
            if (GiNaC::abs(integral_constants(i, 0)) < GiNaC::pow(10, -20))
                continue;
            for (int j = 0; j < N_cols; j++)
                temp(j, 0) = (general_solution_zero(j, i) * integral_constants(i, 0)).expand()
                            .subs(eq_zero.x == 0, GiNaC::subs_options::algebraic);
            final_result = final_result.add(temp);
        }
    }
    return final_result;
}


GiNaC::matrix run_solve::link(const std::vector<boundary_condition>& bc) {
    auto points = get_eval_points();
    std::cerr << "RunSolve: evaluation points are " << points << "\n";
    if (points.nops() == 1) {
        auto first_step = infinity_solve(bc, GiNaC::ex_to<GiNaC::numeric>(points[0]));
        return zero_solve(GiNaC::ex_to<GiNaC::numeric>(points[0]), first_step);
    } else {
        auto first_step = infinity_solve(bc, GiNaC::ex_to<GiNaC::numeric>(points[0]));
        int npoints = points.nops();
        for (int i = 1; i < npoints; i++)
            first_step = finite_solve(GiNaC::ex_to<GiNaC::numeric>(points[i - 1]), GiNaC::ex_to<GiNaC::numeric>(points[i]), first_step);
        return zero_solve(GiNaC::ex_to<GiNaC::numeric>(points[npoints - 1]), first_step);
    }
}

