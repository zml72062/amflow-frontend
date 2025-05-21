#include "border.hpp"
#include "utils.hpp"
#include "ratsolver.h"


#define SUBS(m, r)  (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).subs((r), GiNaC::subs_options::algebraic)))
#define NORMAL(m)   (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).normal()))
#define EXPAND(m)   (GiNaC::ex_to<GiNaC::matrix>(((GiNaC::ex)(m)).expand()))


static GiNaC::matrix append_row(const GiNaC::matrix& _matrix, const GiNaC::matrix& _row) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r + 1, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
    for (int j = 0; j < c; j++)
        newmat(r, j) = _row(0, j);
    return newmat;
}


border_solve::border_solve(const GiNaC::matrix& _coeff, const GiNaC::symbol& _var, const GiNaC::symbol& _eps, const YAML::Node& _config)
    : symdiffeq(NORMAL(SUBS(_coeff, _var == 1 / _var).mul_scalar(-GiNaC::pow(_var, -2))), _var, _eps) {
    GiNaC::lst rules{_var == x, _eps == eps};
    reciprocal_coeff = SUBS(_coeff, rules);

    try_order = 20;
    eps_value = GiNaC::ex(1) / 10000;
    try {
        if (has_non_null_key(_config, "ExpansionOption")) {
            auto expansion_opt = _config["ExpansionOption"];
            if (has_non_null_key(expansion_opt, "TryOrder")) {
                try_order = expansion_opt["TryOrder"].as<int>();
            }
            if (has_non_null_key(expansion_opt, "TryEpsilon")) {
                eps_value = GiNaC::parser()(expansion_opt["TryEpsilon"].as<std::string>());
            }
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "BOrderSolve: failed while reading expansion options\n";
        exit(1);
    }
    std::cerr << "BOrderSolve: set try_order = " << try_order << "\n";
    std::cerr << "BOrderSolve: set eps_value = " << eps_value << "\n";

    raw_coeff        = SUBS(raw_coeff,        eps == eps_value);
    coeff            = SUBS(coeff,            eps == eps_value);
    transform        = SUBS(transform,        eps == eps_value);
    reciprocal_coeff = SUBS(reciprocal_coeff, eps == eps_value);
}


border_solve border_solve::from_string(const std::string& _str, const YAML::Node& _config) {
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
    return border_solve(coeff, eta, eps, _config);
}


std::vector<std::pair<GiNaC::ex, std::vector<int>>> 
border_solve::get_boundary_form() {
    clear_all_reduction();
    moser_reduction(0);

    if (!(bool)(get_moser_struct(0).moser_rank() == 1))
        throw std::runtime_error("not implemented: infinity is an irregular singularity");

    GiNaC::matrix A0(N(), N()), A_ana(N(), N());
    for (int i = 0; i < N(); i++)
        for (int j = 0; j < N(); j++)
            A0(i, j) = coeff(i, j).series(x, 1).coeff(x, -1).normal();

    reg_struct.initialize(A0);
    regular_reduction(0);
    auto J = reg_struct.J;

    // find all exponents (we do not care logarithms)
    GiNaC::lst exponent_list;
    std::vector<std::vector<int>> columns;
    for (int j = 0; j < N(); j++) {
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
    std::cerr << "BOrderSolve: allowed exponents are " << exponent_list << "\n";

    // now we can remove spurious singular part of coefficients
    for (int i = 0; i < N(); i++) {
        for (int j = 0; j < N(); j++) {
            std::cerr << "BOrderSolve: setting up coefficient matrix (" << i << ", " << j << ")          \r";
            A_ana(i, j) = (coeff(i, j) - J(i, j) / x).normal();
        }
    }
    std::cerr << "\nBOrderSolve: done setting up coefficient matrix\n";

    ratsolver_regular solver(A_ana, x, try_order, reg_struct);
    GiNaC::symbol l("l");
    auto out = SUBS(EXPAND(NORMAL(transform.mul(solver.solution()))), GiNaC::log(x) == l);
    int num_behaviors = columns.size();
    std::vector<std::pair<GiNaC::ex, std::vector<int>>> result;
    for (int i = 0; i < num_behaviors; i++) {
        auto exponent = exponent_list[i];
        auto span     = columns[i];
        int  crank    = 0, 
             lspan    = span.size(),
             order    = 0;
        GiNaC::matrix coeff(0, lspan);
        std::vector<int> master_orders(N(), -1);
        while (crank < lspan) {
            for (int j = 0; j < N(); j++) {
                GiNaC::matrix next_row(1, lspan);
                for (int q = 0; q < lspan; q++)
                    next_row(0, q) = (out(j, span[q]) / GiNaC::pow(x, exponent))
                        .expand().collect(l).coeff(l, 0)
                                 .collect(x).coeff(x, order);
                
                auto new_coeff = append_row(coeff, next_row);
                if ((int)new_coeff.rank() > crank) {
                    master_orders[j] = order;
                    coeff = new_coeff;
                    crank++;
                }
                if (crank == lspan)
                    break;
            }
            if (crank == lspan)
                break;
                
            order++;
            if (order >= try_order - 3) {
                std::cerr << "BOrderSolve: try_order too small, try increasing TryOrder\n";
                exit(1);
            }
        }
        result.push_back(std::make_pair(exponent, master_orders));
    }

    return result;
}


