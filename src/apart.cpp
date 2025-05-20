#include "apart.hpp"


/**
 * Add all factors of `poly` to `lst`.
 */
static void set_factors(const GiNaC::ex& poly, GiNaC::lst& lst) {
    if (!GiNaC::is_exactly_a<GiNaC::mul>(poly)) {
        lst.append(poly);
        return;
    }

    for (auto iter = poly.begin(); iter != poly.end(); ++iter)
        set_factors(*iter, lst);
}


/**
 * Check whether `factor` is a linear function of `vars`.
 */
static bool is_a_correct_factor(const GiNaC::ex& factor, const GiNaC::lst& vars) {
    for (auto& var: vars) {
        if (!GiNaC::is_exactly_a<GiNaC::numeric>(factor.diff(GiNaC::ex_to<GiNaC::symbol>(var)).expand()))
            return false;
    }
    return true;
}


GiNaC::lst generalized_div(const GiNaC::ex& numer, const GiNaC::ex& denom, const GiNaC::symbol& var) {
    if (numer.degree(var) < denom.degree(var)) {
        return GiNaC::lst{0, numer};
    }

    if (denom.degree(var) == 0) {
        return GiNaC::lst{numer / denom, 0};
    }

    int quo_deg = numer.degree(var) - denom.degree(var),
        rem_deg = denom.degree(var) - 1;
    
    GiNaC::ex  quo_form = 0, rem_form = 0;
    GiNaC::lst param_list;
    for (int i = 0; i <= quo_deg; i++) {
        GiNaC::symbol quoi("quo" + std::to_string(i));
        quo_form += (quoi * GiNaC::pow(var, i));
        param_list.append(quoi);
    }
    for (int i = 0; i <= rem_deg; i++) {
        GiNaC::symbol remi("rem" + std::to_string(i));
        rem_form += (remi * GiNaC::pow(var, i));
        param_list.append(remi);
    }
    
    int total = numer.degree(var);
    GiNaC::ex  lhs_form = (quo_form * denom + rem_form).expand().collect(var);
    GiNaC::lst eqns;
    for (int j = 0; j <= total; j++) {
        eqns.append(lhs_form.coeff(var, j) == numer.coeff(var, j));
    }
    auto result = GiNaC::lsolve(eqns, param_list);
    
    return GiNaC::lst{quo_form.subs(result, GiNaC::subs_options::algebraic),
                      rem_form.subs(result, GiNaC::subs_options::algebraic)};
}


std::vector<parfrac> parfrac_expansion(const GiNaC::ex& frac, const GiNaC::symbol& var) {
    auto numer_denom = frac.normal().numer_denom();
    auto numer = numer_denom[0], denom = numer_denom[1];
    auto quo_rem = generalized_div(numer, denom, var);
    auto quo = quo_rem[0], rem = quo_rem[1];

    if ((bool)(rem.normal() == 0)) {
        if ((bool)(quo.normal() == 0))
            return std::vector<parfrac>();
        
        std::vector<parfrac> final_result;
        auto quo_numer_denom = quo.normal().numer_denom();
        auto quo_numer = quo_numer_denom[0].expand().collect(var),
             quo_denom = quo_numer_denom[1];
        int  quo_deg = quo_numer.degree(var);
        for (int i = 0; i <= quo_deg; i++) {
            auto numer_coeff = quo_numer.coeff(var, i);
            if ((bool)(numer_coeff == 0))
                continue;
            parfrac par;
            par.coeff = numer_coeff / quo_denom;
            par.denoms.append(var);
            par.npowers.push_back(-i);
            final_result.push_back(par);
        }
        return final_result;
    }
    
    auto factorized = GiNaC::factor(denom);
    GiNaC::lst factor_list;
    set_factors(factorized, factor_list);
    GiNaC::lst param_list, elem_list, info_list;
    int symcnt = 0;
    for (auto& factor: factor_list) {
        if ((bool)(factor.diff(var) == 0))
            continue;
        
        if (GiNaC::is_exactly_a<GiNaC::power>(factor)) {
            auto base = factor.begin()->expand();
            if (!is_a_correct_factor(base, {var})) {
                std::cerr << "Illegal denominators detected!\n";
                exit(1);
            }
            auto iter = factor.begin();
            ++iter;
            int exp = GiNaC::ex_to<GiNaC::numeric>(iter->expand()).to_int();
            for (int i = 1; i <= exp; i++) {
                elem_list.append(GiNaC::pow(base, -i));
                param_list.append(GiNaC::symbol("c" + std::to_string(symcnt++)));
                info_list.append(GiNaC::lst{base, i});
            }
        } else {
            if (!is_a_correct_factor(factor, {var})) {
                std::cerr << "Illegal denominators detected!\n";
                exit(1);
            }
            elem_list.append(1 / factor);
            param_list.append(GiNaC::symbol("c" + std::to_string(symcnt++)));
            info_list.append(GiNaC::lst{factor, 1});
        }
    }

    GiNaC::ex try_ex = 0;
    int N_elem = elem_list.nops();
    for (int i = 0; i < N_elem; i++)
        try_ex += (elem_list[i] * param_list[i]);
    try_ex = (try_ex * factorized).normal().expand().collect(var);
    int deg = try_ex.degree(var);
    
    auto rem_numer_denom = rem.normal().numer_denom();
    auto rem_numer = rem_numer_denom[0].expand().collect(var),
         rem_denom = rem_numer_denom[1];
    GiNaC::lst eqs;
    for (int i = 0; i <= deg; i++)
        eqs.append(try_ex.coeff(var, i) == rem_numer.coeff(var, i));
    auto result = GiNaC::ex_to<GiNaC::lst>(GiNaC::lsolve(eqs, param_list));
    
    std::vector<parfrac> final_result;
    for (int i = 0; i < N_elem; i++) {
        parfrac par;
        par.coeff = param_list[i].subs(result, GiNaC::subs_options::algebraic);
        if (!(bool)(par.coeff == 0)) {
            par.coeff = par.coeff / rem_denom;
            par.denoms.append(info_list[i][0]);
            par.npowers.push_back(GiNaC::ex_to<GiNaC::numeric>(info_list[i][1]).to_int());
            final_result.push_back(par);
        }
    }

    if (!(bool)(quo.normal() == 0)) {
        auto quo_numer_denom = quo.normal().numer_denom();
        auto quo_numer = quo_numer_denom[0].expand().collect(var),
             quo_denom = quo_numer_denom[1];
        int  quo_deg = quo_numer.degree(var);
        for (int i = 0; i <= quo_deg; i++) {
            auto numer_coeff = quo_numer.coeff(var, i);
            if ((bool)(numer_coeff == 0))
                continue;
            parfrac par;
            par.coeff = numer_coeff / quo_denom;
            par.denoms.append(var);
            par.npowers.push_back(-i);
            final_result.push_back(par);
        }
    }

    return final_result;
}


std::vector<parfrac> parfrac_expansion(const GiNaC::ex& frac, const GiNaC::lst& vars) {
    std::vector<parfrac> init_list;
    parfrac init;
    init.coeff = frac;
    init_list.push_back(init);
    for (auto& var: vars) {
        std::vector<parfrac> new_list;
        for (auto& par: init_list) {
            auto new_result = parfrac_expansion(par.coeff, GiNaC::ex_to<GiNaC::symbol>(var));
            for (auto& app: new_result) {
                parfrac to_append;
                to_append.coeff   = app.coeff;
                to_append.denoms  = par.denoms;
                to_append.denoms.append(app.denoms[0]);
                to_append.npowers = par.npowers;
                to_append.npowers.push_back(app.npowers[0]);
                new_list.push_back(to_append);
            }
        }
        init_list = new_list;
    }
    return init_list;
}


