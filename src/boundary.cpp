#include "boundary.hpp"
#include "apart.hpp"
#include "kira.hpp"
#include "utils.hpp"


static std::vector<GiNaC::lst> all_combinations(const GiNaC::lst& haystack, int r) {
    if (r == 0) {
        std::vector<GiNaC::lst> result;
        result.push_back(GiNaC::lst());
        return result;
    }

    int n = haystack.nops();
    if (n < r)
        return std::vector<GiNaC::lst>();
    
    std::vector<GiNaC::lst> result;
    for (int i = 0; i < n - r + 1; i++) {
        GiNaC::lst new_haystack;
        for (int j = i + 1; j < n; j++)
            new_haystack.append(haystack[j]);
        auto partial_combinations = all_combinations(new_haystack, r - 1);
        for (auto& p_comb: partial_combinations) {
            GiNaC::lst lneedle;
            lneedle.append(haystack[i]);
            for (auto& needle: p_comb)
                lneedle.append(needle);
            result.push_back(lneedle);
        }
    }
    return result;
}


/**
 * Determine whether a group of new loop momenta `new_loops` are 
 * independent. If so, representing the old loop momenta `original_loops`
 * in terms of a linear combination of `new_loops`, and return the
 * coefficient matrix.
 */
static std::pair<bool, GiNaC::matrix> original_loops_by_new_loops(const GiNaC::lst& new_loops, const GiNaC::lst& original_loops) {
    int L = new_loops.nops();
    GiNaC::matrix coeff(L, L);
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            coeff(i, j) = new_loops[i].diff(GiNaC::ex_to<GiNaC::symbol>(original_loops[j])).expand();

    if ((bool)(coeff.determinant().expand() == 0))
        return std::make_pair(false, GiNaC::matrix(0, 0));
    else
        return std::make_pair(true, coeff.inverse());
}


static std::vector<std::vector<int>> all_bitset(int N) {
    std::vector<std::vector<int>> base;
    if (N == 0) {
        base.push_back(std::vector<int>());
        return base;
    }

    auto one_less_zero = all_bitset(N - 1), 
         one_less_one  = all_bitset(N - 1);
    
    for (auto& v: one_less_zero)
        v.push_back(0);
    
    for (auto& v: one_less_one)
        v.push_back(1);

    std::size_t sz = one_less_zero.size();
    for (std::size_t i = 0; i != sz; i++) {
        base.push_back(one_less_zero[i]);
        base.push_back(one_less_one[i]);
    }
    return base;
}


static bool is_large_momentum(const GiNaC::matrix& expansion_coeff, const std::vector<int>& large_rule) {
    int L = expansion_coeff.cols();
    bool large = false;
    for (int i = 0; i < L; i++) {
        if (large_rule[i] && !(bool)(expansion_coeff(0, i).expand() == 0)) {
            large = true;
            break;
        }
    }
    return large;
}


std::vector<boundary_region> boundary_region::all_boundary_regions(integralfamily& family, const std::vector<int>& top_sector) {
    if (top_sector.size() != family.propagators.nops()) {
        std::cerr << "BoundaryRegion: wrong top sector input\n";
        exit(1);
    }

    /* loop momenta flowing on each propagator */
    GiNaC::lst all_loop_momenta = family.loop_momenta();
    
    GiNaC::lst interested_loop_momenta,     /* loop momenta on each top-level propagator */
               interested_branches;         /* loop momenta on each top-level branch */
    int N = top_sector.size();
    for (int i = 0; i < N; i++) {
        if (top_sector[i]) {
            interested_loop_momenta.append(all_loop_momenta[i]);
            bool already_exist = false;
            for (auto& q_exist: interested_branches) {
                if ((bool)((q_exist - all_loop_momenta[i]).expand() == 0)) {
                    already_exist = true;
                    break;
                }
            }
            if (!already_exist)
                interested_branches.append(all_loop_momenta[i]);
        }
    }

    /* ways to select L branches as potentially independent loop momenta */
    auto choose_loops = all_combinations(interested_branches, family.loops.nops());
    /* transformation from original loop momenta to new loop momenta */
    std::vector<GiNaC::matrix> transforms;
    for (auto& choice: choose_loops) {
        auto result = original_loops_by_new_loops(choice, family.loops);
        if (result.first)
            transforms.push_back(result.second);
    }

    /* ways to assign large/small to each new loop momentum */
    auto all_large_rules = all_bitset(family.loops.nops());
    
    /* index into "interested_branches" where the branch is cut */
    std::set<int> cut_position;
    if (family.cut.size() != 0) {
        GiNaC::lst cut_momenta;
        int N_props = all_loop_momenta.nops();
        for (int i = 0; i < N_props; i++) {
            if (family.cut[i])
                cut_momenta.append(all_loop_momenta[i]);
        }

        int N_branches = interested_branches.nops();
        for (int j = 0; j < N_branches; j++) {
            auto this_branch_momentum = interested_branches[j];
            for (auto& cut_q: cut_momenta) {
                if ((bool)((this_branch_momentum - cut_q).expand() == 0)) {
                    cut_position.insert(j);
                    break;
                }
            }
        }
    }

    std::map<std::vector<int>, std::pair<GiNaC::matrix, std::vector<int>>> all_regions;
    for (auto& trans: transforms) {
        std::vector<GiNaC::matrix> expansion_matrices;
        int NL = family.loops.nops();
        GiNaC::matrix original_coeff(1, NL);
        for (auto& branch: interested_branches) {
            for (int i = 0; i < NL; i++)
                original_coeff(0, i) = branch.diff(GiNaC::ex_to<GiNaC::symbol>(family.loops[i]));            
            expansion_matrices.push_back(original_coeff.mul(trans));
        } 

        for (auto& rule: all_large_rules) {
            int NB = expansion_matrices.size();
            std::vector<int> large_situation;
            bool cut_invalid = false;
            for (int i = 0; i < NB; i++) {
                large_situation.push_back((int)is_large_momentum(expansion_matrices[i], rule));
                if (cut_position.find(i) != cut_position.end() && large_situation[i]) {
                    cut_invalid = true;
                    break;
                }
            }
            if (!cut_invalid && all_regions.find(large_situation) == all_regions.end())
                all_regions[large_situation] = std::make_pair(trans, rule);
        }
    }

    std::vector<boundary_region> final_regions;
    for (auto& pair_: all_regions) {
        boundary_region reg(family);
        reg.transform = pair_.second.first;
        reg.is_large  = pair_.second.second;
        final_regions.push_back(reg);
    }
    return final_regions;
}


bool boundary_region::is_trivial_region(const std::vector<int>& top_sector,
                                        const char* workdir, const YAML::Node& config) {
    integralfamily presubf = preliminary_subfamily(top_sector);
    GiNaC::lst props_from_original;
    int nprop = top_sector.size();
    for (int i = 0; i < nprop; i++)
        if (top_sector[i])
            props_from_original.append(new_leading_propagators[0][i]);
    
    unsigned long new_top_sector = 0;
    for (int i = 0; i < nprop; i++) {
        auto propi = presubf.propagators[i];
        bool exist = false;
        for (auto& original_prop: props_from_original) {
            if ((bool)((original_prop - propi).expand() == 0)) {
                exist = true;
                break;
            }
        }
        if (exist)
            new_top_sector |= (1ul << i);
    }

    std::cerr << "BoundaryRegion: boundary region has top sector " << new_top_sector << "\n";
    
    kira_agent agent(config, presubf);
    auto trivial_sectors = agent.trivial_sectors(workdir, new_top_sector);
    for (auto& trivsect: trivial_sectors)
        if (trivsect == new_top_sector)
            return true;
    
    return false;
}


void boundary_region::compute_new_propagators() {
    if (new_propagators.nops() != pfamily->propagators.nops()) {
        new_propagators = GiNaC::lst();
        int L = pfamily->loops.nops();
        GiNaC::matrix loops_matrix(L, 1);
        for (int i = 0; i < L; i++)
            loops_matrix(i, 0) = pfamily->loops[i];
        GiNaC::matrix new_loops_matrix = transform.mul(loops_matrix);
        GiNaC::lst loops_transform_rules;
        for (int i = 0; i < L; i++)
            loops_transform_rules.append(pfamily->loops[i] == new_loops_matrix(i, 0).expand());
        
        for (auto& prop: pfamily->propagators)
            new_propagators.append(prop.subs(loops_transform_rules, GiNaC::subs_options::algebraic));
    }
}


void boundary_region::compute_new_leading_propagators() {
    compute_new_propagators();
    if (new_leading_propagators[0].nops() != new_propagators.nops()) {
        new_leading_propagators[0] = GiNaC::lst();
        new_leading_propagators[1] = GiNaC::lst();
        new_leading_propagators[2] = GiNaC::lst();
        leading_factors            = GiNaC::lst();

        int L = pfamily->loops.nops();
        for (auto& prop: new_propagators) {
            GiNaC::symbol sqrteta("sqrteta");
            GiNaC::lst    rules;
            
            for (int i = 0; i < L; i++)
                if (is_large[i])
                    rules.append(pfamily->loops[i] == pfamily->loops[i] * sqrteta);
            rules.append((*pfamily->psymbols)["eta"] == GiNaC::pow(sqrteta, 2));

            auto prop_after_sub = prop.subs(rules, GiNaC::subs_options::algebraic).expand().collect(sqrteta);
            
            auto square_coeff = prop_after_sub.coeff(sqrteta, 2).expand();
            if ((bool)(square_coeff == 0)) {
                new_leading_propagators[0].append(prop);
                new_leading_propagators[1].append(0);
                new_leading_propagators[2].append(0);
                leading_factors.append(1);
            } else {
                new_leading_propagators[0].append(square_coeff);
                new_leading_propagators[1].append(prop_after_sub.coeff(sqrteta, 1).expand());
                new_leading_propagators[2].append(prop_after_sub.coeff(sqrteta, 0).expand());
                leading_factors.append((*pfamily->psymbols)["eta"]);
            }
        }
    }
}


GiNaC::ex boundary_region::compute_expansion_in_new_leading_propagators(const GiNaC::ex& prop) {
    compute_new_leading_propagators();
    if (new_propagator_symbols.nops() != new_propagators.nops()) {
        new_propagator_symbols = GiNaC::lst();
        int n = new_propagators.nops();
        for (int i = 0; i < n; i++)
            new_propagator_symbols.append(GiNaC::symbol("D" + std::to_string(i)));
    }

    // expand in scalar products
    int n_loops = pfamily->loops.nops(), n_legs = pfamily->indeplegs.nops();
    int S = n_loops * (n_loops + 1) / 2 + n_loops * n_legs;
    GiNaC::matrix coeff_sps(1, S);
    GiNaC::ex     bias_sps;
    int q = 0;
    for (int i = 0; i < n_loops; i++) {
        GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->loops[i])).expand();
        coeff_sps(0, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->loops[i])).expand() / 2;
        for (int j = i + 1; j < n_loops; j++)
            coeff_sps(0, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->loops[j])).expand();
    }
    for (int i = 0; i < n_loops; i++) {
        GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->loops[i])).expand();
        for (int j = 0; j < n_legs; j++)
            coeff_sps(0, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(pfamily->indeplegs[j])).expand();
    }

    GiNaC::lst zero_rules;
    for (auto& lp: pfamily->loops)
        zero_rules.append(lp == 0);
    bias_sps = prop.subs(zero_rules, GiNaC::subs_options::algebraic)
                   .subs(pfamily->sps_numeric_rules, GiNaC::subs_options::algebraic)
                   .subs(pfamily->numeric_rules, GiNaC::subs_options::algebraic);
    
    auto coeff_new = coeff_sps.mul(coeff_sps_to_new_props);
    auto bias_new  = coeff_sps.mul(bias_sps_to_new_props)(0, 0) + bias_sps;
    for (int i = 0; i < S; i++)
        bias_new += coeff_new(0, i) * new_propagator_symbols[i];
    return bias_new.expand();
}


integralfamily boundary_region::preliminary_subfamily(const std::vector<int>& top_sector) {
    if (top_sector.size() != pfamily->propagators.nops()) {
        std::cerr << "BoundaryRegion: wrong top sector input\n";
        exit(1);
    }
    
    integralfamily newfamily;
    compute_new_leading_propagators();

    // select leading propagators corresponding to the top sector
    GiNaC::lst top_level_leading_props;
    int nprop = top_sector.size();
    for (int i = 0; i < nprop; i++)
        if (top_sector[i])
            top_level_leading_props.append(new_leading_propagators[0][i]);

    std::vector<int> newcut;
    if (pfamily->cut.size() != 0)
        for (int i = 0; i < nprop; i++)
            if (top_sector[i])
                newcut.push_back(pfamily->cut[i]);
    
    newfamily.name         = pfamily->name;
    newfamily.psymbols     = pfamily->psymbols;
    newfamily.loops        = pfamily->loops;
    newfamily.indeplegs    = pfamily->indeplegs;
    newfamily.propagators  = top_level_leading_props;
    newfamily.cut          = newcut;
    newfamily.prescription = pfamily->prescription;

    newfamily.conservation_rules = pfamily->conservation_rules;
    newfamily.numeric_rules      = pfamily->numeric_rules;
    newfamily.sps_numeric_rules  = pfamily->sps_numeric_rules;

    auto complete_info = newfamily.sps_in_complete_propagator();
    newfamily.propagators  = complete_info.prop;
    newfamily.cut          = complete_info.cut;
    memo_sps               = complete_info.sps;
    coeff_sps_to_new_props = complete_info.coeff;
    bias_sps_to_new_props  = complete_info.bias;
    
    return newfamily;
}


GiNaC::ex boundary_region::exponent(const integral& intgr) {
    compute_new_leading_propagators();
    GiNaC::ex f = 0;

    int total_large = 0;
    for (auto& largeq: is_large)
        total_large += largeq;
    f += (total_large * (*pfamily->psymbols)["d"] / 2);

    int N_ind = intgr.indices.size();
    for (int i = 0; i < N_ind; i++) {
        if ((bool)((leading_factors[i] - 1).expand() == 0))
            continue;
        f -= intgr.indices[i];
    }

    return f;
}


GiNaC::ex boundary_region::factor(const integral& intgr) {
    compute_new_leading_propagators();
    GiNaC::ex f = GiNaC::pow(GiNaC::abs(transform.determinant()), (*pfamily->psymbols)["d"]);
    
    int total_large = 0;
    for (auto& largeq: is_large)
        total_large += largeq;
    f *= GiNaC::pow((*pfamily->psymbols)["eta"], total_large * (*pfamily->psymbols)["d"] / 2);
    
    int N_ind = intgr.indices.size();
    for (int i = 0; i < N_ind; i++)
        f /= GiNaC::pow(leading_factors[i], intgr.indices[i]);

    return f;
}


static void process_an_integrand(const GiNaC::ex& integrand, const GiNaC::lst& denoms, 
                                 std::vector<std::pair<GiNaC::lst, std::vector<std::vector<int>>>>& raw_integrals,
                                 std::map<std::pair<int, int>, GiNaC::ex>& current_result) {
    auto expansion = parfrac_expansion(integrand, denoms);
    for (auto& par: expansion) {
        int current_nfamily = raw_integrals.size();
        int exist_pos = -1;
        for (int i = 0; i < current_nfamily; i++) {
            if ((bool)(raw_integrals[i].first == par.denoms)) {
                exist_pos = i;
                break;
            }
        }
        if (exist_pos == -1) {
            raw_integrals.push_back({par.denoms, {}});
            exist_pos = raw_integrals.size() - 1;
        }

        int current_nintegral = raw_integrals[exist_pos].second.size();
        int exist_intgr_pos = -1;
        for (int i = 0; i < current_nintegral; i++) {
            if (raw_integrals[exist_pos].second[i] == par.npowers) {
                exist_intgr_pos = i;
                break;
            }
        }
        if (exist_intgr_pos == -1) {
            raw_integrals[exist_pos].second.push_back(par.npowers);
            exist_intgr_pos = raw_integrals[exist_pos].second.size() - 1;
        }

        if (current_result.find({exist_pos, exist_intgr_pos}) == current_result.end()) {
            current_result[{exist_pos, exist_intgr_pos}] = par.coeff;
        } else {
            current_result[{exist_pos, exist_intgr_pos}] += par.coeff;
        }
    }
}


std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
          std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>
boundary_region::subintegrals(const std::vector<int>& top_sector, const integral& intgr, int border) {
    if (border < 0)
        return {{}, {}};

    integralfamily presubf = preliminary_subfamily(top_sector);
    GiNaC::lst expanded_leading_propagators[3];
    int n = new_leading_propagators[0].nops();
    for (int q = 0; q < 3; q++)
        for (int i = 0; i < n; i++)
            expanded_leading_propagators[q].append(compute_expansion_in_new_leading_propagators(new_leading_propagators[q][i]));
    
    GiNaC::symbol invsqrteta("invsqrteta");
    GiNaC::ex     res = 1;
    for (int i = 0; i < n; i++) {
        res *= GiNaC::pow(
            expanded_leading_propagators[0][i] 
              + invsqrteta * expanded_leading_propagators[1][i] 
              + GiNaC::pow(invsqrteta, 2) * expanded_leading_propagators[2][i], 
        -intgr.indices[i]);
    }
    res = GiNaC::series_to_poly(res.series(invsqrteta, border * 2 + 1));
    
    std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
              std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>  final_result;
    std::vector<std::pair<GiNaC::lst, std::vector<std::vector<int>>>> raw_integrals;
    for (int i = 0; i <= border; i++) {
        GiNaC::ex integrand = res.coeff(invsqrteta, i * 2).expand();
        std::map<std::pair<int, int>, GiNaC::ex> current_order_result;
        if (GiNaC::is_exactly_a<GiNaC::add>(integrand)) {
            for (auto iter = integrand.begin(); iter != integrand.end(); ++iter) {
                process_an_integrand(*iter, new_propagator_symbols, raw_integrals, current_order_result);
            }
        } else {
            process_an_integrand(integrand, new_propagator_symbols, raw_integrals, current_order_result);
        }
        final_result.second.push_back(current_order_result);
    }
    
    for (auto& family_integral: raw_integrals) {
        auto raw_denoms = family_integral.first;
        std::shared_ptr<integralfamily> subf = std::make_shared<integralfamily>();
        subf->name      = presubf.name;
        subf->psymbols  = presubf.psymbols;
        subf->loops     = presubf.loops;
        subf->indeplegs = presubf.indeplegs;

        GiNaC::lst inv_subs_rules;
        int n_propagators = presubf.propagators.nops();
        for (int i = 0; i < n_propagators; i++)
            inv_subs_rules.append(new_propagator_symbols[i] == presubf.propagators[i]);
        subf->propagators = GiNaC::ex_to<GiNaC::lst>(
            ((GiNaC::ex)raw_denoms).subs(inv_subs_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.conservation_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.numeric_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.sps_numeric_rules, GiNaC::subs_options::algebraic).expand()
        );
        subf->cut                = presubf.cut;
        subf->prescription       = presubf.prescription;
        subf->conservation_rules = presubf.conservation_rules;
        subf->numeric_rules      = presubf.numeric_rules;
        subf->sps_numeric_rules  = presubf.sps_numeric_rules;

        final_result.first.push_back({subf, {}});
        for (auto& intgr_idx: family_integral.second)
            final_result.first.back().second.push_back(integral(*subf, intgr_idx));
    }

    return final_result;
}


std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
          std::vector<std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>>
boundary_region::subintegrals(const std::vector<int>& top_sector, const integral_list& intgrs, const std::vector<int>& border) {
    if (intgrs.size() != border.size()) {
        std::cerr << "BoundaryRegion: master integral input and border input do not match length\n";
        exit(1);
    }
    
    integralfamily presubf = preliminary_subfamily(top_sector);
    GiNaC::lst expanded_leading_propagators[3];
    int n = new_leading_propagators[0].nops();
    for (int q = 0; q < 3; q++)
        for (int i = 0; i < n; i++)
            expanded_leading_propagators[q].append(compute_expansion_in_new_leading_propagators(new_leading_propagators[q][i]));

    std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
              std::vector<std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>> final_result;
    std::vector<std::pair<GiNaC::lst, std::vector<std::vector<int>>>>             raw_integrals;
    int nmasters = intgrs.size();
    for (int master_idx = 0; master_idx < nmasters; master_idx++) {
        final_result.second.push_back({});
        if (border[master_idx] < 0)
            continue;
        
        GiNaC::symbol invsqrteta("invsqrteta");
        GiNaC::ex     res = 1;
        for (int i = 0; i < n; i++) {
            res *= GiNaC::pow(
                expanded_leading_propagators[0][i] 
                + invsqrteta * expanded_leading_propagators[1][i] 
                + GiNaC::pow(invsqrteta, 2) * expanded_leading_propagators[2][i], 
            -intgrs[master_idx].indices[i]);
        }
        res = GiNaC::series_to_poly(res.series(invsqrteta, border[master_idx] * 2 + 1));

        for (int i = 0; i <= border[master_idx]; i++) {
            GiNaC::ex integrand = res.coeff(invsqrteta, i * 2).expand();
            std::map<std::pair<int, int>, GiNaC::ex> current_order_result;
            if (GiNaC::is_exactly_a<GiNaC::add>(integrand)) {
                for (auto iter = integrand.begin(); iter != integrand.end(); ++iter) {
                    process_an_integrand(*iter, new_propagator_symbols, raw_integrals, current_order_result);
                }
            } else {
                process_an_integrand(integrand, new_propagator_symbols, raw_integrals, current_order_result);
            }
            final_result.second.back().push_back(current_order_result);
        }
    }
    
    for (auto& family_integral: raw_integrals) {
        auto raw_denoms = family_integral.first;
        std::shared_ptr<integralfamily> subf = std::make_shared<integralfamily>();
        subf->name      = presubf.name;
        subf->psymbols  = presubf.psymbols;
        subf->loops     = presubf.loops;
        subf->indeplegs = presubf.indeplegs;

        GiNaC::lst inv_subs_rules;
        int n_propagators = presubf.propagators.nops();
        for (int i = 0; i < n_propagators; i++)
            inv_subs_rules.append(new_propagator_symbols[i] == presubf.propagators[i]);
        subf->propagators = GiNaC::ex_to<GiNaC::lst>(
            ((GiNaC::ex)raw_denoms).subs(inv_subs_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.conservation_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.numeric_rules, GiNaC::subs_options::algebraic).expand()
                                   .subs(presubf.sps_numeric_rules, GiNaC::subs_options::algebraic).expand()
        );
        subf->cut                = presubf.cut;
        subf->prescription       = presubf.prescription;
        subf->conservation_rules = presubf.conservation_rules;
        subf->numeric_rules      = presubf.numeric_rules;
        subf->sps_numeric_rules  = presubf.sps_numeric_rules;

        final_result.first.push_back({subf, {}});
        for (auto& intgr_idx: family_integral.second)
            final_result.first.back().second.push_back(integral(*subf, intgr_idx));
    }

    return final_result;
}


std::string boundary_region::to_string() const {
    std::ostringstream out;
    int nloops = pfamily->loops.nops();
    GiNaC::matrix loop_mat(nloops, 1);
    GiNaC::symbol sqrteta("sqrteta");
    for (int i = 0; i < nloops; i++) {
        loop_mat(i, 0) = pfamily->loops[i];
        if (is_large[i])
            loop_mat(i, 0) *= sqrteta;
    }
    auto rhs = transform.mul(loop_mat);
    GiNaC::lst rules;
    for (int i = 0; i < nloops; i++)
        rules.append(pfamily->loops[i] == rhs(i, 0));
    out << rules;
    return out.str();
}

