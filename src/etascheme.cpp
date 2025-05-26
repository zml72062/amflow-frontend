#include "family.hpp"
#include "utils.hpp"
#include <algorithm>


static bool is_cut(const component& comp, int prop_idx) {
    return comp.pfamily->cut.size() != 0 && comp.pfamily->cut[prop_idx];
}


static int is_subset(const GiNaC::lst& sublist, const std::vector<int>& fulllist) {
    std::set<int> fullset(fulllist.begin(), fulllist.end());
    for (auto& idx: sublist) {
        auto int_idx = GiNaC::ex_to<GiNaC::numeric>(idx).to_int();
        if (fullset.find(int_idx) == fullset.end())
            return -1;
        fullset.erase(int_idx);
    }
    return *fullset.begin();
}


static std::vector<int> eta_position(const component& comp, eta_scheme scheme) {
    switch (scheme) {
        case Prescription: {
            std::cerr << "EtaScheme: currently, we do not support Prescription\n";
            return {};
        }
        case Mass: {
            GiNaC::lst mass_list = comp.pfamily->masses();
            // (mass, corresponding list of massive propagators)
            std::vector<std::pair<GiNaC::ex, std::vector<int>>> grouped_propagators;
            for (auto& idx: comp.propagator_indices) {
                if (is_cut(comp, idx))
                    continue;

                auto current_mass = mass_list[idx].expand();
                // mass must be non-zero
                if ((bool)(current_mass == 0))
                    continue;
                
                int existing = -1, cur_len = grouped_propagators.size();
                for (int i = 0; i < cur_len; i++) {
                    if ((bool)((grouped_propagators[i].first - current_mass).expand() == 0)) {
                        existing = i;
                        grouped_propagators[i].second.push_back(idx);
                        break;
                    }
                }
                if (existing == -1)
                    grouped_propagators.push_back({current_mass, {idx}});
            }

            // sort the propagator groups such that groups with fewer
            // propagators are at the front
            std::sort(grouped_propagators.begin(), grouped_propagators.end(),
                [](const auto& x, const auto& y) {
                    return x.second.size() < y.second.size();
                }
            );

            if (grouped_propagators.size() == 0)
                return {};
            else
                return grouped_propagators[0].second;
        }
        case Propagator: {
            // (branch length, propagator index)
            std::vector<std::pair<int, int>> branch_length;
            for (auto& idx: comp.propagator_indices) {
                if (is_cut(comp, idx))
                    continue;

                auto collected_U = comp.U.expand().collect(comp.pfamily->feynman_params[idx]);
                auto coeff_idx   = collected_U.coeff(comp.pfamily->feynman_params[idx], 1);
                std::vector<int> branch_indices;
                for (auto& bidx: comp.propagator_indices) {
                    if (!coeff_idx.has(comp.pfamily->feynman_params[bidx]))
                        branch_indices.push_back(bidx);
                }
                branch_length.push_back({branch_indices.size(), idx});
            }

            // sort the propagators such that propagators with longer
            // branch length are at the front
            std::sort(branch_length.begin(), branch_length.end(),
                [](const auto& x, const auto& y) {
                    return x.first > y.first;
                }
            );

            if (branch_length.size() == 0)
                return {};
            else
                return {branch_length[0].second};
        }
        case Branch: {
            int nprops = comp.propagator_indices.size();
            std::vector<int> classified(nprops, 0);
            std::vector<std::vector<int>> branches;
            for (int i = 0; i < nprops; i++) {
                if (is_cut(comp, comp.propagator_indices[i])) {
                    classified[i] = 1;
                    continue;
                }

                if (!classified[i]) {
                    classified[i] = 1;
                    branches.push_back({comp.propagator_indices[i]});

                    // find all propagators that are in the same branch
                    // as "comp.propagator_indices[i]"
                    auto collected_U = comp.U.expand().collect(comp.pfamily->feynman_params[comp.propagator_indices[i]]);
                    auto coeff_idx   = collected_U.coeff(comp.pfamily->feynman_params[comp.propagator_indices[i]], 1);
                    for (int j = i + 1; j < nprops; j++) {
                        if (!classified[j] && !coeff_idx.has(comp.pfamily->feynman_params[comp.propagator_indices[j]])) {
                            classified[j] = 1;
                            branches.back().push_back(comp.propagator_indices[j]);
                        }
                    }
                }
            }

            // sort the branches such that shorter branches are at the front
            std::sort(branches.begin(), branches.end(),
                [](const auto& x, const auto& y) {
                    return x.size() < y.size();
                }
            );
            
            if (branches.size() == 0)
                return {};
            else
                return branches[0];
        }
        case Loop: {
            std::vector<std::vector<int>> U_terms;
            auto U_expanded = comp.U.expand();
            if (GiNaC::is_exactly_a<GiNaC::add>(U_expanded)) {
                for (auto iter = U_expanded.begin(); iter != U_expanded.end(); ++iter) {
                    std::vector<int> current_term;
                    for (auto& idx: comp.propagator_indices) {
                        if (iter->has(comp.pfamily->feynman_params[idx])) {
                            current_term.push_back(idx);
                        }
                    }
                    U_terms.push_back(current_term);
                }
            } else {
                std::vector<int> current_term;
                for (auto& idx: comp.propagator_indices) {
                    if (U_expanded.has(comp.pfamily->feynman_params[idx])) {
                        current_term.push_back(idx);
                    }
                }
                U_terms.push_back(current_term);
            }

            int nloops = comp.num_loops();
            GiNaC::lst comp_feynman_params;
            for (auto& idx: comp.propagator_indices) {
                comp_feynman_params.append(idx);
            }
            auto sub_products = all_combinations(comp_feynman_params, nloops - 1);
            
            std::vector<std::vector<int>> loops;
            for (auto& sub_prod: sub_products) {
                std::set<int> prop_set;
                for (auto& term: U_terms) {
                    auto rest = is_subset(sub_prod, term);
                    if (rest == -1)
                        continue;
                    // "rest" should not be cut
                    if (is_cut(comp, rest))
                        continue;
                    prop_set.insert(rest);
                }
                if (prop_set.size() == 0)
                    continue;

                std::vector<int> prop_vec(prop_set.begin(), prop_set.end());
                bool exist = false;
                for (auto& loop: loops) {
                    if (loop == prop_vec) {
                        exist = true;
                        break;
                    }
                }
                if (!exist)
                    loops.push_back(prop_vec);
            }

            // sort the loops such that shorter loops are at the front
            std::sort(loops.begin(), loops.end(),
                [](const auto& x, const auto& y) {
                    return x.size() < y.size();
                }
            );
            
            if (loops.size() == 0)
                return {};
            else
                return loops[0];
        }
        case All: {
            std::vector<int> noncut;
            for (auto& idx: comp.propagator_indices) {
                if (is_cut(comp, idx))
                    continue;
                noncut.push_back(idx);
            }
            return noncut;
        }
    }
    std::cerr << "EtaScheme: unknown scheme to insert eta\n";
    return {};
}


static std::vector<int> eta_position(const std::vector<component>& comps, eta_scheme scheme) {
    if (scheme == Prescription || scheme == All) {
        // add eta for every component
        std::vector<int> final_candidates;
        for (auto& comp: comps) {
            std::vector<int> candidates;
            if (!comp.is_vacuum())
                candidates = eta_position(comp, scheme);
            else if (!comp.is_singlemass_vacuum())
                candidates = eta_position(comp, Branch);
            else
                candidates = {};

            for (auto& cand: candidates) {
                final_candidates.push_back(cand);
            }
        }
        return final_candidates;
    } else {
        // choose a component to add eta
        for (auto& comp: comps) {
            std::vector<int> candidates;
            if (!comp.is_vacuum())
                candidates = eta_position(comp, scheme);
            else if (!comp.is_singlemass_vacuum())
                candidates = eta_position(comp, Branch);
            else
                candidates = {};

            if (candidates.size() != 0)
                return candidates;
        }
        return {};
    }
}


std::shared_ptr<integralfamily> integralfamily::insert_eta(const std::vector<int>& top_sector, const std::vector<eta_scheme>& schemes) {
    generate_components(top_sector);
    std::vector<int> positions;
    for (auto scheme: schemes) {
        positions = eta_position(all_components, scheme);
        if (positions.size() != 0)
            break;
    }
    if (positions.size() == 0) { // no place to insert "eta"
        return nullptr;
    }

    std::shared_ptr<integralfamily> newfamily = std::make_shared<integralfamily>();
    newfamily->name        = name;
    newfamily->psymbols    = psymbols;
    newfamily->loops       = loops;
    newfamily->indeplegs   = indeplegs;
    newfamily->propagators = propagators;
    for (auto& idx: positions) {
        newfamily->propagators[idx] -= (*psymbols)["eta"];
    }
    newfamily->cut                = cut;
    newfamily->prescription       = prescription;
    newfamily->conservation_rules = conservation_rules;
    newfamily->numeric_rules      = numeric_rules;
    newfamily->sps_numeric_rules  = sps_numeric_rules;
    return newfamily;
}


