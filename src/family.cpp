#include "family.hpp"
#include "utils.hpp"
#include <sstream>


template <class T>
static void add_symbols(T&& seq_strings, GiNaC::symtab& table) {
    for (const std::string& str: seq_strings)
        add_symbol(str, table);
}

template <class T>
static void add_symbols(T&& seq_strings, const GiNaC::symtab& table, GiNaC::lst& list) {
    for (const std::string& str: seq_strings) {
        list.append(table.at(str));
    }
}

static GiNaC::matrix adjugate(const GiNaC::matrix& M) {
    auto n = M.rows(), m = M.cols();
    if (n != m)
        throw std::runtime_error("adjugate(): non-square matrix");
    
    GiNaC::matrix adj(n, n);
    if (n == 1) {
        adj(0, 0) = 1;
        return adj;
    }

    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = 0; j < n; j++) {
            adj(i, j) = GiNaC::pow(-1, i + j) * 
                GiNaC::ex_to<GiNaC::matrix>(
                    GiNaC::reduced_matrix(M, j, i)
                ).determinant();
        }
    }

    return adj;
}


integralfamily integralfamily::from_yaml(const YAML::Node& _node, GiNaC::symtab* _psymbols) {
    integralfamily family;
    family.psymbols = _psymbols;
    
    try {
        // read family name
        family.name = _node["AMFlowInfo"]["Family"].as<std::string>();

        // read loop momenta
        if (has_non_null_key(_node["AMFlowInfo"], "Loop")) {
            auto loops_ = _node["AMFlowInfo"]["Loop"].as<std::vector<std::string>>();
            add_symbols(loops_, *family.psymbols);
            add_symbols(loops_, *family.psymbols, family.loops);
        }
        
        // read external legs
        std::set<std::string> set_legs;
        if (has_non_null_key(_node["AMFlowInfo"], "Leg")) {
            auto legs_ = _node["AMFlowInfo"]["Leg"].as<std::vector<std::string>>();
            add_symbols(legs_, *family.psymbols);
            for (const std::string& leg: legs_) {
                set_legs.insert(leg);
            }
        }

        // read momentum conservation
        if (has_non_null_key(_node["AMFlowInfo"], "Conservation")) {
            auto conservation_ = _node["AMFlowInfo"]["Conservation"].as<std::vector<std::string>>();
            if (conservation_.size() != 2) {
                std::cerr << "IntegralFamily: wrong conservation rule format\n";
                exit(1);
            }
            set_legs.erase(conservation_[0]);

            GiNaC::parser pleg(*family.psymbols, true);
            family.conservation_rules.append(pleg(conservation_[0]) == pleg(conservation_[1]));
        }
        add_symbols(set_legs, *family.psymbols, family.indeplegs);

        // read numeric rules
        if (has_non_null_key(_node["AMFlowInfo"], "Numeric")) {
            auto numerics_ = _node["AMFlowInfo"]["Numeric"].as<std::vector<YAML::Node>>();
            GiNaC::parser pvalue;
            for (auto& numeric: numerics_) {
                auto numeric_pair = numeric.as<std::vector<std::string>>();
                if (numeric_pair.size() != 2) {
                    std::cerr << "IntegralFamily: wrong numeric format\n";
                    exit(1);
                }
                add_symbol(numeric_pair[0], *family.psymbols);

                family.numeric_rules.append((*family.psymbols)[numeric_pair[0]]
                                             == pvalue(numeric_pair[1]));
            }
        }

        // read propagators
        if (has_non_null_key(_node["AMFlowInfo"], "Propagator")) {
            auto propagators_ = _node["AMFlowInfo"]["Propagator"].as<std::vector<std::string>>();
            GiNaC::parser pprop(*family.psymbols, true);
            for (auto& propagator: propagators_)
                family.propagators.append(pprop(propagator)
                    .subs(family.conservation_rules, GiNaC::subs_options::algebraic)
                    .subs(family.numeric_rules, GiNaC::subs_options::algebraic)
                );
        }

        // read cuts
        if (has_non_null_key(_node["AMFlowInfo"], "Cut")) {
            family.cut = _node["AMFlowInfo"]["Cut"].as<std::vector<int>>();
            if (family.cut.size() != family.propagators.nops()) {
                std::cerr << "IntegralFamily: wrong cut input\n";
                exit(1);
            }
        }

        // read prescriptions
        if (has_non_null_key(_node["AMFlowInfo"], "Prescription")) {
            family.prescription = _node["AMFlowInfo"]["Prescription"].as<std::vector<int>>();
            if (family.prescription.size() != family.loops.nops()) {
                std::cerr << "IntegralFamily: wrong prescription input\n";
                exit(1);
            }
        }

        // read and compute scalar product numerics
        if (has_non_null_key(_node["AMFlowInfo"], "Replacement")) {
            auto replacements_ = _node["AMFlowInfo"]["Replacement"].as<std::vector<YAML::Node>>();
            GiNaC::parser preplace(*family.psymbols, true);
            int num_indeplegs = family.indeplegs.nops(),
                num_repls     = replacements_.size(),
                num_spss      = num_indeplegs * (num_indeplegs + 1) / 2;
            GiNaC::matrix coeff    (num_repls, num_spss),
                          augmented(num_repls, num_spss + 1),
                          bias     (num_repls, 1);
            GiNaC::lst    sps;  // scalar products
            
            // loop over all replacements
            for (int r = 0; r < num_repls; r++) {
                auto replacement_pair = replacements_[r].as<std::vector<std::string>>();
                if (replacement_pair.size() != 2) {
                    std::cerr << "IntegralFamily: wrong replacement format\n";
                    exit(1);
                }
                GiNaC::ex lhs = preplace(replacement_pair[0]).subs(family.conservation_rules, GiNaC::subs_options::algebraic)
                                                             .subs(family.numeric_rules, GiNaC::subs_options::algebraic),
                          rhs = preplace(replacement_pair[1]).subs(family.numeric_rules, GiNaC::subs_options::algebraic);
                
                // fill LHS into coefficient matrix and augmented matrix
                int n = 0;
                for (int i = 0; i < num_indeplegs; i++) {
                    GiNaC::ex deriv = lhs.diff(GiNaC::ex_to<GiNaC::symbol>(family.indeplegs[i]));
                    for (int j = i; j < num_indeplegs; j++) {
                        coeff(r, n)     = deriv.diff(GiNaC::ex_to<GiNaC::symbol>(family.indeplegs[j])) / (i == j ? 2 : 1);
                        augmented(r, n) = coeff(r, n);
                        n++;
                        if (r == 0)
                            sps.append(family.indeplegs[i] * family.indeplegs[j]);
                    }
                }

                // fill RHS into augmented matrix
                augmented(r, num_spss) = rhs;
                bias(r, 0)             = rhs;
            }

            if ((int)coeff.rank() < num_spss) {
                std::cerr << "IntegralFamily: insufficient replacement rules for all independent external scalar products\n";
                exit(1);
            }
            if ((int)augmented.rank() > num_spss) {
                std::cerr << "IntegralFamily: inconsistent replacement rules\n";
                exit(1);
            }
            
            // solve scalar products
            GiNaC::matrix temp(num_repls, 1);
            for (int r = 0; r < num_repls; r++)
                temp(r, 0) = GiNaC::symbol("var" + std::to_string(r));
            GiNaC::matrix result = coeff.solve(temp, bias);
            for (int n = 0; n < num_spss; n++)
                family.sps_numeric_rules.append(sps[n] == result(n, 0));
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "IntegralFamily: failed to load integral family information\n";
        exit(1);
    }

    if (!family.is_complete_propagator()) {
        std::cerr << "IntegralFamily: propagators not complete\n";
        std::cerr << "Try using the following set of propagators:\n";
        sps_struct suggestion = family.sps_in_complete_propagator();
        std::cerr << suggestion.prop << "\n";
        exit(1);
    }

    return family;
}


void integralfamily::compute_symanzik(const std::vector<int>& top_sector) {
    if (top_sector.size() != propagators.nops()) {
        std::cerr << "IntegralFamily: wrong top sector input\n";
        exit(1);
    }

    feynman_params = GiNaC::lst();
    int n_props = top_sector.size();
    for (int i = 0; i < n_props; i++)
        feynman_params.append(GiNaC::symbol("x" + std::to_string(i)));
    
    GiNaC::ex  combined_denominator = 0;
    for (int i = 0; i < n_props; i++)
        if (top_sector[i])
            combined_denominator += (feynman_params[i] * propagators[i]);
    
    GiNaC::lst all_zeros;
    for (auto& l: loops)
        all_zeros.append(l == 0);

    GiNaC::ex J = -combined_denominator.subs(all_zeros, GiNaC::subs_options::algebraic);
    int n_loop = loops.nops();
    GiNaC::matrix V(n_loop, 1), M(n_loop, n_loop);
    for (int i = 0; i < n_loop; i++) {
        GiNaC::ex deriv = combined_denominator.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i]));
        V(i, 0) = -deriv.subs(all_zeros, GiNaC::subs_options::algebraic) / 2;
        for (int j = 0; j < n_loop; j++)
            M(i, j) = deriv.diff(GiNaC::ex_to<GiNaC::symbol>(loops[j])) / 2;
    }

    U = M.determinant().expand();
    F = (V.transpose().mul(adjugate(M)).mul(V)(0, 0) + J * U
         ).expand().subs(conservation_rules, GiNaC::subs_options::algebraic)
                   .subs(numeric_rules, GiNaC::subs_options::algebraic)
                   .subs(sps_numeric_rules, GiNaC::subs_options::algebraic);
    
    GiNaC::lst prop_masses = masses();
    F0 = 0;
    for (int i = 0; i < n_props; i++)
        if (top_sector[i])
            F0 += (feynman_params[i] * prop_masses[i]);

    top_sector_mem = top_sector;
}


void integralfamily::generate_components(const std::vector<int>& top_sector) {
    if (all_components.size() != 0 && top_sector_mem == top_sector)
        // cached
        return;

    compute_symanzik(top_sector);
    GiNaC::ex  U_factored = GiNaC::factor(U);
    GiNaC::lst U_factors;
    set_factors(U_factored, U_factors);

    all_components = std::vector<component>();
    int nprop = feynman_params.nops();
    for (auto& factor: U_factors) {
        component comp;
        comp.propagator_indices = std::vector<int>();
        for (int i = 0; i < nprop; i++)
            if (factor.has(feynman_params[i]))
                comp.propagator_indices.push_back(i);
        comp.U       = factor;
        comp.pfamily = this;
        all_components.push_back(comp);
    }
}


std::string integralfamily::to_string() const {
    std::ostringstream out;
    out << "Name:                   " << name        << "\n"
        << "Loop momenta:           " << loops       << "\n"
        << "Independent legs:       " << indeplegs   << "\n"
        << "Propagators:            " << propagators << "\n"
        << "\n";
    
    if (cut.size() > 0)
        out << "Cut:                    {" << cut          << "}\n";
    if (prescription.size() > 0)
        out << "Prescription:           {" << prescription << "}\n";

    out << "\n";
    out << "Numerics:               " << numeric_rules     << "\n"
        << "Scalar products:        " << sps_numeric_rules << "\n";
    
    return out.str();
}


bool integralfamily::is_complete_propagator() const {
    int n_loops = loops.nops(), n_legs = indeplegs.nops();
    int S = n_loops * (n_loops + 1) / 2 + n_loops * n_legs;
    if ((int)propagators.nops() != S)
        return false;

    GiNaC::matrix coeff(S, S);
    for (int p = 0; p < S; p++) {
        int q = 0;
        GiNaC::ex prop = propagators[p].expand();
        for (int i = 0; i < n_loops; i++) {
            GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand();
            coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand() / 2;
            for (int j = i + 1; j < n_loops; j++)
                coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(loops[j])).expand();
        }
        for (int i = 0; i < n_loops; i++) {
            GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand();
            for (int j = 0; j < n_legs; j++)
                coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(indeplegs[j])).expand();
        }
    }
    
    if ((bool)(coeff.determinant().expand() == 0))
        return false;
    
    return true;
}


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


#define SUBMATRIX(m, r, nr, c, nc) GiNaC::ex_to<GiNaC::matrix>(GiNaC::sub_matrix((m), (r), (nr), (c), (nc)))


sps_struct integralfamily::sps_in_complete_propagator() const {
    int n_loops = loops.nops(), n_legs = indeplegs.nops();
    int S = n_loops * (n_loops + 1) / 2 + n_loops * n_legs;
    GiNaC::lst       full_props;
    std::vector<int> full_cuts = cut;
    
    for (auto& prop: propagators) {
        full_props.append(prop.expand());
    }
    for (int i = 0; i < n_loops; i++) {
        full_props.append(GiNaC::pow(loops[i], 2));
        for (int j = i + 1; j < n_loops; j++)
            full_props.append(GiNaC::pow(loops[i] + loops[j], 2).expand());
    }
    for (int i = 0; i < n_loops; i++)
        for (int j = 0; j < n_legs; j++)
            full_props.append(GiNaC::pow(loops[i] + indeplegs[j], 2).expand());
    
    if (cut.size() != 0) {
        for (int i = 0; i < n_loops; i++)
            for (int j = i; j < n_loops; j++)
                full_cuts.push_back(0);
        for (int i = 0; i < n_loops; i++)
            for (int j = 0; j < n_legs; j++)
                full_cuts.push_back(0);
    }

    GiNaC::matrix coeff(full_props.nops(), S);
    for (int p = 0; p < (int)full_props.nops(); p++) {
        int q = 0;
        GiNaC::ex prop = full_props[p];
        for (int i = 0; i < n_loops; i++) {
            GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand();
            coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand() / 2;
            for (int j = i + 1; j < n_loops; j++)
                coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(loops[j])).expand();
        }
        for (int i = 0; i < n_loops; i++) {
            GiNaC::ex i_deriv = prop.diff(GiNaC::ex_to<GiNaC::symbol>(loops[i])).expand();
            for (int j = 0; j < n_legs; j++)
                coeff(p, q++) = i_deriv.diff(GiNaC::ex_to<GiNaC::symbol>(indeplegs[j])).expand();
        }
    }

    std::vector<int> indep_rows;
    GiNaC::matrix selected_indep_rows(0, S);
    int next_row = 0, accum_rank = 0;
    while (accum_rank != S) {
        auto new_selected_indep_rows = append_row(selected_indep_rows, SUBMATRIX(coeff, next_row++, 1, 0, S));
        if ((int)new_selected_indep_rows.rank() > accum_rank) {
            indep_rows.push_back(next_row - 1);
            selected_indep_rows = new_selected_indep_rows;
            accum_rank++;
        }
    }
    GiNaC::lst       complete_props;
    std::vector<int> newcut;
    for (auto idx: indep_rows)
        complete_props.append(full_props[idx]);
    if (cut.size() != 0)
        for (auto idx: indep_rows)
            newcut.push_back(full_cuts[idx]);
    
    GiNaC::lst sps;
    for (int i = 0; i < n_loops; i++)
        for (int j = i; j < n_loops; j++)
            sps.append(loops[i] * loops[j]);
    for (int i = 0; i < n_loops; i++)
        for (int j = 0; j < n_legs; j++)
            sps.append(loops[i] * indeplegs[j]);

    GiNaC::matrix coefficients = selected_indep_rows.inverse();
    GiNaC::matrix prop_matrix(S, 1), sps_matrix(S, 1);
    for (int i = 0; i < S; i++) {
        prop_matrix(i, 0) = complete_props[i];
        sps_matrix(i, 0)  = sps[i];
    }
    
    GiNaC::matrix bias = GiNaC::ex_to<GiNaC::matrix>(
        ((GiNaC::ex)(sps_matrix.sub(coefficients.mul(prop_matrix))))
            .subs(sps_numeric_rules, GiNaC::subs_options::algebraic)
            .subs(numeric_rules, GiNaC::subs_options::algebraic)
    );

    sps_struct output;
    output.prop  = complete_props;
    output.sps   = sps;
    output.coeff = coefficients;
    output.bias  = bias;
    output.cut   = newcut;
    return output;
}


GiNaC::lst integralfamily::momenta() const {
    GiNaC::lst momenta_list;
    for (auto& prop: propagators) {
        bool success = false;
        GiNaC::ex principal_loop_momentum;
        for (auto& l: loops) {
            auto top_coeff = prop.expand().collect(l).coeff(l, 2);
            if (!(bool)(top_coeff == 0)) {
                principal_loop_momentum = l;
                success = true;
                break;
            }
        }
        if (!success) {
            std::cerr << "IntegralFamily: failed to determine momentum flowing on propagator " << prop << "\n";
            return {};
        }
        auto collected_prop = prop.expand().collect(principal_loop_momentum);
        auto a = collected_prop.coeff(principal_loop_momentum, 2),
             b = collected_prop.coeff(principal_loop_momentum, 1);
        momenta_list.append(principal_loop_momentum + b / (2 * a));
    }
    return momenta_list;
}


GiNaC::lst integralfamily::masses() const {
    GiNaC::lst mass_list;
    for (auto& prop: propagators) {
        bool success = false;
        GiNaC::ex principal_loop_momentum;
        for (auto& l: loops) {
            auto top_coeff = prop.expand().collect(l).coeff(l, 2);
            if (!(bool)(top_coeff == 0)) {
                principal_loop_momentum = l;
                success = true;
                break;
            }
        }
        if (!success) {
            std::cerr << "IntegralFamily: failed to determine momentum flowing on propagator " << prop << "\n";
            return {};
        }
        auto collected_prop = prop.expand().collect(principal_loop_momentum);
        auto a = collected_prop.coeff(principal_loop_momentum, 2),
             b = collected_prop.coeff(principal_loop_momentum, 1),
             c = collected_prop.coeff(principal_loop_momentum, 0);
        mass_list.append(((b * b) / (4 * a) - c).expand()
            .subs(conservation_rules, GiNaC::subs_options::algebraic)
            .subs(numeric_rules, GiNaC::subs_options::algebraic)
            .subs(sps_numeric_rules, GiNaC::subs_options::algebraic));
    }
    return mass_list;
}


GiNaC::lst integralfamily::loop_momenta() const {
    GiNaC::ex momenta_ = momenta();
    GiNaC::lst simplify_rules;
    for (auto& leg: indeplegs)
        simplify_rules.append(leg == 0);
    return GiNaC::ex_to<GiNaC::lst>(momenta_.subs(simplify_rules, GiNaC::subs_options::algebraic));
}


std::vector<int> integralfamily::propagator_prescription() const {
    std::vector<int> prop_pres;
    int n_pres = prescription.size();
    for (auto& prop: propagators) {
        if (n_pres == 0) {
            prop_pres.push_back(1);
            continue;
        }
        
        // now "n_pres" should equal number of loops
        std::vector<int> pre_pos, pre_neg, pre_zero;
        for (int i = 0; i < n_pres; i++) {
            if (prop.has(loops[i])) {
                if (prescription[i] == 1) {
                    pre_pos.push_back(i);
                } else if (prescription[i] == -1) {
                    pre_neg.push_back(i);
                } else if (prescription[i] == 0) {
                    pre_zero.push_back(i);
                }
            }
        }

        bool no_pos = (pre_pos.size() == 0), no_neg = (pre_neg.size() == 0);
        if (no_pos && no_neg) {
            prop_pres.push_back(0);
        } else if (no_pos) {
            prop_pres.push_back(-1);
        } else if (no_neg) {
            prop_pres.push_back(1);
        } else {
            prop_pres.push_back(PRES_FAILED);
        }
    }
    return prop_pres;
}


bool integralfamily::is_zero(const std::vector<int>& top_sector) {
    // a zero top sector must be trivial
    bool all_zero_top = true;
    for (auto& ind: top_sector) {
        if (ind) {
            all_zero_top = false;
            break;
        }
    }
    if (all_zero_top)
        return true;

    compute_symanzik(top_sector);
    auto G = (U + F).expand();
    
    auto eq_generate = G;
    GiNaC::lst vars;
    int nvars = 0, nprop = propagators.nops();
    for (int i = 0; i < nprop; i++) {
        if (top_sector[i]) {
            vars.append(GiNaC::symbol("var" + std::to_string(nvars++)));
            eq_generate -= (G.diff(GiNaC::ex_to<GiNaC::symbol>(feynman_params[i])) * feynman_params[i] * vars[nvars - 1]);
        }
    }

    GiNaC::lst eqns = {eq_generate};
    for (int i = 0; i < nprop; i++) {
        if (top_sector[i]) {
            GiNaC::lst neweqns;
            for (auto& eq: eqns) {
                auto expanded = eq.expand().collect(feynman_params[i]);
                int deg = expanded.degree(feynman_params[i]);
                for (int d = 0; d <= deg; d++)
                    neweqns.append(expanded.coeff(feynman_params[i], d));
            }
            eqns = neweqns;
        }
    }
    
    GiNaC::lst final_eqns;
    for (auto& eq: eqns)
        final_eqns.append(eq == 0);
    auto solution = GiNaC::lsolve(final_eqns, vars);
    if (solution.nops() == 0)
        return false;
    return true;
}

