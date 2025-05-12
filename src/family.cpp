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

    return family;
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

