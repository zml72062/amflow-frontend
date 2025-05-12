#include "reduction.hpp"
#include "utils.hpp"


integral_reducer::integral_reducer(const YAML::Node& _node, const integralfamily& _family)
    : agent(_node, _family), pfamily(&_family) {
    mrank = 3;
    mdot = 0;
    try {
        if (has_non_null_key(_node, "ReduceOption")) {
            auto reduce_options = _node["ReduceOption"];
            if (has_non_null_key(reduce_options, "BlackBoxRank"))
                mrank = reduce_options["BlackBoxRank"].as<int>();
            if (has_non_null_key(reduce_options, "BlackBoxDot"))
                mdot = reduce_options["BlackBoxDot"].as<int>();
        }
    } catch (YAML::BadConversion& ex) {
        std::cerr << "KiraAgent: failed while reading reduction options\n";
        exit(1);
    }
    
    std::cerr << "KiraAgent: Set BlackBoxRank = " << mrank << "\n";
    std::cerr << "KiraAgent: Set BlackBoxDot = "  << mdot  << "\n";
}


static unsigned long top_sector_num(const integral_list& intlist) {
    auto topsect = top_sector(intlist);
    int n = topsect.size();
    unsigned long topsectnum = 0;
    for (int i = 0; i < n; i++)
        topsectnum |= (((unsigned long)topsect[i]) << i);
    return topsectnum;
}


static int max_rank(const integral_list& intlist) {
    int rank = 0;
    for (auto& intgr: intlist) {
        int irank = intgr.rank();
        if (irank > rank)
            rank = irank;
    }
    return rank;
}


static int max_dot(const integral_list& intlist) {
    int dot = 0;
    for (auto& intgr: intlist) {
        int idot = intgr.dot();
        if (idot > dot)
            dot = idot;
    }
    return dot;
}


std::pair<integral_list, GiNaC::matrix>
integral_reducer::reduce(const integral_list& target, const integral_list& preferred, const char* workdir) {
    integral_list all_integrals;
    for (auto& tint: target)
        all_integrals.push_back(tint);
    for (auto& pint: preferred)
        all_integrals.push_back(pint);
    if (all_integrals.size() == 0) {
        std::cerr << "KiraAgent: either target integral list or preferred master integral list should be non-empty\n";
        exit(1);
    }
    
    auto topsectnum = top_sector_num(all_integrals);
    int rank = std::max(mrank, max_rank(all_integrals));
    int dot  = std::max(mdot,  max_dot (all_integrals));

    agent.reduce(workdir, target, preferred, topsectnum, rank, dot, true);
    return agent.reduce(workdir, target, preferred, topsectnum, rank, dot, false);
}


std::pair<integral_list, GiNaC::matrix>
integral_reducer::diffeq(const integral_list& preferred, const char* workdir) {
    if (preferred.size() == 0) {
        std::cerr << "KiraAgent: preferred master integral list should be non-empty\n";
        exit(1);
    }

    int n = preferred[0].indices.size();
    auto topsectnum = top_sector_num(preferred);
    int rank = std::max(mrank, max_rank(preferred));
    int dot  = std::max(mdot,  max_dot(preferred) + 1);
    auto masters = agent.reduce(workdir, {}, preferred, topsectnum, rank, dot, true).first;
    int nmasters = masters.size();

    GiNaC::lst denom_deriv;
    auto eta = GiNaC::ex_to<GiNaC::symbol>((*pfamily->psymbols)["eta"]);
    for (auto& denom: pfamily->propagators) {
        denom_deriv.append(denom.diff(eta).expand());
    }

    std::vector<std::vector<std::pair<integral, GiNaC::ex>>> to_reduce;
    std::map<std::vector<int>, GiNaC::matrix>                reduce_result;
    for (auto& mi: masters) {
        to_reduce.push_back(std::vector<std::pair<integral, GiNaC::ex>>());
        for (int i = 0; i < n; i++) {
            if (!(bool)(denom_deriv[i] == 0) && mi.indices[i] != 0) {
                GiNaC::ex prefactor = -denom_deriv[i] * mi.indices[i];
                std::vector<int> newint(n, 0);
                for (int j = 0; j < n; j++)
                    newint[j] = mi.indices[j] + (j == i ? 1 : 0);
                integral newintegral(*pfamily, newint);
                to_reduce.back().push_back(std::make_pair(newintegral, prefactor));
                reduce_result[newint] = GiNaC::matrix(0, 0);
            }
        }
    }

    integral_list all_to_reduce;
    for (auto& reduce_key: reduce_result)
        all_to_reduce.push_back(integral(*pfamily, reduce_key.first));
    auto result = agent.reduce(workdir, all_to_reduce, preferred, topsectnum, rank, dot, false).second;
    int nreduce = all_to_reduce.size();
    auto reduce_iter = reduce_result.begin();
    for (int i = 0; i < nreduce; i++) {
        reduce_iter->second = GiNaC::ex_to<GiNaC::matrix>(GiNaC::sub_matrix(result, i, 1, 0, nmasters));
        ++reduce_iter;
    }
    
    GiNaC::matrix coeff_matrix(nmasters, nmasters);
    for (int i = 0; i < nmasters; i++) {
        auto pre_reduce_list = to_reduce[i];
        for (auto& pr: pre_reduce_list) {
            auto line = reduce_result[pr.first.indices].mul_scalar(pr.second);
            for (int j = 0; j < nmasters; j++)
                coeff_matrix(i, j) += line(0, j);
        }
    }

    return std::make_pair(masters, coeff_matrix);
}

