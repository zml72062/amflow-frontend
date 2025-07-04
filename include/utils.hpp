#ifndef UTILS_HPP
#define UTILS_HPP


#include <chrono>
#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>


#define START_TIME(description) auto description##_begin = std::chrono::high_resolution_clock::now()
#define END_TIME(description) auto description##_end = std::chrono::high_resolution_clock::now(); \
    std::chrono::duration<double, std::milli> description##_time_ms = description##_end - description##_begin
#define PRINT_TIME(description) std::cerr << "Takes " << description##_time_ms.count() << " ms on " << #description << "()" << std::endl


inline bool has_non_null_key(const YAML::Node& node, const std::string& key) {
    auto node_map = node.as<std::map<std::string, YAML::Node>>();
    return node_map.find(key) != node_map.end() && node[key].Type() != YAML::NodeType::Null;
}


inline void add_symbol(const std::string& str, GiNaC::symtab& table) {
    if (table.find(str) == table.end())
        table[str] = GiNaC::symbol(str);
}


template <class T>
inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& vector) {
    int n = vector.size();
    for (int i = 0; i < n; i++)
        out << (i ? "," : "") << vector[i];
    return out;
}


inline GiNaC::lst matrix_to_lst(const GiNaC::matrix& mat) {
    GiNaC::lst result;
    int r = mat.rows(), c = mat.cols();
    for (int i = 0; i < r; i++) {
        GiNaC::lst arow;
        for (int j = 0; j < c; j++) {
            arow.append(mat(i, j));
        }
        result.append(arow);
    }
    return result;
}


/**
 * Add all factors of `poly` to `lst`.
 */
inline void set_factors(const GiNaC::ex& poly, GiNaC::lst& lst) {
    if (!GiNaC::is_exactly_a<GiNaC::mul>(poly)) {
        lst.append(poly);
        return;
    }

    for (auto iter = poly.begin(); iter != poly.end(); ++iter)
        set_factors(*iter, lst);
}


inline std::vector<GiNaC::lst> all_combinations(const GiNaC::lst& haystack, int r) {
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


inline GiNaC::matrix append_row(const GiNaC::matrix& _matrix, const GiNaC::matrix& _row) {
    int r = _matrix.rows(), c = _matrix.cols();
    GiNaC::matrix newmat(r + 1, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            newmat(i, j) = _matrix(i, j);
    for (int j = 0; j < c; j++)
        newmat(r, j) = _row(0, j);
    return newmat;
}


#endif // UTILS_HPP
