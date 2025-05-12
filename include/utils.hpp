#ifndef UTILS_HPP
#define UTILS_HPP


#include <chrono>
#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>


#define START_TIME(description) auto description##_begin = std::chrono::high_resolution_clock::now()
#define END_TIME(description) auto description##_end = std::chrono::high_resolution_clock::now(); \
    std::chrono::duration<double, std::milli> description##_time_ms = description##_end - description##_begin
#define PRINT_TIME(description) std::cout << "Takes " << description##_time_ms.count() << " ms on " << #description << "()" << std::endl


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


#endif // UTILS_HPP
