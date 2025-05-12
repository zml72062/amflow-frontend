#ifndef INTEGRAL_HPP
#define INTEGRAL_HPP


/**
 * integral.hpp - Define a Feynman integral.
 */


#include "family.hpp"


struct integral {
    // initializers
    integral(const integralfamily& _family)
        : pfamily(&_family), indices(_family.propagators.nops(), 0) { }

    template <class T>
    integral(const integralfamily& _family, T&& index_sequence)
        : pfamily(&_family) {
        for (int idx: index_sequence)
            indices.push_back(idx);
        
        if (indices.size() != _family.propagators.nops()) {
            std::cerr << "AMFlow: integral indices have wrong length\n";
            exit(1);
        }
    }

    integral(const integralfamily& _family, std::initializer_list<int> index_sequence)
        : pfamily(&_family) {
        for (int idx: index_sequence)
            indices.push_back(idx);

        if (indices.size() != _family.propagators.nops()) {
            std::cerr << "AMFlow: integral indices have wrong length\n";
            exit(1);
        }
    }

    std::string    to_string() const;
    
    std::vector<int>  sector() const;
    int               nprops() const;
    int               dot   () const;
    int               rank  () const;

    const integralfamily*   pfamily;
    std::vector<int>        indices;
};


typedef std::vector<integral> integral_list;

std::vector<int>           top_sector   (const integral_list& int_list);
std::vector<int>           top_position (const integral_list& int_list);
std::vector<integral_list> split_targets(const integral_list& int_list);


inline std::ostream& operator<<(std::ostream& out, const integral& integral) {
    out << integral.to_string();
    return out;
}


#endif // INTEGRAL_HPP
