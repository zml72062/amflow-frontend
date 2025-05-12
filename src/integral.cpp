#include "integral.hpp"
#include "utils.hpp"
#include <sstream>


std::string integral::to_string() const {
    std::ostringstream out;
    out << pfamily->name << "[" << indices << "]";
    return out.str();
}


std::vector<int> integral::sector() const {
    int n = indices.size();
    std::vector<int> sect(n, 0);
    for (int i = 0; i < n; i++)
        sect[i] = (int)(indices[i] > 0);
    return sect;
}


int integral::nprops() const {
    int count = 0;
    for (auto& ind: indices)
        if (ind > 0)
            count++;
    return count;
}


int integral::dot() const {
    int count = 0;
    for (auto& ind: indices)
        if (ind > 0)
            count += (ind - 1);
    return count;
}


int integral::rank() const {
    int count = 0;
    for (auto& ind: indices)
        if (ind <= 0)
            count -= ind;
    return count;
}


std::vector<int> top_sector(const integral_list& int_list) {
    if (int_list.size() == 0)
        return std::vector<int>();
    
    int n = int_list[0].indices.size();
    std::vector<int> top_sect(n, 0);
    for (auto& intgr: int_list) {
        auto sect = intgr.sector();
        for (int i = 0; i < n; i++)
            if (sect[i] > top_sect[i])
                top_sect[i] = 1;
    }

    return top_sect;
}


std::vector<int> top_position(const integral_list& int_list) {
    auto top_sect = top_sector(int_list);
    if (top_sect.size() == 0)
        return std::vector<int>();
    
    std::vector<int> top_pos;
    int n = top_sect.size();
    for (int i = 0; i < n; i++)
        if (top_sect[i])
            top_pos.push_back(i);
    
    return top_pos;
}


static bool is_subsector(const std::vector<int>& sect1, const std::vector<int>& sect2) {
    bool is_subsect = true;
    int n = sect1.size();
    for (int i = 0; i < n; i++) {
        if (sect1[i] > sect2[i]) {
            is_subsect = false;
            break;
        }
    }
    return is_subsect;
}


static bool is_true_subsector(const std::vector<int>& sect1, const std::vector<int>& sect2) {
    return is_subsector(sect1, sect2) && !is_subsector(sect2, sect1);
}


std::vector<integral_list> split_targets(const integral_list& int_list) {
    // find sector of each integral
    std::set<std::vector<int>> sectors;
    for (auto& intgr: int_list)
        sectors.insert(intgr.sector());

    // find all top sectors
    std::set<std::vector<int>> top_sectors;
    for (auto& sector: sectors) {
        bool is_a_subsect = false;
        for (auto& another: sectors) {
            if (is_true_subsector(sector, another)) {
                is_a_subsect = true;
                break;
            }
        }
        if (!is_a_subsect)
            top_sectors.insert(sector);
    }

    // classify integrals in the list
    int n = int_list.size(), n_visited = 0;
    auto top_iter = top_sectors.begin();
    std::vector<int> visited(n, 0);
    std::vector<integral_list> result;
    while (n_visited < n && top_iter != top_sectors.end()) {
        result.push_back(integral_list());
        auto current_sector = *top_iter;
        for (int i = 0; i < n; i++) {
            if (!visited[i] && is_subsector(int_list[i].sector(), current_sector)) {
                result.back().push_back(int_list[i]);
                visited[i] = 1;
                n_visited++;
            }
        }
        ++top_iter;
    }

    return result;
}

