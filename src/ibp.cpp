#include "ibp.hpp"
#include "utils.hpp"
#include <fstream>


std::string ibphelper::get_master_id(const std::string& input) {
    std::size_t start = input.find('[') + 1,
                end   = input.find(']');
    return input.substr(start, end - start);
}


integral ibphelper::to_integral(const std::string& input) {
    std::string id = get_master_id(input);
    std::vector<int> nums;
    const char* ptr = id.c_str();
    char* endptr;
    while (true) {
        nums.push_back(std::strtol(ptr, &endptr, 10));
        if (*endptr == '\0')
            return integral(*pfamily, nums);
        ptr = endptr + 1;
    }
}


integral_list ibphelper::get_master_integral_list() {
    integral_list masters;
    for (auto& mi: master_table)
        masters.push_back(to_integral(mi.first));

    return masters;
}


void ibphelper::load_masters(const char* masters_path) {
    std::ifstream masters_file(masters_path);
    std::string master;
    while (true) {
        masters_file >> master;
        if (masters_file.eof())
            break;
        if (master.find(family_name()) != std::string::npos) {
            add_symbol(family_name() + "[" + get_master_id(master) + "]",
                       master_table);
        }
    }
    masters_file.close();
}


void ibphelper::load_ibps(const char* ibp_path) {
    std::ifstream ibp_file(ibp_path);
    std::string ibp, current_key;
    std::size_t _asterisk, counter = 0;
    GiNaC::parser coeff_parser(*psymbols);

    while (true) {
        ibp_file >> ibp;
        if (ibp_file.eof())
            break;
        if (ibp.find(family_name()) != std::string::npos) {
            if ((_asterisk = ibp.find('*')) == std::string::npos) { // an IBP head
                counter++;
                current_key = family_name() + "[" + get_master_id(ibp) + "]";
                ibp_table[current_key] = 0;
                std::cerr << "Processing the " << counter << "-th IBP relation" << "\r";
            } else { // an IBP body
                GiNaC::ex current_integral = master_table.at(
                    family_name() + "[" + get_master_id(ibp) + "]"
                );
                GiNaC::ex coefficient      = coeff_parser(ibp.substr(_asterisk + 1));
                ibp_table[current_key]    += (current_integral * coefficient);
            }
        }
    }

    // add IBP entries for master integrals themselves
    for (auto& master: master_table)
        ibp_table[master.first] = master.second;

    std::cerr << std::endl << "Done!" << std::endl;
    ibp_file.close();
}


void ibphelper::print_masters(std::ostream& out) {
    for (auto& master: master_table)
        out << master.first << "\n";
}


void ibphelper::print_ibps(std::ostream& out) {
    for (auto& ibp: ibp_table) {
        out << ibp.first << " = " << ibp.second << ";\n";
    }
}

