#ifndef SINGLESETUP_HPP
#define SINGLESETUP_HPP

/**
 * singlesetup.hpp - Set up a single system.
 */


#include "reduction.hpp"
#include "boundary.hpp"


/**
 * An integral system means
 * 
 *      - an integral family
 *      - a list of target integrals
 * 
 * An integral system comes along with
 * 
 *      - a working directory
 *      - a configure file
 *      - a global symbol table
 * 
 * In AMFlow, integral systems are fundamental working units.
 */
class integral_system {
public:
    integral_system(integralfamily&                       family, 
                    const std::vector<std::vector<int>>&  integrals, 
                    const std::vector<std::vector<int>>&  preferred, 
                    const YAML::Node&                     config, 
                    const char*                           workdir);

    void reduce_targets       ();
    void build_diffeq         ();
    void determine_boundaries ();
    void determine_border     ();
    void determine_direction  ();
    void setup_subfamilies    ();

private:
    integralfamily*                       pfamily;
    GiNaC::symtab*                        psymbols;
    const YAML::Node*                     pconfig;
    integral_list                         targets;
    integral_list                         preferred_masters;
    integral_list                         masters;
    integral_reducer                      reducer;
    std::string                           workdirstr;
    std::string                           tempdirstr;
    GiNaC::ex                             d0;
    std::vector<int>                      master_sector;
    std::vector<boundary_region>          regions;
    std::vector<std::vector<int>>         borders;
};


#endif // SINGLESETUP_HPP
