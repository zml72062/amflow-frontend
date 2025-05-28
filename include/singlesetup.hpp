#ifndef SINGLESETUP_HPP
#define SINGLESETUP_HPP

/**
 * singlesetup.hpp - Set up a single system.
 */


#include "reduction.hpp"
#include "boundary.hpp"
#include "etascheme.hpp"
#include <memory>


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

    /**
     * Find a set of master integrals, and reduce target integrals to 
     * them.
     */
    void reduce_targets       ();
    
    /**
     * Assuming that a set of master integrals have been found, build
     * differential equation with respect to eta.
     */
    void build_diffeq         ();

    /**
     * Assuming that a set of master integrals have been found, find
     * non-trivial boundary regions around eta=infinity.
     */
    void determine_boundaries ();

    /**
     * Assuming that differential equation of master integrals has been
     * built, and boundary regions have been found, find boundary order
     * for each boundary region.
     */
    void determine_border     ();

    /**
     * Assuming that a set of master integrals have been found, find a
     * direction to integrate eta from infinity to zero. The direction
     * is either "NegIm" or "Im".
     */
    void determine_direction  ();

    /**
     * Assuming that boundary regions have been determined, setup all
     * subfamilies arising from the boundary regions, and reduce all
     * boundary integrals.
     */
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>> 
    setup_subfamilies           ();

    /**
     * Assuming that a set of master integrals have been found, insert
     * eta to the current integral family to construct a new family of
     * integrals, and find master integrals of the new family.
     */
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
    setup_eta                   ();

    /**
     * Assuming that the current integral family is a single-mass vacuum
     * type ending family, and that master integrals have been found for
     * the current integral family, setup all subfamilies corresponding
     * to connected components of the current family, and reduce master
     * integrals of the current family to master integrals of the subfa-
     * milies.
     */
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
    setup_singlemass_vacuum     ();

    /**
     * Assuming that the current integral family is a purely phase-space
     * type ending family, and that master integrals have been found for
     * the current integral family, setup a new family of loop integrals
     * whose imaginary part corresponds to the current family.
     */
    std::vector<std::pair<std::shared_ptr<integralfamily>, integral_system>>
    setup_purely_phasespace     ();

    friend class amflow;
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
    std::vector<eta_scheme>               eta_mode;
    std::vector<int>                      master_sector;
    std::vector<boundary_region>          regions;
    std::vector<std::vector<int>>         borders;
};


#endif // SINGLESETUP_HPP
