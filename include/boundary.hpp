#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP


#include "integral.hpp"
#include <memory>


/**
 * A boundary region refers to a tuple of the following:
 *      (a) a way to choose L independent loop momenta
 *      (b) a way to assign relative scales to those loop momenta,
 *          each of which being "large" (\~sqrt(eta)) or "small" (\~1).
 */
struct boundary_region {
    boundary_region(integralfamily& family): pfamily(&family) { }

    /**
     * Find all boundary regions of a given family depending on "eta".
     */
    static std::vector<boundary_region> all_boundary_regions(integralfamily& family, const std::vector<int>& top_sector);

    /**
     * Determine whether the current region is trivial.
     * 
     * Note: The arguments `workdir` and `config` are for the need of 
     * calling Kira.
     */
    bool is_trivial_region(const std::vector<int>& top_sector,
                           const char* workdir, const YAML::Node& config);

    /**
     * Determine whether the current region is trivial.
     */
    bool is_trivial_region(const std::vector<int>& top_sector);

    /**
     * Compute the overall factor in the leading asymptotic expansion
     * of integral `intgr` in the vicinity of "eta" around infinity.
     */
    GiNaC::ex factor(const integral& intgr);

    /**
     * Compute the overall exponent in the leading asymptotic expansion
     * of integral `intgr` in the vicinity of "eta" around infinity.
     */
    GiNaC::ex exponent(const integral& intgr);
    
    /**
     * Generate boundary integrals corresponding to `intgr` up to expansion
     * order `border` in the vicinity of "eta" around infinity.
     * 
     * @returns a pair consisting of
     *      (a) a list of (family, integral_list) pairs, each entry 
     *          consisting of
     *              - shared pointer to a subfamily
     *              - a list of integrals belonging to this subfamily
     *      (b) a vector of mappings, whose i-th entry (i starting from
     *          0) records information regarding the expansion term
     *          proportional to [eta^i]. The mapping maps a (int, int)
     *          index to a linear combination coefficient. 
     * 
     *          The (int, int) index can be used to index into the
     *          (family, integral_list) table, to refer to a specific
     *          boundary integral.
     */
    std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
              std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>
    subintegrals(const std::vector<int>& top_sector, const integral& intgr, int border);

    /**
     * Generate boundary integrals corresponding to a list of master integrals
     * `intgrs`, the expansion order of each integral determined by the
     * corresponding entry in `border`.
     * 
     * @returns a pair consisting of
     *      (a) a list of (family, integral_list) pairs, each entry 
     *          consisting of
     *              - shared pointer to a subfamily
     *              - a list of integrals belonging to this subfamily
     *      (b) a vector of vector of mappings, each entry storing
     *          information for one master integral. The entry is
     *          a vector of mappings, whose i-th entry (i starting from
     *          0) records information regarding the expansion term
     *          proportional to [eta^i]. The mapping maps a (int, int)
     *          index to a linear combination coefficient. 
     * 
     *          The (int, int) index can be used to index into the
     *          (family, integral_list) table, to refer to a specific
     *          boundary integral.
     */
    std::pair<std::vector<std::pair<std::shared_ptr<integralfamily>, integral_list>>,
              std::vector<std::vector<std::map<std::pair<int, int>, GiNaC::ex>>>>
    subintegrals(const std::vector<int>& top_sector, const integral_list& intgrs, const std::vector<int>& border);

    std::string to_string() const;
    
    integralfamily*  pfamily;
    GiNaC::matrix    transform;
    GiNaC::lst       new_propagators;
    GiNaC::lst       new_leading_propagators[3];
    GiNaC::lst       leading_factors;
    std::vector<int> is_large;

private:
    void      compute_new_propagators();
    void      compute_new_leading_propagators();
    GiNaC::ex compute_expansion_in_new_leading_propagators(const GiNaC::ex& prop);

    GiNaC::lst       memo_sps;
    GiNaC::matrix    coeff_sps_to_new_props;
    GiNaC::matrix    bias_sps_to_new_props;
    GiNaC::lst       new_propagator_symbols;

    /**
     * Generate a preliminary integral family for the boundary integrals 
     * in this region.
     */
    integralfamily preliminary_subfamily(const std::vector<int>& top_sector);
};


#endif // BOUNDARY_HPP
