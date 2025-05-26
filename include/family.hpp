#ifndef FAMILY_HPP
#define FAMILY_HPP


/**
 * family.hpp - Define a Feynman integral family.
 */


#include <yaml-cpp/yaml.h>
#include <ginac/ginac.h>
#include <memory>
#include "etascheme.hpp"


struct sps_struct {
    GiNaC::lst       prop;
    GiNaC::lst       sps;
    GiNaC::matrix    coeff;
    GiNaC::matrix    bias;
    std::vector<int> cut;
};


struct integralfamily;

struct component {
    std::vector<int> propagator_indices;
    GiNaC::ex        U;
    integralfamily*  pfamily;

    /**
     * Determine the number of loops in this connected component.
     */
    int num_loops() const;

    /**
     * Determine whether this component corresponds to a vacuum diagram.
     */
    bool is_vacuum() const;

    /**
     * Determine whether this component corresponds to a single-mass
     * vacuum diagram.
     */
    bool is_singlemass_vacuum() const;

    /**
     * Determine whether this component corresponds to a purely phase-space
     * diagram (i.e., all propagators are legally cut).
     */
    bool is_purely_phasespace() const;

    /**
     * Determine whether this component corresponds to an ending diagram
     * (i.e., single-mass vacuum or purely phase-space).
     */
    bool is_ending() const;
};


struct integralfamily {
    /**
     * Create an integral family from the configure information stored in a
     * YAML node.
     * 
     * @param _node configure YAML node
     * @param _psymbols pointer to an external symbol table
     */
    static integralfamily from_yaml(const YAML::Node& _node, GiNaC::symtab* _psymbols);

    std::string to_string() const;      // for debug

    /******************************************************/
    /*       SECTOR-INDEPENDENT FIELDS AND METHODS        */
    /******************************************************/

    std::string         name;
    GiNaC::symtab*      psymbols;

    GiNaC::lst          loops;
    GiNaC::lst          indeplegs;
    GiNaC::lst          propagators;
    std::vector<int>    cut;
    std::vector<int>    prescription;

    GiNaC::lst      conservation_rules;     // momentum conservation
    GiNaC::lst      numeric_rules;          // numeric kinematics
    GiNaC::lst      sps_numeric_rules;      // numeric scalar products

    /**
     * Decide whether the current list of propagators is complete.
     */
    bool       is_complete_propagator() const;

    /**
     * First, extend the current list of propagators into a complete
     * basis. Then represent all [S = L(L+1)/2 + NL] scalar products 
     * into linear combinations of the basis.
     * 
     * @returns An `sps_struct`, consisting of
     *      - `prop`:  the list of complete propagators
     *      - `sps`:   the list of all scalar products
     *      - `coeff`: coefficient matrix (S by S)
     *      - `bias`:  bias matrix (S by 1)
     *      - `cut`:   whether each of the complete propagators is
     *                 cut, if the original propagators have cuts
     * `coeff` and `bias` satisfy `coeff * prop + bias = sps`.
     */
    sps_struct sps_in_complete_propagator() const;

    /**
     * Determine the momentum that flows on each propagator.
     */
    GiNaC::lst momenta() const;

    /**
     * Determine the mass of each propagator.
     */
    GiNaC::lst masses() const;

    /**
     * Determine the loop momentum that flows on each propagator.
     */
    GiNaC::lst loop_momenta() const;


    /******************************************************/
    /*        SECTOR-DEPENDENT FIELDS AND METHODS         */
    /******************************************************/

    // memory for top sector
    std::vector<int>    top_sector_mem;

    // symanzik polynomials
    GiNaC::lst          feynman_params;
    GiNaC::ex           U;
    GiNaC::ex           F;
    GiNaC::ex           F0;     // mass terms in F

    // connected components
    std::vector<component> all_components;

    /**
     * Compute Symanzik polynomials U, F and F0.
     */
    void compute_symanzik(const std::vector<int>& top_sector);

    /**
     * Generate all connected components for the corresponding Feynman
     * diagram in the given top sector.
     */
    void generate_components(const std::vector<int>& top_sector);

    
    /**
     * Insert "eta" to some propagators, according to the list of schemes
     * to insert "eta".
     * 
     * @param top_sector top sector of the master integrals
     * @param schemes a list of schemes indicating strategies to insert
     * "eta"
     * 
     * @returns when there is a valid way to insert "eta" according to the
     * input list of schemes, return a `std::shared_ptr` to a new integral 
     * family, which represents the one with "eta" inserted; when there is
     * no valid way of inserting "eta", return a `std::shared_ptr` that 
     * holds a `nullptr`.
     */
    std::shared_ptr<integralfamily> insert_eta(const std::vector<int>& top_sector, const std::vector<eta_scheme>& schemes);

    /**
     * Determine whether this integral family is a single-mass vacuum type
     * ending family. 
     * 
     * The definition of a single-mass vacuum type ending family is: all
     * connected components correspond to single-mass vacuum diagrams.
     */
    bool is_singlemass_vacuum_ending(const std::vector<int>& top_sector);

    /**
     * Determine whether this integral family is a purely phase-space type
     * ending family.
     * 
     * The definition of a purely phase-space type ending family is: all
     * except one connected components correspond to single-mass vacuum
     * diagrams, with the only left one corresponding to a purely phase-
     * space diagram.
     */
    bool is_purely_phasespace_ending(const std::vector<int>& top_sector);

    /**
     * Determine whether the given top sector is a zero sector.
     */
    bool is_zero(const std::vector<int>& top_sector);
};


#endif // FAMILY_HPP
