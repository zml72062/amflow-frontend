#ifndef AMFLOW_HPP
#define AMFLOW_HPP


#include "family.hpp"
#include "singlesetup.hpp"
#include <queue>


/**
 * The workflow of AMFlow can be formally characterized as a Finite State 
 * Machine containing two passes: FORWARD PASS and BACKWARD PASS.
 * 
 * Each NODE of the finite state machine corresponds to an integral system,
 * which must be of one of the following six SYSTEM TYPEs.
 * 
 * Each EDGE (i.e., TRANSITION) of the finite state machine must be of one
 * of the following five REDUCTION TYPEs.
 * 
 * In the FORWARD PASS, AMFlow starts from an initial integral system and
 * performs an iterative reduction procedure, and ends up reducing all 
 * integral families to TRIVIAL or ZERO ones.
 * 
 * In the BACKWARD PASS, AMFlow follows the inverse order of the reduction,
 * computing target integrals in each generated system, and eventually 
 * evaluates target integrals in the initial system.
 */


/******   System types   ******/
enum class system_type {
    ZERO,       /* top-level sector is a zero sector      */
    TRIVIAL,    /* no loops, trivially evaluates to 1     */
    NORMAL,     /* not an ending system, eta not inserted */
    ETA,        /* not an ending system, eta inserted     */
    VACUUM,     /* single-mass vacuum type ending system  */
    CUT         /* purely phase-space type ending system  */
};

/******  Reduction types  ******/
enum class reduction_type {
    NO_REDUCE,          /* for ZERO or TRIVIAL systems               */
    INSERT_ETA,         /* NORMAL -> ETA                             */
    BOUNDARY_REDUCE,    /* ETA    -> NORMAL or VACUUM or CUT or ZERO */
    VACUUM_REDUCE,      /* VACUUM -> NORMAL or TRIVIAL               */
    CUT_REDUCE          /* CUT    -> NORMAL                          */
};


struct system_entry {
    system_type         systype;
    reduction_type      transtype;
    integral_system     system;
    std::vector<int>    subsystems;
};


class amflow {
public:
    amflow(const char* config_path, 
           const std::vector<std::vector<int>>& initial_targets,
           const std::vector<std::vector<int>>& initial_preferred,
           const char* work_dir);

    void forward();

    static void usage();
private:
    // message
    const std::string author  = "Xiao Liu and Yan-Qing Ma";
    const std::string version = "1.1";
    const std::string release = "5-Jun-2022";
    void startup();

    YAML::Node  config;
    std::string base_workdir;

    // dependency
    std::string kira;
    std::string fermat;
    void checkdep();
    
    // global symbol table
    GiNaC::symtab globals;

    // integral family
    integralfamily family;

    // list of generated integral systems
    std::vector<system_entry> systems;

    // list of generated integral families
    std::vector<std::shared_ptr<integralfamily>> subfamilies;

    // forward pass
    int             nextnum;
    std::queue<int> enroute;
    void report_forward_status();
};


#endif // AMFLOW_HPP
