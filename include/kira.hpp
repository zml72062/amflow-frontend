#ifndef KIRA_HPP
#define KIRA_HPP


#include "ibp.hpp"


class kira_agent {

public:
    /**
     * Instantiate a Kira agent from configure information stored in a YAML
     * node and a given integral family.
     * 
     * @param _node configure YAML node
     * @param _family an integral family
     */
    kira_agent(const YAML::Node& _node, const integralfamily& _family);
    
    /**
     * Let Kira find master integrals of a given integral family, and
     * optionally reduce a list of target integrals to masters.
     * 
     * @param workdir Kira working directory
     * @param target a list of target integrals
     * @param preferred a list of preferred masters
     * @param top_sector the top sector
     * @param rank the maximum rank
     * @param dot the maximum dot
     * @param masters_only when set `true`, only find master integrals
     * and do not reduce
     * 
     * @returns a pair of objects:
     *    (1) an integral list consisting of all master integrals
     *    (2) a matrix whose (i, j) element is the expansion coefficient
     *        of the i-th target integral into the j-th master integral
     * 
     * When `masters_only` is set `true`, the returned coefficient matrix
     * will be empty.
     */
    std::pair<integral_list, GiNaC::matrix>
    reduce(const char* workdir, const integral_list& target, 
           const integral_list& preferred, unsigned long top_sector, 
           unsigned rank, unsigned dot, bool masters_only);


    /**
     * Let Kira find trivial sectors of a given integral family which
     * are subsectors of a given top sector.
     * 
     * @param workdir Kira working directory
     * @param top_sector the top sector
     * 
     * @returns a list of trivial sectors
     */
    std::vector<unsigned long> trivial_sectors(const char* workdir, unsigned long top_sector);

private:
    std::string           kira;
    std::string           fermat;
    const integralfamily* pfamily;
    ibphelper             helper;
    int                   integral_order;
    enum {
        Masters, 
        Kira, 
        FireFly, 
        Mixed, 
        NoFactorScan
    }                     reduction_mode;
    GiNaC::ex             d0;
    int                   nthreads;

    // exported from "pfamily"
    const std::string&         name()         { return pfamily->name; }
    const GiNaC::symtab*       psymbols()     { return pfamily->psymbols; }
    const GiNaC::lst&          loops()        { return pfamily->loops; }
    const GiNaC::lst&          indeplegs()    { return pfamily->indeplegs; }
    const GiNaC::lst&          propagators()  { return pfamily->propagators; }
    const std::vector<int>&    cut()          { return pfamily->cut; }
    const std::vector<int>&    prescription() { return pfamily->prescription; }
    const GiNaC::lst&  conservation_rules() { return pfamily->conservation_rules; }
    const GiNaC::lst&  numeric_rules()      { return pfamily->numeric_rules; }
    const GiNaC::lst&  sps_numeric_rules()  { return pfamily->sps_numeric_rules; }

    /**
     * Write the content of "integralfamilies.yaml" to an output stream.
     * 
     * @param stream the output stream
     * @param top_sector the top sector
     */
    void write_integralfamilies_yaml(std::ostream& stream, unsigned long top_sector);

    /**
     * Write the content of "kinematics.yaml" to an output stream.
     * 
     * @param stream the output stream
     */
    void write_kinematics_yaml      (std::ostream& stream);

    /**
     * Write the content of "jobs.yaml" to an output stream.
     * 
     * @param stream the output stream
     * @param top_sector the top sector
     * @param rank the maximum rank
     * @param dot the maximum dot
     */
    void write_jobs_yaml            (std::ostream& stream, unsigned long top_sector,
                                     unsigned rank, unsigned dot);

    /**
     * Write preferred integrals to an output stream.
     * 
     * @param stream the output stream
     * @param preferred a list of preferred integrals
     */
    void write_preferred            (std::ostream& stream, const integral_list& preferred);

    /**
     * Write target integrals to an output stream.
     * 
     * @param stream the output stream
     * @param target a list of target integrals
     */
    void write_target               (std::ostream& stream, const integral_list& target);

    /**
     * Run Kira at current working directory.
     * 
     * @returns 0 if Kira normally exited, -1 if an error occurred
     */
    int run_kira();
};


#endif // KIRA_HPP
