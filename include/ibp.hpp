#ifndef IBP_HPP
#define IBP_HPP


#include "integral.hpp"


class ibphelper {
public:
    /**
     * Instantiate an IBP helper from an integral family and a given
     * symbol table.
     * 
     * @param _family the integral family
     * @param _symbols the symbol table
     */
    ibphelper(const integralfamily& _family, const GiNaC::symtab& _symbols)
        : pfamily(&_family), psymbols(&_symbols) { }

    /**
     * Load from a file containing a list of master integrals.
     * 
     * @param masters_path path to the file of master integrals
     */
    void load_masters(const char* masters_path);

    /**
     * Load from a file containing IBP reduction results. One must call
     * `load_masters()` before calling this method.
     * 
     * @param ibp_path path to the file of IBP results
     */
    void load_ibps   (const char* ibp_path);

    // for debug
    void print_masters(std::ostream& out); 
    void print_ibps   (std::ostream& out); 

    std::string   get_master_id            (const std::string& input);
    integral      to_integral              (const std::string& input);
    integral_list get_master_integral_list ();

    friend class kira_agent;
private:
    const integralfamily*   pfamily;
    const GiNaC::symtab*    psymbols;
    GiNaC::symtab           master_table;
    GiNaC::symtab           ibp_table;

    const std::string&      family_name() { return pfamily->name; }
};


#endif // IBP_HPP
