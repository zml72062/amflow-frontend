#include "family.hpp"


int component::num_loops() const {
    auto U_expand = U.expand();
    GiNaC::ex try_term = 
        GiNaC::is_exactly_a<GiNaC::add>(U_expand) ? 
            *U_expand.begin() : U_expand;
    
    // how many variables are multiplied together in "try_term"
    int cnt = 0;
    for (auto& idx: propagator_indices)
        if (try_term.has(pfamily->feynman_params[idx]))
            cnt++;
    return cnt;
}


bool component::is_vacuum() const {
    /**
     * A component corresponds to a vacuum diagram if
     * 
     *      (a) the F polynomial without mass terms, seen as a 
     *          polynomial of Feynman parameters belonging to 
     *          the component, has degree at most equal to the 
     *          number of loops in this component
     * 
     *      (b) all propagators in this component are not cut
     */
    auto F_nomass = (pfamily->F - pfamily->U * pfamily->F0).expand();
    int  nloops   = num_loops();
    
    if (GiNaC::is_exactly_a<GiNaC::add>(F_nomass)) {
        for (auto iter = F_nomass.begin(); iter != F_nomass.end(); ++iter) {
            GiNaC::ex term = *iter;
            int cnt = 0;
            for (auto& idx: propagator_indices)
                if (term.has(pfamily->feynman_params[idx]))
                    cnt++;
            if (cnt >= nloops + 1)
                return false;
        }
    } else {
        int cnt = 0;
        for (auto& idx: propagator_indices)
            if (F_nomass.has(pfamily->feynman_params[idx]))
                cnt++;
        if (cnt >= nloops + 1)
            return false;
    }

    // now we have confirmed that F polynomial without mass terms
    // have degree <= num_loops of this component with respect to
    // Feynman parameters belonging to this component
    
    if (pfamily->cut.size() == 0) // no cut at all
        return true;
    
    for (auto& idx: propagator_indices)
        if (pfamily->cut[idx])
            return false;
    
    return true;
}


bool component::is_singlemass_vacuum() const {
    GiNaC::lst family_masses = pfamily->masses();
    int cnt_nonzero_mass = 0;
    for (auto& idx: propagator_indices)
        if (!(bool)(family_masses[idx].expand() == 0))
            cnt_nonzero_mass++;
    
    return is_vacuum() && (cnt_nonzero_mass == 1);
}


bool component::is_purely_phasespace() const {
    /**
     * A component corresponds to a purely phase-space diagram if
     * 
     *      (a) all propagators in this component are cut
     * 
     *      (b) number of propagators in this component equals number
     *          of loops of this component plus 1
     */
    if (pfamily->cut.size() == 0) // no cut at all
        return false;
    
    for (auto& idx: propagator_indices)
        if (!pfamily->cut[idx])
            return false;
    
    return (int)propagator_indices.size() == (num_loops() + 1);
}


bool component::is_ending() const {
    return is_singlemass_vacuum() || is_purely_phasespace();
}


bool integralfamily::is_singlemass_vacuum_ending(const std::vector<int>& top_sector) {
    generate_components(top_sector);
    for (auto& comp: all_components) {
        if (!comp.is_singlemass_vacuum())
            return false;
    }
    return loops.nops() != 0;
}


bool integralfamily::is_purely_phasespace_ending(const std::vector<int>& top_sector) {
    generate_components(top_sector);
    int cnt_phase = 0;
    for (auto& comp: all_components) {
        if (!comp.is_ending())
            return false;
        if (comp.is_purely_phasespace())
            cnt_phase++;
    }
    return cnt_phase == 1;
}


