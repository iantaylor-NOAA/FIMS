
#ifndef FIMS_POPULATION_DYNAMICS_RECRUITMENT_RECRUITMENT_RCALLBACK_HPP
#define FIMS_POPULATION_DYNAMICS_RECRUITMENT_RECRUITMENT_RCALLBACK_HPP


#include "recruitment_base.hpp"

namespace fims_popdy {

/** @brief BevertonHolt class that returns the Beverton Holt SR
 * from fims_math.
 *
 * @param logit_steep Recruitment relative to unfished recruitment at
 * 20% of unfished spawning biomass. Should be a value between 0.2 and 1.0.
 */
template <typename Type>
struct RecruitmentRCallback : public RecruitmentBase<Type> {
    
    Rcpp::Function f;
    Rcpp::List arguments;
    
    
    virtual const Type evaluate(const Type& spawners, const Type& phi_0) {
        arguments["spawners"] = spawners;
        arguments["phi_0"] = phi_0;
        
       return Rcpp::as<Type>(f(arguments);
        
    }
    
    
    
};



#endif
