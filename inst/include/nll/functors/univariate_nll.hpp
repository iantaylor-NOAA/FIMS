#ifndef UNIVARIATE_NLL_HPP
#define UNIVARIATE_NLL_HPP

#include "../../common/model_object.hpp"
#include "../../interface/interface.hpp"

namespace fims {
    namespace nll {

        /**
         * Normal Likelihood for Data.
         */
        template<typename Type>
        struct NormalDataNLL : public NLLFunctorBase {
            
    
        virtual const Type evaluate(std::shared_ptr<NLLInfo<Type> >& info) {
            std::shared_ptr<DataObject<Type> > observed = info->data;
            for (size_t i = 0; i < this->observed.size(); i++) {
                Type x = observed->at(i);
                Type u = info->expected[i];
                Type sd = fims_math::exp(info->parameters->log_obs_error[i])
                nll -= keep(i) * dnorm(x, u, sd, true);
                //code for osa cdf method
                nll -= keep.cdf_lower(i) * log( pnorm(x, u, sd) );
                nll -= keep.cdf_upper(i) * log( 1.0 - pnorm(x, u, sd) );
                SIMULATE_F(of){
                    observed->at(i) = rnorm(u, sd);
                }
                REPORT_F(observed->at(i), of);
            }
            retun nll;
        }

        /**
         * Normal Likelihood for Parameter.
         */
        template<typename Type>
        struct NormalParmNLL : public NLLFunctorBase {
            
    
        virtual const Type evaluate(std::shared_ptr<NLLInfo<Type> >& info) {
            //how to link in any parameter vector?
            for (size_t i = 0; i < this->observed.size(); i++) {
                Type x = info->process[i];
                Type u = info->parameters->mu; //for RE model, this would always be 0
                Type sd = fims_math::exp(info->parameters->log_proc_error);
                //do not use keep indicator vectors on priors or random effects
                nll -= dnorm(x, u, sd, true);
                SIMULATE_F(of){
                    info->process[i] = rnorm(u, sd);
                }
                //would need to preserve the names somehow
                REPORT_F(info->process, of);
            }
            retun nll;
        }

    } //namespace nll
} //namespace fims