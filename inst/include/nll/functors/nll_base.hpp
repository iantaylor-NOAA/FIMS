
#ifndef NLL_BASE_HPP
#define NLL_BASE_HPP

#include "../../common/model_object.hpp"
#include "../../interface/interface.hpp"

namespace fims {
    namespace nll {

        /**
         * Structure to encapsulate nll relevant information.
         */
        template<typename Type>
        struct UnivariateNLLInfo {
            //separate NLLInfo for RE/Priors?
            // - obs is parameter, not data
            // - do not multiply by keep vector to calculate osa resiudals
            // last meeting we floating the idea of making all data parameters
            //  - these can be fixed via mapping
            //  - several people are using this approach
            std::shared_ptr<DataObject<Type> > data;
            //defined in interface.hpp within #ifdef TMB; requires data is a vector
            data_indicator(keep,data);
            //should priors be internal or external to FIMS (or both?)
            fims::Vector<Type> priors; // weights, keep, or priors
            //how are expected values set?
            fims::Vector<Type> expected;
            fims::Vector<Type> parameters;
            #ifdef TMB_MODEL
                ::objective_function<Type>
                *of;  // :: references global namespace, defined in src/FIMS.cpp,
            #endif
        };
      

        /**
         * Base class for NLL Functors.
         */
        template<typename Type>
        struct NLLFunctorBase {
            virtual const Type evaluate(std::shared_ptr<NLLInfo<Type> >& info) = 0;
            //osa is calculated simultaneous to nll, so can't be calculated separately
            //virtual void compute_OSA(std::shared_ptr<LikelihoodInfo<Type> >& info) = 0;
            virtual void compute_metrics(std::shared_ptr<NLLInfo<Type> >& info) = 0;
        }


        template<typename Type>
        struct NLLComponent {
            std::shared_ptr<NLLInfo<Type> > info;
            std::shared_ptr<NLLFunctorBase<Type> > nll_functor;

            //these are calculated using the TMB::oneStepPredict(obj) function and run after optimization 
            //bool record_residuals;
            //fims::Vector<REAL_T> residuals;

            //bool compute_osa;

            //metrics
            Type chi_square = 0.0;
            Type g_test = 0.0;
            Type rmse = 0.0;
            Type rmsle = 0.0;
            Type r_squared = 0.0;
            Type AIC = 0.0; //Akaike’s Information Criterion.
            Type BIC = 0.0; //Bayesian Information Criterion
            uint32_t k; //number of parameters. Used for AIC and BIC calculation.

            const Type evaluate() {
                return this->nll_functor->evaluate(this->info);
            }
        };



    } //namespace nll

}  // namespace fims

#endif /* NLL_BASE_HPP */