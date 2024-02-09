#ifndef FIMS_COMMON_NLL_HPP
#define FIMS_COMMON_NLL_HPP

#include "fims_vector.hpp"

namespace fims {
    namespace common {

        /**
         * Structure to encapsulate nll relevant information.
         */
        template<typename Type>
        struct LikelihoodInfo {
            std::shared_ptr<DataObject<Type> > data;
            fims::Vector<Type> lambdas; // weights, keep, or priors

        };

        /**
         * Base class for NLL Functors.
         */
        template<typename Type>
        struct LikelihoodFunctorBase {
            virtual const Type Evaluate(std::shared_ptr<LikelihoodInfo<Type> >& info) = 0;
            virtual void ComputeOSA(std::shared_ptr<LikelihoodInfo<Type> >& info) = 0;
            virtual void ComputeMetrics(std::shared_ptr<LikelihoodInfo<Type> >& info) = 0;
        };

        /**
         * Logistic Likelihood.
         */
        template<typename Type>
        struct LogisticNLL : public LikelihoodFunctorBase {

            virtual const Type Evaluate(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                return static_cast<Type> (0);
            }

            virtual void ComputeOSA(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                throw std::runtime_error("Not yet implemented.\n");
            }

            virtual void ComputeMetrics(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                throw std::runtime_error("Not yet implemented.\n");
            }
        };

        template<typename Type>
        struct MultinomialNLL : public LikelihoodFunctorBase {

            virtual const Type Evaluate(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                return static_cast<Type> (0);
            }

            virtual void ComputeOSA(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                throw std::runtime_error("Not yet implemented.\n");
            }

            virtual void ComputeMetrics(std::shared_ptr<LikelihoodInfo<Type> >& info) {
                throw std::runtime_error("Not yet implemented.\n");
            }
        };

        template<typename Type>
        struct NLLComponent {
            std::shared_ptr<LikelihoodInfo<Type> > info;
            std::shared_ptr<LikelihoodFunctorBase<Type> > nll_functor;

            bool record_residuals;
            fims::Vector<REAL_T> residuals;

            bool compute_osa;

            //metrics
            Type chi_square = 0.0;
            Type g_test = 0.0;
            Type rmse = 0.0;
            Type rmsle = 0.0;
            Type r_squared = 0.0;
            Type AIC = 0.0; //Akaike’s Information Criterion.
            Type BIC = 0.0; //Bayesian Information Criterion
            uint32_t k; //number of parameters. Used for AIC and BIC calculation.

            const Type Evaluate() {
                return this->nll_functor->Evaluate(this->info);
            }
        };



    }
}
