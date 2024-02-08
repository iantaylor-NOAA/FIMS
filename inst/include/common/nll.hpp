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

            const Type Evaluate() {
                return this->nll_functor->Evaluate(this->info);
            }
        };



    }
}
