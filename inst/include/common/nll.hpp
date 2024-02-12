#ifndef FIMS_COMMON_NLL_HPP
#define FIMS_COMMON_NLL_HPP

#include "fims_vector.hpp"

#define SQUARE(x) ((x)*(x))

namespace fims {
    namespace common {

        /**
         * Structure to encapsulate nll relevant information.
         */
        template<typename Type>
        struct LikelihoodInfo {
            std::shared_ptr<DataObject<Type> > data;
            fims::Vector<Type> lambdas; // weights, keep, or priors
            fims::Vector<Type> expected;
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
        struct LognormalNLL : public LikelihoodFunctorBase {
            Type sigma = 0.2;
            
            bool use_bias_correction = false;

            virtual const Type Evaluate(std::shared_ptr<LikelihoodInfo<Type> >& info) {

                std::shared_ptr<DataObject<Type> > observed = info->data;
                size_t i, j, k;
                REAL_T obs, se, se2, cv;
                this->years = observed->imax;
                this->seasons = observed->jmax;
                this->ages = observed->kmax;
                this->neff = static_cast<REAL_T> (0.0);
                Type nll = static_cast<REAL_T> (0.0);
                Type expected;

                Type nll1;
                Type nll2;
                switch (observed->dimensions) {
                    case 2:

                        if (this->use_bias_correction) {
                            for (i = 0; i < this->years; i++) {
                                for (j = 0; j < this->seasons; j++) {
                                    size_t index = i * this->seasons + j;
                                    expected = info->expected[index];

                                    REAL_T obs = observed->get(i, j);
                                    
                                    if (obs != observed->missing_value) {
                                        
                                        cv = observed->get_error(i, j);
                                        se2 = std::log(cv * cv + 1.0);
                                        se = std::sqrt(se2) / std::sqrt(std::log(M_E));
                                        nll1 += info->lambda->get(i, j) * std::log(se);
                                        //                                     nll2 += this->lambda->get(i, j) * SQUARE((mas::log(obs) - mas::log(expected)))/se;
                                        nll2 += info->lambda->get(i, j) * SQUARE((mas::log((obs / expected)) / se) + 0.5 * se);
                                    }
                                }

                            }
                        } else {
                            for (i = 0; i < this->years; i++) {
                                for (j = 0; j < this->seasons; j++) {
                                    size_t index = i * this->seasons + j;
                                    expected = info->expected[index];

                                    REAL_T obs = observed->get(i, j);
                                    
                                    if (obs != observed->na_value) {
                                        
                                        cv = observed->get_error(i, j);
                                        se2 = std::log(cv * cv + 1.0);
                                        se = std::sqrt(se2) / std::sqrt(std::log(M_E));
                                        nll1 += info->lambda->at(i, j) * std::log(se);
                                        nll2 += info->lambda->at(i, j) * SQUARE(mas::log((obs / expected))) / se2;
                                    }
                                }

                            }
                        }
                        nll = nll1 + 0.5 * nll2;

                        break;
                    case 3:

                        throw std::invalid_argument("Expected 2 dimensional observation data for log-normal likelihood function.");

                        break;
                }



                return snll;
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
