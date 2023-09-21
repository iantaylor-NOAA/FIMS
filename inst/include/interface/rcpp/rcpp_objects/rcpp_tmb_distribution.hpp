/*
 * File:   rcpp_distributions.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE in the
 * source folder for reuse information.
 *
 */
#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_TMB_DISTRIBUTION_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_TMB_DISTRIBUTION_HPP

#include "../../../distributions/distributions.hpp"
#include "../../interface.hpp"
#include "rcpp_interface_base.hpp"

/**
 * @brief Distributions Rcpp Interface
 *
 */
class DistributionsInterfaceBase : public FIMSRcppInterfaceBase {
public:
    static uint32_t
    id_g; /**< static id of the DistributionsInterfaceBase object */
    uint32_t id; /**< local id of the DistributionsInterfaceBase object */
    static std::map<uint32_t, DistributionsInterfaceBase *> live_objects; /**<
  map relating the ID of the DistributionsInterfaceBase to the
  DistributionsInterfaceBase objects */

    DistributionsInterfaceBase() {
        this->id = DistributionsInterfaceBase::id_g++;
        DistributionsInterfaceBase::live_objects[this->id] = this;
        FIMSRcppInterfaceBase::fims_interface_objects.push_back(this);
    }

    virtual ~DistributionsInterfaceBase() {
    }

    /** @brief get_id method for child distribution interface objects to inherit
     * **/
    virtual uint32_t get_id() = 0;

    /** @brief evaluate method for child distribution interface objects to inherit
     * **/
    virtual double evaluate(bool do_log) = 0;
};

uint32_t DistributionsInterfaceBase::id_g =
        1; /**< static id of the DistributionsInterfaceBase object */
std::map<uint32_t,
DistributionsInterfaceBase
*> /**< local id of the DistributionsInterfaceBase object */
DistributionsInterfaceBase::live_objects; /**<
  map relating the ID of the DistributionsInterfaceBase to the
  DistributionsInterfaceBase objects */

/**
 * @brief Rcpp interface for Dnorm as an S4 object. To instantiate
 * from R:
 * dnorm_ <- new(fims$TMBDnormDistribution)
 *
 */
class DnormDistributionsInterface : public DistributionsInterfaceBase {
public:
    Parameter x; /**< observed data */
    Parameter mean; /**< mean of x for the normal distribution **/
    Parameter sd; /**< sd of x for the normal distribution **/

    DnormDistributionsInterface() : DistributionsInterfaceBase() {
    }

    virtual uint32_t get_id() {
        return this->id;
    }

    virtual ~DnormDistributionsInterface() {
    }

    /**
     * @brief Evaluate normal probability density function, default returns the
     * log of the pdf
     *
     * @tparam T
     * @return log pdf
     */
    virtual double evaluate(bool do_log) {
        fims::Dnorm<double> dnorm;
        dnorm.x = this->x.value_m;
        dnorm.mean = this->mean.value_m;
        dnorm.sd = this->sd.value_m;
        return dnorm.evaluate(do_log);
    }

#ifdef TMB_MODEL

    template<typename Type>
    bool add_to_fims_tmb_internal() {
        std::shared_ptr<fims::Information < Type>> info =
                fims::Information<Type>::GetInstance();

        std::shared_ptr<fims::Dnorm < Type>> distribution =
                std::make_shared<fims::Dnorm < Type >> ();

        // interface to data/parameter value
        distribution->id = this->id;
        distribution->x = this->x.value_m;
        // set relative info
        distribution->mean = this->mean.value_m;
        distribution->sd = this->sd.value_m;

        info->distribution_models[distribution->id] = distribution;

        return true;

    }

    /**
     * @brief adds the dnorm distribution and its parameters to the TMB model
     */
    virtual bool add_to_fims_tmb() {

        this->add_to_fims_tmb_internal<TMB_FIMS_REAL_TYPE>();
        this->add_to_fims_tmb_internal<TMB_FIMS_FIRST_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_SECOND_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_THIRD_ORDER>();

        return true;
    }

#endif
};

/**
 * @brief Rcpp interface for Dlnorm as an S4 object. To instantiate
 * from R:
 * dlnorm_ <- new(fims$TMBDlnormDistribution)
 *
 */
class DlnormDistributionsInterface : public DistributionsInterfaceBase {
public:
    Parameter x; /*!< observation */
    Parameter meanlog; /*!< mean of the distribution of log(x) */
    Parameter sdlog; /*!< standard deviation of the distribution of log(x) */

    DlnormDistributionsInterface() : DistributionsInterfaceBase() {
    }

    virtual ~DlnormDistributionsInterface() {
    }

    /**
     * @brief get the id of the Dlnorm distributions interface class object
     */
    virtual uint32_t get_id() {
        return this->id;
    }

    /**
     * @brief Evaluate lognormal probability density function, default returns the
     * log of the pdf
     *
     * @tparam T
     * @return log pdf
     */
    virtual double evaluate(bool do_log) {
        fims::Dlnorm<double> dlnorm;
        dlnorm.x = this->x.value_m;
        dlnorm.meanlog = this->meanlog.value_m;
        dlnorm.sdlog = this->sdlog.value_m;
        return dlnorm.evaluate(do_log);
    }

#ifdef TMB_MODEL

    template<typename Type>
    bool add_to_fims_tmb_internal() {
        std::shared_ptr<fims::Information < Type> > info =
                fims::Information<Type>::GetInstance();

        std::shared_ptr < fims::Dlnorm <Type> > distribution =
                std::make_shared<fims::Dlnorm <Type> > ();

        // set relative info
        distribution->id = this->id;
        distribution->x = this->x.value_m;
        distribution->meanlog = this->meanlog.value_m;
        distribution->sdlog = this->sdlog.value_m;

        info->distribution_models[distribution->id] = distribution;

        return true;
    }

    /**
     * @brief adds the dlnorm distribution and its parameters to the TMB model
     */
    virtual bool add_to_fims_tmb() {

        this->add_to_fims_tmb_internal<TMB_FIMS_REAL_TYPE>();
        this->add_to_fims_tmb_internal<TMB_FIMS_FIRST_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_SECOND_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_THIRD_ORDER>();

        return true;
    }

#endif
};

/**
 * @brief Rcpp interface for Dmultinom as an S4 object. To instantiate
 * from R:
 * dmultinom_ <- new(fims$TMBDmultinomDistribution)
 *
 */
// template <typename T>

class DmultinomDistributionsInterface : public DistributionsInterfaceBase {
public:
    Rcpp::NumericVector x; /*!< Vector of length K of integers */
    Rcpp::NumericVector p; /*!< Vector of length K, specifying the probability
   for the K classes (note, unlike in R these must sum to 1). */

    DmultinomDistributionsInterface() : DistributionsInterfaceBase() {
    }

    virtual ~DmultinomDistributionsInterface() {
    }

    virtual uint32_t get_id() {
        return this->id;
    }

    /**
     * @brief Evaluate multinom probability density function, default returns the
     * log of the pdf
     *
     * @tparam T
     * @return log pdf
     */
    virtual double evaluate(bool do_log) {
        fims::Dmultinom<double> dmultinom;
        // Decale TMBVector in this scope
        typedef
        typename fims::ModelTraits<TMB_FIMS_REAL_TYPE>::EigenVector TMBVector;
        dmultinom.x = TMBVector(x.size()); // Vector from TMB
        dmultinom.p = TMBVector(p.size()); // Vector from TMB
        for (int i = 0; i < x.size(); i++) {
            dmultinom.x[i] = x[i];
            dmultinom.p[i] = p[i];
        }
        return dmultinom.evaluate(do_log);
    }

#ifdef TMB_MODEL

    template<typename T>
    bool add_to_fims_tmb_internal() {
        typedef typename fims::ModelTraits<T>::EigenVector Vector0;
        std::shared_ptr<fims::Information < T> > info =
                fims::Information<T>::GetInstance();

        std::shared_ptr<fims::Dmultinom < T>> distribution =
                std::make_shared<fims::Dmultinom < T > > ();

        distribution->id = this->id;
        distribution->x = Vector0(x.size());
        distribution->p = Vector0(p.size());

        for (int i = 0; i < x.size(); i++) {
            distribution->x[i] = x[i];
            distribution->p[i] = p[i];
        }

        info->distribution_models[distribution->id] = distribution;

        return true;
    }

    virtual bool add_to_fims_tmb() {

        this->add_to_fims_tmb_internal<TMB_FIMS_REAL_TYPE>();
        this->add_to_fims_tmb_internal<TMB_FIMS_FIRST_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_SECOND_ORDER>();
        this->add_to_fims_tmb_internal<TMB_FIMS_THIRD_ORDER>();

        return true;
    }

#endif
};
#endif
