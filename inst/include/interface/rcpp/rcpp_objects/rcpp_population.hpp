/*
 * File:   rcpp_population.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE file
 * for reuse information.
 */
#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_POPULATION_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_POPULATION_HPP

#include "../../../population_dynamics/population/population.hpp"
#include "rcpp_interface_base.hpp"

/**
 * Population Rcpp interface
 */

/**
 * @brief PopulationInterfaceBase class should be inherited to
 * define different Rcpp interfaces for each possible Population function
 */
class PopulationInterfaceBase : public FIMSRcppInterfaceBase {
public:
    static uint32_t id_g; /**< static id of the population interface base*/
    uint32_t id; /**< id of the population interface base */
    // live objects in C++ are objects that have been created and live in memory
    static std::map<uint32_t, PopulationInterfaceBase*>
    live_objects; /**< map associating the ids of PopulationInterfaceBase to
                         the objects */

    PopulationInterfaceBase() {
        this->id = PopulationInterfaceBase::id_g++;
        /* Create instance of map: key is id and value is pointer to
        PopulationInterfaceBase */
        PopulationInterfaceBase::live_objects[this->id] = this;
        PopulationInterfaceBase::fims_interface_objects.push_back(this);
    }

    virtual ~PopulationInterfaceBase() {
    }

    /** @brief get_id method for child classes to inherit */
    virtual uint32_t get_id() = 0;
};

uint32_t PopulationInterfaceBase::id_g = 1;
std::map<uint32_t, PopulationInterfaceBase*>
PopulationInterfaceBase::live_objects;

/**
 * @brief Rcpp interface for a new Population. To instantiate
 * from R:
 * population <- new(population)
 */
class PopulationInterface : public PopulationInterfaceBase {
public:
    uint32_t nages; /**< number of ages */
    uint32_t nfleets; /**< number of fleets */
    uint32_t nseasons; /**< number of seasons */
    uint32_t nyears; /**< number of years */
    uint32_t maturity_id; /**< id of the maturity function*/
    uint32_t growth_id; /**< id of the growth function*/
    uint32_t recruitment_id; /**< id of the recruitment function*/
    Rcpp::NumericVector log_M; /**< log of the natural mortality of the stock*/
    Rcpp::NumericVector log_init_naa; /**<log of the initial numbers at age*/
    Rcpp::NumericVector ages; /**<vector of ages in the population; length nages*/
    Rcpp::NumericVector proportion_female; /**<doule representing the proportion
                                            of female individuals */
    bool estimate_M; /**<whether parameter should be estimated*/
    bool estimate_initNAA; /**<whether parameter should be estimated*/
    bool estimate_prop_female; /**<whether proportion female should be estimated*/

    enum ForecastType {
        CURRENT = 0,
        MEAN,
        UPPER,
        LOWER
    };

    PopulationInterface() : PopulationInterfaceBase() {
    }

    virtual ~PopulationInterface() {
    }

    virtual uint32_t get_id() {
        return this->id;
    }

    /**
     * @brief Set the unique id for the Maturity object
     *
     * @param maturity_id Unique id for the Maturity object
     */
    void SetMaturity(uint32_t maturity_id) {
        this->maturity_id = maturity_id;
    }

    /**
     * @brief Set the unique id for the growth object
     *
     * @param growth_id Unique id for the growth object
     */
    void SetGrowth(uint32_t growth_id) {
        this->growth_id = growth_id;
    }

    /**
     * @brief Set the unique id for the Maturity object
     *
     * @param recruitment_id Unique id for the Maturity object
     */
    void SetRecruitment(uint32_t recruitment_id) {
        this->recruitment_id = recruitment_id;
    }

    /** @brief evaluate the population function */
    virtual void evaluate() {
        fims_popdy::Population<double> population;
        return population.Evaluate();
    }

    Rcpp::List GetDerivedQuantities() {
        Rcpp::List results;

        std::shared_ptr<fims_info::Information<double> > info =
                fims_info::Information<double>::GetInstance();

        std::shared_ptr<fims_popdy::Population<double> > population =
                info->populations[this->id];

        Rcpp::NumericVector biomass(population->biomass.size());
        Rcpp::NumericVector spawning_biomass(population->spawning_biomass.size());
        Rcpp::NumericVector recruitment(population->expected_recruitment.size());
        for (size_t i = 0; i < biomass.size(); i++) {
            biomass[i] = population->biomass[i];
            spawning_biomass[i] = population->spawning_biomass[i];
            recruitment[i] = population->expected_recruitment[i];
        }


        //collapse catch into a single vector for all fleets
        Rcpp::NumericVector total_catch(population->nyears);

        for (size_t i = 0; i < population->nyears; i++) {
            double sum = 0;
            for (size_t j = 0; j < population->nfleets; j++) {
                size_t index = (i * population->nfleets) + j;
                sum += population->expected_catch[index];
            }
            total_catch[i] = sum;
        }
        Rcpp::NumericMatrix naa(population->nyears, population->nages);
        Rcpp::NumericMatrix F(population->nyears, population->nages);
        Rcpp::NumericMatrix Z(population->nyears, population->nages);
        for (size_t i = 0; i < population->nyears; i++) {
            for (size_t j = 0; j < population->nages; j++) {
                size_t index = (i * population->nages) + j;
                naa(i, j) = population->numbers_at_age[index];
                F(i, j) = population->mortality_F[index];
                Z(i, j) = population->mortality_Z[index];
            }
        }



        results = Rcpp::List::create(
                Rcpp::Named("id") = this->id,
                Rcpp::Named("info_type") = "derived_quantities",
                Rcpp::Named("biomass") = biomass,
                Rcpp::Named("spawning_biomass") = spawning_biomass,
                Rcpp::Named("recruitment") = recruitment,
                Rcpp::Named("numbers_at_age") = naa,
                Rcpp::Named("total_catch") = total_catch,
                Rcpp::Named("Z") = Z,
                Rcpp::Named("F") = F);

        return results;
    }

    /*
     * Returns the forecast values of derived quantities given the 
     * current (last year of model fit), mean, 95 percent upper confidence 
     * bound, 95 percent upper confidence bound of fishing mortality.
     */
    Rcpp::List Forecast(int pyears) {
        Rcpp::List results;
        std::shared_ptr<fims_info::Information<double> > info =
                fims_info::Information<double>::GetInstance();

        std::shared_ptr<fims_popdy::Population<double> > population =
                info->populations[this->id];

        results = Rcpp::List::create(
                Rcpp::Named("current") = this->forecast_internal(pyears, CURRENT),
                Rcpp::Named("mean") = this->forecast_internal(pyears, MEAN),
                Rcpp::Named("upper") = this->forecast_internal(pyears, UPPER),
                Rcpp::Named("lower") = this->forecast_internal(pyears, LOWER));

        return results;

    }

    Rcpp::List forecast_internal(int pyears, ForecastType t) {

        Rcpp::List results;

        std::shared_ptr<fims_info::Information<double> > info =
                fims_info::Information<double>::GetInstance();

        std::shared_ptr<fims_popdy::Population<double> > population =
                info->populations[this->id];





        fims_popdy::Population<double> p;


        ///current fishing mortality
        p.nyears = pyears;
        p.ages = population->ages;
        p.nages = population->nages;
        p.id = population->id;
        //        p.nfleets = 
        //        p.fleets.resize(population->nfleets);
        p.nseasons = population->nseasons;
        p.growth = population->growth;
        p.recruitment = population->recruitment;
        p.maturity = population->maturity;



        for (int i = 0; i < population->fleets.size(); i++) {
            std::shared_ptr<fims_popdy::Fleet<double> > fleet = std::make_shared<fims_popdy::Fleet<double> >();
            fleet->Initialize(pyears, p.nages);
            fleet->selectivity = population->fleets[i]->selectivity;
            fleet->log_obs_error = population->fleets[i]->log_obs_error;

            std::pair<double, double> meanStdev1 = this->mean_and_stdev(population->fleets[i]->Fmort);
            std::pair<double, double> meanStdev2 = this->mean_and_stdev(population->fleets[i]->log_Fmort);

            Rcpp::Rcout << meanStdev1.first << " --- " << meanStdev1.second << std::endl;

            if (t == CURRENT) {
                std::fill(fleet->Fmort.begin(), fleet->Fmort.end(), population->fleets[i]->Fmort[population->fleets[i]->Fmort.size() - 1]);
                std::fill(fleet->log_Fmort.begin(), fleet->log_Fmort.end(), population->fleets[i]->log_Fmort[population->fleets[i]->log_Fmort.size() - 1]);
            } else if (t == MEAN) {
                std::fill(fleet->Fmort.begin(), fleet->Fmort.end(), meanStdev1.first);
                std::fill(fleet->log_Fmort.begin(), fleet->log_Fmort.end(), meanStdev2.first);
            } else if (t == UPPER) {
                std::fill(fleet->Fmort.begin(), fleet->Fmort.end(), meanStdev1.first + 2 * meanStdev1.second);
                std::fill(fleet->log_Fmort.begin(), fleet->log_Fmort.end(), meanStdev2.first + 2 * meanStdev2.second);
            } else if (t == LOWER) {
                std::fill(fleet->Fmort.begin(), fleet->Fmort.end(), meanStdev1.first - 2 * meanStdev1.second);
                std::fill(fleet->log_Fmort.begin(), fleet->log_Fmort.end(), meanStdev2.first - 2 * meanStdev2.second);
            }


            //            fleet->log_Fmort = population->fleets[i]->log_Fmort;
            fleet->log_q = population->fleets[i]->log_q;
            fleet->is_survey = population->fleets[i]->is_survey;
            fleet->id = population->fleets[i]->id;
            //            fleet->Prepare();
            p.fleets.push_back(fleet);
        }

        p.Initialize(pyears, population->nseasons, population->nages);

        for (int a = 0; a < p.nages; a++) {
            size_t j = (population->nyears - 2) * population->nages + a;
            p.log_init_naa[a] = std::log(population->numbers_at_age[j]);
        }



        p.Prepare();
        p.Evaluate();

        Rcpp::NumericVector biomass(p.biomass.size());
        Rcpp::NumericVector spawning_biomass(p.spawning_biomass.size());
        Rcpp::NumericVector recruitment(p.expected_recruitment.size());
        for (size_t i = 0; i < biomass.size(); i++) {
            biomass[i] = p.biomass[i];
            spawning_biomass[i] = p.spawning_biomass[i];
            recruitment[i] = p.expected_recruitment[i];
        }

        //collapse catch into a single vector for all fleets
        Rcpp::NumericVector total_catch(p.nyears);

        for (size_t i = 0; i < p.nyears; i++) {
            double sum = 0;
            for (size_t j = 0; j < p.nfleets; j++) {
                size_t index = (i * p.nfleets) + j;
                sum += p.expected_catch[index];
            }
            total_catch[i] = sum;
        }

        Rcpp::NumericMatrix naa(p.nyears, p.nages);
        Rcpp::NumericMatrix F(p.nyears, p.nages);
        Rcpp::NumericMatrix Z(p.nyears, p.nages);
        for (size_t i = 0; i < p.nyears; i++) {
            for (size_t j = 0; j < p.nages; j++) {
                size_t index = (i * p.nages) + j;
                naa(i, j) = p.numbers_at_age[index];
                F(i, j) = p.mortality_F[index];
                Z(i, j) = p.mortality_Z[index];
            }
        }



        results = Rcpp::List::create(
                Rcpp::Named("id") = this->id,
                Rcpp::Named("info_type") = "forecast_quantities",
                Rcpp::Named("biomass") = biomass,
                Rcpp::Named("spawning_biomass") = spawning_biomass,
                Rcpp::Named("recruitment") = recruitment,
                Rcpp::Named("numbers_at_age") = naa,
                Rcpp::Named("total_catch") = total_catch,
                Rcpp::Named("Z") = Z,
                Rcpp::Named("F") = F);

        return results;

    }

    std::pair<double, double> mean_and_stdev(fims::Vector<double>& a) {
        size_t n = a.size();
        std::pair<double, double> ret;
        if (n == 0) {
            return ret;
        }

        double sum = 0;
        double sq_sum = 0;
        for (int i = 0; i < n; ++i) {
            sum += a[i];
            sq_sum += a[i] * a[i];
        }
        double mean = sum / n;
        double variance = sq_sum / n - mean * mean;
        ret.first = mean;
        ret.second = std::sqrt(variance);
        return ret;
    }

#ifdef TMB_MODEL

    template <typename Type>
    bool add_to_fims_tmb_internal() {
        std::shared_ptr<fims_info::Information<Type> > info =
                fims_info::Information<Type>::GetInstance();

        std::shared_ptr<fims_popdy::Population<Type> > population =
                std::make_shared<fims_popdy::Population<Type> >();

        // set relative info
        population->id = this->id;
        population->nyears = this->nyears;
        population->nfleets = this->nfleets;
        population->nseasons = this->nseasons;
        population->nages = this->nages;
        if (this->nages == this->ages.size()) {
            population->ages.resize(this->nages);
        } else {
            warning("The ages vector is not of size nages.");
        }

        population->growth_id = this->growth_id;
        population->recruitment_id = this->recruitment_id;
        population->maturity_id = this->maturity_id;
        population->log_M.resize(this->log_M.size());
        population->log_init_naa.resize(this->log_init_naa.size());
        for (int i = 0; i < log_M.size(); i++) {
            population->log_M[i] = this->log_M[i];
            if (estimate_M) {
                info->RegisterParameter(population->log_M[i]);
            }
        }

        for (int i = 0; i < log_init_naa.size(); i++) {
            population->log_init_naa[i] = this->log_init_naa[i];
            if (estimate_initNAA) {
                info->RegisterParameter(population->log_init_naa[i]);
            }
        }
        for (int i = 0; i < ages.size(); i++) {
            population->ages[i] = this->ages[i];
        }
        for (int i = 0; i < proportion_female.size(); i++) {
            population->proportion_female[i] = this->proportion_female[i];
            if (estimate_prop_female) {
                info->RegisterParameter(population->proportion_female[i]);
            }
        }

        // add to Information
        info->populations[population->id] = population;

        return true;
    }

    /** @brief this adds the parameter values and derivatives to the TMB model
     * object */
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
