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

    /**
     * @brief Extracts a list of derived quantities from the population with this id
     * from the information object.
     * 
     * @return Rcpp::List
     */
    Rcpp::List GetDerivedQuantities() {
        //        Rcpp::List results;

        std::shared_ptr<fims_info::Information<double> > info =
                fims_info::Information<double>::GetInstance();

        std::shared_ptr<fims_popdy::Population<double> > population =
                info->populations[this->id];

        return this->GetDerivedQuantitiesInternal(*population);
    }

    /**
     * @brief Forecast to "fyears" given current model parametrization. 
     * Extracts a list of derived quantities from the population with this id
     * from the information object.
     * 
     * @return Rcpp::List
     */
    Rcpp::List Forecast(int fyears) {

        //        Rcpp::List results;

        std::shared_ptr<fims_info::Information<double> > info =
                fims_info::Information<double>::GetInstance();

        std::shared_ptr<fims_popdy::Population<double> > population =
                info->populations[this->id];



        fims_popdy::Population<double> p;
        p.nyears = fyears;
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
            fleet->Initialize(fyears, p.nages);
            fleet->selectivity = population->fleets[i]->selectivity;
            fleet->log_obs_error = population->fleets[i]->log_obs_error;
            std::fill(fleet->Fmort.begin(), fleet->Fmort.end(), population->fleets[i]->Fmort[population->fleets[i]->Fmort.size() - 1]);
            std::fill(fleet->log_Fmort.begin(), fleet->log_Fmort.end(), population->fleets[i]->log_Fmort[population->fleets[i]->log_Fmort.size() - 1]);
            //            fleet->log_Fmort = population->fleets[i]->log_Fmort;
            fleet->log_q = population->fleets[i]->log_q;
            fleet->is_survey = population->fleets[i]->is_survey;
            fleet->id = population->fleets[i]->id;
            //            fleet->Prepare();
            p.fleets.push_back(fleet);
        }

        p.Initialize(fyears, population->nseasons, population->nages);

        for (int a = 0; a < p.nages; a++) {
            size_t j = (population->nyears - 2) * population->nages + a;
            p.log_init_naa[a] = std::log(population->numbers_at_age[j]);
        }



        p.Prepare();
        p.Evaluate();

        return this->GetDerivedQuantitiesInternal(p);

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

private:

    Rcpp::List GetDerivedQuantitiesInternal(const fims_popdy::Population<double>& p) {
        Rcpp::List results;

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



        results["population_id"] = this->id;
        results["info_type"] = "forecast_quantities";
        results["biomass"] = biomass;
        results["spawning_biomass"] = spawning_biomass;
        results["recruitment"] = recruitment;
        results["numbers_at_age"] = naa;
        results["total_catch"] = total_catch;
        results["Z"] = Z;
        results["F"] = F;

        return results;
    }

};

#endif
