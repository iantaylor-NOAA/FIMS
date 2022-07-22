/*! \file fleet.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project.
 * Refer to the LICENSE file for reuse information.
 * 
 * The purpose of this file is to declare the growth functor class
 * which is the base class for all growth functors.
 */
#ifndef FIMS_POPULATION_DYNAMICS_FLEET_HPP
#define FIMS_POPULATION_DYNAMICS_FLEET_HPP

#include "../../common/model_object.hpp"
#include "../../common/data_object.hpp"
#include "../../distributions/distributions.hpp"
#include "../selectivity/selectivity.hpp"

namespace fims {

    /* @brief Base class for all fleets.
     *
     * @tparam T The type of the fleet object.
     * */
    template<typename T>
    struct Fleet : public FIMSObject<T> {
        static uint32_t id_g; /*!< reference id for fleet object*/


        //data objects
        int observed_index_data_id = -999;
        std::shared_ptr<fims::DataObject<double> > observed_index_data;

        int observed_agecomp_data_id = -999;
        std::shared_ptr<fims::DataObject<double> > observed_agecomp_data;

        //likelihood components
        int index_likelihood_id = -999;
        std::shared_ptr<fims::DistributionsBase<T> > index_likelihood;

        int agecomp_likelihood_id = -999;
        std::shared_ptr<fims::DistributionsBase<T> > agecomp_likelihood;

        //selectivity
        int selectivity_id = -999;
        std::shared_ptr<fims::SelectivityBase<T> > selectivity;

        //derived quantities
        std::vector<T> catch_at_age;
        std::vector<T> catch_index;
        std::vector<T> age_composition;

        /** 
         * @brief Constructor.
         */
        Fleet() {
            this->id = Fleet::id_g++;
        }
        //likelihood is a log likelihood. To do: figure out if these should be
        // negative log likelihood here or in the population loop...Andrea will think about this.
        const T likelihood() {
            return this->index_likelihood->evaluate(do_log = true)
                    + this->agecomp_likelihood->evaluate(do_log = true);
        }

    };
    template <class T>
    uint32_t Fleet<T>::id_g = 0;

} // namespace fims

#endif /* FIMS_POPULATION_DYNAMICS_FLEET_HPP */
