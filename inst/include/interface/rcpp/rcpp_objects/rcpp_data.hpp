/*
 * File:   rcpp_fleet.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE file
 * for reuse information.
 */
#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_DATA_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_DATA_HPP

#include "../../../common/information.hpp"
#include "rcpp_interface_base.hpp"

/**
 * @brief Rcpp interface for Data as an S4 object. To instantiate
 * from R:
 * fleet <- new(fims$Data)
 *
 */
class DataInterface : public FIMSRcppInterfaceBase {
public:
    Rcpp::NumericVector observed_data; /*!< The data */
    static uint32_t id_g; /**< static id of the DataInterface object */
    uint32_t id; /**< local id of the DataInterface object */
    static std::map<uint32_t, DataInterface*>
    live_objects; /**< map associating the ids of DataInterface to
      the objects */

    /** @brief constructor
     */
    DataInterface() {
        this->id = DataInterface::id_g++;
        DataInterface::live_objects[this->id] = this;
        FIMSRcppInterfaceBase::fims_interface_objects.push_back(this);
    }

    /** @brief destructor
     */
    virtual ~DataInterface() {
    }

    /** @brief get the ID of the interface base object
     **/
    virtual uint32_t get_id() {
        return this->id;
    }

    /**@brief add_to_fims_tmb dummy method
     *
     */
    virtual bool add_to_fims_tmb() {
        return true;
    };
};
uint32_t DataInterface::id_g = 1;
std::map<uint32_t, DataInterface*> DataInterface::live_objects;

/**
 * @brief Rcpp interface for age comp data as an S4 object. To instantiate
 * from R:
 * acomp <- new(fims$AgeComp)
 */
class AgeCompDataInterface : public DataInterface {
public:
    int amax; /*!< first dimension of the data */
    int ymax; /*!< second dimension of the data */
    Rcpp::NumericVector age_comp_data; /*!<the age composition data*/

    /**
     * @brief constructor
     */
    AgeCompDataInterface(int ymax = 0, int amax = 0) : DataInterface() {
        this->amax = amax;
        this->ymax = ymax;
    }

    /**
     * @brief destructor
     */
    virtual ~AgeCompDataInterface() {
    }

    /** @brief get the ID of the interface base object
     **/
    virtual uint32_t get_id() {
        return this->id;
    }


#ifdef TMB_MODEL

    template<typename T>
    bool add_to_fims_tmb_internal() {
        std::shared_ptr<fims::DataObject < T>> age_comp_data =
                std::make_shared<fims::DataObject < T >> (this->ymax,
                this->amax);

        age_comp_data->id = this->id;
        for (int y = 0; y < ymax; y++) {
            for (int a = 0; a < amax; a++) {
                int index_ya = y * amax + a;
                age_comp_data->at(y, a) = this->age_comp_data[index_ya];
            }
        }

        std::shared_ptr<fims::Information < T>> info =
                fims::Information<T>::GetInstance();

        info->data_objects[this->id] = age_comp_data;
        
         return true;
    }

    /**
     * @brief adds parameters to the model
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
 * @brief Rcpp interface for data as an S4 object. To instantiate
 * from R:
 * fleet <- new(fims$Index)
 */
class IndexDataInterface : public DataInterface {
public:
    int ymax; /*!< second dimension of the data */
    Rcpp::NumericVector index_data; /*!<the age composition data*/

    /**
     * @brief constructor
     */
    IndexDataInterface(int ymax = 0) : DataInterface() {
        this->ymax = ymax;
    }

    /**
     * @brief destructor
     */
    virtual ~IndexDataInterface() {
    }

    /** @brief get the ID of the interface base object
     **/
    virtual uint32_t get_id() {
        return this->id;
    }


#ifdef TMB_MODEL

    template<typename T>
    bool add_to_fims_tmb_internal() {
        std::shared_ptr<fims::DataObject < T>> data =
                std::make_shared<fims::DataObject < T >> (this->ymax);

        data->id = this->id;

        for (int y = 0; y < ymax; y++) {
            data->at(y) = this->index_data[y];
        }

        std::shared_ptr<fims::Information < T>> info =
                fims::Information<T>::GetInstance();

        info->data_objects[this->id] = data;
        return true;
    }

    /**
     *@brief function to add to TMB
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

#endif
