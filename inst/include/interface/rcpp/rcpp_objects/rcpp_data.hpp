/*
 * File:   rcpp_fleet.hpp
 *
 * This File is part of the NOAA, National Marine Fisheries Service
 * Fisheries Integrated Modeling System project. See LICENSE file
 * for reuse information.
 */
#ifndef FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_DATA_HPP
#define FIMS_INTERFACE_RCPP_RCPP_OBJECTS_RCPP_DATA_HPP

#include "../../../common/Information.hpp"
#include "rcpp_interface_base.hpp"

/**
 * @brief Rcpp interface for Fleet as an S4 object. To instantiate
 * from R:
 * fleet <- new(fims$Fleet)
 *
 */
class DataInterface : public FIMSRcppInterfaceBase
{
public:
  Rcpp::NumericVector observed_data; /*!< The data */
  static uint32_t id_g; /**< static id of the DataInterface object */
  uint32_t id;          /**< local id of the DataInterface object */

  /** @brief constructor
  */
  DataInterface(Rcpp::NumericVector data_) { 
    this->id = DataInterface::id_g++; 
  }

  /** @brief destructor
  */
  virtual ~DataInterface() {}
};

uint32_t DataInterface::id_g = 1;

class AgeCompDataInterface : public DataInterface
{
public:
  size_t amax; /*!< first dimension of the data */
  size_t ymax; /*!< second dimension of the data */
  Rcpp::NumericVector age_comp_data; /*!<the age composition data*/

  AgeCompDataInterface(size_t amax = 0, size_t ymax = 0, 
  Rcpp::NumericVector data_ = R_NilValue) : DataInterface(age_comp_data),
   age_comp_data(data_) { }

  virtual ~AgeCompDataInterface() {}

  virtual bool add_to_fims_tmb()
  {
    std::shared_ptr<fims::DataObject<double>> age_comp_data =
        std::make_shared<fims::DataObject<double>>(this->amax, this->ymax);

    age_comp_data->id = this->id;

    for (size_t y = 0; y < ymax; y++)
    {
      for (size_t a = 0; a < amax; a++)
      {
        size_t index_ya = y * amax + a;
        age_comp_data->at(y, a) = this->observed_data[index_ya];
      }
    }

    std::shared_ptr<fims::Information<TMB_FIMS_REAL_TYPE>> d0 =
        fims::Information<TMB_FIMS_REAL_TYPE>::GetInstance();

    d0->data_objects[this->id] = age_comp_data;

    std::shared_ptr<fims::Information<TMB_FIMS_FIRST_ORDER>> d1 =
        fims::Information<TMB_FIMS_FIRST_ORDER>::GetInstance();

    d1->data_objects[this->id] = age_comp_data;
    std::shared_ptr<fims::Information<TMB_FIMS_SECOND_ORDER>> d2 =
        fims::Information<TMB_FIMS_SECOND_ORDER>::GetInstance();

    d2->data_objects[this->id] = age_comp_data;
    std::shared_ptr<fims::Information<TMB_FIMS_THIRD_ORDER>> d3 =
        fims::Information<TMB_FIMS_THIRD_ORDER>::GetInstance();

    
    d3->data_objects[this->id] = age_comp_data;

    return true;
  }
};

class IndexDataInterface : public DataInterface
  {
  public:
    size_t ymax; /*!< second dimension of the data */
    Rcpp::NumericVector index_data; /*!<the age composition data*/

    IndexDataInterface(size_t ymax = 0, 
    Rcpp::NumericVector data_ = R_NilValue) : DataInterface(index_data),
    index_data(data_)  { }

    virtual ~IndexDataInterface() {}

    virtual bool add_to_fims_tmb()
    {

      std::shared_ptr<fims::DataObject<double>> index_data =
        std::make_shared <fims::DataObject<double>>(this->ymax);

    index_data->id = this->id;

    for (size_t y = 0; y < ymax; y++)
    {
        index_data->at(y) = this->observed_data[y];
      }

    std::shared_ptr<fims::Information<TMB_FIMS_REAL_TYPE>> d0 =
        fims::Information<TMB_FIMS_REAL_TYPE>::GetInstance();

    d0->data_objects[this->id] = index_data;

    std::shared_ptr<fims::Information<TMB_FIMS_FIRST_ORDER>> d1 =
        fims::Information<TMB_FIMS_FIRST_ORDER>::GetInstance();

    d1->data_objects[this->id] = index_data;
    std::shared_ptr<fims::Information<TMB_FIMS_SECOND_ORDER>> d2 =
        fims::Information<TMB_FIMS_SECOND_ORDER>::GetInstance();

    d2->data_objects[this->id] = index_data;
    std::shared_ptr<fims::Information<TMB_FIMS_THIRD_ORDER>> d3 =
        fims::Information<TMB_FIMS_THIRD_ORDER>::GetInstance();

    
    d3->data_objects[this->id] = index_data;
    return true;
    }
  };

#endif