#pragma once
#include "DiscretizationBase.hpp"

namespace discretization
{

class DiscretizationDFM : public DiscretizationBase
{
 public:
  DiscretizationEDFM(const std::vector<DiscreteFractureConfig> & dfm_fractures,
                     const std::vector<discretization::ControlVolumeData> & mixed_cv_data,
                     const std::vector<discretization::ConnectionData> & mixed_connection_data,
                     const size_t n_dfm_faces, // const size_t n_cells,
                     gprs_data::SimData & data);
  virtual void build() override;

 protected:
  void extract_control_volume_data_();
  void extract_connection_data_();
  // ---------------------------- Variables --------------------- //
};

}  // end namespace discretization
