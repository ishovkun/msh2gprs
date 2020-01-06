#pragma once
#include "DiscretizationBase.hpp"
#include "ConnectionMap.hpp"

namespace discretization
{

/* This class implements flow discretization for EDFM.
 * Currently it can only interpret discretization data
 * produced by DiscretizationDFM class that was applied
 * after cutting the cells containing EDFM fractures.
 * Thus, this class only computes transmissibilities from the data
 * obtained previously. */
class DiscretizationEDFM : public DiscretizationBase
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
  void build_connection_data_();
  void build_matrix_fracture_(const ConnectionData & con);
  void build_edfm_edfm_(const ConnectionData & con);
  void build_edfm_dfm_(const ConnectionData & con);
  size_t calculate_edfm_faces_() const;
  std::vector<size_t> find_edfm_elements_(const ConnectionData & con);
  void convert_flow_map_to_vector_();
  // ---------------------------- Variables --------------------- //
  // references to combiined mixed external props
  const std::vector<ConnectionData> & m_con;
  const std::vector<ControlVolumeData> & m_cv;
  hash_algorithms::ConnectionMap<ConnectionData> m_con_map;
  const size_t m_min_edfm_index;
  const size_t m_n_edfm_faces;
};

}  // end namespace discretization
