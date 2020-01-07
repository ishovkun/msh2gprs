#pragma once

#include "DiscretizationBase.hpp"
#include "ConnectionMap.hpp"

namespace discretization
{

class DiscretizationDFM : public DiscretizationBase
{
 public:
  DiscretizationDFM(const std::vector<DiscreteFractureConfig> & dfm_fractures,
                    gprs_data::SimData & data);

  virtual void build() override;

  void merge_from_matrix_discretization(const std::vector<discretization::ControlVolumeData> & matrix_cv_data,
                                        const std::vector<discretization::ConnectionData> & matrix_connection_data);

 protected:
  // build data like volumes, depth, poro, etc.
  virtual void build_cell_data_() override;
  // build M-F connection list
  void build_fracture_matrix_connections();
  // build F-F connection list
  void build_fracture_fracture_connections();
  // build map edge to neighboring faces
  hash_algorithms::ConnectionMap<std::vector<size_t>> map_edge_to_faces();
  // build matrix-fracture connection
  void build_matrix_fracture(ConnectionData & con);

  // storage for the properties of dfm fractures
  // //  numbering shift of matrix CVs
  // const size_t shift_matrix;
  // // numbering shift of dfm CVs
  // const size_t shift_dfm;
  // map cell index -> control volume
  // std::vector<std::size_t> m_cell_to_cv;
  // vector of fracture apertures
  std::vector<double> m_apertures;
};


} // end namespace
