#pragma once

#include "DiscretizationBase.hpp"
#include "ConnectionMap.hpp"

namespace discretization
{

class DiscretizationDFM : public DiscretizationBase
{
 public:
  DiscretizationDFM(const DoFNumbering & dof_numbering,
                    gprs_data::SimData & data,
                    std::vector<ControlVolumeData> & cv_data,
                    std::vector<ConnectionData> & connection_data);

  virtual void build() override;
  // reused by edfm
  static void build_matrix_fracture(ConnectionData & con,
                                    const ControlVolumeData & cv_frac,
                                    const ControlVolumeData & cv_cell);

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
  void build_matrix_fracture_(ConnectionData & con);
};


} // end namespace
