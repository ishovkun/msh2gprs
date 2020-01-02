#pragma once

#include "ControlVolumeData.hpp"
#include "ConnectionData.hpp"
#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

/* This is an abstract base class for
 * all discretization classes out there. */
class DiscretizationBase
{
 public:
  DiscretizationBase(const std::vector<DiscreteFractureConfig> & dfm_fractures,
                     gprs_data::SimData & data);

  // get a reference to the face_data vector
  std::vector<ConnectionData> & get_face_data();
  // get a reference to the cell_data vector
  std::vector<ControlVolumeData> & get_cell_data();
  // main method. build the discretization
  virtual void build() = 0;

 protected:
  // build control volumes data
  virtual void build_cell_data_();
  // is a face a dfm fracture
  bool is_fracture(const int marker) const;
  // bould a set of dfm face markers
  void build_dfm_markers_();
  // count the number of dfm faces
  size_t count_dfm_faces_() const;
  // find the maximum control volume index
  size_t find_max_cv_index_() const;
  
  // computed properties
  std::vector<ConnectionData> con_data;
  std::vector<ControlVolumeData> cv_data;

 protected:
  //  input
  const mesh::Mesh & m_grid;
  gprs_data::SimData & m_data;
  const std::vector<DiscreteFractureConfig> & m_dfm_config;
  // stores dfm face markers
  std::set<int> m_dfm_markers;
};

}
