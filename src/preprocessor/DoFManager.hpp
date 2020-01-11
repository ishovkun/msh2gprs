#pragma once

#include "PreprocessorConfig.hpp"  // EDFMMethod
#include "Mesh.hpp"
#include "discretization/DoFNumbering.hpp"
#include <unordered_set>  //provides unordered_set

namespace gprs_data {

using discretization::DoFNumbering;

class DoFManager
{
 public:
  DoFManager(mesh::Mesh & grid,
             const std::vector<int> dfm_markers,
             const std::vector<int> edfm_markers);
  DoFNumbering distribute_unsplit_dofs();
  DoFNumbering distribute_dofs();

 protected:
  inline bool is_dfm_(const int face_marker) const { return m_set_dfm_markers.find(face_marker) != m_set_dfm_markers.end(); }
  inline bool is_edfm_(const int face_marker) const { return m_set_edfm_markers.find(face_marker) != m_set_edfm_markers.end(); }

 private:
  const mesh::Mesh & m_grid;
  const std::unordered_set<int> m_set_dfm_markers, m_set_edfm_markers;
};

}  // end namespace gprs_data
