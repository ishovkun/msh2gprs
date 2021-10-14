#pragma once

#include "PreprocessorConfig.hpp"  // EDFMMethod
#include "Mesh.hpp"
#include "discretization/flow/DoFNumbering.hpp"
#include "discretization/flow/ControlVolumeData.hpp"
#include "discretization/flow/ConnectionData.hpp"
#include <unordered_set>  //provides unordered_set

namespace gprs_data {

using discretization::DoFNumbering;

/**
 * This class builds dof configurations for flow FVM.
 * The key functions are distribute_dofs and distribute_unsplit_dofs,
 * which generate DoFNumbering for matrix and fracture CVs for compartmental EDFM
 * and regualar/projection EDFM.
 */
class DoFManager
{
 public:
  DoFManager(mesh::Mesh & grid,
             const std::vector<int> dfm_markers,
             const std::vector<int> edfm_markers);
  // map cell/face index -> (p)edfm control volume
  // takes into account the fact that even though we split the grid with fractures,
  // edfm fractures don't really split cells (embedded), and, therefore, two split cells
  // constitute a single control volume
  std::shared_ptr<DoFNumbering> distribute_unsplit_dofs();
  // map cell/face index -> cedfm control volume
  std::shared_ptr<DoFNumbering> distribute_dofs();
  // map vertex index -> reservoir dof
  std::shared_ptr<DoFNumbering> distribute_dofs_insim(std::vector<std::vector<size_t>> const & well_vertices) const;
  // make map : vertex index to well index
  std::shared_ptr<DoFNumbering> distribute_vertex_to_well_dofs(std::vector<std::vector<size_t>> const & well_vertices) const;

  static void remap(std::vector<discretization::ControlVolumeData> & cv_data,
                    std::vector<discretization::ConnectionData>    & connection_data,
                    const DoFNumbering                             & old_dofs,
                    const DoFNumbering                             & new_dofs);

 protected:
  inline bool is_dfm_(const int face_marker) const { return m_set_dfm_markers.find(face_marker) != m_set_dfm_markers.end(); }
  inline bool is_edfm_(const int face_marker) const { return m_set_edfm_markers.find(face_marker) != m_set_edfm_markers.end(); }

 private:
  const mesh::Mesh & m_grid;
  const std::unordered_set<int> m_set_dfm_markers, m_set_edfm_markers;
};

}  // end namespace gprs_data
