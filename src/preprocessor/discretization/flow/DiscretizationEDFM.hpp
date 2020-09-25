#pragma once
#include "DiscretizationBase.hpp"
#include "ConnectionMap.hpp"
#include "PreprocessorConfig.hpp"  // edfm_method

namespace discretization
{

/* This class implements full flow TPFA discretization for grid with DFM and EDFM.
 * It supports three treatment methods for DFM: simple, pEEDFM, and cEEDFM.
 * Neither DFM nor EDFM fractures are combined into non-split elements.
 * Split cells are combined if cedfm is not used. */
class DiscretizationEDFM : public DiscretizationBase
{
 public:
  DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                     const DoFNumbering & combined_dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data,
                     const std::vector<int> & edfm_markers);

  virtual ~DiscretizationEDFM() = default;
  void build() override;

 protected:
  void build_control_volume_data_();
  void build_connection_data_();
  void build_matrix_edfm_(ConnectionData & con);
  void build_matrix_dfm_(ConnectionData & con);
  void identify_edfm_faces_();
  std::vector<size_t> find_edfm_elements_(const ConnectionData & con);
  std::vector<std::size_t> create_connections_();
  void build_pedfm_();
  std::vector<const mesh::Face*> pedfm_select_faces_(const mesh::Face & frac_face) const;
  size_t pedfm_find_other_cell_(const mesh::Face & frac, const mesh::Face & other) const;
  // pedfm for a single connection
  // returns true if need to kill the connection
  bool build_pedfm_(ConnectionData & mm_con, ConnectionData & fm_con);
  std::pair<size_t,size_t> find_fracture_cv_and_nonfracture_cv_(const ConnectionData & mm_con, const ConnectionData & fm_con) const;
  // ---------------------------- Variables --------------------- //
  const DoFNumbering & m_split_dofs;
  // internal structures to compute dfm discretization after edfm cell splitting
  std::vector<ControlVolumeData> m_split_cv;
  std::vector<ConnectionData> m_split_con;
  std::unordered_set<size_t> m_edfm_faces;
  std::unordered_set<int> m_edfm_markers;
  hash_algorithms::ConnectionMap<ConnectionData> m_con_map;
};

}  // end namespace discretization
