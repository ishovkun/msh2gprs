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
  DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                     const DoFNumbering & combined_dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data,
                     const std::vector<int> & edfm_markers);
  virtual void build() override;

 protected:
  void build_control_volume_data_();
  void build_connection_data_();
  void build_matrix_edfm_(ConnectionData & con);
  void build_matrix_dfm_(ConnectionData & con);
  void identify_edfm_faces_();
  std::vector<size_t> find_edfm_elements_(const ConnectionData & con);
  void create_connections_();
  // ---------------------------- Variables --------------------- //
  const DoFNumbering & m_split_dofs;
  // internal structures to compute dfm discretization after edfm cell splitting
  std::vector<ControlVolumeData> m_split_cv;
  std::vector<ConnectionData> m_split_con;
  std::unordered_set<size_t> m_edfm_faces;
  std::unordered_set<int> m_edfm_markers;
};

}  // end namespace discretization
