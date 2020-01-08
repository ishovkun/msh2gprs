#pragma once

#include "ControlVolumeData.hpp"
#include "ConnectionData.hpp"
#include "DoFNumbering.hpp"
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
  DiscretizationBase(const DoFNumbering & dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data);

  // main method. build the discretization
  virtual void build() = 0;

 protected:
  // build control volumes data
  virtual void build_cell_data_();

 protected:
  //  input
  const mesh::Mesh & m_grid;
  gprs_data::SimData & m_data;
  // stores dfm face markers
  const DoFNumbering & m_dofs;
  // computed properties
  std::vector<ControlVolumeData> & m_cv_data;
  std::vector<ConnectionData> & m_con_data;
};

}
