#pragma once

#include "DiscretizationBase.hpp"

namespace discretization
{

enum tpfa_method
{
  mo = 0,
  kirill = 1
};

class DiscretizationTPFA : public DiscretizationBase
{
 public:
  DiscretizationTPFA(const DoFNumbering & dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data);

  virtual ~DiscretizationTPFA() = default;

  virtual void build() override;
  static void build_mo(ConnectionData & connection,
                       const ControlVolumeData & cell1,
                       const ControlVolumeData & cell2);

 protected:
  void build_kirill(const mesh::Face & face,
                    ConnectionData                  & data);

  const int m_method;
};

}
