#pragma once
#include "DiscretizationBase.hpp"

namespace discretization {

class DiscretizationINSIM : public DiscretizationBase {
 public:
  // constructor
  DiscretizationINSIM(DoFNumbering const & dof_numbering,
                     gprs_data::SimData & data,
                     std::vector<ControlVolumeData> & cv_data,
                     std::vector<ConnectionData> & connection_data);

  // main method. build the discretization
  void build() override;

  // destructor
  virtual ~DiscretizationINSIM() = default;

 private:
  void build_vertex_data_(size_t vertex);
};

}  // end namespace discretization
