#include "DiscretizationINSIM.hpp"

namespace discretization {

DiscretizationINSIM::DiscretizationINSIM(DoFNumbering const & dof_numbering,
                                         gprs_data::SimData & data,
                                         std::vector<ControlVolumeData> & cv_data,
                                         std::vector<ConnectionData> & connection_data)
    : DiscretizationBase(dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationINSIM::build()
{

}

}  // end namespace discretization
