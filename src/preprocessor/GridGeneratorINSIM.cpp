#include "GridGeneratorINSIM.hpp"

namespace gprs_data {

GridGeneratorINSIM::GridGeneratorINSIM(std::vector<WellConfig> const & wells)
    : _wells(wells)
{}

GridGeneratorINSIM::operator mesh::Mesh() const
{
  return mesh::Mesh();
}

}  // end namespace gprs_data
