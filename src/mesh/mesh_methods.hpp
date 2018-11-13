#pragma once

#include <angem/Point.hpp>
#include <angem/PointSet.hpp>


namespace mesh
{

using Point = angem::Point<3,double>;

Point get_element_center(const angem::PointSet<3,double> & vertices,
                         const std::vector<std::size_t>  & ivertices);

}
