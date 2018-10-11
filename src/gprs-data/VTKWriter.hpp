#pragma once

#include <GElement.hpp>
#include <Point.hpp>
using Point = angem::Point<3,double>;

namespace IO
{

class VTKWriter
{
 public:
  static void write_vtk(const std::vector<Point>                    & vertices,
                        const std::vector<std::vector<std::size_t>> & cells,
                        const std::string                           & fname);

  static void write_vtk(const std::vector<Point>    & vertices,
                        const std::vector<Gelement> & elements,
                        const std::string           & fname);
};

}  // end namespace
