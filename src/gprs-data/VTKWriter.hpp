#pragma once

#include <GElement.hpp>
#include "angem/Point.hpp"
using Point = angem::Point<3,double>;

namespace IO
{

class VTKWriter
{
 public:
  // for fractures
  static void write_surface_geometry(const std::vector<Point>                    & vertices,
                                     const std::vector<std::vector<std::size_t>> & cells,
                                     const std::string                           & fname);
  // for fractures
  static void write_surface_geometry(const std::vector<Point>                    & vertices,
                                     const std::vector<std::vector<std::size_t>> & cells,
                                     std::ofstream                               & out);

  static void write_geometry(const std::vector<Point>                    & vertices,
                        const std::vector<std::vector<std::size_t>> & cells,
                        const std::vector<int>                      & vtk_indices,
                        const std::string                           & fname);

  static void write_geometry(const std::vector<Point>                    & vertices,
                        const std::vector<std::vector<std::size_t>> & cells,
                        const std::vector<int>                      & vtk_indices,
                        std::ofstream                               & out);

  // wicked old timur's Gelement format for reservoir
  static void write_geometry(const std::vector<Point>    & vertices,
                        const std::vector<Gelement> & elements,
                        std::ofstream               & out);

  static void write_geometry(const std::vector<Point>    & vertices,
                             const std::vector<Gelement> & elements,
                             const std::string           & fname);

  static void write_well_trajectory(const std::vector<Point>                              & vertices,
                                    const std::vector<std::pair<std::size_t,std::size_t>> & indices,
                                    const std::string                                     & fname);
  // add cell data to vtk file
  static void add_cell_data(const std::vector<double> & property,
                            const std::string           keyword,
                            std::ofstream             & out);

  static void enter_section_cell_data(const std::size_t n_cells,
                                      std::ofstream & out);

 private:
  VTKWriter();
};

}  // end namespace
