#pragma once

#include "simdata.hpp"
#include "mesh/Mesh.hpp"
#include "VTKWriter.hpp"

namespace gprs_data
{

class OutputDataVTK
{
 public:
  OutputDataVTK(const SimData & sim_data, const mesh::Mesh & grid);
  void write_output(const std::string & output_path);

 private:
  void save_reservoir_data(const std::string & fname);

  const SimData & data;
  const mesh::Mesh & grid;
};

}
