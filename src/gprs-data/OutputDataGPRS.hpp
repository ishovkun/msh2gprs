#pragma once

#include "simdata.hpp"
#include "transes.hpp"

namespace gprs_data
{

class OutputDataGPRS
{
public:
  OutputDataGPRS(SimData & sim_data, mesh::Mesh & grid);

  void write_output(const std::string & output_path);
 private:
  void saveGeometry(const std::string & output_path);
  void saveGeomechDataNewKeywords(const std::string file_name);
  void saveEmbeddedFractureProperties(const std::string file_name);
  void saveBoundaryConditions(const std::string file_name);
  void saveDiscreteFractureProperties(const std::string file_name);
  void saveWells(const std::string file_name);
  void saveFlowMultiScaleData(const std::string file_name);
  void saveMechMultiScaleData(const std::string file_name);

protected:
  SimData & data;
  mesh::Mesh & grid;
  std::vector<mesh::face_iterator> ordered_faces;
};

}
