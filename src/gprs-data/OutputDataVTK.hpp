#pragma once

namespace gprs_data
{

class OutputDataVTK
{
 public:
  OutputDataVTK(const SimData & sim_data, const mesh::Mesh & grid);
  void write_output(const std::string & output_path);

 private:
  const SimData & data;
  const mesh::Mesh & grid;
};

}
