#pragma once

#include "PreprocessorConfig.hpp"
#include "SimData.hpp"
#include "muparser/muParser.h" // parser for user-defined expressions for reservoir data

namespace gprs_data {

class CellPropertyManager
{
 public:
  CellPropertyManager(const CellPropertyConfig & cell_properties,
                      const std::vector<DomainConfig> & domain_configs,
                      SimData & data);
  void generate_properties();

 private:
  void print_setup_message_();
  void assign_expressions_(const DomainConfig& domain,
                           std::vector<mu::Parser> & parsers,
                           std::vector<double> & vars);
  void evaluate_expressions_(const DomainConfig& domain,
                           std::vector<mu::Parser> & parsers,
                           std::vector<double> & vars);

  const CellPropertyConfig        & config;
  const std::vector<DomainConfig> & domains;
  SimData & data;
  // number of default variable in config
  // these variables are not output, so I don't save them
  // should be 3=x+y+z
  const size_t m_shift;
};

}  // end namespace gprs_data
