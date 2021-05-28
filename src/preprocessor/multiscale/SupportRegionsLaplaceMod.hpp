#pragma once
#include "SupportRegionsBase.hpp"
#include "SimData.hpp"

namespace multiscale {

class SupportRegionsLaplaceMod :  public SupportRegionsBase  {
 public:
  SupportRegionsLaplaceMod(std::vector<size_t> const &partition, gprs_data::SimData &data);
  virtual ~SupportRegionsLaplaceMod() = default;
  std::vector<double> const & get(size_t block) const {return _support[block];}

 protected:
  size_t find_center_(std::vector<size_t> const  & region) const;
  gprs_data::SimData const &_data;
  std::vector<std::vector<double>> _support;
};

}  // end namespace multiscale
