#include "CellPropertyManager.hpp"

namespace gprs_data {

CellPropertyManager::CellPropertyManager(const CellPropertyConfig & cell_properties,
                                         SimData & data)
: config(cell_properties), data(data)
{}


}  // end namespace gprs_data
