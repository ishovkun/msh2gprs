#include "BoundaryConditionManager.hpp"

namespace gprs_data {

BoundaryConditionManager::BoundaryConditionManager(const std::vector<BCConfig> & face_config,
                                                   const std::vector<BCNodeConfig> & node_config,
                                                   SimData & data)
    :m_face_config(face_config), m_node_config() , m_data(data)
{}

void BoundaryConditionManager::create_properties()
{

}


}  // end namespace gprs_data
