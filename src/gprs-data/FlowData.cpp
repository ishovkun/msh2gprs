#include <FlowData.hpp>


FlowData::FlowData(const std::size_t max_connections)
    :
    max_connections(max_connections)
{}


void FlowData::merge_elements(const std::size_t updated_element,
                              const std::size_t merged_element)
{
  const std::size_t
      u = updated_element,
      d = merged_element;

  // update cell data
  const double v0 = volumes[u];
  const double v1 = volumes[d];
  volumes[u] += v1;
  poro[u] += (v0*poro[u] + v1*poro[d]) / (v0 + v1);
  depth[u] += (v0*depth[u] + v1*depth[d]) / (v0 + v1);

  for (std::size_t i=0; i<custom_data[u].size(); ++i)
    custom_data[u][i] = (v0*custom_data[u][i] + v1*custom_data[d][i]) / (v0 + v1);

  // delete merged cell data
  volumes.erase(volumes.begin() + d);
  poro.erase(poro.begin() + d);
  depth.erase(depth.begin() + d);
  custom_data.erase(custom_data.begin() + d);

  // update face data
  // delete merged-updated connection
  const std::size_t hash = hash_value(u, d);
  map_connection.erase(map_connection.find(hash));
  // update connection with other elements

}
