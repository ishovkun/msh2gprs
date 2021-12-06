#pragma once

#include "config/FiniteElementConfig.hpp"
#include "config/MeshConfig.hpp"
#include "config/MultiscaleConfig.hpp"
#include "config/PropertyConfig.hpp"
#include "config/WellConfig.hpp"
#include "config/OutputConfig.hpp"
#include "angem/Polygon.hpp"

#include <map>
#include <memory> // shared / unique_ptr

enum class BoundaryConditionType : int
{
  dirichlet = 1, neumann = 2, constraint = 3  // penalty type forces all disp to be the same
};

// Structure that holds info on mechanical user-defined boundary conditions
struct BCConfig
{
  int label;
  BoundaryConditionType type;
  std::string location_expression;
  std::array<std::string,3> values_expressions;
  // undefined boundary value (do not impose)
  static constexpr double nan = std::numeric_limits<double>::max();
};


// enum class EDFMMethod
// {
//   simple,      // old and frankly shitty but default by Li and Lee 2008
//   projection,  // pEDFM by Tene 2017
//   compartmental, // cEDFM by Chai 2018
// };

enum class FlowDiscretizationType
{
  tpfa_edfm,  // tpfa + standard edfm [Li and Lee 2008]
  tpfa_projection,  // tpfa + projection EDFM [Tene, 2017]
  tpfa_compartmental,  // tpfa + compartmental edfm [Chai, 2018]
  insim,               // wells connected by 1D tubes [Guo & Reynolds, 2017]
};


struct EmbeddedFractureConfig
{
  std::shared_ptr<angem::Polygon<double>> body;  // embedded fractures
  std::size_t n1 = 0, n2 = 0;                    // remeshing params
  double cohesion = 0;
  double friction_angle = 30;
  double dilation_angle = 0;
  double aperture = 1e-3;   // hydraulic aperture of the fracture [m]
  double conductivity = 10; // hydraulic conductivity of dfm fracture [m·md]
  bool coupled = true;      // whether to couple with geomechanics
  size_t region = 0;
};


struct DiscreteFractureConfig
{
  int label;
  double conductivity = 10;
  double aperture = 1e-3;
  bool coupled = true;  // whether to couple with geomechanics
  size_t region = 0;    // property table index
};

enum class FracturePlacement
{
  move_fracture, move_grid
};

enum class EDFMSearchAlgorithm
{
  robust, fast
};

struct EDFMSettings
{
  double min_dist_to_node = 1e-4; // minimum distance to grid vertices relative to cell size
  double vertex_split_tolerance = 2e-4;  // relative tolerance for merging diplicate vertices during split
  FracturePlacement placement = FracturePlacement::move_fracture;
  EDFMSearchAlgorithm algorithm = EDFMSearchAlgorithm::robust;
};

struct DFMSettings
{
  bool split_mech_vertices = true;
};

struct PreprocessorConfig
{
  std::vector<EmbeddedFractureConfig>  embedded_fractures;  //  embedded  fractures
  std::vector<DiscreteFractureConfig>  discrete_fractures;
  std::vector<BCConfig>                bc_faces;
  std::vector<BCConfig>                bc_nodes;
  std::vector<WellConfig>              wells;

  FlowDiscretizationType flow_discretization = FlowDiscretizationType::tpfa_edfm;

  FiniteElementConfig fem;
  EDFMSettings edfm_settings;
  DFMSettings dfm_settings;
  // double edfm_min_dist_to_node = 1e-4;               // minimum distance to grid vertices relative to cell size
  // global container for all cell properties
  CellPropertyConfig cell_properties;

  /* global properties that can only specified as files
   * this names can be used as variables in each subdomain.
   * Don't assign gprs keywords to these props as they won't be outputted.
   */
  std::vector<std::string>             input_property_file_names;
  std::vector<std::string>             input_global_varialbes;
  // expressions and variables for each subdomain
  double geometry_search_tolerance = 1e-10;
  double frac_cell_elinination_factor = 0.2;

  // multiscale
  MultiscaleConfig ms_flow;
  MultiscaleConfig ms_mech;

  // output format
  OutputConfig output;
  // the name of gmsh grid file
  MeshConfig mesh;
};
