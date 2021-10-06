#pragma once

#include "angem/Polygon.hpp"
#include "config/FiniteElementConfig.hpp"
#include "config/MeshConfig.hpp"
#include "config/MultiscaleConfig.hpp"
#include "config/PropertyConfig.hpp"
#include "config/WellConfig.hpp"

#include <map>
#include <memory> // shared / unique_ptr


enum class OutputFormat
{
  gprs, vtk, postprocessor
};

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


enum class EDFMMethod
{
  simple,      // old and frankly shitty but default by Li and Lee 2008
  projection,  // pEDFM by Tene 2017
  compartmental, // cEDFM by Chai 2018
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


struct VTKOutputConfig
{
  std::string flow_reservoir_grid_file      = "flow_reservoir_mesh.vtk";
  std::string mechanics_reservoir_grid_file = "mechanics_reservoir_mesh.vtk";
  std::string fracture_grid_file            = "fractures.vtk";
  std::string wells_file                    = "wells.vtk";
};

struct GPRSOutputConfig
{
  std::string geometry_file          = "gm_geometry.txt";
  std::string mechanics_kwd_file     = "gm_keywords.txt";
  std::string efrac_file             = "gm_SDA.txt";
  std::string discrete_frac_file     = "gm_DFM.txt";
  std::string bcond_file             = "gm_bcond.txt";
  std::string wells_file             = "wells.txt";
  std::string mech_ms_file           = "ms_mech.txt";
  std::string flow_ms_file           = "ms_flow.txt";
  std::string flow_cv_file           = "fl_cell_data.txt";
  std::string flow_connection_file   = "fl_face_data.txt";
  std::string mech_trans_update_file = "gm_update_trans.txt";
  std::string fem_file               = "gm_fem.txt";
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
  EDFMMethod method = EDFMMethod::simple;       // method to simulate flow in embedded fracs
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
  std::vector<OutputFormat> output_formats = {OutputFormat::gprs,
                                              OutputFormat::vtk,
                                              OutputFormat::postprocessor};

  // the name of gmsh grid file
  MeshConfig mesh_config;
  // std::string mesh_file;
  // output file names
  std::string output_dir            = "output";
  // GPRS format
  GPRSOutputConfig gprs_output;
  // VTK format
  VTKOutputConfig vtk_config;
  // postprocessor output file
  std::string postprocessor_file = "postprocessor_config.yaml";
  // postprocessor output directory
  std::string postprocessor_output_dir = "postprocessing";
};
