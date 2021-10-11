#include "YamlParser.hpp"
#include "yaml-cpp/node/parse.h"  // loadfile
#include "angem/Rectangle.hpp"
#include "logger/Logger.hpp"
#include <algorithm>  // std::foreach
#include <string>
#include <cctype>  // tolower


namespace Parsers
{

using std::string;

YamlParser::YamlParser()
{}

void YamlParser::parse_file(const std::string & fname)
{
  YAML::Node main_node = YAML::LoadFile(fname);

  if (fname.substr(fname.size() - 4, fname.size()) != "yaml")
    throw std::invalid_argument("wrong file type");

  for(YAML::const_iterator it=main_node.begin();it != main_node.end();++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "Reading section: " << key  << std::endl;

    if (key == "Mesh file")
    {
      config.mesh.type = MeshType::file;
      config.mesh.file = it->second.as<std::string>();
    }
    else if (key == "Mesh")
      section_mesh(it->second);
    else if (key == "Domain Flow Properties")
      section_domain_props(it->second, VariableType::flow);
    else if (key == "Domain Mechanical Properties")
      section_domain_props(it->second, VariableType::mechanics);
    else if (key == "Embedded Fractures")
      embedded_fracs(it->second);
    else if (key == "Discrete Fractures")
      discrete_fracs(it->second);
    else if (key == "Boundary Conditions")
      boundary_conditions(it->second);
    else if (key == "Wells")
      section_wells(it->second);
    else if (key == "Multiscale")
      section_multiscale(it->second);
    else
      logging::warning() << "Unknown key: " << key << " skipping" << std::endl;
  }
}

void YamlParser::embedded_fracs(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "\treading entry " << key << std::endl;

    if (key == "file")
      config.gprs_output.efrac_file = it->second.as<std::string>();
    // else if (key == "method")
    // {
    //   const std::string str_method = it->second.as<std::string>();
    //   if (str_method == "simple")
    //     config.edfm_settings.method = EDFMMethod::simple;
    //   else if (str_method == "projection")
    //     config.edfm_settings.method = EDFMMethod::projection;
    //   else if (str_method == "compartmental")
    //     config.edfm_settings.method = EDFMMethod::compartmental;
    //   else
    //     throw std::invalid_argument(str_method);
    // }
    else if (key == "mech-method")
    {
      const std::string str_method = it->second.as<std::string>();
      if (str_method == "strong discontinuity")
        config.fem.method = FEMMethod::strong_discontinuity;
      else if (str_method == "discrete" || str_method == "pfem")
        config.fem.method = FEMMethod::polyhedral_finite_element;
      else if (str_method == "mixed")
        config.fem.method = FEMMethod::mixed;
      else throw std::invalid_argument("Unknown EDFM mechanics method " + str_method);
    }
    else if (key == "subdivision")
    {
      const auto values = it->second.as<std::pair<std::string,size_t>>();
      if (values.first == "gmsh")
        config.fem.subdivision_method = PolyhedralFEMSubdivision::gmsh_generate;
      else if (values.first == "refinement")
        config.fem.subdivision_method = PolyhedralFEMSubdivision::refinement;
      else throw std::invalid_argument("Subdivision method unknown" + values.first);
      config.fem.order = values.second;
    }
    else if (key == "integration_rule")
    {
      const auto str_rule = it->second.as<std::string>();
      if (str_rule == "full")
        config.fem.integration_rule = PolyhedronIntegrationRule::Full;
      else if (str_rule == "faces_average")
        config.fem.integration_rule = PolyhedronIntegrationRule::FacesAverage;
      else if (str_rule == "vertices_average")
        config.fem.integration_rule = PolyhedronIntegrationRule::VerticesAverage;
      else if (str_rule == "faces_pointwise")
        config.fem.integration_rule = PolyhedronIntegrationRule::FacesPointwise;
      else if (str_rule == "vertices_pointwise")
        config.fem.integration_rule = PolyhedronIntegrationRule::VerticesPointwise;
    }
    else if (key == "fracture")
    {
      config.embedded_fractures.emplace_back();
      embedded_fracture(it->second, config.embedded_fractures.back());
    }
    else if (key == "min_distance_to_vertices")
    {
      config.edfm_settings.min_dist_to_node = it->second.as<double>();
    }
    else if (key == "relative_vertex_tolerance")
    {
      config.edfm_settings.vertex_split_tolerance = it->second.as<double>();
    }
    else if (key == "placement")
    {
      const auto value = it->second.as<std::string>();
      if ( value == "move_fracture" )
        config.edfm_settings.placement = FracturePlacement::move_fracture;
      else if ( value == "move_grid" )
        config.edfm_settings.placement = FracturePlacement::move_grid;
      else throw std::invalid_argument("wrong edfm placing strategy");
    }
    else if (key == "solver")
    {
      const auto values = it->second.as<std::pair<std::string,double>>();
      if (values.first == "direct")
        config.fem.solver = SolverType::direct;
      else if (values.first == "cg")
        config.fem.solver = SolverType::cg;
      else if (values.first == "msrsb")
        config.fem.solver = SolverType::msrsb;
      config.fem.solver_tolerance = values.second;
    }
    else
      throw std::invalid_argument( "\t\tattribute " + key + " unknown: skipping");
  }
}

void YamlParser::discrete_fracs(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "\treading entry " << key << std::endl;

    if (key == "fracture")
    {
      config.discrete_fractures.emplace_back();
      // check if frac has a label
      auto & conf = config.discrete_fractures.back();
      extract_subnode_value("label", it, conf.label);
      if (conf.label == 0)
        throw std::invalid_argument("Cannot have label 0");
      discrete_fracture(it->second, conf);
    }
    else if (key == "file")
      config.gprs_output.discrete_frac_file = it->second.as<std::string>();
    else if (key == "split_vertices")
      config.dfm_settings.split_mech_vertices = it->second.as<bool>();
    else
      logging::warning() << "\t\tattribute " << key << " unknown: skipping" << std::endl;
  }

}

void YamlParser::embedded_fracture(const YAML::Node       & node,
                                   EmbeddedFractureConfig & conf)
{
  // default parameters
  std::string shape = "Rectangle";
  double height = 1;  // for rectangles and ellipses
  double length = 1;  // for rectangles and ellipses
  double dip = 90;
  double strike = 0;
  double cohesion = 0;
  angem::Point<3,double> center = {0, 0, 0};
  double friction_angle = conf.friction_angle;
  double dilation_angle = conf.dilation_angle;
  double aperture = conf.aperture;
  double conductivity = conf.conductivity;
  std::size_t n1 = conf.n1;
  std::size_t n2 = conf.n2;
  size_t region = conf.region;

  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "\t\treading entry " << key << std::endl;

    if (key == "type")
    {
      shape = it->second.as<std::string>();
      if (shape != "Rectangle")
        throw std::invalid_argument("non-rectangular fractures not supported");
    }
    else if (key == "length")
      length = it->second.as<double>();
    else if (key == "width" or key == "height")
      height = it->second.as<double>();
    else if (key == "strike")
      strike = it->second.as<double>();
    else if (key == "dip")
      dip = it->second.as<double>();
    else if (key == "cohesion")
      cohesion = it->second.as<double>();
    else if (key == "friction angle")
      friction_angle = it->second.as<double>();
    else if (key == "dilation angle")
      dilation_angle = it->second.as<double>();
    else if (key == "aperture")
      aperture = it->second.as<double>();
    else if (key == "conductivity")
      conductivity = it->second.as<double>();
    else if (key == "region")
      region = it->second.as<size_t>();
    else if (key == "remesh")
    {
      std::vector<std::size_t> remesh_pars = it->second.as<std::vector<std::size_t>>();

      if (remesh_pars.size() != 2)
        logging::warning() << "wrong entry in remesh. Aborting" << std::endl << std::flush;

      n1 = remesh_pars[0];
      n2 = remesh_pars[1];

      if ((n1 > 0 and n2 == 0)  or (n2 > 0 and n1 == 0))
        throw std::invalid_argument("wrong entry in remesh");
    }
    else if (key == "center")
      center = it->second.as<std::vector<double>>();
    else throw std::invalid_argument("attribute unknown");
  }

  logging::log() << "Making embedded fracture" << std::endl;
  conf.body = std::make_shared<angem::Rectangle<double>>
      (center, length, height, dip, strike);

  conf.cohesion = cohesion;
  conf.friction_angle = friction_angle;
  conf.dilation_angle = dilation_angle;
  conf.aperture = aperture;
  conf.conductivity = conductivity;
  conf.n1 = n1; conf.n2 = n2;
  conf.region = region;
}

void YamlParser::discrete_fracture(const YAML::Node       & node,
                                   DiscreteFractureConfig & conf)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "\t\treading entry " << key << std::endl;

    if (key == "aperture")
      conf.aperture = it->second.as<double>();
    else if (key == "conductivity")
      conf.conductivity = it->second.as<double>();
    else if (key == "label")
      conf.label = it->second.as<int>();
    else if (key == "region")
      conf.region = it->second.as<size_t>();
  }
}

void YamlParser::section_domain_props(const YAML::Node & node, const VariableType var_type)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    logging::log() << "\tReading key " << key << std::endl;

    if (key == "files")  // these are common for all flow domains
    {
      domain(it->second, var_type, config.cell_properties.files);
    }
    else if ( key == "discretization" && var_type == VariableType::flow ) {
      config.flow_discretization = parse_flow_discretization_( it->second.as<std::string>() );
    }
    else if (key == "domain")
    {
      int label;
      try {
        label = it->second["label"].as<int>();
      }
      catch (YAML::TypedBadConversion<int> & error) {
        throw std::invalid_argument("domain label must be an integer");
      }
      catch (YAML::InvalidNode & error) {
        throw std::invalid_argument("domain label must be specified");
      }

      // DomainConfig & conf = get_domain_config(label);
      config.cell_properties.domains.emplace_back();
      auto & conf = config.cell_properties.domains.back();
      conf.label = label;

      logging::log() << "\tReading domain " << conf.label << std::endl;
      domain(it->second, var_type, conf);
    }
    else throw std::invalid_argument("attribute " + key + " is unknown");
  }
}



// find an item in a vector
template<typename T>
std::size_t find(const T & item, const std::vector<T> & vec)
{
  for (std::size_t i=0; i<vec.size(); ++i)
    if (vec[i] == item)
      return i;
  return vec.size();
}


// more generic implement of the previous func
template<typename iterable>
std::size_t find(const typename iterable::value_type & item,
                 const iterable & container)
{
  std::size_t counter = 0;
  for (const auto & iter_item : container)
    if (iter_item == item)
      return counter;
    else
      counter++;
  return counter;
}

void YamlParser::domain(const YAML::Node & node, const VariableType var_type, DomainConfig & conf)
{
  conf.type = var_type;
  std::size_t exp_counter = conf.expressions.size();

  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    const std::string value = it->second.as<std::string>();

    logging::log() << "\t\treading entry " << key << std::endl;

    if (key == "Coupled") // coupling with mechanics
    {
      conf.coupled = it->second.as<bool>();
      if (!conf.coupled) logging::warning() << "domain " << conf.label << " is decoupled !!!" << std::endl;
    }
    else if (key == "label")  // already parsed at the top level
      continue;
    else if (key == "files")
    {
      // these are to read files that contain properties.
      // these properties can only be used as variables in expressions.
      // they will not be output.
      for (auto it_var = it->second.begin(); it_var != it->second.end(); ++it_var) {
        config.input_global_varialbes.push_back(it_var->first.as<std::string>());
        config.input_property_file_names.push_back(it_var->second.as<std::string>());
      }
    }
    else
    {
      // save property name and expression
      conf.variables.push_back(key);
      conf.expressions.push_back(value);
      // const std::size_t ind = find(key, config.cell_properties.all_vars);
      // if (ind == config.cell_properties.all_vars.size())  // new variable
      // {
      //   config.cell_properties.all_vars.push_back(key);
      //   // special case - service variable (not outputted)
      //   if ( find(key, config.cell_properties.special_keywords) < config.cell_properties.special_keywords.size())
      //     config.cell_properties.expression_type.push_back(ExpressionDomainType::service);
      //   else
      //     config.cell_properties.expression_type.push_back(var_type);
      // }

      // // save positions in the global list
      // conf.local_to_global_vars[exp_counter] = ind;
      // conf.global_to_local_vars[ind]         = exp_counter;
      // exp_counter++;
    }
  }
}

void YamlParser::boundary_conditions(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading entry " << key << std::endl;

    if (key == "file")
      config.gprs_output.bcond_file = it->second.as<std::string>();
    else if (key == "Faces")
      boundary_conditions_faces(it->second);
    else if (key == "Dirichlet nodes")
      boundary_conditions_nodes(it->second);
    else
    {
      throw std::invalid_argument("Attribute " + key + " is unknown");
    }
  }

}

void YamlParser::boundary_conditions_faces(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\treading entry " << key << std::endl;

    if (key == "face")
    {
      config.bc_faces.emplace_back();
      auto & conf = config.bc_faces.back();
      bc_face(it->second, conf);
    }
    else
    {
      throw std::invalid_argument("\t\tAttribute " + key + " is unknown");
    }
  }
}

void YamlParser::bc_face(const YAML::Node & node,
                         BCConfig         & conf)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\t\treading entry " << key << std::endl;

    if (key == "type")
    {
      string strkey = it->second.as<string>();
      std::for_each(strkey.begin(), strkey.end(), [](char & c){c = ::tolower(c);});
      if ( strkey == "1" || strkey == "dirichlet" )
        conf.type = BoundaryConditionType::dirichlet;
      else if ( strkey == "2" || strkey == "neumann" )
        conf.type = BoundaryConditionType::neumann;
      else if ( strkey == "3" || strkey == "constraint" )
        conf.type = BoundaryConditionType::constraint;
      else throw std::invalid_argument("Wrong BC type");
    }
    else if (key == "label")
      conf.label = it->second.as<int>();
    else if (key == "location")
      conf.location_expression = it->second.as<std::string>();
    else if (key == "value")
    {
      const std::vector<std::string> str_values = it->second.as<std::vector<std::string>>();
      if (str_values.size() != 3)
        throw std::invalid_argument("value should have exactly 3 entries");

      for (std::size_t i=0; i<3; ++i)
      {
        conf.values_expressions[i] = str_values[i];
      }
    }
    else
      throw std::invalid_argument("attribute " + key + " is unknown");
  }
}

void YamlParser::boundary_conditions_nodes(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\t\treading entry " << key << std::endl;

    if (key == "node")
    {
      config.bc_nodes.emplace_back();
      auto & conf = config.bc_nodes.back();
      bc_node(it->second, conf);
    }
    else
    {
      throw std::invalid_argument("Attribute " + key + " is unknown");
    }
  }
}

void YamlParser::bc_node(const YAML::Node & node, BCConfig & conf)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\t\treading entry " << key << std::endl;

    if (key == "location")
      conf.location_expression = it->second.as<std::string>();
    else if (key == "value")
    {
      conf.values_expressions = it->second.as<std::array<std::string,3>>();
    }
    else
      throw std::invalid_argument("\t\tAttribute " + key + " is unknown");
  }
}

// DomainConfig & YamlParser::get_domain_config(const int label)
// {
  // find current domain id in existing config
  // int counter = 0;
  // for (const auto & domain : config.domains)
  // {
  //   if (label == domain.label)
  //     break;
  //   else
  //     counter++;
  // }

  // DomainConfig * p_conf;
  // if (counter == config.domains.size())
  // {
  //   config.domains.emplace_back();
  //   config.domains.back().label = label;
  //   return config.domains.back();
  // }
  // else
  //   return config.domains[counter];
// }

void YamlParser::section_wells(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;

    if (key == "file")
      config.gprs_output.wells_file = it->second.as<std::string>();
    else if (key == "well")
    {
      config.wells.emplace_back();
      read_well(it->second, config.wells.back());
    }
    else
      std::cout << "\tSkipping unknown keyword" << std::endl;
  }
}

void YamlParser::read_well(const YAML::Node & node,
                           WellConfig       & well)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\treading entry " << key << std::endl;

    if (key == "name")
      well.name = it->second.as<std::string>();
    else if (key == "radius")
      well.radius = it->second.as<double>();
    else if (key == "force_connect_fractures")
      well.force_connect_fractures = it->second.as<bool>();
    else if (key == "nodes")
    {
      const std::string line = it->second.as<std::string>();
      std::istringstream iss(line);
      std::vector<std::string> tokens{std::istream_iterator<std::string>{iss},
                                      std::istream_iterator<std::string>{}};
      if (tokens.size() % 3 != 0)
      {
        std::cout << "\t\tinvalid entry" << std::endl;
        abort();
      }

      const std::size_t n_coord = tokens.size() / 3;
      well.coordinates.resize(n_coord);
      for (std::size_t i=0; i<n_coord; ++i)
        for (int d=0; d<3; ++d)
          well.coordinates[i][d] = std::atof(tokens[3*i+d].c_str());
    }
    else if (key == "perforations")
    {
      std::vector<int> tokens = it->second.as<std::vector<int>>();
      for (const auto & token : tokens)
        well.perforated.push_back(static_cast<bool>(token));
    }
    else
      std::cout << "\t\tunknown key; skipping" << std::endl;
  }

  if (well.coordinates.empty())
  {
    std::cout << "invalid entry" << std::endl;
    abort();
  }

  if (well.coordinates.size() > 1)
  {
    if (well.perforated.empty())
      well.perforated.resize(well.coordinates.size() - 1, true);
    else if (well.perforated.size() != well.coordinates.size() - 1)
      throw std::invalid_argument("Invalid perforations size for well " + well.name);
  }
}

void YamlParser::section_multiscale(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;

    if (key == "flow") {
      subsection_multiscale(it->second, config.ms_flow);
    }
    else if ( key == "mechanics" ) {
      subsection_multiscale(it->second, config.ms_mech);
      assert( config.ms_mech.support_type == MSSupportType::mechanics );
    }
    else
    {
      std::cout << "\tunknown keyword aborting" << std::endl;
      abort();
    }
  }
}

void YamlParser::section_mesh(const YAML::Node & node)
{
  auto & conf = config.mesh;
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;
    if (key == "file")
    {
      conf.type = MeshType::file;
      conf.file = it->second.as<std::string>();
    }
    else if (key == "cartesian")
    {
      conf.type = MeshType::cartesian;
      subsection_grid_cartesian(it->second);
    }
    else if (key == "refinement")
    {
      subsection_grid_refinement(it->second);
    }
    else if (key == "insim") {
      conf.type = MeshType::insim;
      subsection_grid_insim(it->second);
    }
    else throw std::invalid_argument("Unknown key " + key);
  }
}

void YamlParser::subsection_grid_insim(const YAML::Node & node)
{
  auto & conf = config.mesh.insim;
  for (auto it = node.begin(); it!=node.end(); ++it) {
    const std::string key = it->first.as<std::string>();
    if ( key == "padding_fraction" )
      conf.padding_fraction = it->second.as<double>();
    else if ( key == "cell_label" )
      conf.cell_label = it->second.as<int>();
  }
}

void YamlParser::subsection_grid_cartesian(const YAML::Node & node)
{
  auto & conf = config.mesh.cartesian;
  std::array<size_t,3> dimens = {1, 1, 1};
  angem::Point<3,double> origin = {0,0,0};
  angem::Point<3,double> corner = {1,1,1};
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    if (key == "dimens")
      dimens = it->second.as<std::array<size_t,3>>();
    else if (key == "origin")
      origin = it->second.as<std::vector<double>>();
    else if (key == "corner")
      corner = it->second.as<std::vector<double>>();
  }

  const auto diff = corner - origin;
  for (auto value : dimens)
    if (value == 0) throw std::invalid_argument("Cartesian grid dimension must be > 0");
  conf.dx.assign(dimens[0], diff[0] / dimens[0]);
  conf.dy.assign(dimens[1], diff[1] / dimens[1]);
  conf.dz.assign(dimens[2], diff[2] / dimens[2]);
  conf.origin = origin;
}

void YamlParser::subsection_grid_refinement(const YAML::Node & node)
{
  auto & conf  = config.mesh.refinement;
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\treading key: " << key << std::endl;
    if (key == "type")
    {
      const std::string value = it->second.as<std::string>();
      if (value == "aspect_ratio")
        conf.type = RefinementType::aspect_ratio;
      else throw std::invalid_argument("invalid refinement type");
    }
    else if (key == "aspect_ratio")
      conf.aspect_ratio = it->second.as<double>();
    else if (key == "max_level")
      conf.max_level = it->second.as<size_t>();
    else throw std::invalid_argument("invalid refinement key");
  }

}

void YamlParser::subsection_multiscale(const YAML::Node & node, MultiscaleConfig & conf) {
  conf.part_type = MSPartitioning::no_partitioning;  // default

  for (auto it = node.begin(); it!=node.end(); ++it) {

    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;
    if ( key == "support" )
    {
      const auto value = it->second.as<std::string>();
      if (value == "msrsb") conf.support_type = MSSupportType::msrsb;
      else if (value == "graph") conf.support_type = MSSupportType::graph;
      else if (value == "mechanics") conf.support_type = MSSupportType::mechanics;
    }
    else if ( key == "blocks" )
    {
      try {
        conf.n_blocks = it->second.as<std::array<size_t,3>>();
        conf.part_type = MSPartitioning::geometric;
      }
      catch (YAML::BadConversion const & e)
      {
        conf.n_blocks = {it->second.as<size_t>(), 1, 1};
        conf.part_type = MSPartitioning::metis;
      }
    }
    else
    {
      std::cout << "\tunknown keyword aborting" << std::endl;
      abort();
    }
  }
}

FlowDiscretizationType YamlParser::parse_flow_discretization_(std::string const & value) const
{
  if (value == "tpfa edfm")
    return FlowDiscretizationType::tpfa_edfm;
  else if (value == "tpfa projection")
    return FlowDiscretizationType::tpfa_projection;
  else if (value == "tpfa compartmental")
    return FlowDiscretizationType::tpfa_compartmental;
  else if (value == "insim")
    return FlowDiscretizationType::insim;
  else throw std::invalid_argument("Invalid flow discretization type: " + value);
  return FlowDiscretizationType::tpfa_edfm;  // shut compiler up
}

}  // end namespace
