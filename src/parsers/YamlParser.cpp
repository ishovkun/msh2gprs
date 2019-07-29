#include <YamlParser.hpp>
#include "yaml-cpp/node/parse.h"  // loadfile
#include <angem/Rectangle.hpp>


namespace Parsers
{

YamlParser::YamlParser()
{}


void YamlParser::parse_file(const std::string & fname)
{
  YAML::Node main_node = YAML::LoadFile(fname);

  for(YAML::const_iterator it=main_node.begin();it != main_node.end();++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "Reading section: " << key  << std::endl;

    if (key == "Mesh file")
      config.mesh_file = it->second.as<std::string>();
    else if (key == "Domain Flow Properties")
      section_domain_props(it->second, 0);
    else if (key == "Domain Mechanical Properties")
      section_domain_props(it->second, 1);
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
      std::cout << "Unknown key: " << key << " skipping" << std::endl;
  }

}


void YamlParser::embedded_fracs(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading entry " << key << std::endl;

    if (key == "file")
      config.efrac_file = it->second.as<std::string>();
    else if (key == "method")
    {
      const std::string str_method = it->second.as<std::string>();
      if (str_method == "simple")
        config.edfm_method = EDFMMethod::simple;
      else if (str_method == "projection")
        config.edfm_method = EDFMMethod::projection;
      else
        throw std::invalid_argument(str_method);
    }
    else if (key == "fracture")
    {
      config.fractures.emplace_back();
      embedded_fracture(it->second, config.fractures.back());
    }
    else
      std::cout << "\t\tattribute " << key << " unknown: skipping" << std::endl;
  }
}


void YamlParser::discrete_fracs(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading entry " << key << std::endl;

    if (key == "fracture")
    {
      config.discrete_fractures.emplace_back();
      // check if frac has a label
      auto & conf = config.discrete_fractures.back();
      extract_subnode_value("label", it, conf.label);
      if (conf.label == 0)
      {
        std::cout << "cannot have label 0" << std::endl;
        exit(-1);
      }
      discrete_fracture(it->second, conf);
    }
    else if (key == "file")
      config.discrete_frac_file = it->second.as<std::string>();
    else
      std::cout << "\t\tattribute " << key << " unknown: skipping" << std::endl;
  }

}


void YamlParser::embedded_fracture(const YAML::Node       & node,
                                   EmbeddedFractureConfig & conf)
{
  // default parameters
  std::string shape = "Rectangle";
  double height = 1;  // for rectangles and ellipses
  double length = 1;  // for rectangles and ellipses
  double dip = 0;
  double strike = 0;
  double cohesion = 0;
  angem::Point<3,double> center = {0, 0, 0};
  double friction_angle = conf.friction_angle;
  double dilation_angle = conf.dilation_angle;
  double aperture = conf.aperture;
  double conductivity = conf.conductivity;
  std::size_t n1 = conf.n1;
  std::size_t n2 = conf.n2;

  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\treading entry " << key << std::endl;

    if (key == "type")
    {
      shape = it->second.as<std::string>();
      assert (shape == "Rectangle");
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
    else if (key == "remesh")
    {
      std::vector<std::size_t> remesh_pars =
          it->second.as<std::vector<std::size_t>>();

      if (remesh_pars.size() != 2)
        std::cout << "wrong entry in remesh. Aborting" << std::endl << std::flush;

      n1 = remesh_pars[0];
      n2 = remesh_pars[1];

      if ((n1 > 0 and n2 == 0)  or (n2 > 0 and n1 == 0))
        std::cout << "wrong entry in remesh. Aborting" << std::endl << std::flush;
    }
    else if (key == "center")
      center = it->second.as<std::vector<double>>();
    else
      std::cout << "\t\tattribute " << key
                << " unknown: skipping" << std::endl;
  }

  std::cout << "Making embedded fracture" << std::endl;
  conf.body = std::make_shared<angem::Rectangle<double>>
      (center, length, height, dip, strike);

  conf.cohesion = cohesion;
  conf.friction_angle = friction_angle;
  conf.dilation_angle = dilation_angle;
  conf.aperture = aperture;
  conf.conductivity = conductivity;
  conf.n1 = n1; conf.n2 = n2;
}


void YamlParser::discrete_fracture(const YAML::Node       & node,
                                   DiscreteFractureConfig & conf)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\treading entry " << key << std::endl;

    if (key == "aperture")
      conf.aperture = it->second.as<double>();
    else if (key == "conductivity")
      conf.conductivity = it->second.as<double>();
    else if (key == "label")
      conf.label = it->second.as<int>();
  }
}


void YamlParser::section_domain_props(const YAML::Node & node,
                                      const int          var_type)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    // std::cout << "\treading entry " << key << std::endl;

    if (key == "file")  // where to ouput properties
    {
      if (var_type == 0)
        config.domain_file = it->second.as<std::string>();
      else if (var_type == 1)
        config.mechanics_domain_file = it->second.as<std::string>();
      continue;
    }
    else if (key == "domain")
    {
      int label;
      try {
        label = it->second["label"].as<int>();
        if (label == 0)
        {
          std::cout << "domain label cannot be 0" << std::endl;
          exit(-1);
        }
      }
      catch (YAML::TypedBadConversion<int> & error)
      {
        std::cout << "domain label must be an integer" << std::endl;
        std::cout << "aborting" << std::endl;
        exit(-1);
      }
      catch (YAML::InvalidNode & error)
      {
        std::cout << "domain label must be specified" << std::endl;
        exit(-1);
      }

      DomainConfig & conf = get_domain_config(label);
      std::cout << "\tReading domain " << conf.label << std::endl;
      domain(it->second, var_type, conf);
    }
    else
      std::cout << "attribute " << key << " unknown: skipping" << std::endl;
  }

}


void YamlParser::domain(const YAML::Node & node,
                        const int          var_type,
                        DomainConfig     & conf)
{
  std::size_t exp_counter = conf.expressions.size();

  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    const std::string value = it->second.as<std::string>();

    std::cout << "\t\treading entry " << key << std::endl;

    if (key == "Coupled") // coupling with mechanics
    {
      conf.coupled = it->second.as<bool>();
      if (!conf.coupled)
        std::cout << "domain " <<conf.label << " is decoupled !!!" << std::endl;
    }
    else if (key == "label") // coupling with mechanics
      continue;
    else
    {
      // save property name and expression
      conf.variables.push_back(key);
      conf.expressions.push_back(value);
      const std::size_t ind = find(key, config.all_vars);
      if (ind == config.all_vars.size())  // new variable
      {
        config.all_vars.push_back(key);
        // special case - service variable (not outputted)
        if ( find(key, config.special_keywords) <
             config.special_keywords.size())
          config.expression_type.push_back(-1);
        else
          config.expression_type.push_back(var_type);
      }

      // save positions in the global list
      conf.local_to_global_vars[exp_counter] = ind;
      conf.global_to_local_vars[ind]         = exp_counter;
      exp_counter++;
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
      config.bcond_file = it->second.as<std::string>();
    else if (key == "Faces")
      boundary_conditions_faces(it->second);
    else if (key == "Dirichlet nodes")
      boundary_conditions_nodes(it->second);
    else
      std::cout << "skipping entry " << key << std::endl;
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
      std::cout << "\t\tattribute " << key << " unknown: skipping" << std::endl;
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
      conf.type = it->second.as<int>();
    else if (key == "label")
      conf.label = it->second.as<int>();
    else if (key == "value")
    {
      const std::vector<std::string> str_values =
          it->second.as<std::vector<std::string>>();

      for (std::size_t i=0; i<3; ++i)
      {
        if (str_values[i] == "nan")
          conf.value[i] = config.nan;
        else
          conf.value[i] = std::atof(str_values[i].c_str());
      }
    }
    else
      std::cout << "attribute " << key << " unknown: skipping" << std::endl;
  }
}


void YamlParser::boundary_conditions_nodes(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\t\treading entry " << key << std::endl;

    if (key == "search tolerance")
    {
      config.node_search_tolerance = it->second.as<double>();
      continue;
    }
    else if (key == "node")
    {
      config.bc_nodes.emplace_back();
      auto & conf = config.bc_nodes.back();
      bc_node(it->second, conf);
    }
    else
      std::cout << "attribute " << key << " unknown: skipping" << std::endl;
  }
}


void YamlParser::bc_node(const YAML::Node & node, BCNodeConfig & conf)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\t\t\treading entry " << key << std::endl;

    if (key == "coord")
    {
      const std::vector<double> v_coord = it->second.as<std::vector<double>>();
      conf.coord = v_coord;
    }
    else if (key == "value")
    {
      const std::vector<std::string> str_values =
          it->second.as<std::vector<std::string>>();
      for (std::size_t i=0; i<3; ++i)
      {
        if (str_values[i] == "nan")
          conf.value[i] = config.nan;
        else
          conf.value[i] = std::atof(str_values[i].c_str());
      }
    }
    else
      std::cout << "attribute " << key << " unknown: skipping" << std::endl;
  }
}


DomainConfig & YamlParser::get_domain_config(const int label)
{
  // find current domain id in existing config
  int counter = 0;
  for (const auto & domain : config.domains)
  {
    if (label == domain.label)
      break;
    else
      counter++;
  }

  DomainConfig * p_conf;
  if (counter == config.domains.size())
  {
    config.domains.emplace_back();
    config.domains.back().label = label;
    return config.domains.back();
  }
  else
    return config.domains[counter];
}


void YamlParser::section_wells(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;

    if (key == "file")
      config.wells_file = it->second.as<std::string>();
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

      // for (std::size_t i=0; i<n_coord; ++i)
      //   std::cout << well.coordinates[i] << std::endl;
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
      well.perforated.resize(well.coordinates.size(), true);
    else if (well.perforated.size() != well.coordinates.size() - 1)
    {
      std::cout << "invalid entry" << std::endl;
      abort();
    }
  }
}


void YamlParser::section_multiscale(const YAML::Node & node)
{
  for (auto it = node.begin(); it!=node.end(); ++it)
  {
    const std::string key = it->first.as<std::string>();
    std::cout << "\treading key " << key << std::endl;

    if (key == "Flow file")
      config.flow_ms_file = it->second.as<std::string>();
    if (key == "Mech file")
      config.mech_ms_file = it->second.as<std::string>();
    else if (key == "blocks")
      config.n_multiscale_blocks = it->second.as<std::size_t>();
    else if (key == "flow")
    {
      const auto value = it->second.as<std::string>();
      if (value == "no")
        config.multiscale_flow = MSPartitioning::no_partitioning;
      else if (value == "msrsb")
        config.multiscale_flow = MSPartitioning::method_msrsb;
      else if (value == "mrst")
        config.multiscale_flow = MSPartitioning::method_mrst_flow;
      else
      {
        std::cout << "\tunknown keyword aborting" << std::endl;
        abort();
      }
    }
    else if (key == "mechanics")
    {
      const auto value = it->second.as<std::string>();
      if (value == "no")
        config.multiscale_mechanics = MSPartitioning::no_partitioning;
      // else if (value == "msrsb")
      //   config.multiscale_flow = MSPartitioning::method_msrsb;
      else if (value == "srfem")
        config.multiscale_mechanics = MSPartitioning::method_mechanics;
      else
      {
        std::cout << "\tunknown keyword aborting" << std::endl;
        abort();
      }

    }
    else
      std::cout << "\tSkipping unknown keyword" << std::endl;
  }
}

}  // end namespace
