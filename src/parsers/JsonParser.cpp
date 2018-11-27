#include <JsonParser.hpp>

#include <string>
#include <fstream>  // std::ifstream
#include <iostream>  // debug
#include <angem/Rectangle.hpp>


namespace Parsers
{

JsonParser::JsonParser()
{}


void
JsonParser::parse_file(const std::string & fname)
{
  const std::size_t str_len = fname.size();
  if (fname.substr(str_len - 4, str_len) == "json")
    parse(fname);
  else
  {
    std::cout << "file type not supported" << std::endl;
    abort();
  }
}


SimdataConfig & JsonParser::get_config()
{
  return config;
}


void
JsonParser::parse(const std::string & fname)
{
  // read a JSON file
  std::ifstream input(fname);
  nlohmann::json jparser;
  input >> jparser;

  // iterate sections
  for (nlohmann::json::iterator section_it = jparser.begin();
       section_it != jparser.end(); ++section_it)
  {
    std::cout << "Entering section " << section_it.key() << '\n';

    if (section_it.key() == "Domain Flow Properties")
      domain_props(section_it, 0);
    else if (section_it.key() == "Domain Mechanics Properties")
      domain_props(section_it, 1);
    else if (section_it.key() == "Embedded Fractures")
      embedded_fracs(section_it);
    else if (section_it.key() == "Discrete Fractures")
      discrete_fracs(section_it);
    else if (section_it.key() == "Boundary conditions")
      boundary_conditions(section_it);
    else if (section_it.key() == "Mesh file")
      config.mesh_file = (*section_it).get<std::string>();
    else
      std::cout << "Skipping section " << section_it.key() << std::endl;
  }  // end section loop
}


void
JsonParser::domain_props(const nlohmann::json::iterator & section_it,
                         const int                        var_type)
{
  nlohmann::json::iterator
      domain_it = (*section_it).begin(),
      domain_end = (*section_it).end();

  // iterate domains
  for (;domain_it != domain_end; ++domain_it)
  {
    if (domain_it.key() == comment)
      continue;
    if (domain_it.key() == "file")
    {
      if (var_type == 0)
        config.domain_file = (*domain_it).get<std::string>();
      else if (var_type == 1)
        config.mechanics_domain_file = (*domain_it).get<std::string>();
      continue;
    }

    // check whether the label has been used already
    const int label = std::atoi(domain_it.key().c_str());
    // find current domain id in existing config
    int counter = 0;
    for (const auto & domain : config.domains)
    {
      if (label == domain.label)
        break;
      else
        counter++;
    }

    // get pointer to the existing config or create if none
    DomainConfig * p_conf;
    if (counter == config.domains.size())
    {
      config.domains.emplace_back();
      p_conf = &config.domains.back();
    }
    else
      p_conf = &(config.domains[counter]);

    DomainConfig & conf = *p_conf;

    // set props
    conf.label = label;
    std::cout << "Reading domain " << conf.label << std::endl;

    nlohmann::json::iterator
        attrib_it = (*domain_it).begin(),
        attrib_end = (*domain_it).end();

    std::size_t exp_counter = conf.expressions.size();
    for (;attrib_it != attrib_end; ++attrib_it)
    {
      const std::pair<std::string, std::string> key_value =
          get_pair((*attrib_it).begin());

      if (key_value.first == comment)
        continue;
      else if (key_value.first == "Coupled")
      {
        conf.coupled = static_cast<bool>(std::stoi(key_value.second));
        if (!conf.coupled)
          std::cout << "domain " <<conf.label << " is decoupled !!!" << std::endl;
        continue;
      }
      else
      {
        std::cout << "\treading expression " << key_value.first << std::endl;
        // save property name and expression
        conf.variables.push_back(key_value.first);
        conf.expressions.push_back(key_value.second);
        // check if in global list and append if so
        const std::size_t ind = find(key_value.first, config.all_vars);
        if (ind == config.all_vars.size())  // new variable
        {
          config.all_vars.push_back(key_value.first);
          // special case - service variable (not outputted)
          if ( find(key_value.first, config.special_keywords) <
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
}


std::pair<std::string,std::string>
JsonParser::get_pair(const nlohmann::json::iterator & pair_it)
{
  return std::make_pair(pair_it.key(), (*pair_it).get<std::string>());
}


void
JsonParser::embedded_fracs(const nlohmann::json::iterator & section_it)
{
  nlohmann::json::iterator
      frac_it = (*section_it).begin(),
      frac_end = (*section_it).end();

  for (;frac_it != frac_end; ++frac_it)
  {
    if (frac_it.key() == comment)
      continue;
    else if (frac_it.key() == "file")
    {
      config.efrac_file = (*frac_it).get<std::string>();
      continue;
    }
    else
    {
      config.fractures.emplace_back();
      embedded_fracture((*frac_it).begin(), (*frac_it).end(),
                        config.fractures.back());
    }
  }

}


void
JsonParser::boundary_conditions(const nlohmann::json::iterator & it)
{
  nlohmann::json::iterator
      bc_it = (*it).begin(),
      bc_end = (*it).end();

  for (;bc_it != bc_end; ++bc_it)
  {
    std::cout << "parsing " << bc_it.key() << std::endl;

    if (bc_it.key() == comment)
      continue;
    else if (bc_it.key() == "file")
    {
      config.bcond_file = (*bc_it).get<std::string>();
      continue;
    }
    else if (bc_it.key() == "Faces")
    {
      boundary_conditions_faces((*bc_it).begin(), (*bc_it).end());
      continue;
    }
    else if (bc_it.key() == "Dirichlet nodes")
    {
      boundary_conditions_nodes((*bc_it).begin(), (*bc_it).end());
      continue;
    }
    else
    {
      std::cout << "skipping " << bc_it.key() << std::endl;
    }

  }
}


void
JsonParser::boundary_conditions_faces(nlohmann::json::iterator it,
                                      const nlohmann::json::iterator & end)
{
  for (; it != end; ++it)
  {
    std::cout << "parsing entry " << it.key() << std::endl;
    if (it.key() == comment)
      continue;

    config.bc_faces.emplace_back();
    auto & conf = config.bc_faces.back();
    conf.label = std::atoi(it.key().c_str());
    // if (conf.label >= 0)
    //   throw std::invalid_argument("boundary labels should be negative");

    boundary_conditions_face((*it).begin(), (*it).end(), conf);
  }
}


void
JsonParser::boundary_conditions_face(nlohmann::json::iterator         it,
                                     const nlohmann::json::iterator & end,
                                     BCConfig & conf)
{
  for (; it != end; ++it)
  {
    const auto key = it.key();
    if (key == comment)
      continue;
    else if (key == "type")
      conf.type = (*it).get<int>();
    else if (key == "value")
    {
      const std::vector<std::string> str_values =
          (*it).get<std::vector<std::string>>();
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


void
JsonParser::boundary_conditions_nodes(nlohmann::json::iterator         it,
                                      const nlohmann::json::iterator & end)
{
  for (; it != end; ++it)
  {
    std::cout << "parsing entry "<< it.key() << std::endl;

    if (it.key() == comment)
      continue;
    else if (it.key() == "search tolerance")
    {
      config.node_search_tolerance = (*it).get<double>();
      continue;
    }

    config.bc_nodes.emplace_back();
    auto & conf = config.bc_nodes.back();
    boundary_conditions_node((*it).begin(), (*it).end(), conf);
  }
}


void
JsonParser::boundary_conditions_node(nlohmann::json::iterator it,
                                     const nlohmann::json::iterator & end,
                                     BCNodeConfig & conf)
{
  for (; it != end; ++it)
  {
    if (it.key() == comment)
      continue;
    else if (it.key() == "coord")
    {
      const std::vector<double> v_coord =
          (*it).get<std::vector<double>>();
      conf.coord = v_coord;
    }
    else if (it.key() == "value")
    {
      const std::vector<std::string> str_values =
          (*it).get<std::vector<std::string>>();
      for (std::size_t i=0; i<3; ++i)
      {
        if (str_values[i] == "nan")
          conf.value[i] = config.nan;
        else
          conf.value[i] = std::atof(str_values[i].c_str());
      }
    }
    else
      std::cout << "attribute " << it.key() << " unknown: skipping" << std::endl;
  }
}

void JsonParser::embedded_fracture(nlohmann::json::iterator it,
                               const nlohmann::json::iterator & end,
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

  for (; it != end; ++it)
  {
    const auto key = it.key();
    std::cout << "\tparsing entry " << key << std::endl;
    if (key == comment)
      continue;
    else if (key == "type")
    {
      shape = (*it).get<std::string>();
      assert (shape == "Rectangle");
    }
    else if (key == "height")
      height = (*it).get<double>();
    else if (key == "length")
      length = (*it).get<double>();
    else if (key == "strike")
      strike = (*it).get<double>();
    else if (key == "dip")
      dip = (*it).get<double>();
    else if (key == "cohesion")
      cohesion = (*it).get<double>();
    else if (key == "friction angle")
      friction_angle = (*it).get<double>();
    else if (key == "dilation angle")
      dilation_angle = (*it).get<double>();
    else if (key == "aperture")
      aperture = (*it).get<double>();
    else if (key == "conductivity")
      conductivity = (*it).get<double>();
    else if (key == "remesh")
    {
      std::vector<std::size_t> remesh_pars = (*it).get<std::vector<std::size_t>>();
      if (remesh_pars.size() != 2)
        std::cout << "wrong entry in remesh. Aborting" << std::endl << std::flush;

      n1 = remesh_pars[0];
      n2 = remesh_pars[1];

      if ((n1 > 0 and n2 == 0)  or (n2 > 0 and n1 == 0))
        std::cout << "wrong entry in remesh. Aborting" << std::endl << std::flush;
    }
    else if (key == "center")
      center = (*it).get<std::vector<double>>();
    else
      std::cout << "attribute " << key
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


void
JsonParser::discrete_fracs(const nlohmann::json::iterator & section_it)
{
  nlohmann::json::iterator
      frac_it = (*section_it).begin(),
      frac_end = (*section_it).end();

  for (;frac_it != frac_end; ++frac_it)
  {
    if (frac_it.key() == comment)
      continue;
    else
    {
      const int label = std::atoi(frac_it.key().c_str());
      if (label == 0)
        continue;
      // if (label < 0)
      //   throw std::invalid_argument("fracture labels should be positive");

      std::cout << "parsing fracture " << label << std::endl;
      config.discrete_fractures.emplace_back();
      config.discrete_fractures.back().label = label;
      discrete_fracture((*frac_it).begin(), (*frac_it).end(),
                        config.discrete_fractures.back());
    }

  }
}


void JsonParser::discrete_fracture(nlohmann::json::iterator it,
                                   const nlohmann::json::iterator & end,
                                   DiscreteFractureConfig & conf)
{
  for (; it != end; ++it)
  {
    const auto key = it.key();
    std::cout << "\tparsing entry " << key << std::endl;
    if (key == comment)
      continue;
    else if (key == "aperture")
      conf.aperture = (*it).get<double>();
    else if (key == "conductivity")
      conf.conductivity = (*it).get<double>();

  }
}

}  // end namespace
