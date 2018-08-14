#include <Parser.hpp>

#include <string>
#include <fstream>  // std::ifstream
#include <iostream>  // debug
#include <Rectangle.hpp>

namespace Parsers
{

Parser::Parser()
{}


void
Parser::parse_file(const std::string & fname)
{
  std::size_t str_len = fname.size();
  if (fname.substr(str_len - 4, str_len) == "json")
    parse_json(fname);
  else
  {
    std::cout << "file type not supported" << std::endl;
    abort();
  }
}


SimdataConfig & Parser::get_config()
{
  return config;
}


void
Parser::parse_json(const std::string & fname)
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

    if (section_it.key() == "Domain Properties")
      domain_props_json(section_it);
    else if (section_it.key() == "Embedded Fractures")
      embedded_fracs_json(section_it);
    else if (section_it.key() == "Boundary conditions")
      boundary_conditions_json(section_it);
    else if (section_it.key() == "Mesh file")
      config.mesh_file = (*section_it).get<std::string>();
    else
      std::cout << "Skipping section " << section_it.key() << std::endl;
  }  // end section loop
}


void
Parser::domain_props_json(const nlohmann::json::iterator & section_it)
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
      config.domain_file = (*domain_it).get<std::string>();
      continue;
    }

    // define local config
    config.domains.emplace_back();
    auto & conf = config.domains.back();

    // set props
    conf.label             = std::atoi(domain_it.key().c_str());
    std::cout << "Reading domain " << conf.label << std::endl;

    nlohmann::json::iterator
        attrib_it = (*domain_it).begin(),
        attrib_end = (*domain_it).end();

    std::size_t exp_counter = 0;
    for (;attrib_it != attrib_end; ++attrib_it)
    {
      const std::pair<std::string, std::string> key_value =
          get_pair_json((*attrib_it).begin());

      if (key_value.first == comment)
        continue;
      else
      {
        std::cout << "reading expression " << key_value.first << "\t";
        // save property name and expression
        conf.variables.push_back(key_value.first);
        conf.expressions.push_back(key_value.second);
        // check if in global list and append if so
        const std::size_t ind = find(key_value.first, config.all_vars);
        if (ind == config.all_vars.size())
          config.all_vars.push_back(key_value.first);
        // save positions in the global list
        conf.local_to_global_vars[exp_counter] = ind;
        conf.global_to_local_vars[ind] = exp_counter;
        exp_counter++;
      }

    }
 }
}


std::pair<std::string,std::string>
Parser::get_pair_json(const nlohmann::json::iterator & pair_it)
{
  return std::make_pair(pair_it.key(), (*pair_it).get<std::string>());
}


void
Parser::embedded_fracs_json(const nlohmann::json::iterator & section_it)
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
Parser::boundary_conditions_json(const nlohmann::json::iterator & it)
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
Parser::boundary_conditions_faces(nlohmann::json::iterator it,
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
    boundary_conditions_face((*it).begin(), (*it).end(), conf);
  }
}


void
Parser::boundary_conditions_face(nlohmann::json::iterator         it,
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
Parser::boundary_conditions_nodes(nlohmann::json::iterator         it,
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
Parser::boundary_conditions_node(nlohmann::json::iterator it,
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
      std::cout << "attribute " << it.key()
                << " unknown: skipping" << std::endl;
  }
}

void Parser::embedded_fracture(nlohmann::json::iterator it,
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
  double friction_angle = 30;
  double dilation_angle = 0;

  for (; it != end; ++it)
  {
    const auto key = it.key();
    std::cout << "parsing entry " << key << std::endl;
    if (key == comment)
      continue;
    else if (key == "type")
    {
      assert (shape == (*it).get<std::string>());
      shape = (*it).get<std::string>();
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
    else if (key == "center")
      center = (*it).get<std::vector<double>>();
    else
      std::cout << "attribute " << key
                << " unknown: skipping" << std::endl;
  }

  conf.body = std::make_shared<angem::Rectangle<double>>
      (center, length, height, dip, strike);

  conf.cohesion = cohesion;
  conf.friction_angle = friction_angle;
  conf.dilation_angle = dilation_angle;
}

}  // end namespace
