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
      // std::cout << "   " <<  *((*section_it).begin()) << std::endl;
    if (section_it.key() == "Rock Properties")
      rock_props_json(section_it);
    else if (section_it.key() == "Embedded Fractures")
      embedded_fracs_json(section_it);
    else if (section_it.key() == "Boundary conditions")
      boundary_conditions_json(section_it);
    else
      std::cout << "Skipping section " << section_it.key() << std::endl;
  }  // end section loop
}


void
Parser::rock_props_json(const nlohmann::json::iterator & section_it)
{
  nlohmann::json::iterator
      domain_it = (*section_it).begin(),
      domain_end = (*section_it).end();

  // iterate domains
  for (;domain_it != domain_end; ++domain_it)
  {
    // define local config
    config.domains.emplace_back();
    auto & conf = config.domains.back();

    // set props
    conf.label             = std::atoi(domain_it.key().c_str());
    conf.model             = (*domain_it)["model"].get<int>();
    conf.biot              = (*domain_it)["biot"].get<double>();
    conf.porosity          = (*domain_it)["porosity"].get<double>();
    conf.permeability      = (*domain_it)["permeability"].get<double>();
    conf.density           = (*domain_it)["density"].get<double>();
    conf.young_modulus     = (*domain_it)["young modulus"].get<double>();
    conf.poisson_ratio     = (*domain_it)["poisson ratio"].get<double>();
    conf.thermal_expansion = (*domain_it)["thermal expansion"].get<double>();
    conf.heat_capacity     = (*domain_it)["heat capacity"].get<double>();
    conf.temperature       = (*domain_it)["temperature"].get<double>();
    conf.pressure          = (*domain_it)["pressure"].get<double>();
    conf.ref_temperature   = (*domain_it)["reference temperature"].get<double>();
    conf.ref_pressure      = (*domain_it)["reference pressure"].get<double>();
 }
}

void
Parser::embedded_fracs_json(const nlohmann::json::iterator & section_it)
{
  nlohmann::json::iterator
      frac_it = (*section_it).begin(),
      frac_end = (*section_it).end();

  for (;frac_it != frac_end; ++frac_it)
  {

    config.fractures.emplace_back();
    auto & frac = config.fractures.back();

    if ((*frac_it)["type"] != "Rectangle")
    {
      std::cout << "Shape "<< (*frac_it)["type"]
                <<" not implemented: skipping" << std::endl;
      config.fractures.pop_back();
      continue;
    }

    nlohmann::json::iterator
        attrib_it = (*frac_it).begin(),
        attrib_end = (*frac_it).end();

    double height = 1;  // for rectangles and ellipses
    double length = 1;  // for rectangles and ellipses
    double dip = 0;
    double strike = 0;
    angem::Point<3,double> center;
    for (;attrib_it != attrib_end; ++attrib_it)
    {
      const auto key = attrib_it.key();
      if (key == "height")
        height = (*attrib_it).get<double>();
      else if (key == "length")
          length = (*attrib_it).get<double>();
      else if (key == "dip")
        dip = (*attrib_it).get<double>();
      else if (key == "strike")
        strike = (*attrib_it).get<double>();
      else if (key == "center")
        center = (*attrib_it).get<std::vector<double>>();
      else if (key == "cohesion")
        frac.cohesion = (*attrib_it).get<double>();
      else if (key == "friction angle")
        frac.friction_angle = (*attrib_it).get<double>();
      else if (key == "dilation angle")
        frac.dilation_angle = (*attrib_it).get<double>();
      else if (key == "_comment_")
        continue;
      else
        std::cout << "attribute " << key
                  << " unknown: skipping" << std::endl;
    }

    frac.body = angem::Rectangle<double>(center, length,
                                         height, dip, strike);
  }

}


void
Parser::boundary_conditions_json(const nlohmann::json::iterator & section_it)
{
  nlohmann::json::iterator
      bc_it = (*section_it).begin(),
      bc_end = (*section_it).end();

  for (;bc_it != bc_end; ++bc_it)
  {
    config.bcs.emplace_back();
    auto & conf = config.bcs.back();

    // std::cout << *bc_it << std::endl;
    conf.label = std::atoi(bc_it.key().c_str());

    nlohmann::json::iterator
        attrib_it = (*bc_it).begin(),
        attrib_end = (*bc_it).end();

    for (;attrib_it != attrib_end; ++attrib_it)
    {
      const auto key = attrib_it.key();
      if (key == "type")
        conf.type  = (*attrib_it).get<int>();
      else if (key == "value")
      {
        const std::vector<std::string> str_values =
            (*attrib_it).get<std::vector<std::string>>();
        for (std::size_t i=0; i<3; ++i)
        {
          if (str_values[i] == "nan")
            conf.value[i] = conf.nan;
          else
            conf.value[i] = std::atof(str_values[i].c_str());
        }

      }
      else if (key == "_comment_")
        continue;
      else
        std::cout << "attribute " << key << " unknown: skipping" << std::endl;
    }
  }
}

}  // end namespace
