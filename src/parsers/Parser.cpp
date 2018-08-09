#include <Parser.hpp>

#include <string>
#include <fstream>  // std::ifstream
#include <iostream>  // debug

namespace Parsers
{

Parser::Parser()
{}


void
Parser::parse_file(const std::string & fname)
{
  // std::cout << "parser parses " << fname << std::endl;
  std::size_t str_len = fname.size();
  std::cout << fname.substr(str_len-4, str_len) << std::endl;
  if (fname.substr(str_len-4, str_len) == "json")
    parse_json(fname);
  else
  {
    std::cout << "file type not supported" << std::endl;
    abort();
  }
}


void
Parser::parse_json(const std::string & fname)
{
  // read a JSON file
  std::ifstream input(fname);
  nlohmann::json jparser;
  input >> jparser;
}

}  // end namespace
