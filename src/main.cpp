#include <gprs-data/simdata.hpp>
#include <gprs-data/femout.hpp>
// #include <gprs_data/transes.hpp>
#include <parsers/JsonParser.hpp>
#include <parsers/YamlParser.hpp>
#include <uint256/uint256_t.h>
#include <parsers/GmshReader.hpp>
#include <mesh/Mesh.hpp>

#include <string>
#include <experimental/filesystem>

namespace filesystem = std::experimental::filesystem;
using Path = filesystem::path;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    std::cout << "please specify a config file."
              << std::endl
              << "Example: "
              << "msh2gprs config.json"
              << std::endl
              << "an example config is distributed with this code"
              << std::endl;
    return 0;
  }
  if (argc > 2)
  {
    std::cout << "what the hell did you pass?" << std::endl;
    return 1;
  }

  const std::string fname_config = argv[1];  // config file
  const Path path_config(fname_config);
  if (!filesystem::exists(path_config))
  {
    std::cout << "config file does not exist:" << std::endl;
    std::cout << filesystem::absolute(path_config) << std::endl;
    return 0;
  }

  std::cout << "parsing " << fname_config << std::endl;
  const std::size_t str_len = fname_config.size();
  SimdataConfig config;
  if (fname_config.substr(str_len - 4, str_len) == "json")
  {
    Parsers::JsonParser parser;
    parser.parse_file(fname_config);
    config = parser.get_config();
  }
  else if (fname_config.substr(str_len - 4, str_len) == "yaml")
  {
    std::cout << "warning: yaml parser is still experimental" << std::endl;
    Parsers::YamlParser parser;
    parser.parse_file(fname_config);
    config = parser.get_config();
  }

  // get path of the config file -- mesh is searched for in relative path
  const Path config_dir_path = path_config.parent_path();
  Path path_gmsh = config_dir_path / config.mesh_file;

  // check if mesh file exists
  if (!filesystem::exists(path_gmsh))
  {
    std::cout << "msh file does not exist:" << std::endl;
    std::cout << filesystem::absolute(path_gmsh) << std::endl;
    return 0;
  }

  std::cout << "reading ";
  std::cout << filesystem::absolute(path_gmsh) << std::endl;
  mesh::Mesh msh;
  Parsers::GmshReader::read_input(filesystem::absolute(path_gmsh), msh);

  SimData * pSimData;
  pSimData = new SimData(msh, config);

  // cout << "Read gmsh data" << endl;
  // pSimData->readGmshFile();

  // cout << "Extract all polygons (slow)" << endl;
  // pSimData->extractInternalFaces();

  // cout << "Convert GMSH FEM mesh into SIM data" << endl;
  // pSimData->convertGmsh2Sim();

  cout << "Fill 3D rock properties" << endl;
  pSimData->defineRockProperties();

  cout << "Make SDA properties" << endl;
  pSimData->defineEmbeddedFractureProperties();

  cout << "Create physical facets" << endl;
  cout << " ( bnd & frac faces )" << endl;
  pSimData->definePhysicalFacets();

  // cout << endl << "Convert FEM mesh into FVM mesh" << endl;
  // pSimData->handleConnections();

  std::cout << "computing reservoir transes" << std::endl;
  pSimData->computeReservoirTransmissibilities();

  // cout << "Compute cell clipping and EDFM transmissibilities" << endl;
  // pSimData->computeCellClipping();

  // // cout << "Merge small edfm cells" << endl;
  // // pSimData->mergeSmallFracCells();

  // // std::cout << "mesh fractures" << std::endl;
  // // pSimData->meshFractures();

  // std::cout << "Compute transmissibilities between edfm fracs" << std::endl;
  // pSimData->computeTransBetweenDifferentEfracs();

  // // // cout << "Create simple wells" << endl;
  // // // pSimData->createSimpleWells();

  if (pSimData->dfm_faces.size() > 0)
  {
    cout << "Split FEM mesh on internal surfaces" << endl;
    pSimData->splitInternalFaces();
  }

  // // cout << "Write FEM mesh data\n";
  // // OutputData output_data(pSimData);

  // // const std::string output_dir = std::string(filesystem::absolute(config_dir_path)) + "/";
  // // std::cout << "output directory: " << output_dir << std::endl;
  // // output_data.write_output(output_dir);

  // // // if no frac remove vtk files
  // // if (pSimData->vEfrac.empty())
  // // {
  // //   const std::string efrac_vtk_file = output_dir + "efrac.vtk";
  // //   if (filesystem::exists(efrac_vtk_file))
  // //   {
  // //     std::cout << "cleanup old efrac.vtk file" << std::endl;
  // //     filesystem::remove(efrac_vtk_file);
  // //   }
  // // }
  // // if (pSimData->nDFMFracs == 0)
  // // {
  // //   const std::string dfm_vtk_file = output_dir + "dfm.vtk";
  // //   if (filesystem::exists(dfm_vtk_file))
  // //   {
  // //     std::cout << "cleanup old dfm.vtk file" << std::endl;
  // //     filesystem::remove(dfm_vtk_file);
  // //   }
  // // }

  // delete pSimData;
  // return 0;
}
