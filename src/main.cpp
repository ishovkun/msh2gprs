#include "simdata.hpp"
#include "femout.hpp"
#include "transes.hpp"
#include <Parser.hpp>

#include <experimental/filesystem>
namespace filesystem = std::experimental::filesystem;

class SimData;
class CalcTranses;
class OutputData;
class tetgenio;
class tetrahedralize;

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

  Parsers::Parser parser;
  std::cout << "parsing " << fname_config << std::endl;
  parser.parse_file(fname_config);
  const SimdataConfig config = parser.get_config();

  // get path of the config file -- mesh is searched for in relative path
  const filesystem::path path_config(fname_config);
  const std::string config_dir = path_config.parent_path().string() + "/";
  const std::string fname_gmsh = config_dir + config.mesh_file;  // msh file

  std::string outstream;
  SimData * pSimData;
  pSimData = new SimData(fname_gmsh, config);

  cout << "Read gmsh data" << endl;
  pSimData->readGmshFile();

  cout << "Extract all polygons (slow)" << endl;
  pSimData->extractInternalFaces();

  cout << "Convert GMSH FEM mesh into SIM data" << endl;
  pSimData->convertGmsh2Sim();

  cout << "Fill 3D rock properties" << endl;
  pSimData->defineRockProperties();

  cout << "Make SDA properties" << endl;
  pSimData->defineEmbeddedFractureProperties();

  cout << "Compute cell clipping" << endl;
  pSimData->computeCellClipping();

  cout << "Create physical facets" << endl;
  cout << " ( bnd & frac faces )" << endl;
  pSimData->definePhysicalFacets();

  cout << "Create bc stress & disp" << endl;
  pSimData->defineStressAndDispOnBoundary();

  cout << endl << "Convert FEM mesh into FVM mesh" << endl;
  pSimData->handleConnections();
  CalcTranses * pTranses;
  pTranses = new CalcTranses(pSimData);
  pTranses->createKarimiData();
  cout << "Extract  transes from FVM mesh" << endl;
  pTranses->createKarimiApproximation();

  cout << "Create simple wells" << endl;
  pSimData->createSimpleWells();

  cout << "Split FEM mesh on internal surfaces" << endl;
  pSimData->splitInternalFaces();

  cout << "Write FEM mesh data\n";
  OutputData * pOut;
  pOut = new OutputData(pSimData);

  const auto output_path = path_config.parent_path().string() + "/";
  pOut->writeGeomechDataNewKeywords(output_path);
  pTranses->outputKarimi(output_path);
}
