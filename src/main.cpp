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
  if (argc < 3)
  {
    std::cout << "please specify msh and config files."
              << std::endl
              << "Example: "
              << "msh2gprs domain.msh config.json"
              << std::endl
              << "for more information type: msh2gprs --help"
              << std::endl;
    return 0;
  }
  if (argc > 3)
  {
    std::cout << "what the hell did you pass?" << std::endl;
    return 1;
  }

  std::string fname_gmsh = argv[1];  // msh file
  std::string fname_config = argv[2];  // config file

  Parsers::Parser parser;
  std::cout << "parsing " << fname_config << std::endl;
  parser.parse_file(fname_config);
  SimdataConfig config = parser.get_config();

  std::string outstream;
  SimData * pSimData;
  pSimData = new SimData(fname_gmsh, config);

  cout << "Read setup data" << endl;
  pSimData->readSetupValues();

  // cout << "Reserve boundary conditions" << endl;
  // pSimData->initilizeBoundaryConditions();

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

  // // cout << "Split FEM mesh on internal surfaces" << endl;
  // // pSimData->splitInternalFaces();

  cout << "Write FEM mesh data\n";
  OutputData * pOut;
  filesystem::path path_config(fname_config);
  pOut = new OutputData(pSimData);

  std::cout << path_config.parent_path().string() << std::endl;
  pOut->writeGeomechDataNewKeywords(path_config.parent_path().string() + "/");
  // pTranses->outputKarimi();
}
