#include "simdata.hpp"
#include "femout.hpp"
#include "transes.hpp"
#include <Parser.hpp>

class SimData;
class CalcTranses;
class OutputData;
class tetgenio;
class tetrahedralize;

int main(int argc, char *argv[])
{
  string instream = argv[1];
  std::cout << "parsing " << instream << std::endl;

  Parsers::Parser parser;
  parser.parse_file(instream);


  // string outstream;
  // SimData * pSimData;
  // pSimData = new SimData(instream);

  // cout << "Read setup data" << endl;
  // pSimData->readSetupValues();

  // cout << "Reserve boundary conditions" << endl;
  // pSimData->initilizeBoundaryConditions();

  // cout << "Read gmsh data" << endl;
  // pSimData->readGmshFile();

  // cout << "Extract all polygons (slow)" << endl;
  // pSimData->extractInternalFaces();

  // cout << "Convert GMSH FEM mesh into SIM data" << endl;
  // pSimData->convertGmsh2Sim();

  // cout << "Fill 3D rock properties" << endl;
  // pSimData->defineRockProperties();

  // cout << "Make SDA properties" << endl;
  // pSimData->defineEmbeddedFractureProperties();

  // cout << "Create physical facets" << endl;
  // cout << " ( bnd & frac faces )" << endl;
  // pSimData->definePhysicalFacets();

  // cout << "Create bc stress & disp" << endl;
  // pSimData->defineStressAndDispOnBoundary();

  // cout << endl << "Convert FEM mesh into FVM mesh" << endl;
  // pSimData->handleConnections();
  // CalcTranses * pTranses;
  // pTranses = new CalcTranses(pSimData);
  // pTranses->createKarimiData();
  // cout << "Extract  transes from FVM mesh" << endl;
  // pTranses->createKarimiApproximation();

  // cout << "Create simple wells" << endl;
  // pSimData->createSimpleWells();

  // // cout << "Split FEM mesh on internal surfaces" << endl;
  // // pSimData->splitInternalFaces();

  // cout << "Write FEM mesh data\n";
  // OutputData * pOut;
  // pOut = new OutputData(pSimData);

  // pOut->writeGeomechDataNewKeywords();
  // pTranses->outputKarimi();
}
