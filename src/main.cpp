// #include <preprocessor/simdata.hpp>
// #include <preprocessor/OutputDataGPRS.hpp>
// #include "preprocessor/OutputDataVTK.hpp"
// #include <parsers/JsonParser.hpp>
// #include <parsers/YamlParser.hpp>
// #include <uint256/uint256_t.h>
// #include <parsers/GmshReader.hpp>
// #include <mesh/Mesh.hpp>

#include "preprocessor/Preprocessor.hpp"
#include <iostream>
#include <string>
#include <experimental/filesystem>

namespace filesystem = std::experimental::filesystem;
using Path = filesystem::path;

int main(int argc, char *argv[])
{
  // process cmd arguments
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
    std::cout << "Please provide only a single input argument" << std::endl;
    return 1;
  }

  // config file
  const std::string fname_config = argv[1];
  const Path path_config(fname_config);
  try {
    // read stuff
    gprs_data::Preprocessor preprocessor(path_config);
    // let the fun begin
    preprocessor.run();
  }
  catch (const std::exception& error)
  {
    std::cout << "Error: " << error.what() << std::endl;
    return 1;
  }

  return 0;

  // // MESH
  // // get path of the config file -- mesh is searched for in relative path
  // const Path config_dir_path = path_config.parent_path();
  // Path path_gmsh = config_dir_path / config.mesh_file;

  // // check if mesh file exists
  // if (!filesystem::exists(path_gmsh))
  // {
  //   std::cout << "msh file does not exist:" << std::endl;
  //   std::cout << filesystem::absolute(path_gmsh) << std::endl;
  //   return 0;
  // }

  // std::cout << "reading ";
  // std::cout << filesystem::absolute(path_gmsh) << std::endl;
  // mesh::Mesh msh;
  // try
  // {
  //   Parsers::GmshReader::read_input(filesystem::absolute(path_gmsh), msh);
  // }
  // catch (std::exception & e)
  // {
  //   std::cout << "error while reading gmsh file:" << std::endl;
  //   std::cout << e.what() << std::endl;
  //   exit(1);
  // }

  // if (msh.n_cells() == 0)
  // {
  //   std::cout << "mesh has not cells. aborting" << std::endl;
  //   return 0;
  // }

  // // do preprocessing
  // gprs_data::SimData preprocessor = gprs_data::SimData(msh, config);

  // cout << "Fill 3D rock properties" << endl;
  // preprocessor.defineRockProperties();

  // cout << "Make SDA properties" << endl;
  // preprocessor.defineEmbeddedFractureProperties();

  // cout << "Create physical facets" << endl;
  // preprocessor.definePhysicalFacets();

  // std::cout << "computing reservoir transes" << std::endl;
  // preprocessor.computeReservoirTransmissibilities();

  // std::cout << "Handle flow embedded fractures" << std::endl;
  // preprocessor.handleEmbeddedFractures();

  // // timur's legacy
  // // cout << "Create simple wells" << endl;
  // // pSimData->createSimpleWells();

  // std::cout << "Setup wells" << std::endl;
  // preprocessor.setupWells();

  // if (!preprocessor.dfm_faces.empty())
  // {
  //   cout << "Split FEM mesh on internal surfaces" << endl;
  //   preprocessor.splitInternalFaces();
  // }

  // std::cout << "compute connections between mech and flow elements" << std::endl;
  // preprocessor.handleConnections();

  // // multiscale
  // std::cout << "build multiscale data" << std::endl;
  // preprocessor.build_multiscale_data();

  // const std::string output_dir = std::string(filesystem::absolute(config_dir_path)) + "/";
  // std::cout << "output directory: " << output_dir << std::endl;
  // // if no frac remove vtk files
  // if (preprocessor.vEfrac.empty())
  // {
  //   const std::string efrac_vtk_file = output_dir + "efrac.vtk";
  //   if (filesystem::exists(efrac_vtk_file))
  //   {
  //     std::cout << "cleanup old efrac.vtk file" << std::endl;
  //     filesystem::remove(efrac_vtk_file);
  //   }
  // }
  // if (preprocessor.n_flow_dfm_faces == 0)
  // {
  //   const std::string dfm_vtk_file = output_dir + "dfm.vtk";
  //   if (filesystem::exists(dfm_vtk_file))
  //   {
  //     std::cout << "cleanup old dfm.vtk file" << std::endl;
  //     filesystem::remove(dfm_vtk_file);
  //   }
  // }

  // // OUPUT
  // cout << "Write Output data\n";
  // for (auto format : config.output_formats)
  // {
  //   switch (format) {
  //     case OutputFormat::gprs :
  //       {
  //         std::cout << "Output gprs format" << std::endl;
  //         gprs_data::OutputDataGPRS output_data(preprocessor, msh);
  //         output_data.write_output(output_dir);
  //         break;
  //       }
  //       case OutputFormat::vtk :
  //         {
  //           std::cout << "Output vtk format" << std::endl;
  //           gprs_data::OutputDataVTK output_data(preprocessor, msh);
  //           output_data.write_output(output_dir);
  //           break;
  //         }
  //   }
  // }

}
