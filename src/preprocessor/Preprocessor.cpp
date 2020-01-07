#include "Preprocessor.hpp"
#include "parsers/YamlParser.hpp"
#include "parsers/GmshReader.hpp"
#include "CellPropertyManager.hpp"
#include "EmbeddedFractureManager.hpp"
#include "DiscreteFractureManager.hpp"
#include "discretization/DiscretizationTPFA.hpp"
#include "discretization/DiscretizationDFM.hpp"
#include "discretization/DiscretizationEDFM.hpp"
#include <string>

namespace gprs_data {

Preprocessor::Preprocessor(const Path config_file_path)
{
  // read configuration file
  read_config_file_(config_file_path);
  // infer grid file path
  const Path config_dir_path = config_file_path.parent_path();
  const Path grid_file_path = config_dir_path / config.mesh_file;
  read_mesh_file_(grid_file_path);
}

void Preprocessor::run()
{
  DiscreteFractureManager dfm_mgr(config.discrete_fractures, data);
  const size_t n_dfm_faces = dfm_mgr.count_dfm_faces();
  const size_t n_cells = data.grid.n_cells();

  /* Split cells due to edfm intersection */
  EmbeddedFractureManager edfm_mgr(config.embedded_fractures, config.edfm_method, data);
  edfm_mgr.split_cells();
  /* Since we split edfm faces, we pretend that they are dfm fractures
   * to reuse the discretization code. */
  const std::vector<DiscreteFractureConfig> edfm_faces_conf = edfm_mgr.generate_dfm_config();
  // combine dfm and edfm configs
  const std::vector<DiscreteFractureConfig> combined_fracture_config =
      DiscreteFractureManager::combine_configs(config.discrete_fractures,
                                               edfm_faces_conf);
  // manages dofs of dfm and edfm after splitting
  DiscreteFractureManager fracture_flow_mgr(combined_fracture_config, data);
  fracture_flow_mgr.distribute_properties();

  // property manager for grid with split cells (due to edfm splitting)
  CellPropertyManager property_mgr(config.cell_properties, config.domains, data);
  property_mgr.generate_properties();

  // flow dof numbering
  fracture_flow_mgr.build_flow_cv_numbering();

  // edfm + dfm discretization
  // we will use only the edfm part from it
  discretization::DiscretizationDFM discr_edfm_dfm(combined_fracture_config, data);
  discr_edfm_dfm.build();

  // build edfm discretization from mixed dfm-edfm discretization
  discretization::DiscretizationEDFM discr_edfm(combined_fracture_config,
                                                discr_edfm_dfm.get_cell_data(),
                                                discr_edfm_dfm.get_face_data(),
                                                n_dfm_faces, data);
  discr_edfm.build();
  // now coarsen all cells to remove edfm splits and compute normal
  // flow dfm-matrix discretizations
  data.grid.coarsen_cells();
  property_mgr.generate_properties();
  discretization::DiscretizationTPFA matrix_discr(config.discrete_fractures, data);
  dfm_mgr.distribute_properties();
  discretization::DiscretizationDFM discr_dfm(config.discrete_fractures, data);
  discr_dfm.build();

  // merge dfm and matrix discretizations
  discr_dfm.merge_matrix_discretization(matrix_discr.get_cell_data(),
                                        matrix_discr.get_face_data());

  // finally, merge edfm, dfm, and matrix discretizations

  // TODO: write code for combining flow data
  assert( false && "Write code for combining flow data" );
}

void Preprocessor::read_config_file_(const Path config_file_path)
{
  if (!filesystem::exists(config_file_path))
  {
    const std::string error_msg =
        "config file does not exist: " +
        std::string(filesystem::absolute(config_file_path));
    throw std::invalid_argument(error_msg);
  }

  std::cout << "reading ";
  std::cout << filesystem::absolute(config_file_path) << std::endl;

  const std::string fname = config_file_path.filename();
  const std::size_t str_len = fname.size();
  if (fname.substr(str_len - 4, str_len) == "yaml")
  {
    Parsers::YamlParser parser;
    parser.parse_file(fname);
    config = parser.get_config();
  }
  else
  {
    std::cout << "Only .yaml configuration files are supported" << std::endl;
    throw std::invalid_argument("File type not supported");
  }
}

void Preprocessor::read_mesh_file_(const Path mesh_file_path)
{
  if (!filesystem::exists(mesh_file_path))
  {
    const std::string msg = "grid file does not exist:" +
                            std::string(filesystem::absolute(mesh_file_path));
    throw std::invalid_argument(msg);
  }

  std::cout << "reading ";
  std::cout << filesystem::absolute(mesh_file_path) << std::endl;

  // check filetype
  const std::string fname = mesh_file_path.filename();
  const std::size_t str_len = fname.size();

  if (fname.substr(str_len - 3, str_len) != "msh")
    throw std::invalid_argument("Only .msh files produced by Gmsh are supported");

  Parsers::GmshReader::read_input(filesystem::absolute(mesh_file_path), data.grid);
}

}  // end namespace gprs_data
