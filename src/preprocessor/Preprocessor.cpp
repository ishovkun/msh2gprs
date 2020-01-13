#include "Preprocessor.hpp"
#include "parsers/YamlParser.hpp"
#include "parsers/GmshReader.hpp"
#include "CellPropertyManager.hpp"
#include "EmbeddedFractureManager.hpp"
#include "DiscreteFractureManager.hpp"
#include "BoundaryConditionManager.hpp"
#include "discretization/DiscretizationTPFA.hpp"
#include "discretization/DiscretizationDFM.hpp"
#include "discretization/DiscretizationEDFM.hpp"
#include "DoFManager.hpp"
#include "WellManager.hpp"
#include "OutputDataVTK.hpp"
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
  m_output_dir = config_dir_path / config.output_dir;
}

void Preprocessor::run()
{
  DiscreteFractureManager dfm_mgr(config.discrete_fractures, data);

  /* Split cells due to edfm intersection */
  EmbeddedFractureManager edfm_mgr(config.embedded_fractures, config.edfm_method, data);
  edfm_mgr.split_cells();

  /* Since we split edfm faces, we pretend that they are dfm fractures
   * to reuse the discretization code. */
  const std::vector<DiscreteFractureConfig> edfm_faces_conf = edfm_mgr.generate_dfm_config();
  // combine dfm and edfm configs
  const std::vector<DiscreteFractureConfig> combined_fracture_config =
      DiscreteFractureManager::combine_configs(config.discrete_fractures, edfm_faces_conf);

  // manages properties of dfm and edfm after splitting
  DiscreteFractureManager fracture_flow_mgr(combined_fracture_config, data);
  fracture_flow_mgr.distribute_properties();
 
  // property manager for grid with split cells (due to edfm splitting)
  CellPropertyManager property_mgr(config.cell_properties, config.domains, data);
  property_mgr.generate_properties();

  // // flow dof numbering
  const std::vector<int> edfm_markers = edfm_mgr.get_face_markers();
  DoFManager dof_manager(data.grid, dfm_mgr.get_face_markers(), edfm_mgr.get_face_markers());
  DoFNumbering split_dofs = dof_manager.distribute_dofs();
  DoFNumbering unsplit_dofs = dof_manager.distribute_unsplit_dofs();

  // build edfm discretization from mixed dfm-edfm discretization
  discretization::DiscretizationEDFM discr_edfm(split_dofs, unsplit_dofs,
                                                data, data.cv_data, data.flow_connection_data,
                                                edfm_markers);
  discr_edfm.build();

  edfm_mgr.build_edfm_grid();
  // generate geomechanics sda properties
  edfm_mgr.distribute_mechanical_properties();
  // map mechanics cells to control volumes
  property_mgr.map_mechanics_to_control_volumes(unsplit_dofs);
  // map sda cells to flow dofs
  edfm_mgr.map_mechanics_to_control_volumes(unsplit_dofs);

  // // setup wells
  // if (config.wells.empty())
  // {
  //   std::cout << "setup wells" << std::endl;
  //   WellManager well_mgr(config.wells, data);
  //   well_mgr.setup();
  // }

  // now coarsen all cells to remove edfm splits and compute normal
  // flow dfm-matrix discretizations
  // data.grid.coarsen_cells();
  // property_mgr.generate_properties();
  // discretization::DiscretizationTPFA matrix_discr(config.discrete_fractures, data);

  BoundaryConditionManager bc_mgr(config.bc_faces, config.bc_nodes, data);
  bc_mgr.create_properties();
  write_output_();
}

void Preprocessor::create_output_dir_()
{
  if (filesystem::exists(m_output_dir))
  {
    std::cout << "cleaning directory " << m_output_dir << std::endl;
    filesystem::remove_all(m_output_dir);
  }
  std::cout << "creating directory " << m_output_dir << std::endl;
  filesystem::create_directory(m_output_dir);
}

void Preprocessor::write_output_()
{
  std::cout << "Write Output data\n";
  create_output_dir_();
  for (auto format : config.output_formats)
  {
    switch (format) {
      case OutputFormat::gprs :
        {
          // std::cout << "Output gprs format" << std::endl;
          // gprs_data::OutputDataGPRS output_data(preprocessor, msh);
          // output_data.write_output(output_dir);
          break;
        }
        case OutputFormat::vtk :
          {
            std::cout << "Output vtk format" << std::endl;
            gprs_data::OutputDataVTK output_data(data, config.vtk_config);
            output_data.write_output(m_output_dir);
            break;
          }
    }
  }
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
