#include "Preprocessor.hpp"
#include "parsers/YamlParser.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "mesh/io/VTKReader.hpp"
#include "BoundaryConditionManager.hpp"
#include "discretization/mechanics/DiscretizationFEM.hpp"
#include "discretization/flow/DiscretizationTPFA.hpp"
#include "discretization/flow/DiscretizationEDFM.hpp"
#include "discretization/flow/DiscretizationDFM.hpp"
#include "GridEntityNumberingManager.hpp"
#include "MultiScaleDataMech.hpp"
#include "DoFManager.hpp"
#include "WellManager.hpp"
#include "OutputDataVTK.hpp"
#include "OutputDataGPRS.hpp"
#include "OutputDataPostprocessor.hpp"
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
  // property manager for grid with split cells (due to edfm splitting)
  pm_property_mgr = std::make_shared<CellPropertyManager>(config.cell_properties, config.domains, data);
  std::cout << "Generating properties" << std::endl;
  pm_property_mgr->generate_properties();

  // create discrete fracture manager
  std::cout << "Initializing Fracture managers" << std::endl;
  pm_dfm_mgr = std::make_shared<DiscreteFractureManager>(config.discrete_fractures, data);
  pm_edfm_mgr = std::make_shared<EmbeddedFractureManager>(config.embedded_fractures, config.edfm_method,
                                                          config.edfm_min_dist_to_node, data);

  // copy geomechanics grid since base grid will be split
  data.geomechanics_grid = data.grid;

  /* Split cells due to edfm intersection */
  std::cout << "Splitting cells..." << std::flush;
  pm_edfm_mgr->split_cells();
  std::cout << "OK" << std::endl << std::flush;

  std::cout << "building flow discretization" << std::endl;
  build_flow_discretization_();

  // build edfm grid for vtk output
  // pm_edfm_mgr->build_edfm_grid(*data.flow_numbering);

  // build dfm flow and mechanics grids (for vtk output)
  data.fracture_grid = pm_cedfm_mgr->build_dfm_grid(data.grid);

  if (data.has_mechanics)
    build_geomechanics_discretization_();

  // Coupling
  // map mechanics cells to control volumes
  if (data.has_mechanics)
    pm_property_mgr->map_mechanics_to_control_volumes(*data.flow_numbering, data.geomechanics_grid);

  // map dfm flow grid to flow dofs
  data.dfm_cell_mapping = pm_cedfm_mgr->map_dfm_grid_to_flow_dofs(data.grid, *data.flow_numbering);

  // remove fine cells
  if (config.edfm_method != EDFMMethod::compartmental)
  {
    data.grid.coarsen_cells();
    // pm_property_mgr->coarsen_cells();
  }

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
          std::cout << "Output gprs format" << std::endl;
          gprs_data::OutputDataGPRS output_data(data, config.gprs_output);
          output_data.write_output(m_output_dir);
          break;
        }
        case OutputFormat::vtk :
          {
            std::cout << "Output vtk format" << std::endl;
            gprs_data::OutputDataVTK output_data(data, config.vtk_config);
            output_data.write_output(m_output_dir);
            break;
          }
        case OutputFormat::postprocessor :
          {
            std::cout << "Output postprocessor format" << std::endl;
            gprs_data::OutputDataPostprocessor output_data(data, config,
                                                           *data.flow_numbering,
                                                           m_output_dir);
            output_data.write_output(config.postprocessor_file);
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
    parser.parse_file(filesystem::absolute(config_file_path));
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

  std::cout << "reading " << filesystem::absolute(mesh_file_path) << std::endl;

  // check filetype
  const std::string fname = mesh_file_path.filename();
  const std::size_t str_len = fname.size();
  const std::string extension = fname.substr(str_len - 3, str_len);

  std::cout << "extension:" << extension << std::endl;
  if (extension == "msh")
    GmshInterface::read_msh(filesystem::absolute(mesh_file_path), data.grid);
  else if (extension == "vtk")
  {
    mesh::io::VTKReader reader(filesystem::absolute(mesh_file_path), data.grid);
    const auto & keys = reader.get_cell_data_keys();
    if (!keys.empty())
    {
      assert( keys.size() == 1 );
      assert( keys[0] == "marker" );
      auto & grid = data.grid;
      size_t icell = 0;
      auto & data = reader.get_cell_data();
      for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
        cell->set_marker(static_cast<int>(data[0][icell++]));
    }
  }
  else
    throw std::invalid_argument("Only .msh files produced by Gmsh are supported");
}

void Preprocessor::build_flow_discretization_()
{
  // add properties for refined cells
  pm_property_mgr->downscale_properties();
  /* Since we split edfm faces, we pretend that they are dfm fractures
   * to reuse the discretization code. */
  const std::vector<DiscreteFractureConfig> edfm_faces_conf = pm_edfm_mgr->generate_dfm_config();
  // combine dfm and edfm configs
  const std::vector<DiscreteFractureConfig> combined_fracture_config =
      DiscreteFractureManager::combine_configs(config.discrete_fractures, edfm_faces_conf);
  // manages properties of dfm and edfm after splitting
  pm_cedfm_mgr = std::make_shared<DiscreteFractureManager>(combined_fracture_config, data);
  pm_cedfm_mgr->distribute_properties();

  // flow dof numbering
  const std::vector<int> edfm_markers = pm_edfm_mgr->get_face_markers();
  // generate dofs for split and unsplit flow CVs
  DoFManager dof_manager(data.grid, pm_dfm_mgr->get_face_markers(), edfm_markers);
  std::shared_ptr<DoFNumbering> p_split_dofs = dof_manager.distribute_dofs();
  std::shared_ptr<DoFNumbering> p_unsplit_dofs = dof_manager.distribute_unsplit_dofs();

  // build edfm discretization from mixed dfm-edfm discretization
  discretization::DiscretizationEDFM discr_edfm(*p_split_dofs, *p_unsplit_dofs, data, data.cv_data,
                                                data.flow_connection_data, edfm_markers, config.edfm_method);
  // if we do cedfm use the split matrix dof numbering
  // else use unsplit matrix dofs
  if ( config.edfm_method == EDFMMethod::compartmental )
    data.flow_numbering = p_split_dofs;
  else
    data.flow_numbering = p_unsplit_dofs;

  std::cout << "build discr" << std::endl;
  discr_edfm.build();

  // setup wells
  if (!config.wells.empty())
  {
    std::cout << "setup wells" << std::endl;
    WellManager well_mgr(config.wells, data, *data.flow_numbering);
    well_mgr.setup();
  }

  // used for coupling later on
  if ( config.edfm_method != EDFMMethod::compartmental )
    pm_dfm_mgr->distribute_properties();

  std::cout << "done flow discretization" << std::endl;
}

void Preprocessor::build_geomechanics_discretization_()
{
  if (config.fem.method == strong_discontinuity)
  {
    // generate geomechanics sda properties
    pm_edfm_mgr->distribute_mechanical_properties();
    // map sda cells to flow dofs
    pm_edfm_mgr->map_mechanics_to_control_volumes(*data.flow_numbering);
  }

  if (config.multiscale_mechanics == MSPartitioning::method_mechanics)
  {
    multiscale::MultiScaleDataMech ms_handler(data.geomechanics_grid, config.n_multiscale_blocks);
    ms_handler.build_data();
    ms_handler.fill_output_model(data.ms_mech_data);
  }

  auto p_frac_mgr = pm_dfm_mgr;

  std::vector<int> dfm_markers;
  if (config.fem.method == FEMMethod::polyhedral_finite_element ||
      config.fem.method == FEMMethod::mixed)
  {
    data.geomechanics_grid = data.grid;
    p_frac_mgr = pm_cedfm_mgr;
  }

  GridEntityNumberingManager mech_numbering_mgr(data.geomechanics_grid);
  data.mech_numbering = std::shared_ptr<discretization::DoFNumbering>(mech_numbering_mgr.get_numbering());

  // split geomechanics DFM faces
  std::cout << "splitting faces of DFM fractures" << std::endl;
  p_frac_mgr->split_faces(data.geomechanics_grid);

  // build mechanics boundary conditions
  std::cout << "Building mechanics boundary conditions" << std::endl;
  BoundaryConditionManager bc_mgr(config.bc_faces, config.bc_nodes, data);

  std::cout << "Building FEM discretization" << std::endl;
  dfm_markers = p_frac_mgr->get_face_markers();
  discretization::DiscretizationFEM fem_discr(data.geomechanics_grid, config.fem, dfm_markers,
                                              bc_mgr.get_neumann_face_markers());
  data.fe_cell_data = fem_discr.get_cell_data();
  data.fe_face_data = fem_discr.get_face_data();
  data.fe_frac_data = fem_discr.get_fracture_data();
}


}  // end namespace gprs_data
