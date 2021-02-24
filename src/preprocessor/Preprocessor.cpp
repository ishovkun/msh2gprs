#include "Preprocessor.hpp"
#include "parsers/YamlParser.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "mesh/CartesianMeshBuilder.hpp"
#include "mesh/io/VTKReader.hpp"
#include "mesh/RefinementAspectRatio.hpp"
#include "BoundaryConditionManager.hpp"
#include "discretization/mechanics/DiscretizationFEM.hpp"
#include "discretization/flow/DiscretizationTPFA.hpp"
#include "discretization/flow/DiscretizationDFM.hpp"
#include "discretization/flow/DiscretizationEDFM.hpp"
#include "discretization/flow/DiscretizationPEDFM.hpp"
#include "GridEntityNumberingManager.hpp"
#include "MultiScaleDataMech.hpp"
#include "DoFManager.hpp"
#include "WellManager.hpp"
#include "OutputDataVTK.hpp"
#include "OutputDataGPRS.hpp"
#include "OutputDataPostprocessor.hpp"
#include "logger/Logger.hpp"
#include <string>

namespace gprs_data {

Preprocessor::Preprocessor(const Path config_file_path)
{
  // read configuration file
  logging::log() << "reading config file: " << config_file_path << std::endl;
  read_config_file_(config_file_path);
  logging::log() << "finished reading config" << std::endl;
  // infer grid file path
  const Path config_dir_path = config_file_path.parent_path();
  setup_grid_(config_dir_path);
  m_output_dir = config_dir_path / config.output_dir;
}

void Preprocessor::setup_grid_(const Path config_dir_path)
{
  if (config.mesh_config.type == MeshType::file)
  {
    const Path grid_file_path = config_dir_path / config.mesh_config.file;
    read_mesh_file_(grid_file_path);
  }
  else if (config.mesh_config.type == MeshType::cartesian)
  {
    logging::log() << "Building Cartesian grid";
    data.grid = mesh::CartesianMeshBuilder(config.mesh_config.cartesian);
  }
  else throw std::invalid_argument("Invalid mesh format");
}

void Preprocessor::run()
{
  // property manager for grid with split cells (due to edfm splitting)
  pm_property_mgr = std::make_shared<CellPropertyManager>(config.cell_properties, config.domains, data);
  logging::log() << "Generating properties" << std::endl;
  pm_property_mgr->generate_properties();

  logging::log() << "Initializing grid searcher" << std::endl;
  data.grid_searcher = std::make_unique<GridIntersectionSearcher>(data.grid);

  // create discrete fracture manager
  logging::log() << "Initializing Fracture managers" << std::endl;
  pm_dfm_mgr = std::make_shared<DiscreteFractureManager>(config.discrete_fractures, data);
  pm_edfm_mgr = std::make_shared<EmbeddedFractureManager>(config.embedded_fractures,
                                                          config.edfm_settings, data);

  // copy geomechanics grid since base grid will be split
  data.geomechanics_grid = data.grid;

  /* Split cells due to edfm intersection */
  logging::important() << "Splitting cells..." << std::flush;
  pm_edfm_mgr->split_cells();
  logging::important() << "Finished splitting cells" << std::endl << std::flush;

  if (config.mesh_config.refinement.type != RefinementType::none)
  {
    logging::important() << "Peforming Grid Refinement" << std::endl;
    if (config.mesh_config.refinement.type == RefinementType::aspect_ratio)
      mesh::RefinementAspectRatio refinement(data.grid, config.mesh_config.refinement.aspect_ratio,
                                             config.mesh_config.refinement.max_level);
    logging::important() << "Finished Grid Refinement" << std::endl;
  }

  // leave this line commented out while not debugging,
  // this is helpful when I debug :-)
  // mesh::IO::VTKWriter::write_geometry(data.grid, "debug.vtk");

  logging::important() << "Building flow discretization" << "\n";
  build_flow_discretization_();
  logging::important() << "Finished flow discretization" << std::endl;


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
  if (config.edfm_settings.method != EDFMMethod::compartmental)
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
    logging::debug() << "cleaning directory " << m_output_dir << std::endl;
    filesystem::remove_all(m_output_dir);
  }
  logging::debug() << "creating directory " << m_output_dir << std::endl;
  filesystem::create_directory(m_output_dir);
}

void Preprocessor::write_output_()
{
  logging::log() << "Write Output data\n";
  create_output_dir_();
  for (auto format : config.output_formats)
  {
    switch (format) {
      case OutputFormat::gprs :
        {
          logging::log() << "Output gprs format" << std::endl;
          gprs_data::OutputDataGPRS output_data(data, config.gprs_output);
          output_data.write_output(m_output_dir);
          break;
        }
        case OutputFormat::vtk :
          {
            logging::log() << "Output vtk format" << std::endl;
            gprs_data::OutputDataVTK output_data(data, config.vtk_config);
            output_data.write_output(m_output_dir);
            break;
          }
        case OutputFormat::postprocessor :
          {
            logging::log() << "Output postprocessor format" << std::endl;
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
    const std::string error_msg = "config file does not exist: " +
                                  std::string(filesystem::absolute(config_file_path));
    throw std::invalid_argument(error_msg);
  }

  logging::important() << "Reading ";
  logging::important() << filesystem::absolute(config_file_path) << std::endl;

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
    logging::warning() << "Only .yaml configuration files are supported" << std::endl;
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

  logging::important() << "Reading mesh file: " << filesystem::absolute(mesh_file_path) << std::endl;

  // check filetype
  const std::string fname = mesh_file_path.filename();
  const std::size_t str_len = fname.size();
  const std::string extension = fname.substr(str_len - 3, str_len);

  logging::debug() << "extension: " << extension << std::endl;
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

  using namespace discretization;
  std::unique_ptr<DiscretizationBase> flow_discr;
  switch (config.edfm_settings.method)
  {
    case (EDFMMethod::simple):
      flow_discr = std::make_unique<DiscretizationEDFM>
          (*p_split_dofs, *p_unsplit_dofs, data, data.cv_data,
           data.flow_connection_data, edfm_markers);
      break;
    case (EDFMMethod::projection):
      flow_discr = std::make_unique<DiscretizationPEDFM>
          (*p_split_dofs, *p_unsplit_dofs, data, data.cv_data,
           data.flow_connection_data, *pm_edfm_mgr);
    break;
    case (EDFMMethod::compartmental):
      flow_discr = std::make_unique<DiscretizationDFM>
          (*p_split_dofs, data, data.cv_data, data.flow_connection_data);
    break;
  }
  // if we do cedfm use the split matrix dof numbering
  // else use unsplit matrix dofs
  if ( config.edfm_settings.method == EDFMMethod::compartmental )
    data.flow_numbering = p_split_dofs;
  else
    data.flow_numbering = p_unsplit_dofs;

  logging::debug() << "invoke discretization class" << std::endl;
  flow_discr->build();

  // setup wells
  if (!config.wells.empty())
  {
    logging::log() << "setup wells" << std::endl;
    WellManager well_mgr(config.wells, data, *data.flow_numbering, config.edfm_settings.method);
    well_mgr.setup();
  }

  // used for coupling later on
  // replace cedfm properties by pure DFM properties
  if ( config.edfm_settings.method != EDFMMethod::compartmental )
    pm_dfm_mgr->distribute_properties();
}

void Preprocessor::build_geomechanics_discretization_()
{
  if (config.fem.method == FEMMethod::strong_discontinuity)
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

  // renumber nodes of some grid faces
  // this is needed for AD-GPRS DFM fractures
  for (auto face = data.geomechanics_grid.begin_active_faces(); face != data.geomechanics_grid.end_active_faces(); ++face)
  {
    if ( face->normal().dot(face->center() - face->neighbors()[0]->center()) > 0 )
    {
      mesh::Face & nc_face = const_cast<mesh::Face&>(*face);
      auto & verts = nc_face.vertices();
      std::reverse(verts.begin(), verts.end());
    }
  }

  GridEntityNumberingManager mech_numbering_mgr(data.geomechanics_grid);
  data.mech_numbering = std::shared_ptr<discretization::DoFNumbering>(mech_numbering_mgr.get_numbering());

  // split geomechanics DFM faces
  if (config.dfm_settings.split_mech_vertices)
  {
    logging::log() << "Splitting faces of DFM fractures" << std::endl;
    p_frac_mgr->split_faces(data.geomechanics_grid);
  }

  // build mechanics boundary conditions
  logging::log() << "Building mechanics boundary conditions" << std::endl;
  BoundaryConditionManager bc_mgr(config.bc_faces, config.bc_nodes, data);

  logging::log() << "Building FEM discretization" << std::endl;
  dfm_markers = p_frac_mgr->get_face_markers();
  discretization::DiscretizationFEM fem_discr(data.geomechanics_grid, config.fem, dfm_markers,
                                              bc_mgr.get_neumann_face_markers());
  data.fe_cell_data = fem_discr.get_cell_data();
  data.fe_face_data = fem_discr.get_face_data();
  data.fe_frac_data = fem_discr.get_fracture_data();
}


}  // end namespace gprs_data
