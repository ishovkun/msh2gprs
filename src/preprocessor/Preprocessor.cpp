#include "Preprocessor.hpp"
#include "GlobalOpts.hpp"
#include "parsers/YamlParser.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "mesh/CartesianMeshBuilder.hpp"
#include "mesh/io/VTKReader.hpp"
#include "mesh/RefinementAspectRatio.hpp"
#include "BoundaryConditionManager.hpp"
#include "discretization/mechanics/DiscretizationStandardFEM.hpp"
#include "discretization/mechanics/DiscretizationPolyhedralFEM.hpp"
#include "discretization/mechanics/DiscretizationPolyhedralFEMOptimized.hpp"
#include "discretization/flow/DiscretizationTPFA.hpp"
#include "discretization/flow/DiscretizationDFM.hpp"
#include "discretization/flow/DiscretizationEDFM.hpp"
#include "discretization/flow/DiscretizationPEDFM.hpp"
#include "discretization/flow/DiscretizationINSIM.hpp"
#include "GridEntityNumberingManager.hpp"
#include "multiscale/MSFlow.hpp"
#include "multiscale/MultiScaleDataMech.hpp"
#include "multiscale/Idea.hpp"
#include "DoFManager.hpp"
#include "WellManager.hpp"
#include "INSIMWellManager.hpp"
#include "output/OutputDataVTK.hpp"
#include "output/OutputDataGPRS.hpp"
#include "output/OutputDataPostprocessor.hpp"
#include "output/OutputDataINSIM.hpp"
#include "logger/Logger.hpp"
#include "GridGeneratorINSIM.hpp"
#include <string>
#include <chrono>  // timing

namespace gprs_data {

Preprocessor::Preprocessor(const Path config_file_path)
{
  // read configuration file
  logging::log() << "reading config file: " << config_file_path << std::endl;
  read_config_file_(config_file_path);
  logging::log() << "finished reading config" << std::endl;
  // infer grid file path
  _input_dir = config_file_path.parent_path();
  setup_grid_(_input_dir);
  m_output_dir = _input_dir / _config.output_dir;
}

void Preprocessor::setup_grid_(const Path config_dir_path)
{
  if ( _config.mesh.type == MeshType::file )
  {
    const Path grid_file_path = config_dir_path / _config.mesh.file;
    read_mesh_file_(grid_file_path);
  }
  else if (_config.mesh.type == MeshType::cartesian)
  {
    logging::log() << "Building Cartesian grid...";
    data.grid = mesh::CartesianMeshBuilder(_config.mesh.cartesian);
    logging::log() << "OK" << std::endl;
  }
  else if (_config.mesh.type == MeshType::radial ) {
    throw std::invalid_argument("Radial grids not implemented yet. Contact me if you need those");
  }
  else if ( _config.mesh.type == MeshType::insim ) {
    std::cout << "Building INSIM grid..." << std::endl;
    data.grid = GridGeneratorINSIM(_config.mesh.insim, _config.wells);
    logging::log() << "OK" << std::endl;
  }
  else throw std::invalid_argument("Invalid mesh format");

  for (auto const & transform : _config.mesh.rotations)
    transform.apply(data.grid.vertices().begin(), data.grid.vertices().end());
}

void Preprocessor::run()
{
  // property manager for grid with split cells (due to edfm splitting)
  pm_property_mgr = std::make_shared<CellPropertyManager>(_config.cell_properties, data, _input_dir);
  logging::log() << "Generating properties" << std::endl;
  pm_property_mgr->generate_properties(data.cell_properties);
  data.property_names = pm_property_mgr->get_property_names();
  data.property_types = pm_property_mgr->get_property_types();

  logging::log() << "Initializing grid searcher" << std::endl;
  data.grid_searcher = std::make_unique<GridIntersectionSearcher>(data.grid);

  // create discrete fracture manager
  logging::log() << "Initializing Fracture managers" << std::endl;
  pm_dfm_mgr = std::make_shared<DiscreteFractureManager>(_config.discrete_fractures, data);
  pm_edfm_mgr = std::make_shared<EmbeddedFractureManager>(_config.embedded_fractures,
                                                          _config.edfm_settings, data);

  // copy geomechanics grid since base grid will be split
  data.geomechanics_grid = data.grid;

  /* Split cells due to edfm intersection */
  logging::important() << "Splitting cells..." << std::flush;
  pm_edfm_mgr->split_cells();
  logging::important() << "Finished splitting cells" << std::endl << std::flush;

  if (_config.mesh.refinement.type != RefinementType::none)
  {
    logging::important() << "Peforming Grid Refinement" << std::endl;
    if (_config.mesh.refinement.type == RefinementType::aspect_ratio)
      mesh::RefinementAspectRatio refinement(data.grid, _config.mesh.refinement.aspect_ratio,
                                             _config.mesh.refinement.max_level);
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

  data.has_mechanics = pm_property_mgr->has_mechanics();
  if (data.has_mechanics) {
    logging::important() << "Building geomechanics discretization" << "\n";
    build_geomechanics_discretization_();
    logging::important() << "finished geomechanics discretization" << "\n";
  }

  // Coupling
  // map mechanics cells to control volumes
  if (data.has_mechanics)
    pm_property_mgr->map_mechanics_to_control_volumes(*data.flow_numbering, data.geomechanics_grid);

  // map dfm flow grid to flow dofs
  data.dfm_cell_mapping = pm_cedfm_mgr->map_dfm_grid_to_flow_dofs(data.grid, *data.flow_numbering);

  // remove fine cells
  if ( _config.flow_discretization != FlowDiscretizationType::tpfa_compartmental )
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

  // decide which formats to use for output
  std::vector<OutputFormat> formats;
  if ( _config.flow_discretization == FlowDiscretizationType::insim )
    formats = {OutputFormat::vtk, OutputFormat::insim};
  else
    formats = {OutputFormat::vtk, OutputFormat::gprs, OutputFormat::postprocessor};

  for (auto format : formats)
  {
    switch (format) {
      case OutputFormat::gprs :
        {
          logging::log() << "Output gprs format" << std::endl;
          gprs_data::OutputDataGPRS output_data(data, _config.gprs_output);
          output_data.write_output(m_output_dir);
          break;
        }
        case OutputFormat::vtk :
          {
            logging::log() << "Output vtk format" << std::endl;
            auto & flags = _config.vtk_config.flags;
            if ( !data.fracture_grid.empty() )
              flags |= VTKOutputFlags::save_fractures;
            if ( !data.wells.empty() )
              flags |= VTKOutputFlags::save_wells;
            if ( _config.flow_discretization == FlowDiscretizationType::insim )
              flags |= VTKOutputFlags::save_flow_graph;

            gprs_data::OutputDataVTK output_data(data, _config.vtk_config);
            output_data.write_output(m_output_dir);
            break;
          }
        case OutputFormat::insim :
          {
            gprs_data::OutputDataINSIM output_data(data, _config.insim_output);
            output_data.write_output(m_output_dir);
            break;
          }
        case OutputFormat::postprocessor :
          {
            logging::log() << "Output postprocessor format" << std::endl;
            gprs_data::OutputDataPostprocessor output_data(data, _config,
                                                           *data.flow_numbering,
                                                           m_output_dir);
            output_data.write_output(_config.postprocessor_file);
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
    _config = parser.get_config();
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
  else throw std::invalid_argument("Only .msh files produced by Gmsh are supported");
}

void Preprocessor::build_flow_discretization_()
{
  data.flow.permeability_idx = pm_property_mgr->get_permeability_keys();
  data.flow.porosity_idx = pm_property_mgr->get_porosity_key();
  data.flow.vmult_idx = pm_property_mgr->get_volume_mult_key();
  data.flow.custom_idx = pm_property_mgr->get_custom_flow_keys();

  // add properties for refined cells
  pm_property_mgr->downscale_properties();
  /* Since we split edfm faces, we pretend that they are dfm fractures
   * to reuse the discretization code. */
  const std::vector<DiscreteFractureConfig> edfm_faces_conf = pm_edfm_mgr->generate_dfm_config();
  // combine dfm and edfm configs
  const std::vector<DiscreteFractureConfig> combined_fracture_config =
      DiscreteFractureManager::combine_configs(_config.discrete_fractures, edfm_faces_conf);
  // manages properties of dfm and edfm after splitting
  pm_cedfm_mgr = std::make_shared<DiscreteFractureManager>(combined_fracture_config, data);
  pm_cedfm_mgr->distribute_properties();

  // flow dof numbering
  const std::vector<int> edfm_markers = pm_edfm_mgr->get_face_markers();
  // generate dofs for split and unsplit flow CVs
  std::shared_ptr<DoFNumbering> p_split_dofs = nullptr;
  std::shared_ptr<DoFNumbering> p_unsplit_dofs = nullptr;
  DoFManager dof_manager(data.grid, pm_dfm_mgr->get_face_markers(), edfm_markers);

  std::unique_ptr<INSIMWellManager> insim_mgr = nullptr;
  if ( _config.flow_discretization == FlowDiscretizationType::insim ) {
    insim_mgr = std::make_unique<INSIMWellManager>( _config.wells, data.grid, *data.grid_searcher );
    auto const well_vertices = insim_mgr->get_well_vertices();
    p_split_dofs = dof_manager.distribute_dofs_insim( well_vertices );
    p_unsplit_dofs = dof_manager.distribute_vertex_to_well_dofs( well_vertices );
  }
  else {
    p_split_dofs = dof_manager.distribute_dofs();
    p_unsplit_dofs = dof_manager.distribute_unsplit_dofs();
  }

  using namespace discretization;
  std::unique_ptr<DiscretizationBase> flow_discr;
  switch ( _config.flow_discretization )
  {
    case ( FlowDiscretizationType::tpfa_edfm ):
      logging::important() << "Simple EDFM is chosen" << std::endl;
      flow_discr = std::make_unique<DiscretizationEDFM>(*p_split_dofs, *p_unsplit_dofs, data, data.flow.cv,
                                                        data.flow.con, edfm_markers);
      break;
    case ( FlowDiscretizationType::tpfa_projection ):
      logging::important() << "Projection EDFM is chosen" << std::endl;
      flow_discr = std::make_unique<DiscretizationPEDFM>(*p_split_dofs, *p_unsplit_dofs, data, data.flow.cv,
                                                         data.flow.con, *pm_edfm_mgr);
    break;
    case ( FlowDiscretizationType::tpfa_compartmental ):
      logging::important() << "Compartmental EDFM is chosen" << std::endl;
      flow_discr = std::make_unique<DiscretizationDFM>(*p_split_dofs, data, data.flow.cv, data.flow.con);
    break;
    case ( FlowDiscretizationType::insim ):
      logging::important() << "INSIM discretization is chosen" << std::endl;
      flow_discr = std::make_unique<DiscretizationINSIM>(*p_split_dofs, *p_unsplit_dofs, data, data.flow.cv, data.flow.con);
    break;
  }
  // if we do cedfm use the split matrix dof numbering
  // else use unsplit matrix dofs
  if ( _config.flow_discretization == FlowDiscretizationType::tpfa_compartmental )
    data.flow_numbering = p_split_dofs;
  else
    data.flow_numbering = p_unsplit_dofs;

  flow_discr->build();

  // setup wells
  if ( !_config.wells.empty() && _config.flow_discretization !=  FlowDiscretizationType::insim )
  {
    logging::log() << "setup wells" << std::endl;
    WellManager well_mgr(_config.wells, data, *data.flow_numbering, _config.flow_discretization);
    well_mgr.setup();
  }
  else if ( _config.flow_discretization == FlowDiscretizationType::insim ) {
    logging::log() << "compute well productivity" << std::endl;
    insim_mgr->assign_dofs( *p_split_dofs );
    insim_mgr->compute_well_indices(data.flow.cv);
    data.well_vtk = insim_mgr->get_well_vtk_data();
    data.wells = insim_mgr->get_wells();
  }

  // multiscale idea
  if ( _config.ms_flow.part_type != MSPartitioning::no_partitioning )
  {
    // multiscale::Idea idea(data.grid, data);
    multiscale::MSFlow ms(data.grid, data, _config.ms_flow);
    ms.fill_output_model(data.ms_flow_data);
  }

  // used for coupling later on
  // replace cedfm properties by pure DFM properties
  if ( _config.flow_discretization != FlowDiscretizationType::tpfa_compartmental )
    pm_dfm_mgr->distribute_properties();
}

void Preprocessor::build_geomechanics_discretization_()
{
  data.output_mech_properties = pm_property_mgr->get_custom_mech_keys();

  if (_config.fem.method == FEMMethod::strong_discontinuity)
  {
    // generate geomechanics sda properties
    pm_edfm_mgr->distribute_mechanical_properties();
    // map sda cells to flow dofs
    pm_edfm_mgr->map_mechanics_to_control_volumes(*data.flow_numbering);
  }

  if (_config.ms_mech.support_type == MSSupportType::mechanics)
  {
    multiscale::MultiScaleDataMech ms_handler(data.geomechanics_grid, _config.ms_mech.n_blocks[0]);
    ms_handler.build_data();
    ms_handler.fill_output_model(data.ms_mech_data);
  }

  auto p_frac_mgr = pm_dfm_mgr;

  std::vector<int> dfm_markers;
  if (_config.fem.method == FEMMethod::polyhedral_finite_element ||
      _config.fem.method == FEMMethod::mixed)
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

  // build mechanics boundary conditions
  logging::log() << "Building mechanics boundary conditions" << std::endl;
  BoundaryConditionManager bc_mgr(_config.bc_faces, _config.bc_nodes, data);

  logging::log() << "Building FEM discretization" << std::endl;
  dfm_markers = p_frac_mgr->get_face_markers();

  // invoke fem discretization class
  using namespace discretization;
  std::unique_ptr<DiscretizationFEMBase> p_discr;
  // WARNING: this part of code can reorder grid vertices
  if (_config.fem.method == FEMMethod::strong_discontinuity)
    p_discr = std::make_unique<DiscretizationStandardFEM>(data.geomechanics_grid, _config.fem,
                                                          dfm_markers, data.neumann_face_indices);
  else if (_config.fem.method == FEMMethod::polyhedral_finite_element) {
    auto const & opts = GlobalOpts::ref();
#ifdef WITH_EIGEN
    if (opts.enable_experimental)
      p_discr = std::make_unique<DiscretizationPolyhedralFEMOptimized>(data.geomechanics_grid, config.fem,
                                                                       dfm_markers, data.neumann_face_indices);
    else
      p_discr = std::make_unique<DiscretizationPolyhedralFEM>(data.geomechanics_grid, config.fem,
                                                              dfm_markers, data.neumann_face_indices);
#elseif
    throw std::invalid_argument("PFEM is not available without Eigen");
#endif
  }
  else throw std::invalid_argument("mechanics discretization is unknown");

  using namespace std::chrono;
  auto const start_time = high_resolution_clock::now();
  p_discr->build();
  auto const end_time = high_resolution_clock::now();
  logging::debug() << "Time to build geomech discretization: "
                   << (duration_cast<milliseconds>(end_time - start_time)).count()
                   << " [ms]" << std::endl;

  // split geomechanics DFM faces
  // This part must go after building discretization since
  // the discretization might do vertex renumbering
  if (_config.dfm_settings.split_mech_vertices)
  {
    logging::log() << "Splitting faces of DFM fractures" << std::endl;
    p_frac_mgr->split_faces(data.geomechanics_grid);
  }

  data.fe_cell_data = p_discr->get_cell_data();
  data.fe_face_data = p_discr->get_face_data();
  data.fe_frac_data = p_discr->get_fracture_data();
  data.isomorphic_groups = p_discr->get_cell_isomorphic_groups();
}


}  // end namespace gprs_data
