// #include "DiscretizationMechanics.hpp"
// #include "DiscretizationPolyhedralFEM.hpp"
// #include "DiscretizationStandardFEM.hpp"
// #include "gmsh_interface/GmshInterface.hpp"
// #include <stdexcept>
// #include "MeshStatsComputer.hpp"
// #include "logger/ProgressBar.hpp"  // provides ProgressBar
// #include "mesh/io/VTKWriter.hpp"   // debugging, provides io::VTKWriter
// // #include "Isomorphism.hpp"         // provides isomorphism


// namespace discretization
// {
// using Point = angem::Point<3,double>;

// DiscretizationMechanics::DiscretizationMechanics(mesh::Mesh & grid, const FiniteElementConfig & config,
//                                                  const std::vector<int> & fracture_face_markers,
//                                                  const std::vector<int> & neumann_face_markers)
//     : _grid(grid), _config( config )
// {
//   #ifndef WITH_EIGEN
//   throw std::runtime_error("Cannot use DFEM method without linking to Eigen");
//   #endif

//   for (const int marker : fracture_face_markers)
//     _fracture_face_orientation[marker] = {angem::Basis<3, double>(), false};
//   for (const int marker : neumann_face_markers)
//     _neumann_face_orientation[marker] = {angem::Basis<3, double>(), false};

//   _face_data.resize( _grid.n_faces() );
//   _cell_data.resize( _grid.n_cells_total() );
//   _frac_data.resize( _grid.n_faces() );

//   std::unique_ptr<DiscretizationFEMBase> p_discr = factory();

//   logging::ProgressBar progress("Build Finite Elements", _grid.n_active_cells());

//   size_t item = 0;
//   std::vector<size_t> order;
//   for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
//   {
//     progress.set_progress(item++);
//     // FiniteElementDataTopology const * p_master = nullptr;
//     p_discr->build(*cell);
//     // if ( known_element_(*cell, order, p_master) )
//     // {
//     //   cell->reorder_vertices(order);
//     //   scale_cell_fem_data_(*cell, *p_master);
//     //   throw std::runtime_error("Awesome");
//     // }
//     // else
//     // {
//     //   std::unique_ptr<FiniteElementBase> p_discr;
//     //   try {
//     //     p_discr = build_element(*cell);
//     //   } catch (std::runtime_error &error) {
//     //     mesh::IO::VTKWriter::write_geometry(
//     //         _grid, *cell, "geometry-" + std::to_string(cell->index()) + ".vtk");
//     //     throw std::runtime_error(error.what());
//     //   }

//       // build_cell_data_(*p_discr, *cell);
//       // build_face_data_(*p_discr, *cell);
//     // }
//   }

//   progress.finalize();
// }

// void DiscretizationMechanics::scale_cell_fem_data_(mesh::Cell const & cell,
//                                              FiniteElementDataTopology const & data)
// {

//   std::cout << "scaling shit" << std::endl;
//   exit(0);
// }

// // bool DiscretizationMechanics::known_element_(mesh::Cell const & cell,
// //                                        std::vector<size_t> & order,
// //                                        FiniteElementDataTopology const *& p_master) const
// // {
// //   size_t const hsh = cell.n_vertices();
// //   if (_cell_data_compressed.count(hsh))
// //   {
// //     auto it = _cell_data_compressed.find(hsh);
// //     for (auto const & master : it->second)
// //     {
// //       auto result = Isomorphism::check(*master.topology, *cell.polyhedron());
// //       if (result.first)
// //       {
// //         order = result.second;
// //         p_master = &master;
// //         return true;
// //       }
// //     }
// //     return false;
// //   }
// //   else return false;
// // }

// void DiscretizationMechanics::build_cell_data_(FiniteElementBase & discr,
//                                                mesh::Cell const & cell)
// {
//   FiniteElementData cell_fem_data = discr.get_cell_data();
//   cell_fem_data.element_index = cell.index();
//   _cell_data[cell.index()] = std::move(cell_fem_data);
//   auto & it = _cell_data_compressed[cell.n_vertices()];
//   it.emplace_back();
//   auto & saved = it.back();
//   saved.data = _cell_data[cell.index()];
//   saved.topology = cell.polyhedron();
// }

// void DiscretizationMechanics::build_face_data_(FiniteElementBase & discr,
//                                                mesh::Cell const & cell)
// {
//   size_t iface = 0;
//   for ( const mesh::Face * face : cell.faces() )
//   {
//     const size_t face_index = face->index();
//     const bool is_fracture = _fracture_face_orientation.find( face->marker() ) !=
//         _fracture_face_orientation.end();
//     const bool is_neumann = _neumann_face_orientation.find( face->marker() ) !=
//         _neumann_face_orientation.end();
//     if (is_fracture)
//       _face_basis = get_basis_(*face, _fracture_face_orientation.find(face->marker())->second);
//     else if (is_neumann)
//       _face_basis = get_basis_(*face, _neumann_face_orientation.find(face->marker())->second);

//     if (is_neumann || is_fracture)
//       if ( _face_data[face->index()].points.empty() )
//       {
//         _face_data[face_index] = discr.get_face_data(iface, _face_basis);
//         _face_data[face_index].element_index = face_index;
//       }

//     if (is_fracture)
//     {
//       _frac_data[face_index].push_back(discr.get_fracture_data(iface, _face_basis));
//       _frac_data[face_index].back().element_index = cell.index();
//     }

//     iface++;
//   }
// }

// void DiscretizationMechanics::analyze_cell_(const mesh::Cell & cell)
// {
// //   std::cout << "vertices: ";
// //   for(auto v : cell.vertices())
// //     std::cout << v << " ";
// //   std::cout << std::endl;

// //   std::cout << "faces:" << std::endl;
// //   for (auto * f : cell.faces())
// //   {
// //     for (auto v : f->vertices())
// //       std::cout << v << " ";
// //     std::cout << std::endl;
// //   }

// //   mesh::IO::VTKWriter::write_geometry(_grid, cell, "output/geometry-" + std::to_string(cell.index()) + ".vtk");

// // #ifdef WITH_EIGEN
// //   PolyhedralElementDirect de(cell, _grid, _config);
// //   // PolyhedralElementMSRSB de(cell, _config);
// //   de.save_shape_functions("output/shape_functions-" + std::to_string(cell.index())+ ".vtk");
// //   de.debug_save_boundary_face_solution("output/face_solutions.vtk");

// //   auto verts = cell.vertices();

// //   auto data =  de.get_cell_data();

// //   std::cout << "============================" << std::endl;
// //   std::cout << "Weights:" << std::endl;
// //   double wsum = 0;
// //   for (size_t q=0; q<data.points.size(); ++q)
// //   {
// //     std::cout << data.points[q].weight << " ";
// //     wsum += data.points[q].weight;
// //   }
// //   std::cout << std::endl;

// //   std::cout << "weights sum = " << wsum << "  (should be " << data.center.weight << ")" << std::endl;


// //   std::cout << "============================" << std::endl;
// //   std::cout << "Values:" << std::endl;
// //   for (size_t q=0; q<data.points.size(); ++q)
// //   {
// //     const double sum = std::accumulate(data.points[q].values.begin(),
// //                                        data.points[q].values.end(), 0.0);
// //     if (std::fabs(sum - 1.0) > 1e-10)
// //     {
// //       std::cout << "Error: q = " << q << " sum values = " << sum << std::endl;
// //     }
// //   }

// //   std::cout << "============================" << std::endl;
// //   std::cout << "Gradients:" << std::endl;
// //   for (size_t q=0; q<data.points.size(); ++q)
// //   {
// //     angem::Point<3,double> sum;
// //     for (auto & p : data.points[q].grads)
// //       sum += p;
// //     std::cout << q << " " << sum << std::endl;
// //   }

// //   exit(0);

// //   // std::cout << "======= shape functions ==========" << std::endl;
// //   // std::cout << "analytic" << std::endl;
// //   // data =  fe.get_cell_data();
// //   // for (size_t q=0; q<data.points.size(); ++q)
// //   // {
// //   //   std::cout << q << ": ";
// //   //   for (size_t v=0; v<verts.size(); ++v)
// //   //     std::cout << data.points[q].values[v] << " ";
// //   //   std::cout << std::endl;
// //   // }
// //   // std::cout << "numeric" << std::endl;
// //   // data =  de.get_cell_data();
// //   // for (size_t q=0; q<data.points.size(); ++q)
// //   // {
// //   //   std::cout << q << ": ";
// //   //   for (size_t v=0; v<verts.size(); ++v)
// //   //     std::cout << data.points[q].values[v] << " ";
// //   //   std::cout << std::endl;
// //   // }

// //   // This proves partition of unity
// //   std::cout << "======= partition of unity ===========" << std::endl;
// //   data =  de.get_cell_data();
// //   {
// //     double maxsum = 0;
// //     for (auto & qp : data.points)
// //     {
// //       double sum = 0;
// //       for (size_t v=0; v<qp.values.size(); ++v)
// //       {
// //         sum += qp.values[v];
// //       }
// //       maxsum = std::max(sum, maxsum);
// //     }
// //     std::cout << "numeric sum = " << std::scientific << std::fabs(maxsum - 1.0) << std::endl << std::defaultfloat;
// //   }

// //   std::cout << "======= patch test ========" << std::endl;
// //   {
// //     Point maxsum = {0,0,0};
// //     for (auto & qp : data.points)
// //     {
// //       for (size_t i=0; i<3; ++i)
// //       {
// //         double sum = 0;
// //         for (size_t v=0; v<qp.values.size(); ++v)
// //           sum += qp.grads[v][i];
// //         maxsum[i] = std::max(maxsum[i], std::fabs(sum));
// //       }

// //     }
// //     std::cout << "grad maxsum = " << maxsum << std::endl;
// //   }

// //   exit(0);
// // #endif
// }

// std::unique_ptr<DiscretizationFEMBase> DiscretizationMechanics::factory()
// {
//   std::unique_ptr<DiscretizationFEMBase> p_discr = nullptr;
//   if (_config.method == FEMMethod::strong_discontinuity)
//     p_discr = std::make_unique<DiscretizationStandardFEM>();
//   else if (_config.method == FEMMethod::polyhedral_finite_element) {
// #ifdef WITH_EIGEN
//     p_discr = std::make_unique<DiscretizationPolyhedralFEM>(_grid, _config);
// #else
//     throw std::runtime_error("Cannot use PDFEM method without linking to Eigen");
// #endif
//   }
//   else throw std::invalid_argument("Invalid mechanics discretization type");

//   return p_discr;
// }


// // std::unique_ptr<FiniteElementBase> DiscretizationMechanics::build_element(const mesh::Cell & cell)
// // {
// //   std::unique_ptr<FiniteElementBase> p_discr;
// //   if (_config.method == FEMMethod::polyhedral_finite_element)
// //   {
// // #ifdef WITH_EIGEN
// //     if (_config.solver == SolverType::direct || _config.solver == SolverType::cg)
// //       p_discr = std::make_unique<PolyhedralElementDirect>(cell, _grid, _config);
// //     else if (_config.solver == SolverType::msrsb)
// //       p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _grid, _config);
// // #else
// //     throw std::runtime_error("Cannot use PDFEM method without linking to Eigen");
// // #endif

// //   }
// //   else if (_config.method == FEMMethod::strong_discontinuity)
// //     p_discr = std::make_unique<StandardFiniteElement>(cell);
// //   else if (_config.method == FEMMethod::mixed)
// //   {
// //     if (cell.vtk_id() == angem::GeneralPolyhedronID)
// //     {
// // #ifdef WITH_EIGEN
// //       if (_config.solver == SolverType::direct || _config.solver == SolverType::cg)
// //         p_discr = std::make_unique<PolyhedralElementDirect>(cell, _grid, _config);
// //       else if (_config.solver == SolverType::msrsb)
// //         p_discr = std::make_unique<PolyhedralElementMSRSB>(cell, _grid, _config);
// // #else
// //     throw std::runtime_error("Cannot use PFEM method without linking to Eigen");
// // #endif
// //     }
// //     else
// //       p_discr = std::make_unique<StandardFiniteElement>(cell);
// //   }

// //   return p_discr;
// // }

// const angem::Basis<3, double> &
// DiscretizationMechanics::get_basis_(const mesh::Face & face,
//                               FaceOrientation &orientation) noexcept
// {
//   if (!orientation.assigned)
//   {
//     orientation.basis = face.polygon().plane().get_basis();
//     orientation.assigned = true;
//   }

//   return orientation.basis;
// }


// }  // end namepsace discretization
