#include "DiscretizationFEM.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include <stdexcept>
#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "PolyhedralElementDirect.hpp"
#endif
#include "StandardFiniteElement.hpp"
#include "ProgressBar.hpp"  // provides ProgressBar


namespace discretization
{
using Point = angem::Point<3,double>;
using api = gprs_data::GmshInterface;

DiscretizationFEM::DiscretizationFEM(const mesh::Mesh & grid, const FiniteElementConfig & config)
    : _grid(grid), _config( config )
{
  #ifndef WITH_GMSH
  throw std::runtime_error("Cannot use DFEM method without linking to GMsh");
  #endif
  #ifndef WITH_EIGEN
  throw std::runtime_error("Cannot use DFEM method without linking to Eigen");
  #endif
}

#ifdef WITH_GMSH

void DiscretizationFEM::build()
{
  // analyze_cell_(_grid.cell(0));
  // analyze_cell_(_grid.cell(1));

  _face_data.resize( _grid.n_faces() );
  _cell_data.resize( _grid.n_cells() );
  utils::ProgressBar progress("Build Finite Elements", _grid.n_active_cells());
  size_t item = 0;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    progress.set_progress(item++);
    api::initialize_gmsh();
    const std::unique_ptr<FiniteElementBase> p_discr = build_element(*cell);
    api::finalize_gmsh();

   FiniteElementData cell_fem_data = p_discr->get_cell_data();
   cell_fem_data.element_index = cell->index();
   _cell_data[cell->index()] = std::move(cell_fem_data);

   std::vector<FiniteElementData> face_data = p_discr->get_face_data();
   size_t iface = 0;
   for ( const mesh::Face * face : cell->faces() )
   {
     if ( _face_data[face->index()].points.empty() )
     {
       face_data[iface].element_index = face->index();
       _face_data[face->index()] = face_data[iface++];
     }
   }
  }
  progress.finalize();
}

void DiscretizationFEM::analyze_cell_(const mesh::Cell & cell)
{
  api::initialize_gmsh();
  PolyhedralElementDirect de(cell, _config);
  StandardFiniteElement fe(cell);
  // DFEMElement discr_element(cell, _msrsb_tol);
  api::finalize_gmsh();
  de.debug_save_shape_functions_("output/shape_functions" + std::to_string(cell.index())+ ".vtk");

  FiniteElementData an_data =  fe.get_cell_data();
  auto verts = cell.vertices();
  std::cout << "+======COORDINATES" << std::endl;
  // std::cout << "analytic" << std::endl;
  // for (size_t q=0; q<an_data.points.size(); ++q)
  // {
  //   Point p (0,00,0);
  //   for (size_t v=0; v<verts.size(); ++v)
  //   {
  //     p += _grid.vertex(verts[v]) * an_data.points[q].values[v];
  //   }
  //   std::cout << "p = " << p << std::endl;
  // }

  auto data =  de.get_cell_data();
  auto qpoints = de.get_cell_integration_points();
  std::cout << "============================" << std::endl;
  std::cout << "Recovery integration points:" << std::endl;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    Point p (0,00,0);
    for (size_t v=0; v<verts.size(); ++v)
    {
      p += _grid.vertex(verts[v]) * data.points[q].values[v];
    }

    std::cout << "diff = " << (p - qpoints[q]).norm()<< std::endl;
  }

  std::cout << "============================" << std::endl;
  std::cout << "Weights:" << std::endl;
  double wsum = 0;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    std::cout << data.points[q].weight << " ";
    wsum += data.points[q].weight;
  }
  std::cout << std::endl;

  std::cout << "weights sum = " << wsum << "  (should be " << data.center.weight << ")" << std::endl;


  // std::cout << "============================" << std::endl;
  // std::cout << "Values:" << std::endl;
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << "\n";
  //   for (size_t v=0; v<verts.size(); ++v)
  //   {
  //     std::cout << data.points[q].values[v] << " "
  //               << an_data.points[q].values[v] << " "
  //               << std::fabs( (data.points[q].values[v] - an_data.points[q].values[v]) / an_data.points[q].values[v]  )
  //               << std::endl;
  //   }
  // }

  std::cout << "============================" << std::endl;
  std::cout << "Gradients:" << std::endl;
  for (size_t q=0; q<data.points.size(); ++q)
  {
    std::cout << q << "\n";
    for (size_t v=0; v<verts.size(); ++v)
    {
      std::cout <<  data.points[q].grads[v] << "\t";
    }
    std::cout << std::endl;
  }


  // std::cout << "======= shape functions ==========" << std::endl;
  // std::cout << "analytic" << std::endl;
  // data =  fe.get_cell_data();
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << ": ";
  //   for (size_t v=0; v<verts.size(); ++v)
  //     std::cout << data.points[q].values[v] << " ";
  //   std::cout << std::endl;
  // }
  // std::cout << "numeric" << std::endl;
  // data =  de.get_cell_data();
  // for (size_t q=0; q<data.points.size(); ++q)
  // {
  //   std::cout << q << ": ";
  //   for (size_t v=0; v<verts.size(); ++v)
  //     std::cout << data.points[q].values[v] << " ";
  //   std::cout << std::endl;
  // }

  // This proves partition of unity
  std::cout << "======= partition of unity ===========" << std::endl;
  data =  de.get_cell_data();
  {
    double maxsum = 0;
    for (auto & qp : data.points)
    {
      double sum = 0;
      for (size_t v=0; v<qp.values.size(); ++v)
      {
        sum += qp.values[v];
      }
      maxsum = std::max(sum, maxsum);
    }
    std::cout << "numeric sum = " << std::scientific << std::fabs(maxsum - 1.0) << std::endl << std::defaultfloat;
  }

  std::cout << "======= patch test ========" << std::endl;
  {
    Point maxsum = {0,0,0};
    for (auto & qp : data.points)
    {
      for (size_t i=0; i<3; ++i)
      {
        double sum = 0;
        for (size_t v=0; v<qp.values.size(); ++v)
          sum += qp.grads[v][i];
        maxsum[i] = std::max(maxsum[i], std::fabs(sum));
      }

    }
    std::cout << "grad maxsum = " << maxsum << std::endl;
  }

  exit(0);
}

std::unique_ptr<FiniteElementBase> DiscretizationFEM::build_element(const mesh::Cell & cell)
{
  std::unique_ptr<FiniteElementBase> p_discr;
  if (_config.method == FEMMethod::polyhedral_finite_element)
  {
    if (_config.solver == direct || _config.solver == cg)
      p_discr = std::make_unique<PolyhedralElementDirect>(cell, _config);
    else if (_config.solver == msrsb)
    {
      throw std::invalid_argument("regression");
      /* p_discr = std::make_unique<DFEMElement>(*cell, _msrsb_tol); */
    }
  }
  else if (_config.method == strong_discontinuity)
    p_discr = std::make_unique<StandardFiniteElement>(cell);
  else if (_config.method == mixed)
  {
    if (cell.vtk_id() == angem::GeneralPolygonID)
    {
      if (_config.solver == direct || _config.solver == cg)
        p_discr = std::make_unique<PolyhedralElementDirect>(cell, _config);
      else if (_config.solver == msrsb)
      {
        throw std::invalid_argument("regression");
        /* p_discr = std::make_unique<DFEMElement>(*cell, _msrsb_tol); */
      }
    }
    else
      p_discr = std::make_unique<StandardFiniteElement>(cell);
  }

  return p_discr;
}


#else
void DiscretizationFEM::build() {}
#endif

}  // end namepsace discretization
