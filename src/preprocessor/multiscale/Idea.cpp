#include "Idea.hpp"
#include "mesh/io/VTKWriter.hpp"    // debugging, provides io::VTKWriter
#include <Eigen/SparseLU>
#include <string>
#include <fstream>
#include <algorithm>  // std::swap

namespace multiscale {

using discretization::ConnectionData;
using discretization::ControlVolumeData;

Idea::Idea(mesh::Mesh const & grid, gprs_data::SimData & data)
    : _grid(grid), _data(data),
      _source(find_center_cell_()),
      _boundary_cells(find_boundary_cells_())
{
  auto [A, rhs] = build_system_();
  A.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);
  solver.analyzePattern(A);
  solver.factorize(A);

  // if (_grid.n_active_cells() < 100)
  // {
  //   std::cout << "matrix: " << std::endl;
  //   std::cout << A << std::endl;
  // }

  if (solver.info() ==  Eigen::ComputationInfo::Success)
  {
    auto sol = solver.solve(rhs);
    size_t n = sol.size();
    std::vector<double> soln(n);
    for (size_t i = 0; i < n; ++i)
      soln[i] = sol[i];

    std::string fname = "output/sf.vtk";
    std::cout << "saving " << fname << std::endl;
    std::ofstream out;
    out.open(fname.c_str());
    mesh::IO::VTKWriter::write_geometry(_grid, out);
    mesh::IO::VTKWriter::enter_section_cell_data(_grid.n_active_cells(), out);
    mesh::IO::VTKWriter::add_data(soln, "solution", out);
    out.close();
  } else {
    std::cout << "matrix: " << std::endl;
    std::cout << A << std::endl;
    throw std::runtime_error("Could not invert your fucking matrix");
  }
}

size_t Idea::find_center_cell_() const
{
  angem::Point<3,double> c;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
    c += cell->center();
  c /= _grid.n_active_cells();

  size_t s = 0;
  double mindist = std::numeric_limits<double>::max();
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    double const dist = cell->center().distance(c);
    if (dist < mindist)
    {
      mindist = dist;
      s = cell->index();
    }
  }
  return s;
}

std::vector<size_t> Idea::find_boundary_cells_() const
{
  std::vector<size_t> ans;
  for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell)
  {
    if ( cell->neighbors().size() < 4 )  // cartesian only
      ans.push_back(cell->index());
  }
  return ans;
}


std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> Idea::build_system_() const
{
  size_t n = _grid.n_active_cells();
    auto A = Eigen::SparseMatrix<double,Eigen::RowMajor>(n, n);
  Eigen::VectorXd R = Eigen::VectorXd::Zero(n);

  build_laplace_terms_(A, R);
  build_special_terms_(A, R);
  impose_bc_(A, R);

  return std::tie(A, R);
}

void Idea::build_laplace_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd & res) const
{
  auto const & cons =  _data.flow_connection_data;

  for (auto const & con : cons)
  {
    size_t const i = con.elements[0];
    size_t const j = con.elements[1];
    A.coeffRef(i, i) += con.coefficients[0];
    A.coeffRef(i, j) -= con.coefficients[0];
    A.coeffRef(j, i) += con.coefficients[1];
    A.coeffRef(j, j) -= con.coefficients[1];
  }
}

double compute(ConnectionData const & con, std::vector<ControlVolumeData> const & cv, angem::Point<3,double> const &c,
               size_t i)
{
  size_t const dofi = con.elements[i];

  // spherical coordinates radius variable
  auto r = cv[dofi].center - c;
  double const r_value = r.norm();
  if (!std::isnormal(r_value)) return 0.f;
  r /= r_value;

  // compute coefficient
  auto const d = con.center - cv[dofi].center;
  // return -2.f / r_value * std::fabs(con.coefficients[i]) * d.dot(r);
  double ans = +2.f / r_value * std::fabs(con.coefficients[i]) * d.dot(r);

    // if (std::fabs(entry1) < 1e-5)
  // std::cout << std::endl;
  // std::cout << "entry1 = " << ans << " (" << con.elements[0] << ")" << std::endl;
  // std::cout << "d*r = " << d.dot(r) << std::endl;
  // std::cout << "d = " << d << std::endl;
  // std::cout << "r = " << r << std::endl;
  return ans;
}

void Idea::build_special_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd & res) const
{
  auto const & cons =  _data.flow_connection_data;
  auto const & cv = _data.cv_data;
  auto const c = _grid.cell(_source).center();

  for (auto const & con : cons)
  {
    double entry1 = compute(con, _data.cv_data, c, 0);
    // if (std::fabs(entry1) < 1e-5)
    //   std::cout << "entry1 " << entry1 << " (" << con.elements[0] << ")" << std::endl;
    A.coeffRef(con.elements[0], con.elements[0]) += entry1;
    A.coeffRef(con.elements[0], con.elements[1]) -= entry1;

    double entry2 = compute(con, _data.cv_data, c, 1);
    A.coeffRef(con.elements[1], con.elements[1]) += entry2;
    A.coeffRef(con.elements[1], con.elements[0]) -= entry2;
  }

}

void Idea::impose_bc_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const
{
  // boundary
  for (auto c : _boundary_cells)
  {
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, c); it; ++it)
      it.valueRef() = (it.row() == it.col()) ? 1.f : 0.f;

    rhs[c] = 0;
  }

  // source
  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, _source); it; ++it)
    it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
  rhs[_source] = 1;
}

}  // end namespace multiscale
