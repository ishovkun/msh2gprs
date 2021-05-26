#include "ShapeFunctionSolver.hpp"

namespace multiscale {

using discretization::ConnectionData;
using discretization::ControlVolumeData;

ShapeFunctionSolver::ShapeFunctionSolver(size_t source, std::vector<size_t> const & bnd_cells,
                                         mesh::Mesh const & grid, gprs_data::SimData & data)
    :
    _source(source), _boundary_cells(bnd_cells),
    _grid(grid), _data(data)
{
  auto [A, rhs] = build_system_();
  A.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);

  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Lower|Eigen::Upper> solver;
  // solver.setMaxIterations(200);
  // solver.setTolerance(1e-10);

  solver.analyzePattern(A);
  solver.factorize(A);

  auto sol = solver.solve(rhs);
  if (solver.info() ==  Eigen::ComputationInfo::Success)
    _solution.resize(sol.size(), 0.f);
  else throw std::runtime_error("Could not solve system");

  for (size_t i = 0; i < sol.size(); ++i)
    _solution[i] = sol[i];
}

std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> ShapeFunctionSolver::build_system_() const
{
  size_t n = _grid.n_active_cells();
    auto A = Eigen::SparseMatrix<double,Eigen::RowMajor>(n, n);
  Eigen::VectorXd R = Eigen::VectorXd::Zero(n);

  build_laplace_terms_(A, R);
  build_special_terms_(A, R);
  impose_bc_(A, R);

  return std::tie(A, R);
}

void ShapeFunctionSolver::build_laplace_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd & res) const
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

double compute(ConnectionData const & con, std::vector<ControlVolumeData> const & cv,
               angem::Point<3,double> const &c, size_t i)
{
  // global index
  size_t const dofi = con.elements[i];

  // spherical coordinates radius variable
  auto r = cv[dofi].center - c;
  double const r_value = r.norm();
  if (!std::isnormal(r_value)) return 0.f;
  r /= r_value;

  // compute coefficient
  auto const d = con.center - cv[dofi].center;
  return 2.f / r_value * std::fabs(con.coefficients[i]) * d.dot(r);
}

void ShapeFunctionSolver::build_special_terms_(Eigen::SparseMatrix<double,Eigen::RowMajor>& A, Eigen::VectorXd & res) const
{
  auto const & cons =  _data.flow_connection_data;
  auto const & cv = _data.cv_data;
  auto const c = _grid.cell(_source).center();

  for (auto const & con : cons)
  {
    double entry1 = compute(con, _data.cv_data, c, 0);
    A.coeffRef(con.elements[0], con.elements[0]) += entry1;
    A.coeffRef(con.elements[0], con.elements[1]) -= entry1;

    double entry2 = compute(con, _data.cv_data, c, 1);
    A.coeffRef(con.elements[1], con.elements[1]) += entry2;
    A.coeffRef(con.elements[1], con.elements[0]) -= entry2;
  }

}

void impose_bc(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs, size_t dof, double value)
{
  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it) {
    // it.valueRef() = (it.row() == it.col()) ? 1.f : 0.f;

      if (it.row() == it.col()) it.valueRef() = 1.f;
      else {
        it.valueRef() = 0.f;
        double const entry = mat.coeffRef(it.col(), dof);
        mat.coeffRef(it.col(), dof) = 0.f;
        rhs[it.col()] -= entry * value;
      }

    rhs[dof] = value;
  }
}

void ShapeFunctionSolver::impose_bc_(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs) const
{
  for (auto c : _boundary_cells)
    impose_bc(mat, rhs, c, 0.f);

  impose_bc(mat, rhs, _source, 1.f);
}


}  // end namespace multiscale
