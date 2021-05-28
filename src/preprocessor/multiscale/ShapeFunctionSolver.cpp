#include "ShapeFunctionSolver.hpp"

namespace multiscale {

using discretization::ConnectionData;
using discretization::ControlVolumeData;

ShapeFunctionSolver::ShapeFunctionSolver(size_t source,
                                         std::vector<size_t> const & sub_region,
                                         std::vector<size_t> const & bnd,
                                         std::vector<double> const & bnd_values,
                                         gprs_data::SimData const & data)
    : _source(source), _region(sub_region), _bnd(bnd), _bnd_values(bnd_values), _data(data),
      _mapping(_data.cv_data.size(), _region.size())
{
  for (size_t i = 0; i < _region.size(); ++i)
    _mapping[_region[i]] = i;

  solve_();
}

ShapeFunctionSolver::ShapeFunctionSolver(size_t source, std::vector<size_t> const & bnd_cells,
                                         gprs_data::SimData const & data)
    : _source(source), _bnd(bnd_cells), _data(data), _region(_data.cv_data.size()),
      _mapping(_data.cv_data.size(), _region.size())
{
  std::iota( _region.begin(), _region.end(), 0 );
  std::iota( _mapping.begin(), _mapping.end(), 0 );

  _bnd_values.assign(_bnd.size(), 0.f);
  _bnd.push_back(source);
  _bnd_values.push_back(1.f);

  solve_();
}

std::tuple<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::VectorXd> ShapeFunctionSolver::build_system_() const
{
  size_t n = _region.size();
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
    if (_mapping[con.elements[0]] < _region.size() && _mapping[con.elements[1]] < _region.size())
  {
    size_t const i = _mapping[ con.elements[0] ];
    size_t const j = _mapping[ con.elements[1] ];
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
  auto const c = cv[_source].center;

  for (auto const & con : cons)
    if (_mapping[con.elements[0]] < _region.size() && _mapping[con.elements[1]] < _region.size())
  {
    double const entry1 = compute(con, _data.cv_data, c, 0);
    size_t const i = _mapping[con.elements[0]];
    size_t const j = _mapping[con.elements[1]];
    A.coeffRef(i, i) += entry1;
    A.coeffRef(i, j) -= entry1;

    double const entry2 = compute(con, _data.cv_data, c, 1);
    A.coeffRef(j, j) += entry2;
    A.coeffRef(j, i) -= entry2;
  }

}

void impose_bc(Eigen::SparseMatrix<double,Eigen::RowMajor>& mat, Eigen::VectorXd & rhs, size_t dof, double value)
{
  for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, dof); it; ++it) {

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
  for (size_t i = 0; i < _bnd.size(); ++i)
    impose_bc(mat, rhs, _mapping[ _bnd[i] ], _bnd_values[i]);
}

void ShapeFunctionSolver::solve_()
{
  auto [A, rhs] = build_system_();
  A.makeCompressed();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(A);

  // Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Lower|Eigen::Upper> solver;
  // solver.setMaxIterations(200);
  // solver.setTolerance(1e-10);

  // solver.analyzePattern(A);
  // solver.factorize(A);

  auto sol = solver.solve(rhs);
  if (solver.info() ==  Eigen::ComputationInfo::Success)
    _solution.resize(sol.size(), 0.f);
  else throw std::runtime_error("Could not solve system");

  for (size_t i = 0; i < sol.size(); ++i)
    _solution[i] = sol[i];
}


}  // end namespace multiscale
