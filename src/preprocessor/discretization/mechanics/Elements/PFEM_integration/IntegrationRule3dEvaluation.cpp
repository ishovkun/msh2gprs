#include "IntegrationRule3dEvaluation.hpp"
#include "../FeValues.hpp"

auto myfunc = [](angem::Point<3,double> x) -> double {
//   // return x[0] + x[1] + x[2];
//   // return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
//   // return x[0]*x[0]*x[0] + x[1]*x[1]*x[1] + x[2]*x[2]*x[2];
//   // return x[0]*x[0]*x[0]*x[0] + x[1]*x[1]*x[1]*x[1] + x[2]*x[2]*x[2]*x[2];
//   // trilinear
//   // return 0.f + x[0] + x[1] + x[2] + x[0]*x[1] + x[0]*x[2] + x[1]*x[2] + x[0]*x[1]*x[2];
//   // bicubic
  return x[0] + x[1]   // x + y
      + x[0]*x[0] + x[1]*x[1] + x[0]*x[1] // x² + y² + xy
      + x[0]*x[1]*x[2] + x[0]*x[0]*x[1]  + x[0]*x[0]*x[0] + x[1]*x[1]*x[0] + x[1]*x[1]*x[1];  // xyz + x²y + x³ + xy² + y³
};
auto f0 = [](angem::Point<3,double> x) -> double { return 1.0;};
auto f1 = [](angem::Point<3,double> x) -> double { return x[0];};
auto f2 = [](angem::Point<3,double> x) -> double { return x[1];};
auto f3 = [](angem::Point<3,double> x) -> double { return x[2];};
auto f4 = [](angem::Point<3,double> x) -> double { return x[0]*x[0];};
auto f5 = [](angem::Point<3,double> x) -> double { return x[0]*x[1];};
auto f6 = [](angem::Point<3,double> x) -> double { return x[1]*x[1];};
auto f7 = [](angem::Point<3,double> x) -> double { return x[0]*x[0]*x[0];};
auto f8 = [](angem::Point<3,double> x) -> double { return x[0]*x[0]*x[1];};
auto f9 = [](angem::Point<3,double> x) -> double { return x[0]*x[1]*x[1];};
auto f10 = [](angem::Point<3,double> x) -> double { return x[1]*x[1]*x[1];};

double eval_fem(mesh::Cell const & cell, std::function<double(angem::Point<3,double>)> f) {
  discretization::FeValues<angem::VTK_ID::HexahedronID> fe_values;
  fe_values.update(cell);
  // q-point coordinates to get function values in them
  std::vector<angem::Point<3,double>> coord(fe_values.n_integration_points());
  for (size_t v = 0; v < cell.n_vertices(); ++v) {
    for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
      coord[q] += cell.vertex_coordinates()[v] * fe_values.value(v, q);
    }
  }

  double integral = 0;
  for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
    integral += fe_values.JxW(q) * f(coord[q]);
  }
  return integral;
}

double one_point(mesh::Cell const & cell, std::function<double(angem::Point<3,double>)> f) {
  return cell.volume() * f(cell.center());
  // return cell.volume() * f(cell.vertex_coordinates()[1]);
  // return cell.volume() * f({1, -1, 0});
}

namespace discretization {

double eval(PolyhedralElementBase const & element, std::function<double(angem::Point<3,double>)> f)
{
  const size_t nv = ElementTraits<angem::TetrahedronID>::n_vertices;
  FeValues<angem::VTK_ID::TetrahedronID> fe_values;
  const auto & grid = element.get_grid();
  double integral = 0;
  for (auto cell = grid.begin_active_cells(); cell != grid.end_active_cells(); ++cell)
  {
    fe_values.update(*cell);
    integral += f(cell->center()) * cell->volume();
    // for (size_t q = 0; q < fe_values.n_integration_points(); ++q) {
    //   integral += myfunc(cell->center()) * fe_values.JxW(q);
    // }
  }
  return integral;
}

IntegrationRule3dEvaluation::
IntegrationRule3dEvaluation(PolyhedralElementBase const & element, TributaryRegion3dBase const  & tributary)
    : IntegrationRule3d(element, tributary)
{
  std::cout.precision(10);

  std::vector<std::function<double(angem::Point<3,double>)>> funcs;
  funcs.push_back(f0); funcs.push_back(f1); funcs.push_back(f2); funcs.push_back(f3);
  funcs.push_back(f4); funcs.push_back(f5); funcs.push_back(f6); funcs.push_back(f7);
  funcs.push_back(f8); funcs.push_back(f9); funcs.push_back(f10);
  // funcs.push_back(myfunc);

  // for (auto const & f : funcs)
  for (int i = 0; i < funcs.size(); ++i)
  {
    auto const & f = funcs[i];
    double integral = eval(element, f);
    // double integral_fem = eval_fem(element.host_cell(), f);
    // double integral_1p = one_point(element.host_cell(), f);
    // std::cout << i << "\t" << integral << "\t" << integral_fem  << "\t" << integral_1p << std::endl;
    std::cout << integral << std::endl;
  }
  std::cout << std::endl;
  for (int i = 0; i < funcs.size(); ++i) {
    auto const & f = funcs[i];
    double integral_fem = eval_fem(element.host_cell(), f);
    double integral_1p = one_point(element.host_cell(), f);
    // std::cout << i << "\t" << integral << "\t" << integral_fem  << "\t" << integral_1p << std::endl;
    std::cout << integral_fem  << " " << integral_1p << std::endl;
  }
  // std::cout << "\nintegral value = "  << integral << std::endl;
  // std::cout << "integral fem = "  << integral_fem << std::endl;
  // std::cout << "integral 1p = "  << integral_1p << std::endl;
  exit(0);
}

}  // end namespace discretization
