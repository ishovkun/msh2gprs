#include <vector>

namespace gprs_data
{

void invert_matrix_3x3(const std::vector<double> &A, std::vector<double> &A_inv);

double determinant_3x3(const std::vector<double> &m);


}  // end namespace gprs_data
