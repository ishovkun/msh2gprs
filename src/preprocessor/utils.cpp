#include "utils.hpp"
#include <cassert>

namespace gprs_data {

void invert_matrix_3x3(const std::vector<double> &A, std::vector<double> &A_inv)
{
  assert(A.size() == 9);
  A_inv.resize(9);
  const double buffer48 = A[4] * A[8];
  const double buffer18 = A[1] * A[8];
  const double buffer61 = A[6] * A[1];
  const double buffer62 = A[6] * A[2];
  const double buffer32 = A[3] * A[2];
  const double buffer05 = A[0] * A[5];
  const double buffer =  (A[5] * buffer61 - buffer62 * A[4] -
                          A[3] * buffer18 + buffer32 * A[7] +
                          A[0] * buffer48 - buffer05 * A[7] );
  A_inv[0] =   (buffer48 - A[5] * A[7]) / buffer;
  A_inv[1] = - (buffer18 - A[2] * A[7]) / buffer;
  A_inv[2] = (A[1] * A[5] - A[2] * A[4]) / buffer;
  A_inv[3] = - (-A[6] * A[5] + A[3] * A[8]) / buffer;
  A_inv[4] = (-buffer62 + A[0] * A[8]) / buffer;
  A_inv[5] = - (-buffer32 + buffer05) / buffer;
  A_inv[6] = (-A[6] * A[4] + A[3] * A[7]) / buffer;
  A_inv[7] = - (-buffer61 + A[0] * A[7]) / buffer;
  A_inv[8] = (-A[3] * A[1] + A[0] * A[4]) / buffer;
}

std::vector<double> transpose3x3(const std::vector<double> &mat)
{
  std::vector<double> result(9, 0.0);
  result[0] = mat[0];
  result[1] = mat[3];
  result[2] = mat[6];
  result[3] = mat[1];
  result[4] = mat[4];
  result[5] = mat[7];
  result[6] = mat[2];
  result[7] = mat[5];
  result[8] = mat[8];
  return result;
}

double determinant_3x3(const std::vector<double> &m)
{
 assert( m.size() == 9 );
 return  -m[0] * m[5] * m[7] +
          m[0] * m[8] * m[4] +
          m[5] * m[1] * m[6] -
          m[8] * m[1] * m[3] -
          m[2] * m[6] * m[4] +
          m[2] * m[3] * m[7];
}

}  // end namespace gprs_data
