#include<PointSet.hpp>

namespace angem
{

template<>
PointSet<3,double>::PointSet(const double tol)
    :
    tol(tol)
{
  __int128 nx = static_cast<__int128>
      (  std::cbrt( std::numeric_limits<__int128>::max() )  );

  upper = {tol*nx / 2, tol*nx / 2, tol*nx / 2};
  lower = -upper;
}

}
