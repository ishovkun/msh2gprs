#pragma once

#include<vector>
using namespace std;
class Vector
{

public:
  double i;
  double j;
  double k;
  template <typename T, size_t N>
      inline
      size_t SizeOfArray( const T(&)[ N ] )
  {
    return N;
  }

  Vector& operator=(const Vector&);
  Vector();
  Vector(double xDir,double yDir,double zDir);
  virtual ~Vector();
  Vector operator +(Vector);
  Vector operator -(Vector);
  double Geti() { return i; }
  void Seti(double val) { i = val; }
  double Getj() { return j; }
  void Setj(double val) { j = val; }
  double Getk() { return k; }
  void Setk(double val) { k = val; }
  static double dot(Vector,Vector);
  static Vector cross(Vector,Vector);
  static Vector maxDot(vector<Vector>,Vector );
  double norm();
protected:

};
