#include "Vector.h"
#include <iostream>
#include <cmath>
using namespace std;

Vector& Vector::operator=(const Vector& rhs)
{
	i=const_cast<Vector&>(rhs).Geti ();
	j=const_cast<Vector&>(rhs).Getj();
	k=const_cast<Vector&>(rhs).Getk();
  return *this;
}
Vector::Vector()
{
	i=0.0;
    j=0.0;
    k=0.0;
}
Vector::Vector(double xDir,double yDir,double zDir)
{
    i=xDir;
    j=yDir;
    k=zDir;
    //ctor
}
Vector Vector::operator +(Vector vec)
{
	Vector temp;
	temp.i=i+vec.i;
	temp.j=j+vec.j;
	temp.k=k+vec.k;
	return temp;
}
Vector Vector::operator -(Vector vec)
{
	Vector temp(0.0,0.0,0.0);
	temp.i=i-vec.i;
	temp.j=j-vec.j;
	temp.k=k-vec.k;
	return temp;
}

double Vector::dot(Vector vec1,Vector vec2)
{
    double result;
    result= vec1.i*vec2.i+vec1.j*vec2.j+vec1.k*vec2.k;
    return result;
}
Vector Vector::cross(Vector vec1,Vector vec2)
{
    Vector result(0.0,0.0,0.0);
    result.i=vec1.j*vec2.k-vec1.k*vec2.j;
    result.j=vec1.k*vec2.i-vec1.i*vec2.k;
    result.k=vec1.i*vec2.j-vec1.j*vec2.i;
    return result;
}

Vector Vector::maxDot(vector<Vector> list,Vector direction)
{
    double max=-1000000.0;
    double dotValue;
    Vector maxVec;
    for(size_t i=0;i<list.size();i++)
    {
        dotValue=dot(list[i],direction);
        if(dotValue>max)
        {
          max=dotValue;
          maxVec=list[i];
        }

    }
    return maxVec;
}
double Vector::norm()
{
	return (sqrt(this->i*this->i+this->j*this->j+this->k*this->k));
}
Vector::~Vector()
{

}
