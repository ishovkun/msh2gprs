#include "Circle.h"
#include <cmath>
#include <iostream>
using namespace std;

Circle::Circle(void)
{
	radius=0.0;
}
Circle::Circle(double r,Vector c,Vector n)
{
	radius=r;
	centre=c;
	normal=n;
}
Vector Circle::getSupport(Vector dir)
{
	//cout<<"Input Circle Dir="<<dir.i<<" "<<dir.j<<" "<<dir.k<<endl;
	Vector crossVec=Vector::cross(normal,Vector::cross(dir,normal));
	if(crossVec.norm()==0)
		return centre;
	double crossMod=crossVec.norm();
	Vector newVec;
	newVec.i=(crossVec.i/crossMod)*radius;
	newVec.j=(crossVec.j/crossMod)*radius;
	newVec.k=(crossVec.k/crossMod)*radius;
	Vector support;
	support=newVec+centre;
	//cout<<"Circle Dir="<<newVec.i<<" "<<newVec.j<<" "<<newVec.k<<endl;
	//cout <<"Circle Support="<<support.i<<" "<<support.j<<" "<<support.k<<endl;
	return support;
	
}
Circle::~Circle(void)
{
}
