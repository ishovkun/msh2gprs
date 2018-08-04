#include <iostream>
#include "Polyhedra.h"


Polyhedra::Polyhedra()
{
}
Polyhedra::Polyhedra(vector<Vector> list)
{
   pointList=list;
    //ctor
}

Polyhedra::~Polyhedra()
{
    //dtor
}
Vector Polyhedra::getSupport(Vector dir)
{
	//std::cout<<"Inside Polyhedra GetSupport"<<std::endl;
	return Vector::maxDot(pointList,dir); 
}