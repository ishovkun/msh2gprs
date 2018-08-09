#pragma once
#include "Polyhedra.h"

class Circle: public Polyhedra
{
	Vector centre;
	Vector normal;
	double radius;
	
public:
	Circle(void);
	Circle(double r,Vector c,Vector n);
	virtual Vector getSupport(Vector);
	~Circle(void);
};
