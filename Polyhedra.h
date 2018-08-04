#ifndef POLYHEDRA_H
#define POLYHEDRA_H

#include "Vector.h"
class Polyhedra
{
    public:
		vector<Vector> pointList;
		Polyhedra();
        Polyhedra(vector<Vector> list);
        virtual ~Polyhedra();
		virtual Vector getSupport(Vector);
    protected:
    private:
};

#endif // POLYHEDRA_H
