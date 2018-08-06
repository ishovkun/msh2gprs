/*
 * This code was shamelessly stolen from
 * https://github.com/ruddradev/gjk_cplusplus
 * and adapted for the local infrastructure
 * Author: Igor Shovkun <ishovkun@stanford.edu>
 */

#pragma once

#include <vector>
#include "Shape.hpp"
#include "Point.hpp"

namespace angem
{

template<typename Scalar>
class CollisionGJK
{
public:
	CollisionGJK();
	bool collide(const Shape<Scalar> & shape1,
               const Shape<Scalar> & shape2);

 private:
	Point<3,Scalar> support(Point<3,Scalar> & direction);
	Point<3,Scalar> updateDirection();

  // attributes
	const Shape<Scalar> * a;
	const Shape<Scalar> * b;
	std::vector<Point<3,Scalar>> simplex;
	Point<3,Scalar> origin;
};


template<typename Scalar>
CollisionGJK<Scalar>::CollisionGJK()
{
	// a = one;
	// b = two;
	simplex.reserve(4);
}


template<typename Scalar>
Point<3,Scalar>
CollisionGJK<Scalar>::support(Point<3,Scalar> & dir)
{
	Point<3,Scalar> support;
	Point<3,Scalar> negDir = origin-dir;
	Point<3,Scalar> aDot = a->support(dir);
	Point<3,Scalar> bDot = b->support(negDir);
	support = aDot - bDot;
	return support;
}


template<typename Scalar>
Point<3,Scalar>
CollisionGJK<Scalar>::updateDirection()
{
	Point<3,Scalar> dir,AB,AC,AD,AO,tPlane;
  std::vector<Point<3,Scalar>> temp;
	temp.reserve(3);
	size_t n = simplex.size();

	switch(n)
	{
	case 2:
		AB = (simplex.at(0)-simplex.at(1));
		AO = origin-simplex.at(1);
		// dir = Vector::cross(Vector::cross(AB, AO), AB);
		dir = (AB.cross(AO)).cross(AB);
		break;

	case 3:
		AB = simplex.at(1) - simplex.at(2);
		AC = simplex.at(0) - simplex.at(2);
		AO = origin - simplex.at(2);
		// tPlane = Vector::cross(AC,AB);
		tPlane = AC.cross(AB);
		// if(Vector::dot(Vector::cross(AC,tPlane), AO) > 0)
		if( (AC.cross(tPlane)).dot(AO) > 0 )
		{
			// dir = Vector::cross(AC,Vector::cross(AO,AC));
			dir = AC.cross(AO.cross(AC));
			simplex.erase(simplex.begin()+1);
		}
		// else if(Vector::dot(Vector::cross(tPlane,AB),AO)<=0 && Vector::dot(Vector::cross(tPlane,AB),AO)<=0)
		else if(dot( cross(tPlane, AB), AB ) <= 0 and
            dot( cross(tPlane, AB), AO ) <= 0)
    {
      // if(Vector::dot(tPlane,AO)>0)
      if(tPlane.dot(AO) > 0)
        dir = tPlane;
      // else if(Vector::dot(tPlane,AO) == 0)
      else if(tPlane.dot(AO) == static_cast<Scalar>(0))
        break;
      else
        dir = origin - tPlane;
    }
		else
		{
			// dir = Vector::cross(AB,Vector::cross(AO,AB));
			dir = AB.cross( AO.cross(AB) );
			simplex.erase(simplex.begin());
		}
		break;

	case 4:
		Point<3,Scalar> A = simplex.at(3);
		Point<3,Scalar> B = simplex.at(2);
		Point<3,Scalar> C = simplex.at(1);
		Point<3,Scalar> D = simplex.at(0);
		Point<3,Scalar> AD = D-A;
		Point<3,Scalar> AC = C-A;
		Point<3,Scalar> AB = B-A;
		// Point<3,Scalar> normADC = Vector::cross(AD,AC);
		Point<3,Scalar> normADC = AD.cross(AC);
		// Point<3,Scalar> normACB = Vector::cross(AC,AB);
		Point<3,Scalar> normACB = AC.cross(AB);
		// Point<3,Scalar> normABD = Vector::cross(AB,AD);
		Point<3,Scalar> normABD = AB.cross(AD);
		Point<3,Scalar> AO = origin - A;

		// if(Vector::dot(normADC, AO)>0)
		if (normADC.dot(AO) > 0)
		{
			dir = normADC;
			simplex.erase(simplex.begin() + 2);
		}
		// else if(Vector::dot(normACB, AO) > 0)
		else if (normACB.dot(AO) > 0)
		{
			dir = normACB;
			simplex.erase(simplex.begin());
		}
		// else if(Vector::dot(normABD,AO) > 0)
		else if (normABD.dot(AO) > 0)
		{
			dir = normABD;
			simplex.erase(simplex.begin() + 1);
		}
		else
			break;
	}
	return dir;
}


template<typename Scalar>
bool
CollisionGJK<Scalar>::collide(const Shape<Scalar> & shape1,
                              const Shape<Scalar> & shape2)
{
  a = & shape1;
  b = & shape2;

	Point<3,Scalar> d(1,-1,-1);
	Point<3,Scalar> s = support(d);
	simplex.push_back(s);
	d = origin - d;

  int iter = 0;
	while(iter < 50)
	{
    iter++;
		s = support(d);
		// if(Vector::dot(d,s) < 0)
		if(d.dot(s) < 0)
		{
			return false;
			break;
		}

		simplex.push_back(s);
		d = updateDirection();
		if(d.x() == 0 and d.y() == 0 and d.z() == 0)
		{
			return true;
			break;
		}
		/*else if(d.k==-1)
      {
			cout<<"No Intersection"<<endl;
			return false;
			break;
      }*/
	}
  return false;
}

}  // end namespace
