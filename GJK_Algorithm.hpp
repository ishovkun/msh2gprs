#pragma once

#include "Point.hpp"
#include "Polygon.hpp"
#include "utils.hpp" // farthest point in the direction

namespace angem
{

template<typename Scalar>
class GJK_Algorithm
{
public:
  GJK_Algorithm();
  ~GJK_Algorithm();

  bool collision(const std::vector<Point<3,Scalar>> &,
                 const std::vector<Point<3,Scalar>> &);

private:

  Point<3,Scalar> support(const std::vector<Point<3,Scalar>> &,
                          const std::vector<Point<3,Scalar>> &,
                          const Point<3,Scalar> &) const;

  bool ContainsOrigin(Point<3,Scalar> &);
  bool line(Point<3,Scalar> &);
  bool triangle(Point<3,Scalar> &);
  bool tetrahedron(Point<3,Scalar> &);

  bool checkTetrahedron(const Point<3,Scalar> &,
                        const Point<3,Scalar> &,
                        const Point<3,Scalar> &,
                        const Point<3,Scalar> &,
                        Point<3,Scalar> & dir);

  Point<3,Scalar> a, b, c, d;
  int nrPointsSimplex = 0;
};


template<typename Scalar>
GJK_Algorithm<Scalar>::GJK_Algorithm()
{}


template<typename Scalar>
GJK_Algorithm<Scalar>::~GJK_Algorithm()
{}


template<typename Scalar>
bool GJK_Algorithm<Scalar>::collision(const std::vector<Point<3,Scalar>> & box1,
                                      const std::vector<Point<3,Scalar>> & box2)
{
  Point<3,Scalar> dir = {1, 1, 1};

  c = support(box1, box2, dir);


  dir = -c; //negative direction

  b = support(box1, box2, dir);

  if (b.dot(dir) < 0)
  {
    return false;
  }

  // dir = vector3::doubleCross(c - b, -b);
  dir = cross_product( cross_product(c-b, -b), c-b );

  nrPointsSimplex = 2; //begin with 2 points in simplex

  int steps = 0;//avoid infinite loop
  while (steps<50)
  {
    a = support(box1, box2, dir);
    if (a.dot(dir) < 0)
    {
      return false;
    }
    else
    {
      if (ContainsOrigin(dir))
      {
        return true;
      }
    }
    steps++;

  }

  return false;
}


template<typename Scalar>
bool GJK_Algorithm<Scalar>::ContainsOrigin(Point<3,Scalar> & dir)
{
  /*if (n == 1)
    {
		return line(a, dir);
    }*/

  if (nrPointsSimplex == 2)
  {
		return triangle(dir);
  }
  else if (nrPointsSimplex == 3)
  {
		return tetrahedron(dir);
  }

  return false;
}


template<typename Scalar>
bool GJK_Algorithm<Scalar>::line(Point<3,Scalar> & dir)
{
  Point<3,Scalar> ab = b - a;
  Point<3,Scalar> ao = -a;//vector3::zero() - a

  //can t be behind b;

  //new direction towards a0
  // dir = vector3::doubleCross(ab, ao);
  dir = ab.cross(ao).cross(ab);

  c = b;
  b = a;
  nrPointsSimplex = 2;

  return false;
}


template<typename Scalar>
bool GJK_Algorithm<Scalar>::triangle(Point<3,Scalar> & dir)
{
  Point<3,Scalar> ao = -a;
  Point<3,Scalar> ab = b - a;
  Point<3,Scalar> ac = c - a;
  Point<3,Scalar> abc = cross_product(ab, ac);

  //point is can't be behind/in the direction of B,C or BC


  Point<3,Scalar> ab_abc = cross_product(ab, abc);
  // is the origin away from ab edge? in the same plane
  //if a0 is in that direction than
  if (ab_abc.dot(ao) > 0)
  {
		//change points
		c = b;
		b = a;

		//dir is not ab_abc because it's not point towards the origin
		// dir = doubleCross(ab,ao);
		dir = ab.cross(ao).cross(ab);

		//direction change; can't build tetrahedron
		return false;
  }


  // Point<3,Scalar> abc_ac = Point<3,Scalar>::cross(abc, ac);
  Point<3,Scalar> abc_ac = abc.cross(ac);

  // is the origin away from ac edge? or it is in abc?
  //if a0 is in that direction than
  if (abc_ac.dot(ao) > 0)
  {
		//keep c the same
		b = a;

		//dir is not abc_ac because it's not point towards the origin
		// dir = Point<3,Scalar>::doubleCross(ac, ao);
		dir = ac.cross(ao).cross(ac);

		//direction change; can't build tetrahedron
		return false;
  }

  //now can build tetrahedron; check if it's above or below
  if (abc.dot(ao) > 0)
  {
    //base of tetrahedron
    d = c;
    c = b;
    b = a;

    //new direction
    dir = abc;
  }
  else
  {
    //upside down tetrahedron
    d = b;
    b = a;
    dir = -abc;
  }

  nrPointsSimplex = 3;

  return false;
}


template<typename Scalar>
bool GJK_Algorithm<Scalar>::tetrahedron(Point<3,Scalar>& dir)
{
  Point<3,Scalar> ao = -a;//0-a
  Point<3,Scalar> ab = b - a;
  Point<3,Scalar> ac = c - a;

  //build abc triangle
  // Point<3,Scalar> abc = Point<3,Scalar>::cross(ab, ac);
  Point<3,Scalar> abc = ab.cross(ac);

  //CASE 1
  if (abc.dot(ao) > 0)
  {
    //in front of triangle ABC
    //we don't have to change the ao,ab,ac,abc meanings
    checkTetrahedron(ao,ab,ac,abc,dir);
  }


  //CASE 2:

  Point<3,Scalar> ad = d - a;

  //build acd triangle
  // Point<3,Scalar> acd = Point<3,Scalar>::cross(ac, ad);
  Point<3,Scalar> acd = ac.cross(ad);

  //same direaction with ao
  if (acd.dot(ao) > 0)
  {

    //in front of triangle ACD
    b = c;
    c = d;
    ab = ac;
    ac = ad;
    abc = acd;

    checkTetrahedron(ao, ab, ac, abc, dir);
  }

  //build adb triangle
  // Point<3,Scalar> adb = Point<3,Scalar>::cross(ad, ab);
  Point<3,Scalar> adb = ad.cross(ab);

  //case 3:

  //same direaction with ao
  if (adb.dot(ao) > 0)
  {

    //in front of triangle ADB

    c = b;
    b = d;

    ac = ab;
    ab = ad;

    abc = adb;
    checkTetrahedron(ao, ab, ac, abc, dir);
  }


  //origin in tetrahedron
  return true;

}

template<typename Scalar>
bool GJK_Algorithm<Scalar>::checkTetrahedron(const Point<3,Scalar>& ao,
                                             const Point<3,Scalar>& ab,
                                             const Point<3,Scalar>& ac,
                                             const Point<3,Scalar>& abc,
                                             Point<3,Scalar>& dir)
{

	//almost the same like triangle checks
	// Point<3,Scalar> ab_abc = Point<3,Scalar>::cross(ab, abc);
	Point<3,Scalar> ab_abc = ab.cross(abc);

  if (ab_abc.dot(ao) > 0)
  {
    c = b;
    b = a;

    //dir is not ab_abc because it's not point towards the origin;
    //ABxA0xAB direction we are looking for
    // dir = Point<3,Scalar>::doubleCross(ab, ao);
    dir = ab.cross(ao).cross(ab);

    //build new triangle
    // d will be lost
    nrPointsSimplex = 2;

    return false;
  }

  Point<3,Scalar> acp = cross_product(abc, ac);

  if (acp.dot(ao) > 0)
  {
    b = a;

    //dir is not abc_ac because it's not point towards the origin;
    //ACxA0xAC direction we are looking for
    // dir = Point<3,Scalar>::doubleCross(ac, ao);
    dir = ac.cross(ao).cross(ac);

    //build new triangle
    // d will be lost
    nrPointsSimplex = 2;

    return false;
  }

  //build new tetrahedron with new base
  d = c;
  c = b;
  b = a;

  dir = abc;

  nrPointsSimplex = 3;

  return false;
}


template<typename Scalar>
Point<3,Scalar> GJK_Algorithm<Scalar>::support(const std::vector<Point<3,Scalar>> & a,
                                               const std::vector<Point<3,Scalar>> & b,
                                               const Point<3,Scalar> & dir) const
{
  Point<3,Scalar> p1 = fathest_point_in_direction(a, dir);
  Point<3,Scalar> p2 = fathest_point_in_direction(b, -dir);

  Point<3,Scalar> p3 = p1 - p2;

	return  p3;
}


}
