#include <iostream>
#include "Collision.h"

using namespace std;

Collision::Collision(Polyhedra* one,Polyhedra* two)
{
	a = one;
	b = two;
	const size_t size = 4;
	simplex.reserve(size);
}


Vector Collision::support(Vector dir)
{
	Vector support;
	Vector negDir = origin - dir;
	Vector aDot = a->getSupport(dir);
	Vector bDot = b->getSupport(negDir);
	support = aDot - bDot;
	return support;
}


Vector Collision::updateDir()
{
	Vector dir,AB,AC,AD,AO,tPlane;
	vector<Vector> temp;
	temp.reserve(3);
	size_t n=simplex.size();

	switch(n)
	{
	case 2:
		AB=(simplex.at(0)-simplex.at(1));
		AO=origin-simplex.at(1);
		dir=Vector::cross(Vector::cross(AB,AO),AB);
		break;
	case 3:
		AB = simplex.at(1) - simplex.at(2);
		AC = simplex.at(0) - simplex.at(2);
		AO = origin - simplex.at(2);
		tPlane = Vector::cross(AC,AB);
		if(Vector::dot(Vector::cross(AC,tPlane), AO) > 0)
		{
			dir = Vector::cross(AC,Vector::cross(AO,AC));
			simplex.erase(simplex.begin()+1);
		}
		else if(Vector::dot(Vector::cross(tPlane,AB),AO)<=0 && Vector::dot(Vector::cross(tPlane,AB),AO)<=0)
				if(Vector::dot(tPlane,AO)>0)
					dir = tPlane;
				else if(Vector::dot(tPlane,AO) == 0)
					break;
				else
					dir = origin - tPlane;
		else
		{
			dir = Vector::cross(AB,Vector::cross(AO,AB));
			simplex.erase(simplex.begin());
		}
		break;
	case 4:
		Vector A = simplex.at(3);
		Vector B = simplex.at(2);
		Vector C = simplex.at(1);
		Vector D = simplex.at(0);
		Vector AD = D-A;
		Vector AC = C-A;
		Vector AB = B-A;
		Vector normADC = Vector::cross(AD,AC);
		Vector normACB = Vector::cross(AC,AB);
		Vector normABD = Vector::cross(AB,AD);
		Vector AO = origin-A;
		if(Vector::dot(normADC, AO)>0)
		{
			dir = normADC;
			simplex.erase(simplex.begin() + 2);
		}
		else if(Vector::dot(normACB,AO) > 0)
		{
			dir = normACB;
			simplex.erase(simplex.begin());
		}
		else if(Vector::dot(normABD,AO) > 0)
		{
			dir = normABD;
			simplex.erase(simplex.begin() + 1);
		}
		else
			break;
	}
	return dir;
}


bool Collision::checkCollision()
{
	Vector d(1,-1,-1);
	Vector s=Collision::support(d);
	simplex.push_back(s);
	d=(origin-d);
	while(true)
	{
		s = Collision::support(d);
		// cout<<" Dir="<<d.i<<" "<<d.j<<" "<<d.k<<endl;
		// cout <<" Support="<<s.i<<" "<<s.j<<" "<<s.k<<endl;
		if(Vector::dot(d,s) < 0)
		{
			// cout<<"No Intersection"<<endl;
			return false;
			break;
		}

		simplex.push_back(s);
		d = Collision::updateDir();
		if(d.i==0 and d.j==0 and d.k==0)
		{
			// cout<<"Intersection"<<endl;
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
}
Collision::~Collision(void)
{
}
