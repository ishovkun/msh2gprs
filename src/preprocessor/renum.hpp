#include <stdio.h>
#include <deque>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <stdlib.h>
using namespace std;

// extern pair< C::iterator, C::iterator > r;
//this will help to sort according to external array
class order_comparator
{
	int * order;
public:
	order_comparator(int * _order) : order(_order) {}
	order_comparator(const order_comparator & b) : order(b.order) {}
	order_comparator & operator =(order_comparator const & b) {order = b.order; return *this;}
	bool operator()(int i, int j) {return order[i] < order[j];}
};

class renum
{
public:
  renum();
  ~renum();
  void convert(const int nc, const int nv, const vector<int>& ia, const vector<int>& ja, vector<int>& rcm);
};