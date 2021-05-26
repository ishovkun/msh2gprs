#include "IndexMinPQ.hpp"
#include <limits>
#include <iostream>
#include <cassert>

namespace algorithms {

static constexpr size_t NOT_EXISTS =  std::numeric_limits<size_t>::max();

IndexMinPQ::IndexMinPQ(const size_t n_keys)
    : _priorities(n_keys, std::numeric_limits<double>::max()),
      _qp(n_keys, NOT_EXISTS)
{}

void IndexMinPQ::enqueue(const size_t i, const double priority)
{
  _priorities[i] = priority;
  const size_t heap_idx = _pq.size();
  _pq.emplace_back(i);
  _qp[i] = heap_idx;
  bubble_up_(heap_idx);
}

bool IndexMinPQ::contains(const size_t i) const
{
  if (i >= _qp.size()) return false;
  else                 return _qp[i] != NOT_EXISTS;
}

void IndexMinPQ::bubble_up_(const size_t heap_idx)
{
  if (has_parent_(heap_idx))
  {
    const size_t p = parent_(heap_idx);
    if ( compare_(heap_idx, p) )
    {
      swap_(heap_idx, p);
      bubble_up_(p);
    }
  }
  else return;
}

void IndexMinPQ::print() const
{
  for (size_t i = 0; i < size(); ++i)
    std::cout << _pq[i] << " (" << _priorities[_pq[i]] << ") ";
  std::cout << std::endl;
}

void IndexMinPQ::setPriority(const size_t i, double new_priority)
{
  assert( contains(i) && "Thid element is not in queue" );
  // const size_t heap_idx = _pq[i];
  const size_t heap_idx = _qp[i];
  const bool down = new_priority > _priorities[i];
  _priorities[i] = new_priority;
  if (down) trickle_down_(heap_idx);
  else bubble_up_(heap_idx);
}

void IndexMinPQ::trickle_down_(const size_t idx)
{
  for (const size_t child : children_(idx))
    if (compare_(child, idx))
    {
      swap_(child, idx);
      trickle_down_(child);
    }
}

std::vector<size_t> IndexMinPQ::children_(const size_t idx) const
{
  std::vector<size_t> ch;
  ch.reserve(2);
  const size_t ch1 = 2*idx + 1;
  const size_t ch2 = 2*idx + 2;
  if (has_index_(ch1)) ch.push_back(ch1);
  if (has_index_(ch2)) ch.push_back(ch2);
  return ch;
}

size_t IndexMinPQ::dequeue()
{
  const size_t root = 0;
  const size_t result = _pq[root];
  _pq[root] = _pq[size() - 1];
  _pq.pop_back();
  _qp[result] = NOT_EXISTS;
  trickle_down_(root);
  return result;
}

}  // end namespace algorithms
