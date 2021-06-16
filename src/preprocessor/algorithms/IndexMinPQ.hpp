#pragma once
#include <vector>
#include <cstddef>  // size_t without std::

namespace algorithms {

/*  Priority queue with capability to decrease key priority
 *  This structure is key to several graph-based alsorithms,
 *  so we use a graph-based terminology here:
 *  index is the number of the vertex
 *  key is the priority (edge weight)
 *  heap-position is the index of the key in the bindary heap structure
 */
class IndexMinPQ {
 public:
  IndexMinPQ(const size_t n_keys);                      // we must know the number of vertices before
  virtual ~IndexMinPQ() = default;

  void enqueue(const size_t i, const double priority);  // associate key with index i
  void setPriority(const size_t i, double new_priority); // change priority of vertex i
  bool contains(size_t i) const;                        // is i an index on the priority queue?
  size_t dequeue();                                     // remove a minimal key and return its associated index
  bool empty() const noexcept {return size() == 0;}     // is the PQ empty
  size_t size() const noexcept {return _pq.size();}     // number of entries in the PQ
  void print() const;

 protected:
  void bubble_up_(const size_t heap_idx);
  void trickle_down_(const size_t heap_idx);
  bool has_parent_(const size_t heap_idx) const noexcept {return heap_idx > 0;}
  bool has_index_(const size_t heap_idx) const noexcept {return heap_idx < size();}
  size_t parent_(const size_t heap_idx) const noexcept {return (heap_idx-1) / 2;}
  std::vector<size_t> children_(const size_t heap_idx) const;
  bool compare_(const size_t hidx1, const size_t hidx2)
  {
    return _priorities[_pq[hidx1]] < _priorities[_pq[hidx2]];
  }

  // input : heap index 1, heap index 2
  void swap_(const size_t i, const size_t j)
  {
    const size_t tmp = _pq[i];
    _pq[i] = _pq[j];
    _pq[j] = tmp;
  }

 private:
  std::vector<double> _priorities;  // priority of index i (vertex)
  std::vector<size_t> _pq;          // index of priority in heap position i (index of vertex in heap position)
  std::vector<size_t> _qp;          // heap position of priority with index i
};

}  // end namespace algorithms
