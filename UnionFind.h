#include <iostream>
#include <cstdlib>
#include <vector>
#include <assert.h>

#ifndef UNION_FIND_DATASTRUCTURE
#define UNION_FIND_DATASTRUCTURE
class UnionFind {
public:
  //FIXME: code optimization?
  //FIXME: add tests
  //FIXME: unsigned int
  /**
   * General constructor.
   *
   * @param[in] n Number of elements.
   * @return Object of type UnionFind.
   */
  UnionFind(int n) {
    assert(n >= 2);

    this->n = n;
    this->nsets = n;
    // vector is zero based, i.e., we add a dummy element here
    // for simplification reasons
    this->root.reserve(n);
    this->size.reserve(n);
    for (int i = 1; i <= n; ++i) {
      this->root[i] = i;
      this->size[i] = 1;
    }
  }

  /**
   * Get root element / representative.
   *
   * @param[in] i Element.
   * @return Root element.
   */
  int getRoot(int i) const {
    assert(i >= 1 & i <= this->n);

    //FIXME: path compression, i.e., if i is deeper in the tree,
    // propagate parent
    while (i != root[i]) {
      i = root[i];
    }
    return i;
  }

  /**
   * Check if two elements are in the same set.
   *
   * @param i, j Elements.
   * @return Boolean indicating whether i and j are in the same set.
   */
  bool find(int i, int j) const {
    return getRoot(i) == getRoot(j);
  }

  /**
   * Set-union operation.
   *
   * @param i, j Set elements.
   */
  void unite(int i, int j) {
    assert(i >= 1 & i <= this->n);
    assert(j >= 1 & j <= this->n);

    int root_i = getRoot(i);
    int root_j = getRoot(j);
    if (root_i != root_j) {
      //FIXME: grow less high tree
      root[root_j] = i;
      size[i] += size[root_j];
      this->nsets -= 1;
    }
  }

  /// Printer for debugging.
  void print() const {
    std::cout << "UnionFind DS: n = " << this->n << std::endl;
    std::vector<bool> output(this->n);
    for (int i = 1; i <= this->n; ++i) {
      if (output[i]) {
        continue;
      }
      for (int j = 1; j <= this->n; ++j) {
        if (getRoot(j) == i) {
          std::cout << j << ", ";
          output[j] = true;
        }
      }
      std::cout << std::endl;
      output[i] = true;
    }
  }

private:
  /// number of elements initially added
  int n;
  /// number of sets
  int nsets;
  /// pointers to root elements
  std::vector<int> root;
  /// sizes of the set of each element
  std::vector<int> size;
};
#endif
