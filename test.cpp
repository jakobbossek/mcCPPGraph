#include "graph.h"
#include "UnionFind.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int main() {
  Graph g(5, 2, false);
  g.print();
  g.addEdge(1, 2, 10, 10);
  g.addEdge(1, 3, 12, 8);
  g.addEdge(1, 4, 2, 3);
  g.addEdge(1, 5, 6, 1);
  g.addEdge(2, 3, 9, 15);
  g.addEdge(4, 5, 20, 25);
  g.print(true);
  // g.removeEdge(1, 2);
  // g.print(true);

  cout << "SPANNING TREE TEST: " << endl;

  Graph mst = g.getMSTKruskal(1);
  mst.print(true);

  // UnionFind UF(10);
  // UF.print();
  // UF.unite(2, 8);
  // UF.print();
  // UF.unite(1, 3);
  // UF.print();
  // UF.unite(1, 8);
  // UF.print();

  exit(0);
}
