#include "graph.h"
#include "UnionFind.h"
#include "debug.h"
#include <iostream>
#include <cstdlib>

using namespace std;

int main() {
  Graph g1(5, 2, false);
  g1.addEdge(1, 2, 10, 10);
  g1.addEdge(1, 3, 12, 8);
  g1.addEdge(1, 4, 2, 3);
  g1.addEdge(1, 5, 6, 1);
  g1.addEdge(2, 3, 9, 15);
  g1.addEdge(4, 5, 20, 25);

  Graph g2(5, 2, false);
  g2.addEdge(1, 2, 10, 10);
  g2.addEdge(1, 3, 12, 8);
  g2.addEdge(3, 4, 10, 3);

  cout << "CONNECTED SUBTREE: " << endl;
  g1.print(true);
  std::vector<int> bfsres = g1.getConnectedSubtree(4, 4);
  for (int i: bfsres) {
    cout << i << ", " << endl;
  }
  Graph inducedg1 = g1.getInducedSubgraph(bfsres);
  std::vector<int> ss = inducedg1.getSumOfEdgeWeights();
  cout << ss[0] << ", " << ss[1] << endl;
  //Graph inducedg1mst = inducedg1.getMSTKruskal(1);
  inducedg1.print(true);

  // cout << "GRAPH INTERSECTION: " << endl;
  // g1.print(true);
  // g2.print(true);

  // Graph gintersection = Graph::getIntersectionGraph(g1, g2);
  // gintersection.print(true);

  cout << "------------------------------------------" << endl;

  cout << "SPANNING TREE TEST: " << endl;

  Graph mst = g1.getMSTKruskal(1);
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
