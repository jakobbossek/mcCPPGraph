#include "graph.h"
#include "UnionFind.h"
#include "debug.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <algorithm>

using namespace std;

double sum(std::vector<double> x) {
  double sum = 0.0;
  for (auto& element : x) {
    sum += element;
  }
  return sum;
}

double mean(std::vector<double> x) {
  return sum(x) / x.size();
}


int main() {

  cout << "DEBUGGING" << endl;

  Graph g = Graph::importFromGrapheratorFile("tests/instances/graph_N50-E164-C0-W2---LHSNG---DEG-STEG-WEG---RWG-RWG_1.graph");
  g.print(false);
  Graph mst = g.getMSTKruskal(1);
  mst.print(false);
  std::vector<std::vector<int>> comps = mst.getConnectedComponents();

  cout << "Is spanning tree: " << mst.isSpanningTree() << " " << comps.size() << endl;

  Graph mst2 = g.getMSTBySubforestMutation(mst, 1);
  mst2.print(false);
  comps = mst2.getConnectedComponents();
  cout << "Is spanning tree: " << mst2.isSpanningTree() << " " << comps.size() << endl;

  Graph mst3 = g.getMSTBySubgraphMutation(mst, 10);
  mst3.print(false);
  comps = mst3.getConnectedComponents();
  cout << "Is spanning tree: " << mst3.isSpanningTree() << " " << comps.size() << endl;

  // cout << "SUBFOREST MUTATION" << endl;
  // cout << "==================" << endl;

  // double timeStart = clock();
  // //Graph g = Graph::importFromGrapheratorFile("tests/instances/graph_N500-E1497-C5-W2---LHSNG-CLUNG---CLDEG-CLSTEG---RWG-RWG_1.graph");
  // Graph g = Graph::importFromGrapheratorFile("tests/instances/graph_N250-E31125-C0-W2---UNG---CEG---RWG-RWG_1.graph");
  // g.print(false);
  // Graph mst = g.getMSTKruskal(0.5);
  // mst.print(false);

  // // small benchmark to check how long it takes to import and measure stuff
  // int nmutations = 250; //g.getV();
  // std::vector<double> runtimes(nmutations);

  // double timeImport = ((double)(clock() - timeStart) / (double)CLOCKS_PER_SEC);
  // cout << "Time for import: " << timeImport << endl;

  // for (int i = 0; i < nmutations; ++i) {
  //   cout << ".";
  //   if (i == nmutations - 1)
  //     cout << endl;
  //   timeStart = clock();

  //   Graph mst2 = mst.getMSTBySubforestMutation((int)(forest.getV() / 2));

  //   runtimes[i] = ((double)(clock() - timeStart) / (double)CLOCKS_PER_SEC);
  // }


  // cout << "Performed " << nmutations << " runs." << endl;
  // cout << "Mean running time   : " << mean(runtimes) << endl;
  // cout << "Overall running time: " << sum(runtimes) << endl;
  // cout << "Overall passed time : " << sum(runtimes) + timeImport << endl;

  // cout << "SUBTREE MUTATION" << endl;
  // cout << "================" << endl;

  // for (int i = 0; i < nmutations; ++i) {
  //   cout << ".";
  //   if (i == nmutations - 1)
  //     cout << endl;
  //   timeStart = clock();

  //   // get random start node
  //   int rndNode = (int)(((double)rand() / (double)RAND_MAX) * g.getV()) + 1;
  //   // maximal size of connected subtree
  //   int maxNodes = 20;
  //   // get subtree reached by extraction from mst
  //   std::vector<int> subtreeNodes = mst.getConnectedSubtree(rndNode, maxNodes);
  //   Graph g2 = g.getInducedSubgraph(subtreeNodes);
  //   // compute locally optimal mst on subgraph
  //   double rndWeight = (double)rand() / (double)RAND_MAX;
  //   Graph subtreeMST = g.getMSTKruskal(rndWeight, g2);

  //   runtimes[i] = ((double)(clock() - timeStart) / (double)CLOCKS_PER_SEC);
  // }

  // cout << "Performed " << nmutations << " runs." << endl;
  // cout << "Mean running time   : " << mean(runtimes) << endl;
  // cout << "Overall running time: " << sum(runtimes) << endl;
  // cout << "Overall passed time : " << sum(runtimes) + timeImport << endl;


  // g.print(false);

  // Graph g1(5, 2, false);
  // g1.addEdge(1, 2, 10, 10);
  // g1.addEdge(1, 3, 12, 8);
  // g1.addEdge(1, 4, 2, 3);
  // g1.addEdge(1, 5, 6, 1);
  // g1.addEdge(2, 3, 9, 15);
  // g1.addEdge(4, 5, 20, 25);

  // Graph g2(5, 2, false);
  // g2.addEdge(1, 2, 10, 10);
  // g2.addEdge(1, 3, 12, 8);
  // g2.addEdge(3, 4, 10, 3);

  // cout << "CONNECTED SUBTREE: " << endl;
  // g1.print(true);
  // std::vector<int> bfsres = g1.getConnectedSubtree(4, 4);
  // for (int i: bfsres) {
  //   cout << i << ", " << endl;
  // }
  // Graph inducedg1 = g1.getInducedSubgraph(bfsres);
  // std::vector<int> ss = inducedg1.getSumOfEdgeWeights();
  // cout << ss[0] << ", " << ss[1] << endl;
  // //Graph inducedg1mst = inducedg1.getMSTKruskal(1);
  // inducedg1.print(true);

  // // cout << "GRAPH INTERSECTION: " << endl;
  // // g1.print(true);
  // // g2.print(true);

  // // Graph gintersection = Graph::getIntersectionGraph(g1, g2);
  // // gintersection.print(true);

  // cout << "------------------------------------------" << endl;

  // cout << "SPANNING TREE TEST: " << endl;

  // Graph mst = g1.getMSTKruskal(1);
  // mst.print(true);

  // // UnionFind UF(10);
  // // UF.print();
  // // UF.unite(2, 8);
  // // UF.print();
  // // UF.unite(1, 3);
  // // UF.print();
  // // UF.unite(1, 8);
  // // UF.print();

  exit(0);
}
