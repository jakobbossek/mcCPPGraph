#include "../graph.h"
#include "gtest/gtest.h"


namespace {

// Tests basic graph construction
TEST(GraphTest, Basic_Methods) {
  Graph g(5, 2, false);
  g.addEdge(1, 2, 10, 10);
  g.addEdge(1, 3, 12, 8);
  g.addEdge(1, 4, 2, 3);
  g.addEdge(1, 5, 6, 1);
  g.addEdge(2, 3, 9, 15);
  g.addEdge(4, 5, 20, 25);
  EXPECT_EQ(5, g.getV());
  EXPECT_EQ(2, g.getW());
  EXPECT_EQ(6, g.getE());
  EXPECT_EQ(4, g.getDegree(1));

  g.removeEdge(1, 2);
  EXPECT_EQ(6, g.getE() + 1);
  EXPECT_EQ(4, g.getDegree(1) + 1);
}

TEST(GraphTest, ConnectedComponents) {
  Graph g(5, 2, false);
  std::vector<std::vector<int>> components = g.getConnectedComponents();
  EXPECT_EQ(5, components.size());
  g.addEdge(1, 2, 1, 1);
  g.addEdge(1, 3, 1, 1);
  g.addEdge(4, 5, 1, 1);
  components = g.getConnectedComponents();
  EXPECT_EQ(2, components.size());

  // check if grapherator graphs are connected
  std::vector<std::string> files {
    "tests/instances/graph_N100-E259-C5-W2---UNG-CLUNG---CLSTEG-CLDEG-DEG---RWG-RWG_1.graph",
    "tests/instances/graph_N250-E729-C0-W2---UNG---DEG---RWG-RWG_1.graph",
    "tests/instances/graph_N50-E153-C5-W2---LHSNG-CLUNG---CLDEG-CLSTEG---RWG-RWG_1.graph"
  };

  for (auto pathToFile: files) {
    Graph g = Graph::importFromGrapheratorFile(pathToFile);
    std::vector<std::vector<int>> components = g.getConnectedComponents();
    EXPECT_EQ(1, components.size());
  }
}

TEST(MSTTest, BasicKruskal) {
  unsigned int V = 5;
  // All edges incident to node 1 are in MST regarding objective 1
  Graph g(V, 2, false);
  g.addEdge(1, 2, 1, 10);
  g.addEdge(1, 3, 1, 10);
  g.addEdge(1, 4, 1, 10);
  g.addEdge(1, 5, 1, 10);
  g.addEdge(2, 3, 10, 10);

  Graph mstg = g.getMSTKruskal(1);
  EXPECT_EQ(mstg.getE(), V - 1);
  EXPECT_EQ(mstg.getDegree(1), V - 1);
  EXPECT_EQ(mstg.getSumOfEdgeWeights()[0], 4);

  // now add the costly edge (2,3) to initial tree
  // I.e., (2, 3) must be in computed spanning tree
  Graph g2(V, 2, false);
  g2.addEdge(2, 3, 10, 10);

  mstg = g.getMSTKruskal(1, g2);
  EXPECT_EQ(mstg.getE(), V - 1);
  EXPECT_EQ(mstg.getDegree(1), V - 2);
  EXPECT_EQ(mstg.getSumOfEdgeWeights()[0], 13); // c(2, 3) + (|V| - 2) * 1 = 13
}

TEST(MSTTest, KruskalOnGrapherator) {
  std::vector<std::string> files {
    "tests/instances/graph_N100-E259-C5-W2---UNG-CLUNG---CLSTEG-CLDEG-DEG---RWG-RWG_1.graph",
    "tests/instances/graph_N250-E729-C0-W2---UNG---DEG---RWG-RWG_1.graph",
    "tests/instances/graph_N50-E153-C5-W2---LHSNG-CLUNG---CLDEG-CLSTEG---RWG-RWG_1.graph"
  };

  for (auto pathToFile: files) {
    Graph g = Graph::importFromGrapheratorFile(pathToFile);
    Graph mst = g.getMSTKruskal(1);
    std::vector<std::vector<int>> components = mst.getConnectedComponents();
    EXPECT_EQ(mst.getV(), g.getV());
    EXPECT_EQ(mst.getE(), g.getV() - 1);
    EXPECT_EQ(1, components.size());
    EXPECT_EQ(g.getV(), components[0].size());

    // now check if subforest mutation works
    Graph mst2 = g.getMSTBySubforestMutation(mst, 10);
    std::vector<std::vector<int>> components2 = mst2.getConnectedComponents();
    EXPECT_EQ(mst2.getV(), g.getV());
    EXPECT_EQ(mst2.getE(), g.getV() - 1);
    EXPECT_EQ(1, components2.size());
    EXPECT_EQ(g.getV(), components2[0].size());

    // now check if subgraph mutation works
    Graph mst3 = g.getMSTBySubgraphMutation(mst, 10);
    std::vector<std::vector<int>> components3 = mst3.getConnectedComponents();
    EXPECT_EQ(mst3.getV(), g.getV());
    EXPECT_EQ(mst3.getE(), g.getV() - 1);
    EXPECT_EQ(1, components3.size());
    EXPECT_EQ(g.getV(), components3[0].size());
  }
}

} // namespace
