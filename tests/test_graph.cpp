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

} // namespace