#include <iostream>
#include <cstdlib>
#include <vector>
#include <assert.h>
#include "debug.h"
#include "UnionFind.h"

#ifndef GRAPH
#define GRAPH

typedef std::pair< int, std::pair <int, int> > Edge;
typedef std::vector< std::vector <Edge> > AdjacencyList;

class Graph {
public:
  Graph(int V, int W, bool directed = false) {
    this->V = V;
    this->E = 0;
    this->W = W;
    this->isDirected = directed;
    for (int i = 1; i <= V; ++i) {
      degrees.push_back(0);
      adjList.push_back(std::vector<Edge>());
    }
  }

  int getV() const {
    return this->V;
  }

  int getE() const {
    return this->E;
  }

  int getW() const {
    return this->W;
  }

  unsigned int getDegree(int v) const {
    assert(v >= 1 && v <= this->V);
    return this->degrees[v];
  }

  void addEdge(int u, int v, int w1, int w2) {
    //FIXME: generalize to multiple edge weights
    //FIXME: use templates to allow integer or double weights
    assert(u >= 1 && u <= this->V);
    assert(v >= 1 && u <= this->V);

    if (!this->hasEdge(u, v)) {
      this->adjList[u].push_back({v, {w1, w2}});
      this->degrees[u] += 1;
      if (!this->isDirected) {
        this->adjList[v].push_back({u, {w1, w2}});
      }
      this->E += 1;
    }
  }

  bool hasEdge(int u, int v) {
    for (int i = 0; i < this->adjList[u].size(); ++i) {
      if (this->adjList[u][i].first == v)
        return true;
    }
    return false;
  }

  std::vector<int> getSumOfEdgeWeights() const {
    std::vector<int> sum(2);
    for (std::vector<Edge> alist: this->adjList) {
      for (Edge edge: alist) {
        sum[1] += edge.second.first;
        sum[2] += edge.second.second;
      }
    }
    if (!isDirected)
      sum[1] /= 2; sum[2] /= 2;
    return sum;
  }

  void removeEdge(const int u, const int v) {
    assert(u >= 1 && u <= this->V);
    assert(v >= 1 && u <= this->V);

    for (int j = 0; j < this->adjList[u].size(); ++j) {
      if (this->adjList[u][j].first == v) {
        // this is ugly!
        this->adjList[u].erase(this->adjList[u].begin() + j);
        this->E -= 1;
        break;
      }
    }

    if (!this->isDirected) {
      for (int j = 0; j < this->adjList[v].size(); ++j) {
        if (this->adjList[v][j].first == u) {
          // this is ugly!
          this->adjList[v].erase(this->adjList[v].begin() + j);
          break;
        }
      }
    }
  }

  void print(bool detailed = false) const {
    std::cout << "Weighted graph: n = "
              << V << ", m = "
              << E << ", p = "
              << W << std::endl;
    if (detailed) {
      // iterate over nodes
      for (int u = 1; u <= this->V; ++u) {
        // iterate over adjacency list
        for (int j = 0; j < this->adjList[u].size(); ++j) {
          std::cout << "c(" << u << ", " << this->adjList[u][j].first << ") = (" << this->adjList[u][j].second.first << ", " << this->adjList[u][j].second.second << ")";
          if (j == (this->adjList[u].size() - 1)) {
            std::cout << std::endl;
          } else {
            std::cout << ", ";
          }
        }
      }
      std::vector<int> sum = this->getSumOfEdgeWeights();
      std::cout << "Sum of edge weights: c(" << sum[1] << ", " << sum[2] << ")" << std::endl;
    }
  }

  //FIXME: weight -> lambda
  //FIXME: weight needs to be double -> edge weights are double
  Graph getMSTKruskal(int weight) {
    assert(weight >= 0 && weight <= 1);

    // represent minimum spanning tree with graph object
    Graph mst(this->getV(), this->getW(), false);

    // init efficient set data structure
    UnionFind UF(this->V);

    // now we need to transform graph to list of edges which can be
    // sorted by weights
    std::vector<std::pair<int, std::pair<int, int>>> edgelist;
    for (int u = 1; u <= this->V; ++u) {
      for (int j = 0; j < this->adjList[u].size(); ++j) {
        // FIXME: ugly as sin! and restricted to two weights!
        int costs = weight * this->adjList[u][j].second.first + (1 - weight) * this->adjList[u][j].second.second;
        // FIXME: how to remove duplicates!?
        edgelist.push_back({costs, {u, this->adjList[u][j].first}});
      }
    }

    // sort edges in increasing order
    sort(edgelist.begin(), edgelist.end());

    // Edge iterator
    std::vector<std::pair<int, std::pair<int, int>>>::iterator it;
    for (it = edgelist.begin(); it != edgelist.end(); it++) {
      // get end nodes
      int u = it->second.first;
      int v = it->second.second;

      if (!UF.find(u, v)) {
        DEBUG("Linking components");
        // link components
        //FIXME: here I need to return the original untransformed weights
        // I.e., "edgelist" should be triple std::vector<std::triple<int, std::pair<int, int>, std::pair<int, int>>>
        mst.addEdge(u, v, it->first, it->first);
        // merge components
        UF.unite(u, v);
      }

      // found spanning tree if number of edge is |V| - 1
      if (mst.getE() == (this->getV() - 1)) {
        DEBUG("Tree has |V| - 1 edges, i.e., it is a spanning tree.");
        break;
      }
    }
    return mst;
  }
private:
  int V;
  int E;
  int W;
  bool isDirected;
  std::vector<unsigned int> degrees;
  AdjacencyList adjList;
};
#endif // GRAPH
