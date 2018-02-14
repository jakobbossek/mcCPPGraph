#include <fstream>
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
    this->directed = directed;
    // Pay attention: we add a dummy element at position 0
    for (int i = 0; i <= V; ++i) {
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

  bool isDirected() const {
    return this->directed;
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
      if (!this->directed) {
        this->adjList[v].push_back({u, {w1, w2}});
      }
      this->E += 1;
    }
  }

  bool hasEdge(int u, int v) const {
    for (int i = 0; i < this->adjList[u].size(); ++i) {
      if (this->adjList[u][i].first == v)
        return true;
    }
    return false;
  }

  std::vector<std::pair<int, int>> getEdges() const {
    std::vector<std::pair<int, int>> edges(this->getE());
    int edgeCounter = 0;
    for (int u = 1; u <= this->getV(); ++u) {
      for (Edge edge : this->adjList[u]) {
        int v = edge.first;
        // if (u,v) with u < v already added, then (v, u) can be skipped
        if (u < v) {
          edges[edgeCounter] = {u, v};
          edgeCounter += 1;
        }
      }
    }
    return edges;
  }

  std::vector<int> getSumOfEdgeWeights() const {
    std::vector<int> sum(2);
    sum[0] = 0;
    sum[1] = 0;
    for (std::vector<Edge> alist: this->adjList) {
      for (Edge edge: alist) {
        sum[0] += edge.second.first;
        sum[1] += edge.second.second;
      }
    }
    if (!this->isDirected())
      sum[0] /= 2; sum[1] /= 2;
    return sum;
  }

  static Graph getIntersectionGraph(const Graph g1, const Graph g2) {
    // sanity checks
    assert(g1.isDirected() == g2.isDirected());
    assert(g1.getV() == g2.getV());
    assert(g1.getW() == g2.getW());

    int V = g1.getV();

    std::cout << "Intersection of graphs with " << V << " nodes" << std::endl;
    // bare skeleton
    Graph g(V, g1.getW(), g1.isDirected());

    std::cout << "so far" << std::endl;

    // now iterate over all edges
    // Should be possible in O(|E|)
    unsigned int u;
    for (u = 1; u <= V; ++u) {
      std::cout << "u = " << u << std::endl;

      for (int j = 0; j < g1.adjList[u].size(); ++j) {
        std::cout << "j = " << j << std::endl;


        unsigned int v = g1.adjList[u][j].first;
        if (g2.hasEdge(u, v)) {
          std::pair <int, int> weight = g1.adjList[u][j].second;
          //g.addEdge(u, v, weight.first, weight.second);
        }
      }
    }
    std::cout << "Finalized" << std::endl;
    return g;
  }

  static Graph importFromGrapheratorFile(std::string pathToFile) {
    std::ifstream infile(pathToFile, std::ios::in);

    char sep = ',';

    // read first line (nnodes, nedges, nclusters, nweights)
    int V, E, C, W;
    infile >> V >> sep >> E >> sep >> C >> sep >> W;

    // init graph
    Graph g(V, W, false);

    // now go to first character and ignore first 4 + |V| lines
    // (meta data + node coordinates)
    //FIXME: have a look ad ignore function
    infile.seekg(0);
    unsigned int linesToSkip = 4 + V;
    if (C > 0)
      linesToSkip += 1; // skip cluster membership
    unsigned int curLine = 1;
    std::string line;
    while (curLine <= linesToSkip) {
      std::getline(infile, line); // discard
      curLine++;
    }

    // now we are ready to read edge costs section
    //FIXME: generalize to >= 2 objectives
    int u, v;
    double c1, c2;

    while (infile.good()) {
      infile >> u >> sep >> v >> sep >> c1 >> sep >> c2;
      g.addEdge(u, v, c1, c2);
    }

    infile.close();
    return g;
  }


  void removeEdge(const int u, const int v) {
    assert(u >= 1 && u <= this->V);
    assert(v >= 1 && u <= this->V);

    for (int j = 0; j < this->adjList[u].size(); ++j) {
      if (this->adjList[u][j].first == v) {
        // this is ugly!
        this->adjList[u].erase(this->adjList[u].begin() + j);
        this->E -= 1;
        this->degrees[u] -= 1;
        break;
      }
    }

    if (!this->directed) {
      for (int j = 0; j < this->adjList[v].size(); ++j) {
        if (this->adjList[v][j].first == u) {
          // this is ugly!
          this->adjList[v].erase(this->adjList[v].begin() + j);
          this->degrees[v] -= 1;
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
      for (unsigned int u = 1; u <= this->getV(); ++u) {
        // iterate over adjacency list
        for (unsigned int j = 0; j < this->adjList[u].size(); ++j) {
          std::cout << "c(" << u << ", " << this->adjList[u][j].first << ") = (" << this->adjList[u][j].second.first << ", " << this->adjList[u][j].second.second << ")";
          if (j == (this->adjList[u].size() - 1)) {
            std::cout << std::endl;
          } else {
            std::cout << ", ";
          }
        }
      }
      // std::vector<int> sum = this->getSumOfEdgeWeights();
      // std::cout << "Sum of edge weights: c(" << sum[1] << ", " << sum[2] << ")" << std::endl;
    }
  }

  std::vector<int> getConnectedSubtree(int startNode, unsigned int maxNodes) {
    assert(startNode >= 1 && startNode <= this->getV());
    assert(maxNodes >= 1 && maxNodes <= this->getV());

    // we need to store node and its predecessor in BFS tree
    std::vector<int> queue;
    queue.push_back(startNode);

    std::vector<bool> done(this->getV());
    for (int i = 0; i <= this->getV(); ++i) {
      done[i] = false;
    }

    std::vector<int> output;

    unsigned int nsel = 0;

    unsigned int curNode = -1;
    while (nsel < maxNodes && queue.size() > 0) {
      // get node
      curNode = queue.back();
      queue.pop_back();
      output.push_back(curNode);
      done[curNode] = true;

      // access adjList and put all neighbours into queue
      for (Edge edge: this->adjList[curNode]) {
        // skip nodes already added
        if (!done[edge.first]) {
          queue.push_back(edge.first);
        }
      }
      nsel += 1;
    }
    return output;
  }

  std::vector<std::vector<int>> getConnectedComponents() const {
    std::vector<std::vector<int>> components;

    unsigned int V = this->getV();
    std::vector<bool> visited(V + 1);

    for (int node = 1; node <= V; ++node) {
      // already visited, i.e., in some component?
      if (visited[node])
        continue;

      // otherwise perform BFS
      std::vector<int> queue = {node};
      std::vector<int> component;

      while (!queue.empty()) {
        // get topmost element from queue
        int curNode = queue.back();
        queue.pop_back();

        // mark as visited
        visited[curNode] = true;

        // add to current component
        component.push_back(curNode);

        // go through adjacency list and eventually add neighbours to queue
        for (auto edge: this->adjList[curNode]) {
          if (!visited[edge.first]) {
            queue.push_back(edge.first);
          }
        }
      }
      components.push_back(component);
    }

    return components;
  }

  Graph getInducedSubgraph(std::vector<int> nodes) {
    // copy constructor
    Graph g(*this);

    // which nodes should be kept
    std::vector<bool> keep(this->getV() + 1);
    for (auto node: nodes) {
      keep[node] = true;
    }

    // now go through all edges and drop, if not both endpoints should be kept
    for (unsigned int u = 1; u <= this->getV(); ++u) {
      for (unsigned int j = 0; j < this->adjList[u].size(); ++j) {
        unsigned int v = this->adjList[u][j].first;
        if (!keep[u] || !keep[v]) {
          g.removeEdge(u, v);
        }
      }
    }

    return g;
  }

  Graph getMSTKruskal(double weight) {
    assert(weight >= 0 && weight <= 1);

    // represent minimum spanning tree with graph object
    Graph initialTree(this->getV(), this->getW(), false);
    UnionFind UF(this->V);

    return this->getMSTKruskal(weight, initialTree, UF);
  }

  Graph getMSTKruskal(double weight, Graph &initialTree) {
    assert(weight >= 0 && weight <= 1);

    std::vector<std::vector<int>> components = initialTree.getConnectedComponents();
    UnionFind UF(this->V, components);

    return this->getMSTKruskal(weight, initialTree, UF);
  }

  Graph getMSTKruskal(double weight, Graph &initialTree, UnionFind &UF) {
    // assert(weight >= 0 && weight <= 1);

    // // represent minimum spanning tree with graph object
    // Graph mst(this->getV(), this->getW(), false);

    // // init efficient set data structure
    // UnionFind UF(this->V);

    // now we need to transform graph to list of edges which can be
    // sorted by weights
    std::vector<std::pair<int, std::pair<int, int>>> edgelist;
    for (int u = 1; u <= this->V; ++u) {
      for (int j = 0; j < this->adjList[u].size(); ++j) {
        // FIXME: ugly as sin! and restricted to two weights!
        double costs = weight * this->adjList[u][j].second.first + (1 - weight) * this->adjList[u][j].second.second;
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
        //DEBUG("Linking components");
        // link components
        //FIXME: here I need to return the original untransformed weights
        // I.e., "edgelist" should be triple std::vector<std::triple<int, std::pair<int, int>, std::pair<int, int>>>
        initialTree.addEdge(u, v, it->first, it->first);
        // merge components
        UF.unite(u, v);
      }

      // found spanning tree if number of edge is |V| - 1
      if (initialTree.getE() == (this->getV() - 1)) {
        //DEBUG("Tree has |V| - 1 edges, i.e., it is a spanning tree.");
        break;
      }
    }
    return initialTree;
  }
private:
  int V;
  int E;
  int W;
  bool directed;
  std::vector<unsigned int> degrees;
  AdjacencyList adjList;
};
#endif // GRAPH
