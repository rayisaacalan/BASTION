// Reinventing the wheel, 1 header file at a time
// TODO: Find a way to do all this using the BOOST C++ libraries
#include <vector>
#include <algorithm>
#include <Rcpp.h>

// For each vertex, we want to store every outbound edge and its weight (int and double)
typedef std::pair<int, double> Pair;

// We also want to store the weight for every edge
struct Edge{
  int parent, child;
  double weight;
};

// The following graph class code is adapted from an assignment submission written for
// Dr. Aakash Tyagi's CSCE 221 class at Texas A&M University
class Graph {
protected:
  // The data structure to represent the graph is an adjacency list
  // https://en.wikipedia.org/wiki/Adjacency_list
  // In short, every vertex is represented by the vector of its edges
  // TODO: assess performance versus adjacency matrix representation
  std::vector<std::vector<Pair>> adjac_list;
public:
  // Construct an empty graph
  Graph() {

  }
  // Construct a graph based on an adjacency list
  Graph(std::vector<std::vector<Pair>> inp_adjac_list) {
    adjac_list = inp_adjac_list;
  }
  // Constructor based on a vector of edges {(parent, child, weight)}
  // and how many total vertices are in the graph (size of exterior vector)
  Graph(std::vector<Edge> &E, int numVertices) {
    // Allocate space for each vertex to be represented
    adjac_list.resize(numVertices);
    // Add each edge to the appropriate subvector
    for(auto elem = E.begin(); elem != E.end(); ++elem) {
      int parent = elem->parent;
      int child = elem->child;
      double weight = elem->weight;
      // Add edge from parent to child
      adjac_list[parent].push_back(std::make_pair(child, weight));
      // Add edge from child to parent [Remove to make graph weighted]
      adjac_list[child].push_back(std::make_pair(parent, weight));
    }
  }

};

// Merging two disconnected subgraphs using an edge
Graph mergeDisconnectedSubGraphs(Graph &g1, Graph &g2, Edge e) {
  // Create merged adjacency list
  std::vector<std::vector<Pair>> new_adjac_list(g1.getGraph());
  // Create space to store result
  new_adjac_list.reserve(g1.nVertex() + g2.nVertex());
  // Append the vertices from g2 onto g1
  std::vector<std::vector<Pair>> g2_adjac_list = g2.getGraph();
  new_adjac_list.insert(new_adjac_list.end(), g2_adjac_list.begin(), g2_adjac_list.end());
  // Create new graph with merged adjacency list
  Graph new_graph(new_adjac_list);
  new_graph.addEdge(e.parent, e.child, e.weight);
  return(new_graph);
}

