

// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <cmath>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>

using namespace Rcpp;

typedef boost::adjacency_list<
  boost::listS, // Store out-edges in an std::list
  boost::vecS, // Store vertices in a std::vector
  boost::undirectedS, // Undirected graph
  boost::no_property,
  boost::property<boost::edge_weight_t, double>
  > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::property_map<Graph, boost::edge_weight_t>::type Weights;

class RSTlearner{
public:
    const Graph* g; // Pointer to 'parent' graph
    Graph currentGraph;
    std::vector<int> membership; // Maps each vertex to a cluster's integer key
    int k; // Number of clusters
    int n; // Number of vertices in the graph
    std::vector<std::vector<int>> clusts;
    // For each cluster (outer vector), the integer key of the vertices belonging to it (inner vector)
    std::vector<int> new_clust_ids;
    // Multi-purpose vector; contains vertex ids of varying contexts
    // After a birth: Vertices belonging to the new cluster
    // After a death: Vertices belonging to the newly unified cluster
    std::vector<int> old_clust_ids;
    // Multi-purpose vector; contains vertex ids of varying contexts
    // After a birth: Vertices from the cluster which was split excluding the ones from new_clust_ids
    // After a death: Vertices belonging only to the cluster being merged which has a higher number


    RSTlearner(const Graph* graph) {
      g = graph; // Initialize g
      boost::copy_graph(*g, currentGraph); // Initialize currentGraph
      n = boost::num_vertices(currentGraph); // Initialize n
      membership.reserve(n);
      minSpanTree(); // Make currentGraph a minimum spanning tree
      updateMembership(); // Initialize k, membership, clusts
      storeGraphState();
    }

    void updateMembership() {
      membership.clear();
      k = boost::connected_components(currentGraph,
                                      boost::make_iterator_property_map(
                                        membership.begin(), boost::get(
                                            boost::vertex_index, currentGraph
                                        )));
      updateClusts();
    }

    void updateClusts() {
      clusts.clear();
      for(auto i = membership.begin(); i != membership.end(); ++i) {
        clusts[*i].push_back(i - membership.begin());
      }
    }

    void minSpanTree() {
      std::vector<Edge> mst_edges;
      // This shouldn't be necessary since the weights will either be given to
      // current graph by way of copying the full graph (on first initialization)
      // or through the hyper step
      //Weights currentWeight = boost::get(boost::edge_weight, currentGraph);
      boost::kruskal_minimum_spanning_tree(currentGraph, std::back_inserter(mst_edges));
      Graph new_graph(mst_edges.begin(), mst_edges.end()/*, currentWeight*/, n);
      boost::copy_graph(new_graph, currentGraph);
    }

    // [TODO] make this function cheaper
    // Maintain a vector of bools that is updated more efficiently at
    // each step instead of a bulk computation
    std::vector<bool> getBetweenEdges() {
      // Initialize result
      std::vector<bool> result;
      // Get list of all edges (iterator pair of beginning and end)
      auto edges_iter_full = boost::edges(*g);
      for(auto i = edges_iter_full.first; i != edges_iter_full.second; ++i) {
        // For each edge, do the source and target vertex belong to the same cluster?
        // If they belong to the same cluster, it is NOT a between edge (false)
        // If they belong to different clusters, it IS a between edge (true)
        result.push_back(membership[(*i).m_source] != membership[(*i).m_target]);
      }
      return(result);
    }

    void birth() {
      storeGraphState();
        Rcpp::NumericVector probs;
        for (int i = 0; i < k; i++) {
            probs.push_back(clusts[i].size() - 1);
        }
        // First, randomly pick cluster to split (based on size, this could change)
        int clust_split = (Rcpp::sample(k, 1, false, probs))(0);
        // What vertices belong to the cluster we chose to split
        std::vector<int> clust_ids = clusts[clust_split];
        // Pick a vertex from that cluster at which the split will occur
        int vertex_split = (Rcpp::sample(clust_ids.size(), 1, false))(0);
        Vertex vertex_out = boost::vertex(vertex_split, currentGraph);
        // Get the edges connected to that vertex
        auto edges_split = boost::out_edges(vertex_out, currentGraph);
        // Pick an edge to remove
        int edge_split = (Rcpp::sample(std::distance(edges_split.first, edges_split.second), 1, false))(0);
        for(int i = 0; i < edge_split; ++i) {
          (edges_split.first)++;
        }
        Edge edge_remove = *(edges_split.first);
        // Remove the edge
        boost::remove_edge(edge_remove, currentGraph);
        // Update the membership and clusts
        updateMembership();
        // Update new_clust_ids and old_clust_ids
        // [TODO] this loop is inefficient; some trickery involving
        // how Boost calculates connected components could do this faster
        // by directly going to which clusts were modified
        new_clust_ids.clear();
        old_clust_ids.clear();
        int new_clust_num = 0;
        int old_clust_num = 0;
        // Determine the ids of the 2 clusters which made up cluster which was split
        for(auto i = clust_ids.begin(); i != clust_ids.end(); ++i) {
          // New clust ids are those belonging to the greater cluster index
          // Old  clust ids are the remaining clust ids not in the new cluster
          // [TODO] a break condition could be added to this loop; once
          // both cluster nums are determined no need to keep going
          if(membership[*i] > new_clust_num) {
            old_clust_num = new_clust_num;
            new_clust_num = membership[*i];
          } else if(membership[*i] > old_clust_num && membership[*i] != new_clust_num) {
            old_clust_num = membership[*i];
          }
        }
        // Get the vertex ids belonging to the old and new cluster
        new_clust_ids = clusts[new_clust_num];
        old_clust_ids = clusts[old_clust_num];
    }

    void death() {
      storeGraphState();
      // First, find which edges are between 2 clusters
      std::vector<bool> betweenness = getBetweenEdges();
      std::vector<int> between_edge_ids;
      for(int i = 0; i < betweenness.size(); ++i) {
        if(betweenness[i]) {
          between_edge_ids.push_back(i);
        }
      }
      // Randomly sample one of the edges to add back to currentGraph
      int returning_edge_id = (Rcpp::sample(between_edge_ids.size(), 1, false))(0);
      // Get a reference to that edge
      auto edge_iter_full = boost::edges(*g).first;
      for(int i = 0; i < returning_edge_id; ++i) {
        ++edge_iter_full;
      }
      Edge returning_edge = *edge_iter_full;
      // Record the clust ids of the vertices on either side of the new edge
      // (old cluster is the one with the larger id)
      int old_clust_num = membership[returning_edge.m_source];
      int new_clust_num = membership[returning_edge.m_target];
      if(new_clust_num > old_clust_num) {
        int temp = new_clust_num;
        new_clust_num = old_clust_num;
        old_clust_num = temp;
      }
      // Update new_clust_ids and old_clust_ids
      old_clust_ids.clear();
      new_clust_ids.clear();
      old_clust_ids = clusts[old_clust_num];
      new_clust_ids = clusts[new_clust_num];
      new_clust_ids.insert(new_clust_ids.end(), old_clust_ids.begin(), old_clust_ids.end());
      // Add that edge back to the current graph
      boost::add_edge(returning_edge.m_source, returning_edge.m_target, currentGraph);
      // Update the membership and clusts
      updateMembership();
    }

    void hyper() {
      storeGraphState();
      // First, find which edges were connecting different clusters
      std::vector<bool> between_edge_ids = getBetweenEdges();
      // Clear the current graph state, reset it to the full graph
      currentGraph.clear();
      boost::copy_graph(*g, currentGraph);
      // For every edge, set the weight as follows:
      auto edge_iter_full = boost::edges(currentGraph).first;
      for(int i = 0; i < between_edge_ids.size(); ++i) {
        // If the edge connects 2 clusters, sample its weight Unif(0.5,1)
        // Otherwise sample its weight Unif(0,0.5)
        if(between_edge_ids[i]) {
          boost::put(boost::edge_weight, currentGraph, currentGraph[*edge_iter_full], R::runif(0.5, 1.0));
        } else {
          boost::put(boost::edge_weight, currentGraph, currentGraph[*edge_iter_full], R::runif(0.0, 0.5));
        }
        edge_iter_full++;
      }
      // Now that we have new edge weights for the full graph, find a new spanning tree
      minSpanTree();
      // Remove any edges connecting previously distinct clusters
      // Get list of current edges (iterator pair of beginning and end)
      auto edges_iter = boost::edges(currentGraph);
      for(auto i = edges_iter.first; i != edges_iter.second; ++i) {
        // If an edge's ends connects vertices of different membership; remove it
        if(membership[(*i).m_source] != membership[(*i).m_target]) {
          boost::remove_edge(i, currentGraph);
        }
      }
    }
    // This only goes back by one step which is fine for the standard
    // set of birth, death, hyper moves following a rejection
    // However, for 'change' moves (death followed by birth) such as implemented
    // in BAST model, the function will have to store the original state and
    // reverse on its own since it would be going back by two steps
    void reverseMove() {
      currentGraph = old_currentGraph;
      membership = old_membership;
      clusts = old_clusts;
      new_clust_ids = old_new_clust_ids;
      old_clust_ids = old_old_clust_ids;
    }

private:
  Graph old_currentGraph;
  std::vector<int> old_membership;
  std::vector<std::vector<int>> old_clusts;
  std::vector<int> old_new_clust_ids;
  std::vector<int> old_old_clust_ids;

  void storeGraphState() {
    old_currentGraph = currentGraph;
    old_membership = membership;
    old_clusts = clusts;
    old_new_clust_ids = new_clust_ids;
    old_old_clust_ids = old_clust_ids;
  }

};

Rcpp::List BASTIONfit(const Rcpp::IntegerMatrix &edges,
                      const Rcpp::NumericVector &weights,
                      const Rcpp::NumericVector &Y,
                      const int MCMC_iter,
                      const int BURNIN,
                      const int THIN,
                      const int n_learners,
                      const int k_max) {
// First, construct the graph from the edge list
  int n_edges = edges.rows();
  int n_verts = Y.size();
  // Create a vector of vertex pairs
  std::vector<std::pair<int, int>> edge_vec;
  for(int i = 0; i < n_edges; i++) {
    edge_vec.push_back(std::pair<int,int>(edges(i, 0), edges(i, 1)));
  }
  // We now have our full graph with initial weights
  Graph g(edge_vec.begin(), edge_vec.end(), weights, n_verts);
  // Initialize a matrix to store each fitted value throughout the iterations
  // Each column represents a weak learner and each row represents a vertex
  NumericMatrix fitted_mus(n_verts, n_learners);
  // Initialize each weak learner ([TODO] further consideration on whether to use
  // list vs vector to store each weak learner is warranted)
  std::vector<RSTlearner> WeakLearners;
  for (int learner = 0; learner < n_learners; learner++) {
    RSTlearner Learner(&g);
    WeakLearners.push_back(Learner);
  }
  // There should be at most k_max <= n_verts clusters
  int max_clusts = (k_max <= n_verts)? k_max : n_verts;
  // 0 is a birth, 1 is a death, 2 is a change, 3 is a hyper
  NumericVector moves = {0, 1, 2, 3};
  for (int iter = 0; iter < MCMC_iter; iter++) {
      for (int learner = 0; learner < n_learners; learner++) {
        int k_m = WeakLearners[learner].k;
        NumericVector LearnerResponse = Y - rowSums(fitted_mus) - fitted_mus( _ , learner);
        NumericVector moveProbability;
        if(k_m == 1) {
          moveProbability = NumericVector::create(0.9, 0, 0, 0.1);
        } else if(k_m == max_clusts) {
          moveProbability = NumericVector::create(0, 0.6, 0.3, 0.1);
        } else {
          moveProbability = NumericVector::create(0.3, 0.3, 0.3, 0.1);
        }
        int move = sample(move, 1, false, moveProbability)(0);
        switch(move) {
        case 0:
          WeakLearners[learner].birth();
          // Calculate acceptance probability acc_prob
          // if(R::runif(0, 1) > acc_prob) // reject
          //    WeakLearners[learner].reverseMove();
          break;
        case 1:
          WeakLearners[learner].death();
          // Calculate acceptance probability acc_prob
          // if(R::runif(0, 1) > acc_prob) // reject
          //    WeakLearners[learner].reverseMove();
          break;
        case 2:
          WeakLearners[learner].death();
          WeakLearners[learner].birth();
          // Calculate acceptance probability acc_prob
          // if(R::runif(0, 1) > acc_prob) // reject
          //    WeakLearners[learner].reverseMove();
          break;
        case 3:
          WeakLearners[learner].hyper();
        }
        // Update the learner column of g w/ new fitted_mus
      }
    // Update sigma squared of y
    // Save the result
  }
  // Return the result
}







