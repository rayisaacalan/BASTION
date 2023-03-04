#ifndef RSTLEARNER_HPP
#define RSTLEARNER_HPP

// No longer using RcppClock
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <cmath>
#include <random>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/copy.hpp>

// This definition enables debug statements to Rcout
// #define RSTLEARNER_DEBUG

#ifdef RSTLEARNER_DEBUG
#define DEBUG_MSG(str) do { Rcpp::Rcout << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

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

class RSTlearner;

class RSTlearner {
public:
  std::mt19937 generator; // PRNG

  const Graph* g; // Pointer to 'parent' graph
  Graph current_graph;
  std::vector<int> membership; // Maps each vertex to a cluster's integer key
  int k; // Number of clusters
  int n; // Number of vertices in the graph
  std::vector<std::vector<int>> clusts;

  /* Storage for proposals */
  Graph prop_graph;
  std::vector<int> prop_membership;
  int prop_k;
  std::vector<std::vector<int>> prop_clusts;
  // For each cluster (outer vector), the integer key of the vertices belonging to it (inner vector)
  std::vector<int> new_clust_ids;
  // Multi-purpose vector; contains vertex ids of varying contexts
  // After a birth: Vertices belonging to the new cluster
  // After a death: Vertices belonging to the newly unified cluster
  std::vector<int> old_clust_ids;
  // Multi-purpose vector; contains vertex ids of varying contexts
  // After a birth: Vertices from the cluster which was split excluding the ones from new_clust_ids
  // After a death: Vertices belonging only to the cluster being merged which has a higher number

  RSTlearner(const Graph* graph,
             const unsigned int seed) :
    n(boost::num_vertices(*graph)),
    generator(seed)
  {
    DEBUG_MSG("RSTlearner: Initializing");
    g = graph; // Initialize g
    boost::copy_graph(*g, prop_graph); // Initialize prop_graph
    // n = boost::num_vertices(*graph); // Initialize n
    DEBUG_MSG("RSTlearner: copied graph");
    minSpanTree(); // Make prop_graph a minimum spanning tree
    DEBUG_MSG("RSTlearner: constructed MST");
    updateMembership(); // Initialize prop_k, prop_membership, prop_clusts
    DEBUG_MSG("RSTlearner: initialized prop_k, prop_membership, prop_clusts");
    storeGraphState();
    accept(); // Move all proposed values to 'current' values
    //printGraph();
    DEBUG_MSG("RSTlearner: Finished Initializing");
  }

  void updateMembership() {
    prop_membership.clear();
    std::vector<int> mem(n);
    prop_k = boost::connected_components(prop_graph,
                                    boost::make_iterator_property_map(
                                      mem.begin(), boost::get(
                                          boost::vertex_index, prop_graph
                                      )));
    /*Rprintf("\n   New Membership:\n");
     for(int i = 0; i < mem.size(); ++i) {
     //Rprintf("   %i ", mem[i]);
     }
     //Rprintf("\n");*/
     prop_membership.assign(mem.begin(), mem.end());
     updateClusts();
  }

  void updateClusts() {
    prop_clusts.clear();
    prop_clusts.reserve(prop_k);
    for(int i = 0; i < prop_k; ++i) {
      prop_clusts.emplace_back(std::vector<int>());
      prop_clusts[i].reserve(n);
    }
    for(int i = 0; i < prop_membership.size(); ++i) {
      prop_clusts[(prop_membership[i])].emplace_back(i);
    }
  }

  void minSpanTree() {
    std::vector<Edge> mst_edges;
    std::vector<std::pair<int,int>> mst_edges_vertices;
    // This shouldn't be necessary since the weights will either be given to
    // current graph by way of copying the full graph (on first initialization)
    // or through the hyper step
    //Weights currentWeight = boost::get(boost::edge_weight, prop_graph);
    boost::kruskal_minimum_spanning_tree(prop_graph, std::back_inserter(mst_edges));
    mst_edges_vertices.reserve(mst_edges.size());
    //Rprintf("\n");
    for(auto i = 0; i < mst_edges.size(); ++i) {
      std::pair<int,int> vertex_pair = std::make_pair(mst_edges[i].m_source, mst_edges[i].m_target);
      mst_edges_vertices.emplace_back(vertex_pair);
      //Rprintf("%i --- %i\n", vertex_pair.first, vertex_pair.second);
    }
    Graph new_graph(mst_edges_vertices.begin(), mst_edges_vertices.end()/*, currentWeight*/, n);
    prop_graph.clear();
    boost::copy_graph(new_graph, prop_graph);
  }

  // [TODO] make this function cheaper
  // Maintain a vector of bools that is updated more efficiently at
  // each step instead of a bulk computation
  std::vector<bool> getBetweenEdges() {
    // Initialize result
    DEBUG_MSG("RSTlearner: Getting between edges");
    std::vector<bool> result;
    // Get list of all edges (iterator pair of beginning and end)
    auto edges_iter_full = boost::edges(*g);
    //Rprintf("         Declared variables\n");
    for(auto i = edges_iter_full.first; i != edges_iter_full.second; ++i) {
      // For each edge, do the source and target vertex belong to the same cluster?
      // If they belong to the same cluster, it is NOT a between edge (false)
      // If they belong to different clusters, it IS a between edge (true)
      //Rprintf("         Checking an edge\n");
      result.emplace_back((membership[(*i).m_source]) != (membership[(*i).m_target]));
      //Rprintf("         Checked an edge\n");
    }
    return(result);
  }

  void birth() {
    if(!can_move) {
      Rcpp::stop("Learner attempted to move before accepting or rejecting proposal.\n");
    }
    //printMembership();
    DEBUG_MSG("RSTlearner: Entering a Birth step\n");
    std::vector<double> probs(k);
    for (int i = 0; i < k; ++i) {
      probs[i] = (clusts[i].size() - 1);
      //Rprintf("        Cluster %i has prob weight %i\n", i, clusts[i].size() - 1);
    }
    // First, randomly pick cluster to split (based on size, this could change)
    int clust_split = samp_idx_discrete(probs);
    //Rprintf("       Splitting cluster %i\n", clust_split);
    // What vertices belong to the cluster we chose to split
    std::vector<int> clust_ids = clusts[clust_split];
    // Pick a vertex from that cluster at which the split will occur
    int vertex_split = clust_ids[samp_idx_discrete(clust_ids.size())];
    //Rprintf("       Splitting at vertex %i\n", vertex_split);
    Vertex vertex_out = boost::vertex(vertex_split, current_graph);
    // Get the edges connected to that vertex
    /*if(boost::degree(vertex_out, current_graph) == 0) {
     printGraph();
     stop("Selected vertex has no edges");
    }*/
    auto edges_split = boost::out_edges(vertex_out, current_graph);
    // Pick an edge to remove
    int num_eligible_edges = std::distance(edges_split.first, edges_split.second);
    for(auto i = edges_split.first; i != edges_split.second; ++i) {
      //Rprintf("       %i -- %i\n", (*i).m_source, (*i).m_target);
    }
    int edge_split = /*(num_eligible_edges == 0) ? 0 :*/ samp_idx_discrete(num_eligible_edges);
    //Rprintf("       Removing the %i out edge\n", edge_split);
    for(int i = 0; i < edge_split; ++i) {
      ++(edges_split.first);
    }
    Edge edge_remove = *(edges_split.first);
    // Remove the edge
    DEBUG_MSG("RSTlearner: Removing edge");
    boost::remove_edge(edge_remove.m_source, edge_remove.m_target, prop_graph);
    // Update the prop_membership, prop_k, and prop_clusts
    updateMembership();
    // Update new_clust_ids and old_clust_ids
      // [TODO] this loop is inefficient; some trickery involving
      // how Boost calculates connected components could do this faster
      // by directly going to which clusts were modified
      //Rprintf("       Calculating old & new clust ids\n");
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
      if(prop_membership[*i] > new_clust_num) {
        old_clust_num = new_clust_num;
        new_clust_num = prop_membership[*i];
      } else if(prop_membership[*i] > old_clust_num && prop_membership[*i] != new_clust_num) {
        old_clust_num = prop_membership[*i];
      }
    }
    DEBUG_MSG("RSTlearner: Fetching clust ids");
    // Get the vertex ids belonging to the old and new cluster
    new_clust_ids = prop_clusts[new_clust_num];
    old_clust_ids = prop_clusts[old_clust_num];
    can_move = false;
  }

  void death() {
    if(!can_move) {
      Rcpp::stop("Learner attempted to move before accepting or rejecting proposal.\n");
    }
    DEBUG_MSG("RSTlearner: Entering a Death step");
    //printMembership();
    //Rprintf("       Current edges:\n");
    //printEdges();
    auto edge_iter_full = boost::edges(*g).first;
    auto edge_iter_full_copy = boost::edges(*g).first;
    // First, find which edges are between 2 clusters
    std::vector<bool> betweenness = getBetweenEdges();
    //Rprintf("       Got betweenness\n");
    std::vector<int> between_edge_ids;
    between_edge_ids.reserve(betweenness.size());
    for(int i = 0; i < betweenness.size(); ++i) {
      if(betweenness[i]) {
        between_edge_ids.emplace_back(i);
      }
      //Rprintf("         Edge %i -- %i between: ",
      //        (*edge_iter_full_copy).m_source,
      //        (*edge_iter_full_copy).m_target);
      //Rcout << betweenness[i] << ", " << std::endl;
      ++edge_iter_full_copy;
    }
    //Rprintf("\n       Found between edge ids: ");
    /*for(int i = 0; i < between_edge_ids.size(); ++i) {
     //Rprintf("%i, ", between_edge_ids[i]);
    }*/
    // Randomly sample one of the edges to add back to current_graph
    int returning_edge_id = between_edge_ids[samp_idx_discrete(between_edge_ids.size())];
    //Rprintf("\n       Returning edge has id %i\n", returning_edge_id);
    // Get a reference to that edge
    for(int i = 0; i < returning_edge_id; ++i) {
      ++edge_iter_full;
    }
    Edge returning_edge = *edge_iter_full;
    DEBUG_MSG("RSTlearner: Found returning edge");
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
    //Rprintf("       Recorded old and new clust_ids\n");
    // Add that edge back to the current graph
    //Rprintf("       Adding edge: %i -- %i\n", returning_edge.m_source, returning_edge.m_target);
    boost::add_edge(returning_edge.m_source, returning_edge.m_target, prop_graph);
    DEBUG_MSG("RSTlearner: Readded edge");
    // Update the membership and clusts
    updateMembership();
    can_move = false;
  }

  void hyper() {
    if(!can_move) {
      Rcpp::stop("Learner attempted to move before accepting or rejecting proposal.\n");
    }
    //storeGraphState();
    //printMembership();
    // First, find which edges were connecting different clusters
    std::vector<bool> between_edge_ids = getBetweenEdges();
    // Clear the current graph state, reset it to the full graph
    prop_graph.clear();
    boost::copy_graph(*g, prop_graph);
    // For every edge, set the weight as follows:
    auto edge_iter_full = boost::edges(prop_graph).first;
    for(int i = 0; i < between_edge_ids.size(); ++i) {
      // If the edge connects 2 clusters, sample its weight Unif(0.5,1)
      // Otherwise sample its weight Unif(0,0.5)
      if(between_edge_ids[i]) {
        boost::put(boost::edge_weight, prop_graph, *edge_iter_full, r_std_unif(0.5, 1.0));
      } else {
        boost::put(boost::edge_weight, prop_graph, *edge_iter_full, r_std_unif(0.0, 0.5));
      }
      ++edge_iter_full;
    }
    // Now that we have new edge weights for the full graph, find a new spanning tree
    minSpanTree();
    // Remove any edges connecting previously distinct clusters
    // Get list of current edges (iterator pair of beginning and end)
    boost::graph_traits<Graph>::edge_iterator ei, ei_end, next;
    tie(ei, ei_end) =boost::edges(prop_graph);
    for(next = ei; ei != ei_end; ei = next) {
      ++next;
      if(membership[(*ei).m_source] != membership[(*ei).m_target])
        boost::remove_edge(*ei, prop_graph);

    }
    // Shouldn't be necessary; membership should not change after hyper step
    //updateMembership();
    //accept();
    //printGraph();
    can_move = false;
  }

  // Utility function to randomly sample from 0 to probs.size() - 1 (integer)
  // (Rcpp::sample(k, 1, false, probs))(0) - 1 ==> samp_idx_discrete(probs)
  int samp_idx_discrete(std::vector<double> probs) {
    return (std::discrete_distribution<int>(probs.begin(), probs.end()))(generator);
  }

  // Utility function to randomly uniform sample from 0 to max_idx - 1 (integer)
  // (Rcpp::sample(vector.size(), 1, false))(0) - 1 ==> samp_idx_discrete(vector.size())
  int samp_idx_discrete(int max_idx) {
    std::vector<double> probs(max_idx, (double) (1.0/max_idx));
    return (std::discrete_distribution<int>(probs.begin(), probs.end()))(generator);
  }

  // Sample single standard uniform random variable
  double r_std_unif() {
    return ((std::uniform_real_distribution<double>(0.0, 1.0))(generator));
  }

  // Sample single non-standard uniform random variable
  double r_std_unif(double a, double b) {
    return ((std::uniform_real_distribution<double>(a, b))(generator));
  }

  void printGraph() {
    Rprintf("\n");
    Rprintf("Graph of size %i with %i clusters:\n", n, k);
    auto edge_iters = boost::edges(current_graph);
    Rprintf("\n");
    for(auto iter = edge_iters.first; iter != edge_iters.second; ++iter) {
      Rprintf("%i -- %i\n", (*iter).m_source, (*iter).m_target);
    }
    Rprintf("\nClusts:\n");
    for(int i = 0; i < clusts.size(); ++i) {
      Rprintf("\n   Clust %i:\n   ", i);
      for(int j = 0; j < clusts[i].size(); ++j) {
        Rprintf("%i, ", (clusts[i])[j]);
      }
    }
    Rprintf("\n Membership:\n");
    for(int i : membership) {
      Rprintf("%i ", i);
    }
    Rprintf("\n");
  }

  void printMembership() {
    Rprintf("\n       Current Membership:\n       ");
    for(int i : membership) {
      Rprintf("%i ", i);
    }
    Rprintf("\n");
  }

  void printEdges() {
    auto edge_iters = boost::edges(current_graph);
    Rprintf("\n");
    for(auto iter = edge_iters.first; iter != edge_iters.second; ++iter) {
      Rprintf("%i -- %i\n", (*iter).m_source, (*iter).m_target);
    }
  }

  void reject() {
    DEBUG_MSG("RSTlearner: Rejecting move");
    prop_graph = current_graph;
    prop_membership = membership;
    prop_clusts = clusts;
    prop_k = k;
    //printMembership();
    can_move = true;
  }

  // Accept the proposed move
  void accept() {
    DEBUG_MSG("RSTlearner: Accepting move");
    current_graph = prop_graph;
    membership = prop_membership;
    clusts = prop_clusts;
    k = prop_k;
    //printMembership();
    can_move = true;
  }

  void storeGraphState() {
    stored_current_graph = current_graph;
    stored_membership = membership;
    stored_clusts = clusts;
    stored_new_clust_ids = new_clust_ids;
    stored_old_clust_ids = old_clust_ids;
    stored_k = k;
  }

  void restoreGraphState() {
    current_graph = stored_current_graph;
    membership = stored_membership;
    clusts = stored_clusts;
    new_clust_ids = stored_new_clust_ids;
    old_clust_ids = stored_old_clust_ids;
    k = stored_k;
  }

protected:
  // Storage for state of the graph (needed for change move)
  Graph stored_current_graph;
  std::vector<int> stored_membership;
  std::vector<std::vector<int>> stored_clusts;
  std::vector<int> stored_new_clust_ids;
  std::vector<int> stored_old_clust_ids;
  int stored_k;
  // Flag for whether the learner is ready to move
  bool can_move;


};


#endif
