
// [[Rcpp::depends(RcppClock)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

#include <RcppClock.h>
#include <thread>
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
      //Rprintf(" Initializing a weak learner");
      g = graph; // Initialize g
      boost::copy_graph(*g, currentGraph); // Initialize currentGraph
      n = boost::num_vertices(currentGraph); // Initialize n
      //Rprintf(".");
      minSpanTree(); // Make currentGraph a minimum spanning tree
      //Rprintf(".");
      updateMembership(); // Initialize k, membership, clusts
      //Rprintf(".");
      storeGraphState();
          //printGraph();
      //Rprintf(" done\n");
    }

    void updateMembership() {
      membership.clear();
      std::vector<int> mem(boost::num_vertices(currentGraph));
      k = boost::connected_components(currentGraph,
                                      boost::make_iterator_property_map(
                                        mem.begin(), boost::get(
                                            boost::vertex_index, currentGraph
                                        )));
      /*Rprintf("\n   New Membership:\n");
      for(int i = 0; i < mem.size(); ++i) {
        //Rprintf("   %i ", mem[i]);
      }
      //Rprintf("\n");*/
      membership.assign(mem.begin(), mem.end());
      updateClusts();
    }

    void updateClusts() {
      clusts.clear();
      clusts.reserve(k);
      for(int i = 0; i < k; ++i) {
        clusts.emplace_back(std::vector<int>());
        clusts[i].reserve(n);
      }
      for(int i = 0; i < membership.size(); ++i) {
        clusts[(membership[i])].emplace_back(i);
      }
    }

    void minSpanTree() {
      std::vector<Edge> mst_edges;
      std::vector<std::pair<int,int>> mst_edges_vertices;
      // This shouldn't be necessary since the weights will either be given to
      // current graph by way of copying the full graph (on first initialization)
      // or through the hyper step
      //Weights currentWeight = boost::get(boost::edge_weight, currentGraph);
      boost::kruskal_minimum_spanning_tree(currentGraph, std::back_inserter(mst_edges));
      mst_edges_vertices.reserve(mst_edges.size());
      //Rprintf("\n");
      for(auto i = 0; i < mst_edges.size(); ++i) {
        std::pair<int,int> vertex_pair = std::make_pair(mst_edges[i].m_source, mst_edges[i].m_target);
        mst_edges_vertices.emplace_back(vertex_pair);
        //Rprintf("%i --- %i\n", vertex_pair.first, vertex_pair.second);
      }
      Graph new_graph(mst_edges_vertices.begin(), mst_edges_vertices.end()/*, currentWeight*/, n);
      currentGraph.clear();
      boost::copy_graph(new_graph, currentGraph);
    }

    // [TODO] make this function cheaper
    // Maintain a vector of bools that is updated more efficiently at
    // each step instead of a bulk computation
    std::vector<bool> getBetweenEdges() {
      // Initialize result
      //Rprintf("         Getting between edges\n");
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
      storeGraphState();
      //printMembership();
      //Rprintf("       Entering a Birth step\n");
        Rcpp::NumericVector probs(k);
        for (int i = 0; i < k; ++i) {
            probs(i) = (clusts[i].size() - 1);
            //Rprintf("        Cluster %i has prob weight %i\n", i, clusts[i].size() - 1);
        }
        // First, randomly pick cluster to split (based on size, this could change)
        int clust_split = (Rcpp::sample(k, 1, false, probs))(0) - 1;
        //Rprintf("       Splitting cluster %i\n", clust_split);
        // What vertices belong to the cluster we chose to split
        std::vector<int> clust_ids = clusts[clust_split];
        // Pick a vertex from that cluster at which the split will occur
        int vertex_split = clust_ids[(Rcpp::sample(clust_ids.size(), 1, false))(0) - 1];
        //Rprintf("       Splitting at vertex %i\n", vertex_split);
        Vertex vertex_out = boost::vertex(vertex_split, currentGraph);
        // Get the edges connected to that vertex
        /*if(boost::degree(vertex_out, currentGraph) == 0) {
          printGraph();
          stop("Selected vertex has no edges");
        }*/
        auto edges_split = boost::out_edges(vertex_out, currentGraph);
        // Pick an edge to remove
        int num_eligible_edges = std::distance(edges_split.first, edges_split.second);
        for(auto i = edges_split.first; i != edges_split.second; ++i) {
          //Rprintf("       %i -- %i\n", (*i).m_source, (*i).m_target);
        }
        int edge_split = /*(num_eligible_edges == 0) ? 0 :*/ ((Rcpp::sample(num_eligible_edges, 1, false))(0) - 1);
        //Rprintf("       Removing the %i out edge\n", edge_split);
        for(int i = 0; i < edge_split; ++i) {
          ++(edges_split.first);
        }
        Edge edge_remove = *(edges_split.first);
        // Remove the edge
        //Rprintf("       Removing edge\n");
        boost::remove_edge(edge_remove, currentGraph);
        // Update the membership and clusts
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
          if(membership[*i] > new_clust_num) {
            old_clust_num = new_clust_num;
            new_clust_num = membership[*i];
          } else if(membership[*i] > old_clust_num && membership[*i] != new_clust_num) {
            old_clust_num = membership[*i];
          }
        }
        //Rprintf("       Fetching clust ids\n");
        // Get the vertex ids belonging to the old and new cluster
        new_clust_ids = clusts[new_clust_num];
        old_clust_ids = clusts[old_clust_num];
    }

    void death() {
      //Rprintf("       Entering a Death step\n");
      //printMembership();
      //Rprintf("       Current edges:\n");
      //printEdges();
      storeGraphState();
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
      // Randomly sample one of the edges to add back to currentGraph
      int returning_edge_id = between_edge_ids[(Rcpp::sample(between_edge_ids.size(), 1, false))(0) - 1];
      //Rprintf("\n       Returning edge has id %i\n", returning_edge_id);
      // Get a reference to that edge
      for(int i = 0; i < returning_edge_id; ++i) {
        ++edge_iter_full;
      }
      Edge returning_edge = *edge_iter_full;
      //Rprintf("       Found returning edge\n");
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
      boost::add_edge(returning_edge.m_source, returning_edge.m_target, currentGraph);
      //Rprintf("       Readded edge\n");
      // Update the membership and clusts
      updateMembership();
    }

    void hyper() {
      storeGraphState();
      //printMembership();
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
          boost::put(boost::edge_weight, currentGraph, *edge_iter_full, R::runif(0.5, 1.0));
        } else {
          boost::put(boost::edge_weight, currentGraph, *edge_iter_full, R::runif(0.0, 0.5));
        }
        ++edge_iter_full;
      }
      // Now that we have new edge weights for the full graph, find a new spanning tree
      minSpanTree();
      // Remove any edges connecting previously distinct clusters
      // Get list of current edges (iterator pair of beginning and end)
      auto edges_iter = boost::edges(currentGraph);
      for(auto i = edges_iter.first; i != edges_iter.second; ++i) {
        // If an edge's ends connects vertices of different membership; remove it
        if(membership[(*i).m_source] != membership[(*i).m_target]) {
          boost::remove_edge(*i, currentGraph);
        }
      }
      // Shouldn't be necessary; membership should not change after hyper step
      //updateMembership();
      //printGraph();
    }

    void printGraph() {
      Rprintf("\n");
      Rprintf("Graph of size %i with %i clusters:\n", n, k);
      auto edge_iters = boost::edges(currentGraph);
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
      auto edge_iters = boost::edges(currentGraph);
      Rprintf("\n");
      for(auto iter = edge_iters.first; iter != edge_iters.second; ++iter) {
        Rprintf("%i -- %i\n", (*iter).m_source, (*iter).m_target);
      }
    }

    // This only goes back by one step which is fine for the standard
    // set of birth, death, hyper moves following a rejection
    // However, for 'change' moves (death followed by birth) such as implemented
    // in BAST model, the function will have to store the original state and
    // reverse on its own since it would be going back by two steps
    void reverseMove() {
      //Rprintf("     Reversing move\n");
      currentGraph = old_currentGraph;
      membership = old_membership;
      clusts = old_clusts;
      new_clust_ids = old_new_clust_ids;
      old_clust_ids = old_old_clust_ids;
      n = old_n;
      k = old_k;
      //printMembership();
    }

private:
  Graph old_currentGraph;
  std::vector<int> old_membership;
  std::vector<std::vector<int>> old_clusts;
  std::vector<int> old_new_clust_ids;
  std::vector<int> old_old_clust_ids;
  int old_n;
  int old_k;

  void storeGraphState() {
    old_currentGraph = currentGraph;
    old_membership = membership;
    old_clusts = clusts;
    old_new_clust_ids = new_clust_ids;
    old_old_clust_ids = old_clust_ids;
    old_n = n;
    old_k = k;
  }

};


// Utility function (similar to unlist)
// Taken from:
// https://stackoverflow.com/a/30201460/14190296

NumericVector combine(const List& list) {
  std::size_t n = list.size();
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  // Allocate the vector
  NumericVector output = no_init(total_length);
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);

    // Update the index
    index += el.size();
  }
  return output;
}

// [[Rcpp::export]]
Rcpp::List BASTIONfit_cpp(const Rcpp::IntegerMatrix &edges,
                      const Rcpp::NumericVector &weights,
                      const Rcpp::NumericVector &Y,
                      const int MCMC_iter,
                      const int BURNIN,
                      const int THIN,
                      const List init_values,
                      const List hyperpars) {
//Rprintf("Entered function body\n");
//Rcpp::Clock clock;
// Fetch hyper parameters
  //clock.tick("init_vals");
  // Fetch the initial values and hyper parameters
    const int k_max = hyperpars["k_max"];
    const int n_learners = hyperpars["n_learners"];
    const double lambda_k = hyperpars["lambda_k"];
    const double lambda_s = hyperpars["lambda_s"];
    const double nu = hyperpars["nu"];
    const double sigmasq_mu = hyperpars["sigmasq_mu"];
// Fetch initial values
    double sigmasq_y = init_values["sigmasq_y"];
    List mu = init_values["mu"];
//Rprintf("Fetched initial values\n");
  //clock.tock("init_vals");
  //clock.tick("graph_build");
// First, construct the graph from the edge list
  int n_edges = edges.rows();
  int n_verts = Y.size();
  // Create a vector of vertex pairs
  std::vector<std::pair<int, int>> edge_vec;
  for(int i = 0; i < n_edges; ++i) {
    edge_vec.emplace_back(std::make_pair(edges(i, 0), edges(i, 1)));
  }
//Rprintf("Created edge vector\n");
  // We now have our full graph with initial weights
  Graph g(edge_vec.begin(), edge_vec.end(), as<std::vector<double>>(weights).begin(), n_verts);
//Rprintf("Constructed graph\n");
  //clock.tock("graph_build");
  //clock.tick("learners_build");
  // Initialize a matrix to store each fitted value throughout the iterations
  // Each column represents a weak learner and each row represents a vertex
  NumericMatrix fitted_mus(n_verts, n_learners);
  // Vector of k's for evaluating log posterior
  NumericVector K(n_learners);
  // Initialize each weak learner ([TODO] further consideration on whether to use
  // list vs vector to store each weak learner is warranted)
  std::vector<RSTlearner> WeakLearners;
  for (int learner = 0; learner < n_learners; ++learner) {
    RSTlearner Learner(&g);
    WeakLearners.emplace_back(Learner);
    K(learner) = (Learner.k);
  }
//Rprintf("Initialized Weak Learners\n");
  // Initialize some values to be used in the MCMC
  // Storage for output values
  int out_length = (((MCMC_iter) - BURNIN) % THIN);
  List mu_out;
  NumericVector sigmasq_y_out(out_length);
  List cluster_out;
  NumericVector log_post_out(out_length);
  // There should be at most k_max <= n_verts clusters
  int max_clusts = (k_max <= n_verts)? k_max : n_verts;
  // 0 is a birth, 1 is a death, 2 is a change, 3 is a hyper
  NumericVector moves = {0, 1, 2, 3};
  //clock.tock("learners_build");
//Rprintf("Initialized MCMC Values\n");
  for (int iter = 0; iter < MCMC_iter; ++iter) {
    //Rprintf("Entering iteration %i\n", iter);
    NumericMatrix MembershipByLearner(n_verts, n_learners);
      for (int learner = 0; learner < n_learners; ++learner) {
        //Rprintf(" Entering Weak Learner %i,", learner);
        int k_m = WeakLearners[learner].k;
        //Rprintf(" Currently %i clusters\n", k_m);
        NumericVector LearnerResponse = Y - (rowSums(fitted_mus) - fitted_mus( _ , learner));
        NumericVector moveProbability;
        if(k_m == 1) {
          moveProbability = NumericVector::create(0.9, 0, 0, 0.1);
        } else if(k_m == max_clusts) {
          moveProbability = NumericVector::create(0, 0.6, 0.3, 0.1);
        } else {
          moveProbability = NumericVector::create(0.3, 0.3, 0.3, 0.1);
        }
        int move = sample(moves, 1, false, moveProbability)(0);
        switch(move) {
          case 0: {
            //Rprintf("   Performing a Birth step\n");
            //clock.tick("birth");
            WeakLearners[learner].birth();
          // Calculate acceptance probability acc_prob
            // Log-prior ratio
            double log_A = log(lambda_k) - log(k_m + 1);
              // If k_m is almost at k_max, rd_new is smaller
            double rd_new = (k_m == (std::min(k_max, n_verts) - 1)) ? 0.6 : 0.3;
            // Log-proposal ratio
            double log_P = log(rd_new) - log(moveProbability(0));
            // Log-likelihood ratio
            double sigma_ratio = sigmasq_y / sigmasq_mu;
            int csize_old = WeakLearners[learner].old_clust_ids.size();
            int csize_new = WeakLearners[learner].new_clust_ids.size();
            NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
            NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
            NumericVector Response_old = LearnerResponse[vid_old];
            NumericVector Response_new = LearnerResponse[vid_new];
            double sum_e_old = sum(Response_old);
            double sum_e_new = sum(Response_new);
            double logdetdiff = -0.5 * (
              log(csize_old + sigma_ratio) +
              log(csize_new + sigma_ratio) -
              log(csize_old + csize_new + sigma_ratio) -
              log(sigma_ratio)
            );
            double quaddiff = (0.5/sigmasq_y) * (
              ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) +
              ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
              ((pow(sum_e_new + sum_e_old, 2))/(csize_old+csize_new+sigma_ratio))
            );
            double log_L = logdetdiff + quaddiff;
            // Calculate acceptance probability
            double log_acc_prob = fmin(0.0, log_A + log_P + log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) > acc_prob) { // reject
              WeakLearners[learner].reverseMove();
            }
            //clock.tock("birth");
            break;
          }
          case 1: {
            //Rprintf("   Performing a Death step\n");
            //clock.tick("death");
            WeakLearners[learner].death();
          // Calculate acceptance probability acc_prob
            // Log-prior ratio
            double log_A = log(k_m) - log(lambda_k);
            // If k_m is almost at k_max, rd_new is smaller
            double rb_new = (k_m == 2) ? 0.9 : 0.3;
            // Log-proposal ratio
            double log_P = log(moveProbability(1)) - log(rb_new);
            // Log-likelihood ratio
            double sigma_ratio = sigmasq_y / sigmasq_mu;
            int csize_old = WeakLearners[learner].old_clust_ids.size();
            int csize_new = WeakLearners[learner].new_clust_ids.size();
            NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
            NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
            NumericVector Response_old = LearnerResponse[vid_old];
            NumericVector Response_new = LearnerResponse[vid_new];
            double sum_e_old = sum(Response_old);
            double sum_e_new = sum(Response_new);
            double logdetdiff = -0.5 * (
                log(csize_new + sigma_ratio) -
                log(csize_old + sigma_ratio) -
                log(csize_new - csize_old + sigma_ratio) +
                log(sigma_ratio)
            );
            double quaddiff = (0.5/sigmasq_y) * (
                ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
                ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) -
                ((pow(sum_e_new - sum_e_old, 2))/(csize_new-csize_old+sigma_ratio))
            );
            double log_L = logdetdiff + quaddiff;
            // Calculate acceptance probability
            double log_acc_prob = fmin(0.0, log_A + log_P + log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) > acc_prob) {
              WeakLearners[learner].reverseMove();
            } // reject
            //clock.tock("death");
            break;
          }
          case 2: {
            //Rprintf("   Performing a Change step\n");
            //clock.tick("change");
            // Store the current graph state (since it will be going forwards
            // by two steps)
            Graph old_currentGraph = WeakLearners[learner].currentGraph;
            std::vector<int> old_membership = WeakLearners[learner].membership;
            std::vector<std::vector<int>> old_clusts = WeakLearners[learner].clusts;
            std::vector<int> old_new_clust_ids = WeakLearners[learner].new_clust_ids;
            std::vector<int> old_old_clust_ids = WeakLearners[learner].old_clust_ids;
            // First, compute log-likelihood ratio for the death move
            WeakLearners[learner].death();

              // Log-likelihood ratio
              double sigma_ratio = sigmasq_y / sigmasq_mu;
              int csize_old = WeakLearners[learner].old_clust_ids.size();
              int csize_new = WeakLearners[learner].new_clust_ids.size();
              NumericVector vid_old = wrap(WeakLearners[learner].old_clust_ids);
              NumericVector vid_new = wrap(WeakLearners[learner].new_clust_ids);
              NumericVector Response_old = LearnerResponse[vid_old];
              NumericVector Response_new = LearnerResponse[vid_new];
              double sum_e_old = sum(Response_old);
              double sum_e_new = sum(Response_new);
              double logdetdiff = -0.5 * (
                log(csize_new + sigma_ratio) -
                  log(csize_old + sigma_ratio) -
                  log(csize_new - csize_old + sigma_ratio) +
                  log(sigma_ratio)
              );
              double quaddiff = (0.5/sigmasq_y) * (
                ((sum_e_new*sum_e_new)/(csize_new+sigma_ratio)) -
                  ((sum_e_old*sum_e_old)/(csize_old+sigma_ratio)) -
                  ((pow(sum_e_new - sum_e_old, 2))/(csize_new-csize_old+sigma_ratio))
              );
              double log_L_death = logdetdiff + quaddiff;

            // Now, compute the log-likelihood ratio for the birth move
            WeakLearners[learner].birth();

              // Log-likelihood ratio
              int csize_old_b = WeakLearners[learner].old_clust_ids.size();
              int csize_new_b = WeakLearners[learner].new_clust_ids.size();
              NumericVector vid_old_b = wrap(WeakLearners[learner].old_clust_ids);
              NumericVector vid_new_b = wrap(WeakLearners[learner].new_clust_ids);
              NumericVector Response_old_b = LearnerResponse[vid_old_b];
              NumericVector Response_new_b = LearnerResponse[vid_new_b];
              double sum_e_old_b = sum(Response_old_b);
              double sum_e_new_b = sum(Response_new_b);
              double logdetdiff_b = -0.5 * (
                log(csize_old_b + sigma_ratio) +
                  log(csize_new_b + sigma_ratio) -
                  log(csize_old_b + csize_new_b + sigma_ratio) -
                  log(sigma_ratio)
              );
              double quaddiff_b = (0.5/sigmasq_y) * (
                ((sum_e_old_b*sum_e_old_b)/(csize_old_b+sigma_ratio)) +
                  ((sum_e_new_b*sum_e_new_b)/(csize_new_b+sigma_ratio)) -
                  ((pow(sum_e_new_b + sum_e_old_b, 2))/(csize_old_b+csize_new_b+sigma_ratio))
              );
              double log_L_birth = logdetdiff_b + quaddiff_b;

            // Add the log likelihood ratios for the two moves
            double log_L = log_L_death + log_L_birth;
            double log_acc_prob = fmin(0.0, log_L);
            double acc_prob = exp(log_acc_prob);
            if(R::runif(0, 1) > acc_prob) { // reject
              // Manually reverse step
              WeakLearners[learner].currentGraph = old_currentGraph;
              WeakLearners[learner].membership = old_membership;
              WeakLearners[learner].clusts = old_clusts;
              WeakLearners[learner].new_clust_ids = old_new_clust_ids;
              WeakLearners[learner].old_clust_ids = old_old_clust_ids;
            }
            //clock.tock("change");
            break;
          }
          case 3: {
            //clock.tick("hyper");
            //Rprintf("   Performing a Hyper step\n");
            WeakLearners[learner].hyper();
            //clock.tock("hyper");
            break;
          }
        }
        //clock.tick("learner_update");
        //Rprintf("   Updating values\n");
        // Update the learner column of fitted_mus w/ new fitted mus
        k_m = WeakLearners[learner].k;
        NumericVector membership = wrap(WeakLearners[learner].membership);
        if((Rcpp::max(membership) + 1) > k_m) {
          WeakLearners[learner].printGraph();
          stop("Something has gone wrong with connected clusters!");
        }
        NumericVector csize(WeakLearners[learner].clusts.size());
        NumericVector ClustResponse(WeakLearners[learner].clusts.size());
        for(int i = 0; i < WeakLearners[learner].clusts.size(); ++i) {
          csize(i) = (WeakLearners[learner].clusts[i].size());
          NumericVector ClustIndices = wrap(WeakLearners[learner].clusts[i]);
          NumericVector ResponseByClust = LearnerResponse[ClustIndices];
          ClustResponse(i) = (sum(ResponseByClust));
        }
        NumericVector Qinv_diag = 1/((csize/sigmasq_y) + 1/(sigmasq_mu));
        NumericVector b = (Qinv_diag * ClustResponse) / sigmasq_y;
        NumericVector NewMus(k_m);
        for(int i = 0; i < k_m; ++i) {
          NewMus(i) = (R::rnorm(b[i], sqrt(Qinv_diag)[i]));
        }
        mu[learner] = NewMus;
        std::vector<double> new_response;
        for(int i = 0; i < n_verts; ++i) {
          new_response.emplace_back(NewMus[WeakLearners[learner].membership[i]]);
        }
        NumericVector NewResponse = wrap(new_response);
        fitted_mus( _ , learner) = NewResponse;
        MembershipByLearner( _ , learner) = membership;
        K[learner] = WeakLearners[learner].k;
        //clock.tock("learner_update");
      }
  // Update sigma squared of y
    // Simplification from R logic:
    // Y_hat = fitted_mus[, n_learners - 1] + Y - (Y - (rowsums(fitted_mus) - fitted_mus[, n_learners - 1]))
    // = fitted_mus[, n_learners - 1] + Y - Y + rowsums(fitted_mus) - fitted_mus[, n_learners - 1]
    // = rowsums(fitted_mus)
    //Rprintf(" Updating Iteration Values\n");
    //clock.tick("iter_update");
    NumericVector Y_hat = rowSums(fitted_mus);
    double scale = 1/(0.5*(nu*lambda_s + sum((Y - Y_hat)*(Y - Y_hat))));
    double shape = (double)(n_verts+nu)/2.0;
    sigmasq_y = 1/(R::rgamma(shape, scale));
    // Save the result
    if(((iter+1) > BURNIN) && ((((iter+1) - BURNIN) % THIN) == 0)) {
      mu_out.push_back(clone(mu));
      sigmasq_y_out.push_back(sigmasq_y);
      cluster_out.push_back(MembershipByLearner);
    // Evaluate log posterior
      // First evaluate log-likelihood
      double log_prior = -((1+(nu/2.0))*log(sigmasq_y))
                         -((nu*lambda_s)/(2.0*sigmasq_y))
                         +sum(-1*(lchoose(n_verts-1,K-1))
                              + K*log(lambda_k)
                              - lfactorial(K))
                         -(sum(pow(combine(mu), 2))/(2.0*sigmasq_mu));
      // Then evaluate log-likelihood
      double log_like = (-n_verts/2.0)*log(sigmasq_y)
                        -(sum((Y - Y_hat)*(Y - Y_hat))/(2.0*sigmasq_y));
      // Add the two
      log_post_out.push_back(log_prior+log_like);
    }
    //clock.tock("iter_update");
  }
//Rprintf("Finished MCMC Loop\n");
  //clock.stop("BASTIONfit_C_timings");
  // Return the result
  List OutputList = List::create(
    Named("mu_out") = mu_out,
    Named("sigmasq_y_out") = sigmasq_y_out,
    Named("cluster_out") = cluster_out,
    Named("log_post_out") = log_post_out
  );
  //Rprintf("Finished computations");
  return(OutputList);
}






