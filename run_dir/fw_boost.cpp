#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>
#include <limits>
#include <unordered_map>
#include <time.h>

// type for weight/distance on each edge
typedef double t_weight;
#define INF 99999.0

// define the graph type
typedef boost::property<boost::edge_weight_t, t_weight> EdgeWeightProperty;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, 
                              boost::no_property, EdgeWeightProperty> Graph;

typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

// Declare a matrix type and its corresponding property map that
// will contain the distances between each pair of vertices.
typedef boost::exterior_vertex_property<Graph, t_weight> DistanceProperty;
typedef DistanceProperty::matrix_type DistanceMatrix;
typedef DistanceProperty::matrix_map_type DistanceMatrixMap;

void printSolution(DistanceMatrix distances, int V)
{
    FILE *outfile = fopen("out_boost.txt", "w");
    for (int i = 0; i < V; i++)
    {
        for (int j = 0; j < V; j++)
        {
            if (distances[i][j] == INF)
                fprintf(outfile, "| %6s ", "INF");
            else
                fprintf(outfile, "| %6.2f ",distances[i][j]);
        }
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

int main(int argc, char* argv[])
{

    struct timespec start, end, floyd_start, floyd_end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return -1;
    }

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Unable to open file";
        return -1;
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line[0] != '%' && !line.empty()) {
            break;
        }
    }

    std::istringstream iss(line);
    int num_vertices, num_edges;
    iss >> num_vertices;
    std::getline(infile, line);
    std::istringstream iss2(line);
    iss2 >> num_edges;


    Graph g(num_vertices);

    int u, v;
    t_weight weight;
    std::unordered_map<int, t_weight> self_loops;

    while (infile >> u >> v >> weight) {
        boost::add_edge(u-1, v-1, weight, g);
        if (u == v) {
            self_loops[u-1] = weight;  // Store self-loop weights
        }
    }

    infile.close();

    WeightMap weight_pmap = boost::get(boost::edge_weight, g);

    // Initialize the distance matrix
    DistanceMatrix distances(num_vertices);
    DistanceMatrixMap dm(distances, g);

    // Set the initial distances to infinity
    for (std::size_t i = 0; i < num_vertices; ++i) {
        for (std::size_t j = 0; j < num_vertices; ++j) {
            distances[i][j] = std::numeric_limits<t_weight>::max();
        }
    }

    // Set the initial distances from each vertex to itself to 0
    for (std::size_t i = 0; i < num_vertices; ++i) {
        distances[i][i] = 0;
    }

    // Update distances with the input edges
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei) {
        int u = boost::source(*ei, g);
        int v = boost::target(*ei, g);
        t_weight w = boost::get(weight_pmap, *ei);
        distances[u][v] = w;
    }
    clock_gettime(CLOCK_MONOTONIC, &floyd_start);
    // Find all pairs shortest paths
    bool valid = floyd_warshall_all_pairs_shortest_paths(g, dm, boost::weight_map(weight_pmap));
    clock_gettime(CLOCK_MONOTONIC, &floyd_end);

    // Check if there are no negative cycles
    if (!valid) {
        std::cerr << "Error - Negative cycle in matrix" << std::endl;
        return -1;
    }

   printSolution(distances, num_vertices);

    clock_gettime(CLOCK_MONOTONIC, &end);

        double total_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        double floyd_time = (floyd_end.tv_sec - floyd_start.tv_sec) + (floyd_end.tv_nsec - floyd_start.tv_nsec) / 1e9;

        printf("Total Time measured: %.3f seconds.\n", total_time);
        printf("Floyd Time measured: %.3f seconds.\n", floyd_time);

    return 0;
}
