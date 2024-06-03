#include <iostream>
#include <fstream>
#include <string>
#include <mpi.h>
#include <vector>

using namespace std;

// Define infinity as a large number (used for initialization)
const double INF = 1e10;

void writeMatrixToFile(const vector<vector<double>> &matrix, int V, const string &filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error opening output file!" << endl;
        return;
    }

    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (matrix[i][j] == INF)
                outfile << "INF ";
            else
                outfile << matrix[i][j] << " ";
        }
        outfile << endl;
    }

    outfile.close();
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) {
            cout << "Check README for usage." << endl;
        }
        MPI_Finalize();
        return -1;
    }

    int V, E;
    ifstream ifile(argv[1]);

    if (!ifile) {
        if (rank == 0) {
            cout << "File not found." << endl;
        }
        MPI_Finalize();
        return -1;
    }

    ifile >> V >> E;

    vector<vector<double>> dist(V, vector<double>(V, INF));
    for (int i = 0; i < V; i++)
        dist[i][i] = 0;

    int u, v;
    double w;
    for (int i = 0; i < E; i++) {
        ifile >> u >> v >> w;
        dist[u - 1][v - 1] = w;
    }
    ifile.close();

    int rows_per_process = V / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank + 1) * rows_per_process;

    // Handle the case where V is not perfectly divisible by size
    if (rank == size - 1) {
        end_row = V;
    }

    // Perform the local computation for this process's portion of the matrix
    for (int k = 0; k < V; k++) {
        // Calculate the root rank for the current k
        int root_rank = k / rows_per_process;
        if (root_rank >= size) {
            root_rank = size - 1;
        }

        // Broadcast the k-th row to all processes
        vector<double> k_row(V);
        if (rank == root_rank) {
            k_row = dist[k];
        }
        MPI_Bcast(&k_row[0], V, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

        for (int i = start_row; i < end_row; i++) {
            for (int j = 0; j < V; j++) {
                if (dist[i][j] > dist[i][k] + k_row[j]) {
                    dist[i][j] = dist[i][k] + k_row[j];
                }
            }
        }
    }

    // Prepare the buffer for gathering the results
    vector<double> local_data;
    for (int i = start_row; i < end_row; i++) {
        local_data.insert(local_data.end(), dist[i].begin(), dist[i].end());
    }

    vector<double> global_data;
    if (rank == 0) {
        global_data.resize(V * V);
    }

    // Gather the results back to the root process
    MPI_Gather(local_data.data(), local_data.size(), MPI_DOUBLE, 
               global_data.data(), local_data.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Convert the gathered 1D data into a 2D matrix
        vector<vector<double>> global_dist(V, vector<double>(V));
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < V; j++) {
                global_dist[i][j] = global_data[i * V + j];
            }
        }

        // Write the result matrix to a file
        writeMatrixToFile(global_dist, V, "output.txt");
    }

    MPI_Finalize();
    return 0;
}
