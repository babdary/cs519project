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

    if (V < size) {
        cout << "Processes cannot be greater than vertices. Vertices: " << V << " Comm size: " << size << endl;
        MPI_Finalize();
        return -1;
    }
    // allocate global data
    vector<vector<double>> dist(V, vector<double>(V, INF));
    for (int i = 0; i < V; i++)
        dist[i][i] = 0; //dist from each node to itself is 0

    int source, target;
    double weight; 
    //store the weights in a matrix
    //this can be read from the input file, but that's inefficent
    for (int i = 0; i < E; i++) {
        ifile >> source >> target >> weight;
        dist[source - 1][target - 1] = min(weight, dist[source - 1][target - 1]); //-1 since input starts from 1 not 0
    }
    ifile.close();
    // 8 + 4 -1 / 4 = 2
    int block_size = (V + size - 1) / size;  // Calculate block size +size -1 is used to handle boundary conditions

    // Allocate local data for each process
    vector<vector<double>> local_dist(block_size, vector<double>(V, INF));
    //
    vector<double> local_data(block_size * (V+1), INF);
    int cc = 0;
    // Distribute rows in a block cyclic manner
    for (int i = 0; i < block_size; i++) {
        int global_row = i * size + rank;
        if (global_row < V) {
            local_dist[i] = dist[global_row];
            for (int j = 0; j < V; j++) {
                local_data[i * V + j+cc] = dist[global_row][j];
            }
            local_data[V+i*V+cc] = global_row;
            cc++;
        }
    }
    
    // Perform the local computation for this process's portion of the matrix
    for (int k = 0; k < V; k++) {
        // Calculate the source rank for the current k
        int root_rank = k % size;
        
        // Broadcast the k-th row to all processes
        vector<double> k_row(V, INF);
        if ((k / size) < block_size && (k % size) == rank) {
            k_row = local_dist[k / size];
        }
        MPI_Bcast(&k_row[0], V, MPI_DOUBLE, root_rank, MPI_COMM_WORLD);

        for (int i = 0; i < block_size; i++) {
            int global_row = i * size + rank;
            if (global_row < V) {
                for (int j = 0; j < V; j++) {
                    if (local_dist[i][j] > local_dist[i][k] + k_row[j]) {
                        local_dist[i][j] = local_dist[i][k] + k_row[j];
                    }
                }
            }
        }
    }

    // Update local_data from local_dist before gathering
    cc = 0;
    for (int i = 0; i < block_size; i++) {
        int global_row = i * size + rank;
        if (global_row < V) {
            for (int j = 0; j < V; j++) {
                local_data[i * V + j+cc] = local_dist[i][j];
            }
            cc++;
        }
    }

    // Prepare the buffer for gathering the results
    vector<double> global_data;
    if (rank == 0) {
        global_data.resize(V * (V+1), INF);
    }

    // Gather the results back to the root process
    MPI_Gather(local_data.data(), block_size * (V+1), MPI_DOUBLE, 
               global_data.data(), block_size * (V+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Convert the gathered 1D data into a 2D matrix
        cc=0;


        vector<vector<double>> global_dist(V, vector<double>(V, INF));
        for (int i = 0; i < V; i++) {
            //cout << "Row index: " << global_data[i*V+V+cc] << "Row: " << global_data[i * V + j] << endl;
            for (int j = 0; j < V; j++) {
                global_dist[global_data[i*V+V+cc]][j] = global_data[i * V + j +cc];
            }
            cc++;
        }

        // Write the result matrix to a file
        writeMatrixToFile(global_dist, V, "output.txt");
    }

    MPI_Finalize();
    return 0;
}
