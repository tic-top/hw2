#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <cassert>
using namespace std;

int IT_NUM = 10;

double f(double x) {
    double y = x * 2.0;
    for (int i = 1; i <= 10; ++i) {
        y += x * cos(y + i) / pow(1.5, i);
    }
    return y;
}

double g(int z) {
    // Calculate min{30, z}
    double min_val = (z < 30) ? z : 30;
    
    // Calculate max{-25, min{30, z}}
    retrun (min_val > -25) ? min_val : -25;
}


void get_first_col(vector<double>& msg, const vector<vector<double>>& A0) {
    for (vector<double> vec: A0) msg.push_back(vec[0]);
}

void get_last_col(vector<double>& msg, const vector<vector<double>>& A0) {
    for (vector<double> vec: A0) msg.push_back(vec[vec.size()-1]);
}

void append_col(vector<vector<double>>& A0) {
    for (vector<double>& vec: A0) vec.push_back(0);
}

void run_parallel(int m, int n, double (*f)(double), int verbose, int P, int ID) {
    int n_of_P = sqrt(P);
    int sub_rows = ceil(m / sqrt(P));
    int sub_cols = ceil(n / sqrt(P));

    int row = floor(ID / sqrt(P));
    int col = ID - row * sqrt(P);

    if (verbose)  {
        cout << "Start process " << ID << " of " << P 
            << ", sub_rows is " << sub_rows 
            << ", sub_cols is " << sub_cols
            << ", i_start = " << sub_rows * row 
            << ", i_end = " << min(sub_rows * (row + 1), m) 
            << ", j_start = " << sub_cols * col 
            << ", j_end = " << min(sub_cols * (col + 1), n) << " ";
    }

    int num_rows = min(sub_rows * (row + 1), m) - sub_rows * row;
    int num_cols = min(sub_cols * (col + 1), n) - sub_cols * col;

    vector<vector<double>> A0(num_row, vector<double> (num_col, 0));

    // initialize A0
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            A0[i][j] = (double)(col*sub_n + j) * sin(row*sub_n + i) + (row*sub_n + i) * cos(col*sub_n + j) + sqrt(row*sub_n + i + col*sub_n + j + 1);
        }
    }

    // check the initialization
    if (verbose) {
        cout << "row = " << row << ", " << "col = " << col << ", ";
        cout << "Contents of " << ID <<  ": ";
        for (int i = 0; i < num_row; i++) {
            for (int j = 0; j < num_col; j++) {
                cout << A0[i][j] << " ";
            }
        }
        cout << "\n";
    }

    // start sending 
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    for (int it = 0; it < IT_NUM; it++) {
        // For each process, send the required columns or rows to neighboring processes
        // Send left column to the left neighbor
        if (col != 0) {  // If not in the first column
            vector<double> msg;
            get_first_col(msg, A0);  // Get the first column of the submatrix
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - 1, 0, MPI_COMM_WORLD);  // Send to the left
        }

        // Send right column to the right neighbor
        if (col != sqrt(P) - 1) {  // If not in the last column
            vector<double> msg;
            get_last_col(msg, A0);  // Get the last column of the submatrix
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + 1, 0, MPI_COMM_WORLD);  // Send to the right
        }

        // Send top row to the top neighbor
        if (row != 0) {  // If not in the first row
            vector<double> msg = A0[0];  // Get the first row of the submatrix
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - sqrt(P), 0, MPI_COMM_WORLD);  // Send to the top
        }

        // Send bottom row to the bottom neighbor
        if (row != sqrt(P) - 1) {  // If not in the last row
            vector<double> msg = A0[sub_rows - 1];  // Get the last row of the submatrix
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + sqrt(P), 0, MPI_COMM_WORLD);  // Send to the bottom
        }

        // Handle diagonal neighbors
        if (row != 0 && col != 0) {  // Top-left diagonal neighbor
            vector<double> msg;
            get_first_col(msg, A0);  // Assuming first column for diagonal communication
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - sqrt(P) - 1, 0, MPI_COMM_WORLD);  // Top-left diagonal
        }

        if (row != 0 && col != sqrt(P) - 1) {  // Top-right diagonal neighbor
            vector<double> msg;
            get_last_col(msg, A0);  // Assuming last column for diagonal communication
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - sqrt(P) + 1, 0, MPI_COMM_WORLD);  // Top-right diagonal
        }

        if (row != sqrt(P) - 1 && col != 0) {  // Bottom-left diagonal neighbor
            vector<double> msg;
            get_first_col(msg, A0);  // Assuming first column for diagonal communication
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + sqrt(P) - 1, 0, MPI_COMM_WORLD);  // Bottom-left diagonal
        }

        if (row != sqrt(P) - 1 && col != sqrt(P) - 1) {  // Bottom-right diagonal neighbor
            vector<double> msg;
            get_last_col(msg, A0);  // Assuming last column for diagonal communication
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + sqrt(P) + 1, 0, MPI_COMM_WORLD);  // Bottom-right diagonal
        }


        // Receiving the row above from the top neighbor
        if (row != 0) {  // If not in the first row
            MPI_Recv(&last_up_row[0], num_col, MPI_DOUBLE, ID - sqrt(P), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Receiving the column to the left from the left neighbor
        if (col != 0) {  // If not in the first column
            MPI_Recv(&last_up_col[0], num_row, MPI_DOUBLE, ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Receiving the row below from the bottom neighbor
        if (row != sqrt(P) - 1) {  // If not in the last row
            MPI_Recv(&last_down_row[0], num_col, MPI_DOUBLE, ID + sqrt(P), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Receiving the column to the right from the right neighbor
        if (col != sqrt(P) - 1) {  // If not in the last column
            MPI_Recv(&last_down_col[0], num_row, MPI_DOUBLE, ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Receiving individual points (corners or adjacent points)
        // l_l (left neighbor point at the same row)
        if (col != 0) {
            MPI_Recv(&l_l, 1, MPI_DOUBLE, ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // l_r (right neighbor point at the same row)
        if (col != sqrt(P) - 1) {
            MPI_Recv(&l_r, 1, MPI_DOUBLE, ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // l_u (top neighbor point at the same column)
        if (row != 0) {
            MPI_Recv(&l_u, 1, MPI_DOUBLE, ID - sqrt(P), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // l_d (bottom neighbor point at the same column)
        if (row != sqrt(P) - 1) {
            MPI_Recv(&l_d, 1, MPI_DOUBLE, ID + sqrt(P), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }


        // calculation
        vector<vector<double>> A(num_row, vector<double> (num_col, 0));
        for (int i = 0; i < num_row; ++i) {
            for (int j = 0; j < num_col; ++j) {
                // Compute the solution based on the received values and local matrix elements
                double f_value = 0.0;

                // Handle boundaries and apply the equation
                if (i == 0 && j == 0) {  // Top-left corner
                    f_value = (f(last_up_row[j]) + f(last_up_row[j+1]) + f(last_down_row[j]) + f(last_down_row[j+1]) + f(A0[i][j])) / 5;
                } 
                else if (i == 0 && j == num_col - 1) {  // Top-right corner
                    f_value = (f(last_up_row[j-1]) + f(last_up_row[j]) + f(last_down_row[j-1]) + f(last_down_row[j]) + f(A0[i][j])) / 5;
                }
                else if (i == num_row - 1 && j == 0) {  // Bottom-left corner
                    f_value = (f(last_up_row[j]) + f(last_up_row[j+1]) + f(last_down_row[j]) + f(last_down_row[j+1]) + f(A0[i][j])) / 5;
                }
                else if (i == num_row - 1 && j == num_col - 1) {  // Bottom-right corner
                    f_value = (f(last_up_row[j-1]) + f(last_up_row[j]) + f(last_down_row[j-1]) + f(last_down_row[j]) + f(A0[i][j])) / 5;
                }
                else if (i == 0) {  // Top edge (not corners)
                    f_value = (f(last_up_row[j-1]) + f(last_up_row[j+1]) + f(last_down_row[j-1]) + f(last_down_row[j+1]) + f(A0[i][j])) / 5;
                }
                else if (i == num_row - 1) {  // Bottom edge (not corners)
                    f_value = (f(last_up_row[j-1]) + f(last_up_row[j+1]) + f(last_down_row[j-1]) + f(last_down_row[j+1]) + f(A0[i][j])) / 5;
                }
                else if (j == 0) {  // Left edge (not corners)
                    f_value = (f(last_up_row[j]) + f(last_up_row[j+1]) + f(last_down_row[j]) + f(last_down_row[j+1]) + f(A0[i][j])) / 5;
                }
                else if (j == num_col - 1) {  // Right edge (not corners)
                    f_value = (f(last_up_row[j-1]) + f(last_up_row[j]) + f(last_down_row[j-1]) + f(last_down_row[j]) + f(A0[i][j])) / 5;
                }
                else {  // Interior points
                    f_value = (f(A0[i-1][j-1]) + f(A0[i-1][j+1]) + f(A0[i+1][j-1]) + f(A0[i+1][j+1]) + f(A0[i][j])) / 5;
                }

                // Store the calculated value in A
                A[i][j] = g(f_value);
            }
        }

        if (verbose && it == IT_NUM-1) {
            cout << "Contents of A of " << ID <<  ": ";
            for (int i = 0; i < num_row; i++) {
                for (int j = 0; j < num_col; j++) {
                    cout << A[i][j] << " ";
                }
                // cout << "\n";
            }
            cout << "\n";
        }

        A0 = A;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Verification
    // sum
    double local_sum = 0;
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            local_sum += A0[i][j];
        }
    }

    if (ID != 0){
        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        double total_sum = local_sum;
        for (int i = 1; i < P; i++) {
            double i_sum = -1;
            MPI_Recv(&i_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += i_sum;
        }
        cout << "Sum is:  " << total_sum << endl;
    }

    // sum of square
    double local_sum = 0;
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            local_sum += A0[i][j] * A0[i][j];
        }
    }
    if (ID != 0){
        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    else {
        double total_sum = local_sum;
        for (int i = 1; i < P; i++) {
            double i_sum = -1;
            MPI_Recv(&i_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum += i_sum;
        }
        cout << "Sum of square is:  " << total_sum << endl;
    }

    if (ID == 0) {
        double end = MPI_Wtime();
        cout << "Time: " << end - start << endl;
    }

    return;

}


int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    // initilize
    vector<vector<double>> A0(m, vector<double> (n, 0));
    
    int P;
    int ID;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    
    // print the input matrix
    int verbose = atoi(argv[3]);
    if (ID == 0) {
        cout << "Root processor " << ID << " is initializing." <<  "\n";
        for (int i = 0; i < m; i++) 
            for (int j = 0; j < n; j++) 
                A0[i][j] = (double)j * sin(i) + i * cos(j) + sqrt(i + j + 1); // double?
        if (verbose) {
            cout << "n is: " << n <<  "\n";
            cout << "A0:\n";
            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) 
                cout << A0[i][j] << " ";
                cout << "\n";
            }
            cout << "\n";
        } 
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    run_parallel(n, &f, verbose, (int)P, (int)ID);

    // Finalize MPI.
    MPI_Finalize();
    return 0;
}