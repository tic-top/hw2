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

// a,b,c,d,e
double h(double a, double b, double c, double d, double e) {
    double z = (f(a) + f(b) + f(c) + f(d) + f(e)) / 5;
    return max ( -25. , min (30. , z));
}

double g(double a, double b, double c, double d, double e) {
    double z = (a + b + c + d + e) / 5;
    return max ( -25. , min (30. , z));
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

void run_serial(int m, int n, int verbose, int P, int ID) {
    // initilize
    vector<vector<double>> A0(m, vector<double> (n, 0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double ai = i;
            double aj = j;
            A0[i][j] = (double)aj * sin(ai) + ai * cos(aj) + sqrt(ai + aj + 1);
        }
    }

    double start = MPI_Wtime();
    for (int it = 0; it < IT_NUM; it++) {
        vector<vector<double>> A(m, vector<double> (n, 0));
        vector<vector<double>> A1(m, vector<double> (n, 0));
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                A1[i][j] = f(A0[i][j]);
            }
        }
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                // A(i, j) = Ao(i, j) if i = 0 or i = m − 1 or j = 0 or j = n − 1 (i.e., it is unchanged along
                if (i == 0 || i == m - 1 || j == 0 || j == n - 1) {
                    A[i][j] = A0[i][j];
                } else {
                    // now is the local level
                    // A[i][j] = h(A0[i][j], A0[i-1][j-1], A0[i+1][j-1], A0[i-1][j+1], A0[i+1][j+1]);
                    A[i][j]= g(A1[i][j], A1[i-1][j-1], A1[i+1][j-1], A1[i-1][j+1], A1[i+1][j+1]);
                }
            }
        }
        A0 = A;
    }
    // sum
    double sum = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            sum += A0[i][j];
        }
    }
    // sum of square
    double sum_square = 0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            sum_square += A0[i][j] * A0[i][j];
        }
    }
    // time
    double end = MPI_Wtime();
    cout << "Sum is:  " << sum << endl;
    cout << "Sum of square is:  " << sum_square << endl;
    cout << "Time: " << end - start << endl;
}

void run_parallel(int m, int n, int verbose, int P, int ID) {
    int n_of_P = sqrt(P);
    double sqrt_P = sqrt(P);
    int sub_rows = ceil((double)m / sqrt_P);
    int sub_cols = ceil((double)n / sqrt_P);
    int row = ID / n_of_P;
    int col = ID % n_of_P;
    // boundary
    int num_row = min(sub_rows * (row + 1), m) - sub_rows * row;
    int num_col = min(sub_cols * (col + 1), n) - sub_cols * col;

    if (verbose)  {
        cout << "Start process " << ID << " of " << P 
            << ", sub_rows is " << sub_rows 
            << ", sub_cols is " << sub_cols
            << ", i_start = " << sub_rows * row 
            << ", i_end = " << min(sub_rows * (row + 1), m) 
            << ", j_start = " << sub_cols * col 
            << ", j_end = " << min(sub_cols * (col + 1), n) 
            << ", num_row = " << num_row
            << ", num_col = " << num_col
            << "\n";  
    }
    vector<vector<double>> A0(num_row, vector<double> (num_col, 0));

    // initialize A0
    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            double ai = sub_rows * row + i;
            double aj = sub_cols * col + j;
            A0[i][j] = (double)aj * sin(ai) + ai * cos(aj) + sqrt(ai + aj + 1);
        }
    }

    // check the initialization
    if (verbose) {
        cout << "row = " << row << ", " << "col = " << col << ", ";
        cout << "Contents of " << ID <<  ": ";
        // for (int i = 0; i < num_row; i++) {
        //     for (int j = 0; j < num_col; j++) {
        //         cout << A0[i][j] << " ";
        //     }
        // }
        cout << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
    double up_row[num_col];
    double down_row[num_col];
    double right_col[num_row];
    double left_col[num_row];
    double t_l, t_r, b_l, b_r;
        
    for (int it = 0; it < IT_NUM; it++) { 
        if (ID == 0 && verbose) {
            cout << "Iteration " << it <<endl;
        }
        // Send diagonal neighbors
        if (row != 0 && col != 0) {  // send to Top-left neighbor
            double msg = A0[0][0];  // Assuming first element for diagonal communication
            MPI_Send(&msg, 1, MPI_DOUBLE, ID - n_of_P - 1, 0, MPI_COMM_WORLD);  // Top-left diagonal
        }
        if (row != 0 && col != n_of_P - 1) {  // send to Top-right neighbor
            double msg = A0[0][num_col - 1];  // Assuming last element for diagonal communication
            MPI_Send(&msg, 1, MPI_DOUBLE, ID - n_of_P + 1, 0, MPI_COMM_WORLD);  // Top-right diagonal
        }
        if (row != n_of_P - 1 && col != 0) {  // send to Bottom-left neighbor
            double msg = A0[num_row - 1][0];  // Assuming first element for diagonal communication
            MPI_Send(&msg, 1, MPI_DOUBLE, ID + n_of_P - 1, 0, MPI_COMM_WORLD);  // Bottom-left diagonal
        }
        if (row != n_of_P - 1 && col != n_of_P - 1) {  // send to Bottom-right neighbor
            double msg = A0[num_row - 1][num_col - 1];  // Assuming last element for diagonal communication
            MPI_Send(&msg, 1, MPI_DOUBLE, ID + n_of_P + 1, 0, MPI_COMM_WORLD);  // Bottom-right diagonal
        }

        // Send l,r,u,b neighbors
        // Send bottom row to the bottom neighbor
        if (row != n_of_P - 1) {  // If not in the last row
            vector<double> msg = A0[sub_rows - 1];  // Get the last row of the submatrix
            if (ID == 0 && verbose) {
                cout << "bottom row msg: " << msg.size() << endl;
            }
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + n_of_P, 0, MPI_COMM_WORLD);  // Send to the bottom
        }
        
        // Send top row to the top neighbor
        if (row != 0) {  // If not in the first row
            vector<double> msg = A0[0];  // Get the first row of the submatrix
            if (ID == 0 && verbose) {
                cout << "top row msg: " << msg.size() << endl;
            }
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - n_of_P, 0, MPI_COMM_WORLD);  // Send to the top
        }

        if (col != 0) {  // If not in the first column
            vector<double> msg;
            get_first_col(msg, A0);  // Get the first column of the submatrix
            if (ID == 0 && verbose) {
                cout << "left col msg: " << msg.size() << endl;
            }
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID - 1, 0, MPI_COMM_WORLD);  // Send to the left
        }

        // Send right column to the right neighbor
        if (col != n_of_P - 1) {  // If not in the last column
            vector<double> msg;
            get_last_col(msg, A0);  // Get the last column of the submatrix
            if (ID == 0 && verbose) {
                // print msg
                cout << "right col msg, ID: "<< ID << " size: " << msg.size() << endl;
                cout << row <<' '<< col << ' ' << n_of_P << endl;
            }
            MPI_Send(&msg[0], msg.size(), MPI_DOUBLE, ID + 1, 0, MPI_COMM_WORLD);  // Send to the right
            if (ID == 0 && verbose) {
                cout << "Finish sending right" << endl;
            }
        }

        // diagonal communication
        // from right
        if (col != n_of_P - 1) {  // If not in the last column
            MPI_Recv(&right_col[0], num_row, MPI_DOUBLE, ID + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from left
        if (col != 0) {  // If not in the first column
            MPI_Recv(&left_col[0], num_row, MPI_DOUBLE, ID - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from top
        if (row != 0) {  // If not in the first row
            MPI_Recv(&up_row[0], num_col, MPI_DOUBLE, ID - n_of_P, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from bottom
        if (row != n_of_P - 1) {  // If not in the last row
            MPI_Recv(&down_row[0], num_col, MPI_DOUBLE, ID + n_of_P, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from top-left
        if (row != 0 && col != 0) {  // Top-left diagonal neighbor
            MPI_Recv(&t_l, 1, MPI_DOUBLE, ID - n_of_P - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from top-right
        if (row != 0 && col != n_of_P - 1) {  // Top-right diagonal neighbor
            MPI_Recv(&t_r, 1, MPI_DOUBLE, ID - n_of_P + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from bottom-left
        if (row != n_of_P - 1 && col != 0) {  // Bottom-left diagonal neighbor
            MPI_Recv(&b_l, 1, MPI_DOUBLE, ID + n_of_P - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // from bottom-right
        if (row != n_of_P - 1 && col != n_of_P - 1) {  // Bottom-right diagonal neighbor
            MPI_Recv(&b_r, 1, MPI_DOUBLE, ID + n_of_P + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        // calculation
        vector<vector<double>> A(num_row, vector<double> (num_col, 0));
        vector<vector<double>> A1(num_row, vector<double> (num_col, 0));

        for (int i = 0; i < num_row; ++i) {
            for (int j = 0; j < num_col; ++j) {
                A1[i][j] = f(A0[i][j]);
            }
        }
        for (int i = 0; i < num_col; ++i){
            up_row[i] = f(up_row[i]);
            down_row[i] = f(down_row[i]);
        }
        for (int i = 0; i < num_row; ++i){
            right_col[i] = f(right_col[i]);
            left_col[i] = f(left_col[i]);
        }
        t_l = f(t_l);
        t_r = f(t_r);
        b_l = f(b_l);
        b_r = f(b_r);

        for (int i = 0; i < num_row; ++i) {
            for (int j = 0; j < num_col; ++j) {
                double global_i = sub_rows * row + i;
                double global_j = sub_cols * col + j;
                if (global_i == 0 || global_i == m - 1 || global_j == 0 || global_j == n - 1) {
                    A[i][j] = A0[i][j];
                } else {
                    // now is the local level
                    // left-up
                    if (i == 0 && j == 0) {
                        A[i][j] = g(A1[i][j], t_l, up_row[j+1], left_col[i+1], A1[i+1][j+1]);
                    }
                    // right-up
                    else if (i == 0 && j == num_col - 1) {
                        A[i][j] = g(A1[i][j], t_r, up_row[j-1], A1[i+1][j-1], right_col[i+1]);
                    }
                    // left-down
                    else if (i == num_row - 1 && j == 0) {
                        A[i][j] = g(A1[i][j], down_row[j+1], b_l, left_col[i-1], A1[i-1][j+1]);
                    }
                    // right-down
                    else if (i == num_row - 1 && j == num_col - 1) {
                        A[i][j] = g(A1[i][j], down_row[j-1], b_r, A1[i-1][j-1], right_col[i-1]);
                    }
                    // top
                    else if (i == 0) {
                        A[i][j] = g(A1[i][j], up_row[j-1], up_row[j+1], A1[i+1][j-1], A1[i+1][j+1]);
                    }
                    // bottom
                    else if (i == num_row - 1) {
                        A[i][j] = g(A1[i][j], down_row[j-1], down_row[j+1], A1[i-1][j-1], A1[i-1][j+1]);
                    }
                    // left
                    else if (j == 0) {
                        A[i][j] = g(A1[i][j], A1[i-1][j+1], A1[i+1][j+1], left_col[i-1], left_col[i+1]);
                    }
                    // right
                    else if (j == num_col - 1) {
                        A[i][j] = g(A1[i][j], A1[i-1][j-1], A1[i+1][j-1], right_col[i-1], right_col[i+1]);
                    }
                    else{
                        A[i][j] = g(A1[i][j], A1[i-1][j-1], A1[i+1][j-1], A1[i-1][j+1], A1[i+1][j+1]);
                    }
                }
            }
        }

        if (verbose && it == IT_NUM-1 && ID == 0) {
            cout << "Contents of A of " << ID <<  ": ";
            for (int i = 0; i < num_row; i++) {
                for (int j = 0; j < num_col; j++) {
                    cout << A[i][j] << " ";
                }
                cout << "\n";
            }
            cout << "\n";
        }

        A0 = A;
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Verification
    // sum
    double local_sum = 0;
    double local_sum_square = 0;

    for (int i = 0; i < num_row; i++) {
        for (int j = 0; j < num_col; j++) {
            local_sum += A0[i][j];
            local_sum_square += A0[i][j] * A0[i][j];
        }
    }

    double global_sum = 0.0 , global_sum2 = 0.0;

    MPI_Reduce(& local_sum , & global_sum , 1, MPI_DOUBLE , MPI_SUM , 0, MPI_COMM_WORLD);
    MPI_Reduce(& local_sum_square , & global_sum2 , 1, MPI_DOUBLE , MPI_SUM , 0, MPI_COMM_WORLD);

    if (ID == 0) {
        double end = MPI_Wtime();
        cout << "Sum is:  " << global_sum << endl;
        cout << "Sum of square is:  " << global_sum2 << endl;
        cout << "Time: " << end - start << endl;
    }

    return;

}


int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int verbose = atoi(argv[3]);
    int parallel = atoi(argv[4]);
    int P;
    int ID;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);
    MPI_Barrier(MPI_COMM_WORLD);
    if (parallel == 1) {
        run_parallel(m, n, verbose, (int)P, (int)ID);
    }
    else {
        run_serial(m, n, verbose, (int)P, (int)ID);
    }

    // Finalize MPI.
    MPI_Finalize();
    return 0;
}