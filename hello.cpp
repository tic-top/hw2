#include <mpi.h>
#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char **argv)
{

    // Initialize the MPI environment.
    MPI_Init(&argc, &argv);

    // Obtain processor id and world size(number of processors)
    int P;
    int ID;
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &ID);

    cout << "Hello world from processor " << ID << " of " << P << endl;

    // Finalize MPI
    MPI_Finalize();
    return 0;
}