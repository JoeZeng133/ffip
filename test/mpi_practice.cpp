#include <mpi.h>
#include <vector>
#include <iostream>
#include <array>
#include <numeric>
#include <cmath>
#include <iomanip>

std::vector<double> recbuf, sendbuf;

//add -use-hwthread-cpus to use hyperthread

int main(int argc, char** argv) {

    int dimx = 10, dimy = 3, dimz = 30;
    int np, rank;
    MPI_Comm cart_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Finalize();
}