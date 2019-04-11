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

    auto division = divide_grid3(dimx, dimy, dimz, np);
    int periods[] = {0, 1, 1};

    MPI_Cart_create(MPI_COMM_WORLD, 3, division.data(), periods, 1, &cart_comm);
    
    int coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, coords);
    // std::cout << "Rank " << rank << " has coordinates " << coords[0] << coords[1] << coords[2] << "\n";

    for(int dir = 0; dir < 3; ++dir) {
        MPI_Status status;
        int output, input = -1;
        int source, dest;
        output = coords[0] * 100 + coords[1] * 10 + coords[2];

        MPI_Cart_shift(cart_comm, dir, 1, &source, &dest);

        // MPI_Cart_coords(cart_comm, prev, 3, coords);
        // std::cout << "Rank " << rank << " has neighbors " << next << " " << prev << " in direction " << dir << "\n";
        MPI_Sendrecv(&output, 1, MPI_INT, dest, dir, &input, 1, MPI_INT, source, dir, cart_comm, &status);
        // std::cout <<  "Rank " << rank << "Receives " << std::setfill('0') << std::setw(3) << input << " in direction" << dir << "\n";
    }

    if (rank == 0)  {
        for(int i = 0; i < argc; ++i) {
            std::cout << argv[i] << "\n";
        }
    }

    MPI_Finalize();
}