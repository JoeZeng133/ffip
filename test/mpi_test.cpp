#include <mpi.h>
#include <vector>
#include <iostream>
#include <array>
#include <numeric>
#include <cmath>
#include <iomanip>

std::vector<double> recbuf, sendbuf;

//divide a grid into np parts with smallest maximum chunk size
auto get_chunk_size(int dim1, int dim2, int dim3, int size1, int size2, int size3) -> size_t {
        return  ((dim1 / size1) + dim1 % size1) * 
                ((dim2 / size2) + dim2 % size2) * 
                ((dim3 / size3) + dim3 % size3);
};


std::array<int, 3> divide_grid3(int dim1, int dim2, int dim3, int np) {
    

    std::array<int, 3> res = {1, 1, np};
    size_t min_size = get_chunk_size(dim1, dim2, dim3, 1, 1, np);

    for(int i = 1; i <= std::sqrt(np); ++i) 
    if(np % i == 0) {

        int rest = np / i;

        for(int j = 1; j <= std::sqrt(rest); ++j)
        if (rest % j == 0) {
            int k = rest / j;
            //(i, j, k)
            if (auto size = get_chunk_size(dim1, dim2, dim3, i, j, k); size < min_size) {
                min_size = size;
                res = {i, j, k};
            }

            //(i, k, j)
            if (auto size = get_chunk_size(dim1, dim2, dim3, i, k, j); size < min_size) {
                min_size = size;
                res = {i, k, j};
            }

        }

        for(int j = 1; j <= std::sqrt(i); ++j) 
        if (i % j == 0) {
            int k = i / j;
            //(rest, j, k)
            if (auto size = get_chunk_size(dim1, dim2, dim3, rest, j, k); size < min_size) {
                min_size = size;
                res = {rest, j, k};
            }

            //(rest, k, j)
            if (auto size = get_chunk_size(dim1, dim2, dim3, rest, k, j); size < min_size) {
                min_size = size;
                res = {rest, k, j};
            }
        }
    }

    return res;
}

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
    std::cout << "Rank " << rank << " has coordinates " << coords[0] << coords[1] << coords[2] << "\n";

    for(int dir = 0; dir < 3; ++dir) {
        MPI_Status status;
        int output, input = -1;
        int source, dest;
        output = coords[0] * 100 + coords[1] * 10 + coords[2];

        MPI_Cart_shift(cart_comm, dir, 1, &source, &dest);

        // MPI_Cart_coords(cart_comm, prev, 3, coords);
        // std::cout << "Rank " << rank << " has neighbors " << next << " " << prev << " in direction " << dir << "\n";
        MPI_Sendrecv(&output, 1, MPI_INT, dest, dir, &input, 1, MPI_INT, source, dir, cart_comm, &status);
        std::cout <<  "Rank " << rank << "Receives " << std::setfill('0') << std::setw(3) << input << " in direction" << dir << "\n";
    }

    MPI_Finalize();
}