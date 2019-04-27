#include <simulation.hpp>
#include <iostream>

using namespace ffip;

int main(int argc, char** argv) {
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if (rank == 0)
    // {
    //     int i = 0;
    //     char hostname[256];
    //     gethostname(hostname, sizeof(hostname));
    //     printf("PID %d on %s ready for attach\n", getpid(), hostname);
    //     fflush(stdout);
    //     while(i==0) sleep(5);
    //     std::cout << "Process Attached\n";
    // }

    std::fstream file{"config.json", std::ios::in};
    std::fstream dbfile{"debug.txt", std::ios::out};
    json config;
    
    file >> config;
    Simulation sim;
    sim.init(config);
    // sim.output_details(dbfile);
    sim.run(config.at("stop condition"), dbfile);
    sim.output();

    MPI_Finalize();
}