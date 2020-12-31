#include <simulation.hpp>
#include <iostream>

using namespace ffip;

int main(int argc, char** argv) {
    int rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::fstream file{"config.json", std::ios::in};
    std::fstream dbfile{"debug.txt", std::ios::out};
    json config;
    
    file >> config;
    Simulation sim;
    sim.init(config);
    // std::cout << "Process " << rank << " start running\n";
    // sim.output_details(dbfile);
    sim.run(config.at("stop condition"), std::cout);
    sim.output();

    MPI_Finalize();
}
