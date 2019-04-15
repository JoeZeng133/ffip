#include <simulation.hpp>
#include <iostream>

using namespace ffip;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    std::fstream file{"config.json", std::ios::in};
    json config;
    
    file >> config;
    Simulation sim;
    sim.init(config);
    sim.run(config.at("stop_condition"));
    sim.output();
}
