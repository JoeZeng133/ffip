/* test plane wave projector
 */
#include <iostream>
#include <simulation.hpp>

using namespace std;
using namespace ffip;

int main(int argc, char const *argv[]) {
	double ur = 1;
	double er = 1;
	double c = c0 / sqrt(ur * er);
	
	double courant = 1 / sqrt(3);
	double dt = 2e-17 / 50;
	double dx = c * dt / courant;
	iVec3 dim{30, 30, 30};
	
	Plane_Wave projector(dx, dt, dim.z);
	Simulation sim(dx, dt, dim);
	int step = 300;
	
	//add PML layers
	double sigma_max = PML::optimal_sigma_max(3, dx, er, ur);
	
	sim.add_PML_layer(new PML(X, High, 6, sigma_max));
	sim.add_PML_layer(new PML(X, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Y, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Y, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Z, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Z, Low,  6, sigma_max));
	
	projector.set_PML(PML(Z, High, 6, sigma_max));
	
	//add dipole source
	int Np = 15;
	ffip::real fp = c / (Np * dx);
	auto ricker_source = Rickerwavelet_Func(fp, 1 / fp);
	auto sin_source = Sinuosuidal_Func(fp);
	
	auto phase = sin_source.get_functor();
	GriddedInterp interp({1, 1, 1}, fVec3{15, 15, 15} * dx, {0, 0, 0}, {1.0});
	projector.set_medium(er, ur);
	projector.set_excitation(phase);
	sim.add_source(new Eigen_Source(projector));
	

	//set background materials
	auto medium1 = make_medium(er, 0, ur, 0);
	sim.set_background_medium(medium1);
	
	//initialization of simulation
//	projector.init();
	sim.init();
	
	fstream fo{"data.out", ios::out};
	if (!fo.is_open())
		throw runtime_error("fail to open data.out");
	
	//run simulation
	for(int i = 0; i < step; ++i) {
		sim.advance(fo, 4);
	}
	//output results
	
	fo.close();
	return 0;
}
