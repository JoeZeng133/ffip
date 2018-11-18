/* far field dipole test
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

	Simulation sim(dx, dt, dim);
	int step = 600;

	//add PML layers
	double sigma_max = PML::optimal_sigma_max(3, dx, er, ur);

	sim.add_PML_layer(new PML(X, High, 6, sigma_max));
	sim.add_PML_layer(new PML(X, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Y, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Y, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Z, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Z, Low,  6, sigma_max));

	//add dipole source
	int Np = 20;
	auto ricker_source = Rickerwavelet_Func(c / (Np * dx), 0);

	auto phase = ricker_source.get_functor();
	GriddedInterp interp({1, 1, 1}, fVec3{15, 15, 15} * dx, {0, 0, 0}, {1.0});
	sim.add_source(new Current_Source(interp, phase, Ez));

	//add probes
	fstream fin{"request.in", ios::in};
	if (!fin.is_open())
		throw runtime_error("fail to open request.in");

	int r;
	fin >> r;
	for (int i = 0; i < r; ++i) {
		double x, y, z, f;
		fin >> x >> y >> z >> f;
		sim.add_probe(new Probe_Frequency{{x, y, z}, f});
	}
	fin.close();

	//set background materials
	sim.set_background_medium(make_medium(er, 0, ur, 0));


	//initialization of simulation
	sim.init();

	fstream fo{"data.out", ios::out};
	if (!fo.is_open())
		throw runtime_error("fail to open data.out");

	//run simulation
	for(int i = 0; i < step; ++i) {
		sim.advance(fo);
	}

	//output results
	sim.output(fo);

	fo.close();
    return 0;
}
