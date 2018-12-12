/* Mie test
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
	iVec3 dim{50, 50, 50};
	
	Plane_Wave projector(dx, dt, dim.z);
	Simulation sim(dx, dt, dim);
	int step = 1500;
	
	//add PML layers
	double sigma_max = PML::optimal_sigma_max(3, dx, er, ur);
	
	sim.add_PML_layer(new PML(X, High, 6, sigma_max));
	sim.add_PML_layer(new PML(X, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Y, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Y, Low,  6, sigma_max));
	sim.add_PML_layer(new PML(Z, High, 6, sigma_max));
	sim.add_PML_layer(new PML(Z, Low,  6, sigma_max));
	
	projector.set_PML(PML(Z, High, 6, sigma_max));
	
	//add plane wave source
	int Np = 30;
	double fp = c / (Np * dx);
	auto ricker_source = Rickerwavelet_Func(fp, 1/fp);
	auto sin_source = Sinuosuidal_Func(fp);
	auto phase = ricker_source.get_functor();
	projector.set_medium(er, ur);
	projector.set_excitation(phase);
	sim.add_source(new Eigen_Source(projector));
	
	//set background medium
	auto bg_medium = make_medium(er, 0, ur, 0);
	sim.set_background_medium(bg_medium);
	
	//set sphere and its medium
	auto sph_medium = make_medium(er, 0, ur, 0);
	sph_medium->add_e_poles(new Lorentz_Pole(0.8, 4e16, 1e16), new Lorentz_Pole(0.5, 6e16, 1e16));
	auto sph1 = make_sphere(dim * (dx / 2), 10 * dx);
	//sim.add_solid(make_solid(sph_medium, sph1));						//geometry model
	sim.add_solid(make_solid(sph_medium, bg_medium,  "objective_geometry.in"));	//inhomgeneous model
	
	//read far field probes
	fstream fin{"probes.in", ios::in};
	if (!fin.is_open())
		throw runtime_error("fail to open request.in");
	
	int r;
	fin >> r;
	for (int i = 0; i < r; ++i) {
		double f, th, phi, rho;
		fin >> th >> phi >> rho >> f;
		sim.add_farfield_probe(f, {th, phi, rho});
	}
	fin.close();
	
	//initializations
	sim.init();
	projector.init();
	
	fstream fo{"output.out", ios::out};
	//run simulation
	for(int i = 0; i < step; ++i) {
		sim.advance(fo);
	}
	
	//output results
	sim.output_farfield(fo);
	
	fo.close();
	
	return 0;
}
