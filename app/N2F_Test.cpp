/* test far field from a Hertzian dipole
 */
#include <iostream>
#include <simulation.hpp>

using namespace std;
using namespace ffip;

void check_fstream(fstream& fs) {
	if (!fs.is_open())
		throw runtime_error("fail to open the file");
}

int main(int argc, char const *argv[]) {
	double ur = 1;
	double er = 1;
	double c = c0 / sqrt(ur * er);
	
	double courant = 1 / sqrt(3);
	double dt = 7.2e-17;
	double dx = c * dt / courant;
	iVec3 dim{10, 10, 10};
	
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
	
	//add Ez dipole source in the center
	int Np = 20;
	double fp = c / (Np * dx);
	auto ricker_source = Rickerwavelet_Func(fp, 1/fp);

	auto phase = ricker_source.get_functor();
	GriddedInterp interp({1, 1, 1}, dim * (dx / 2), {0, 0, 0}, {1.0});
	sim.add_source(new Current_Source(interp, phase, Ez));
	
	//add probes
	fstream fin{"request.in", ios::in};
	check_fstream(fin);
	
	int r;
	fin >> r;
	for (int i = 0; i < r; ++i) {
		double f, th, phi, rho;
		fin >> f >> th >> phi >> rho;
		sim.add_farfield_probe(f, {th, phi, rho});
	}
	fin.close();
	
	//set background materials
	auto medium1 = make_medium(er, 0, ur, 0);
	sim.set_background_medium(medium1);
	
	//initialization of simulation
	sim.init();
	
	fstream fo{"data.out", ios::out};
	check_fstream(fo);
	
	//run simulation
	for(int i = 0; i < step; ++i) {
		sim.advance(fo, 4);
	}
	
	sim.udf_output();
	sim.output_farfield(fo);
	//output results
	
	fo.close();
	return 0;
}
