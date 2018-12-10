/* standard configuration read
 */
#include <iostream>
#include <simulation.hpp>

using namespace std;
using namespace ffip;

int read_basic_config(istream& fin, Simulation& sim) {
	char c;
	fin >> c;
	double dt, dx;
	int dimx, dimy, dimz;
	double er, ur;
	int time_step;
	int PML_thickness;
	double sigma_max;
	
	if (c != '{')
		throw runtime_error("basic format is not right");
	
	fin >> dt >> dx >> dimx >> dimy >> dimz >> time_step >> er >> ur >> PML_thickness;
	sim.setup(dx, dt, {dimx, dimy, dimz});
	sim.set_background_medium(make_medium(er, 0, ur, 0));
	sigma_max = PML::optimal_sigma_max(3, dx, er, ur);
	
	if (PML_thickness) {
		sim.add_PML_layer(new PML(X, High, PML_thickness, sigma_max));
		sim.add_PML_layer(new PML(X, Low,  PML_thickness, sigma_max));
		sim.add_PML_layer(new PML(Y, High, PML_thickness, sigma_max));
		sim.add_PML_layer(new PML(Y, Low,  PML_thickness, sigma_max));
		sim.add_PML_layer(new PML(Z, High, PML_thickness, sigma_max));
		sim.add_PML_layer(new PML(Z, Low,  PML_thickness, sigma_max));
	}
	
	fin >> c;
	if (c != '}')
		throw runtime_error("basic format is not right");
	
	return time_step;
}

Pole_Base* read_pole(istream& fin) {
	string type;
	double rel_perm, freq, damp, inv_relaxation, relaxation;
	char c;
	Pole_Base* res{nullptr};
	
	fin >> c;
	if (c != '{')
		throw runtime_error("pole format is not right");
	
	fin >> type;
	if (type == "Lorentz") {
		fin >> rel_perm >> freq >> damp;
		res =  new Lorentz_Pole(rel_perm, freq, damp);
	}

	if (type == "Drude") {
		fin >> freq >> inv_relaxation;
		res =  new Drude_Pole(freq, inv_relaxation);
	}
	
	if (type == "Deybe") {
		fin >> rel_perm >> relaxation;
		res =  new Deybe_Pole(rel_perm, relaxation);
	}
	
	fin >> c;
	if (c != '}' || res == nullptr)
		throw runtime_error("pole format is not right");
	
	return res;
}

Medium* read_medium(istream& fin) {
	double er, sigma_e, ur, sigma_u;
	Medium* res{nullptr};
	char c;
	
	fin >> c;
	if (c != '{')
		throw runtime_error("medium format is not right");
	
	fin >> er >> sigma_e >> ur >> sigma_u;
	res = make_medium(er, sigma_e, ur, sigma_u);
	
	int n;
	fin >> n;
	for(int i = 0; i < n;++i)
		res->add_e_poles(read_pole(fin));
	
	fin >> c;
	if (c != '}' || res == nullptr)
		throw runtime_error("medium format is not right");
	
	return res;
}

void read_geometry(istream& fin, Simulation& sim, vector<Medium*>& medium_gather) {
	char c;
	string type;
	fin >> c;
	if (c != '{')
		throw runtime_error("Geometry format is not right");
	
	fin >> type;
	
	if (type == "inhom") {
		int idx1, idx2;
		string filename;
		fin >> idx1 >> idx2 >> filename;
		sim.add_solid(make_solid(medium_gather[idx1], medium_gather[idx2], filename));
	}
	
	if (type == "sphere") {
		int idx;
		double radius, x, y, z;
		fin >> idx >> radius >> x >> y >> z;
		sim.add_solid(make_solid(medium_gather[idx], make_sphere(fVec3{x, y, z}, radius)));
	}
	
	fin >> c;
	if (c != '}')
		throw runtime_error("Geometry format is not right");
}

void read_source(istream& fin, Simulation& sim) {
	string type;
	double fp, d;
	int n;
	char c;
	
	fin >> c;
	if (c != '{')
		throw runtime_error("medium format is not right");
	
	fin >> type;
	if (type == "eigen") {
		fin >> n >> fp >> d;
		
		Medium const* bg_medium = sim.get_bg_medium();
		Plane_Wave projector(sim.get_dx(), sim.get_dt(), n);
		auto ricker_source = Rickerwavelet_Func(fp, d);
		
		projector.set_medium(bg_medium->get_e_inf(), bg_medium->get_u_inf());
		projector.set_PML(PML(Direction::Z, High, 10, PML::optimal_sigma_max(3, sim.get_dx(), bg_medium->get_e_inf(), bg_medium->get_u_inf())));
		projector.set_excitation(ricker_source.get_functor());
		
		sim.add_source(new Eigen_Source(projector));
	}
	
	if (type == "dipole") {
		string filename;
		int n;double x, y, z, amp;
		int ctype;
		
		fin >> filename;
		
		fstream dipole_file{filename, ios::in};
		dipole_file >> n;
		for (int i = 0; i < n; ++i) {
			dipole_file >> x >> y >> z >> amp >> fp >> d >> ctype;
			auto ricker_source = Rickerwavelet_Func(fp, d);
			GriddedInterp interp({1, 1, 1}, {x, y, z}, {0, 0, 0}, {amp});
			sim.add_source(new Current_Source(interp, ricker_source.get_functor(), (Coord_Type)ctype));
		}
	}
	
	fin >> c;
	if (c != '}')
		throw runtime_error("medium format is not right");
}

void read_probe(istream& fin, Simulation& sim) {
	int n;
	double x, y, z, freq;
	
	fin >> n;
	for (int i = 0; i < n; ++i) {
		fin >> x >> y >> z >> freq;
		sim.add_probe(new Probe_Frequency({x, y, z}, freq));
	}
}

void read_farfield(istream& fin, Simulation& sim) {
	int n;
	double th, phi, rho, freq;
	
	fin >> n;
	for (int i = 0; i < n; ++i) {
		fin >> th >> phi >> rho >> freq;
		sim.add_farfield_probe(freq, {th, phi, rho});
	}
}

int main(int argc, char const *argv[]) {
	int time_step;
	Simulation sim;
	vector<Medium*> medium_gather;
	fstream fin{"config.in", ios::in};
	fstream fo{"output.out", ios::out};
	fstream probes_output_file;
	fstream probes_input_file;
	fstream farfield_input_file;
	fstream farfield_output_file;
	
	string field;
	fin >> field;
	while(!fin.eof()) {
		bool assigned = false;
		if (field == "basic") {
			time_step = read_basic_config(fin, sim);
			assigned = true;
		}
		
		if (field == "medium") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				medium_gather.push_back(read_medium(fin));
			assigned = true;
			
			fin >> c;
		}
		
		if (field == "geometry") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				read_geometry(fin, sim, medium_gather);
			assigned = true;
			
			fin >> c;
		}
		
		if (field == "source") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				read_source(fin, sim);
			assigned = true;
			
			fin >> c;
		}
		
		if (field == "probe") {
			string filename;
			fin >> filename;
			probes_input_file = fstream{filename, ios::in};
			fin >> filename;
			probes_output_file = fstream{filename, ios::out};
			read_probe(probes_input_file, sim);
			
			assigned = true;
		}
		
		if (field == "farfield") {
			string filename;
			fin >> filename;
			farfield_input_file = fstream{filename, ios::in};
			fin >> filename;
			farfield_output_file = fstream{filename, ios::out};
			read_farfield(farfield_input_file, sim);
			
			assigned = true;
		}
		
		if (!assigned)
			break;
		fin >> field;
	}
	
	sim.init();
	for(int i = 0; i < time_step; ++i) {
		sim.advance(fo);
	}
	
	sim.output(probes_output_file);
	sim.output_farfield(farfield_output_file);

	return 0;
}