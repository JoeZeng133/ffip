/* standard configuration read
 */
#include <iostream>
#include <simulation.hpp>
#include <chrono>

using namespace std;
using namespace ffip;

int read_basic_config(istream& fin, Simulation& sim) {
	char c;
	double dt, dx;
	int dimx, dimy, dimz;
	int sf_thickness;
	int PML_thickness;
	int time_step;
	double er, ur;

	fin >> c;
	if (c != '{')
		throw runtime_error("basic format is not right");
	
	fin >> dt >> dx >> dimx >> dimy >> dimz >> sf_thickness >> PML_thickness >> time_step >> er >> ur;
	sim.setup(dx, dt, {dimx, dimy, dimz});
	sim.set_background_medium(sim.make_medium(er, 0, ur, 0));
	
	if (PML_thickness) {
		sim.add_PML_layer(sim.make_pml(PML_thickness), X, High);
		sim.add_PML_layer(sim.make_pml(PML_thickness), X, Low);
		sim.add_PML_layer(sim.make_pml(PML_thickness), Y, High);
		sim.add_PML_layer(sim.make_pml(PML_thickness), Y, Low);
		sim.add_PML_layer(sim.make_pml(PML_thickness), Z, High);
		sim.add_PML_layer(sim.make_pml(PML_thickness), Z, Low);
	}
	sim.add_sf_layer(sf_thickness);
	sim.chunk_init();
	
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
	
	if (type == "CP") {
		double A, phi, Omega, Gamma;
		fin >> A >> phi >> Omega >> Gamma;
		res = new CP_Pole{A, phi, Omega, Gamma};
	}
	
	fin >> c;
	if (c != '}' || res == nullptr)
		throw runtime_error("pole format is not right");
	
	return res;
}

Medium* read_medium(istream& fin, Simulation& sim) {
	double er, sigma_e, ur, sigma_u;
	Medium* res{nullptr};
	char c;
	
	fin >> c;
	if (c != '{')
		throw runtime_error("medium format is not right");
	
	fin >> er >> sigma_e >> ur >> sigma_u;
	res = sim.make_medium(er, sigma_e, ur, sigma_u);
	
	int n;
	fin >> n;
	for(int i = 0; i < n;++i)
		res->add_e_poles(read_pole(fin));
	
	fin >> c;
	if (c != '}' || res == nullptr)
		throw runtime_error("medium format is not right");
	
	return res;
}

void read_geometry(istream& fin, Simulation& sim, const vector<Medium*>& medium_gather) {
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
		sim.add_solid(sim.make_solid(medium_gather[idx1], medium_gather[idx2], filename));
	}
	
	if (type == "sphere") {
		int idx;
		double radius, x, y, z;
		fin >> idx >> radius >> x >> y >> z;
		sim.add_solid(sim.make_solid(medium_gather[idx], sim.make_sphere(fVec3{x, y, z}, radius)));
	}

	if (type == "box") {
		int idx;
		double x0, y0, z0, x1, y1, z1;
		fin >> idx >> x0 >> y0 >> z0 >> x1 >> y1 >> z1;
		sim.add_solid(sim.make_solid(medium_gather[idx], sim.make_box(fVec3{ x0, y0, z0 }, fVec3{ x1, y1, z1 })));
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
	if (type == "plane") {
		fin >> n >> fp >> d;
		
		Medium const* bg_medium = sim.get_bg_medium();
		Plane_Wave projector(sim.get_dx(), sim.get_dt(), n);
		
		projector.set_medium(bg_medium->get_e_inf(), bg_medium->get_u_inf());
		projector.set_PML(PML(6, 0.8 * 4 / (sim.get_dx() * z0), 1, 0, 3, 1));
		projector.set_excitation(make_ricker_func(fp, d));
		
		sim.add_inc_source(new Inc_Source(projector));
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
			sim.add_dipole(amp, {x, y, z}, (Coord_Type)ctype, make_ricker_func(fp, d));
		}
	}
	
	fin >> c;
	if (c != '}')
		throw runtime_error("medium format is not right");
}

vector<Nearfield_Probe const*> nearfield_probes;
void read_nearfield_probe(istream& fin, Simulation& sim) {
	int n;
	double x, y, z, freq;
	
	fin >> n;
	for (int i = 0; i < n; ++i) {
		fin >> x >> y >> z >> freq;
		nearfield_probes.push_back(sim.add_nearfield_probe(freq, {x, y, z}));
	}
}

vector<fVec3> farfield_pos;
vector<double> farfield_freq;

N2F_Box const* read_farfield_probe(istream& fin, Simulation& sim) {
	int n;
	fVec3 p1, p2;
	double th, phi, rho, freq;
	
	fin >> p1.x >> p1.y >> p1.z;
	fin >> p2.x >> p2.y >> p2.z;
	fin >> n;
	for (int i = 0; i < n; ++i) {
		fin >> th >> phi >> rho >> freq;
		farfield_pos.emplace_back(th, phi, rho);
		farfield_freq.push_back(freq);
	}
	
	return sim.add_n2f_box(p1, p2, farfield_freq);
}

vector<double> flux_freq;
Flux_Box const* read_flux(istream& fin, Simulation& sim) {
	int n;
	fVec3 p1, p2;
	real_arr freq_list;
	
	fin >> p1.x >> p1.y >> p1.z;
	fin >> p2.x >> p2.y >> p2.z;
	fin >> n;
	for(int i = 0; i < n; ++i) {
		double freq;
		fin >> freq;
		flux_freq.push_back(freq);
	}
	
	return sim.add_flux_box(p1, p2, flux_freq);
}

int main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	size_t num_threads = thread::hardware_concurrency();
	set_num_proc(num_threads);

	int time_step;
	Simulation sim;
	sim.set_num_proc(num_threads);

	vector<Medium*> medium_gather;
	fstream fin{"config.in", ios::in};
	fstream fo{"log.out", ios::out};

	if (!fin.is_open())
		throw runtime_error("fail to open configuration file");

	if (!fo.is_open())
		throw runtime_error("fail to open log file");

	fstream nearfield_output_file;
	fstream nearfield_input_file;
	fstream farfield_input_file;
	fstream farfield_output_file;
	fstream flux_output_file;
	fstream flux_input_file;
	
	Flux_Box const* flux_box{nullptr};
	N2F_Box const* n2f_box{nullptr};
	
	string field;
	fin >> field;
	while(!fin.eof()) {
		bool assigned = false;
		if (field == "basic") {
			time_step = read_basic_config(fin, sim);
			assigned = true;
			//cout << field << " Complete \n";
		}
		
		if (field == "medium") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				medium_gather.push_back(read_medium(fin, sim));
			assigned = true;
			fin >> c;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "geometry") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				read_geometry(fin, sim, medium_gather);
			assigned = true;
			
			fin >> c;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "source") {
			int n;
			char c;
			fin >> n >> c;
			
			for(int i = 0; i < n; ++i)
				read_source(fin, sim);
			assigned = true;
			
			fin >> c;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "nearfield") {
			string filename;
			fin >> filename;
			nearfield_input_file = fstream{filename, ios::in};
			fin >> filename;
			nearfield_output_file = fstream{filename, ios::out};
			read_nearfield_probe(nearfield_input_file, sim);
			
			assigned = true;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "farfield") {
			string filename;
			fin >> filename;
			farfield_input_file = fstream{filename, ios::in};
			fin >> filename;
			farfield_output_file = fstream{filename, ios::out};
			n2f_box = read_farfield_probe(farfield_input_file, sim);
			
			assigned = true;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "flux") {
			string filename;
			fin >> filename;
			flux_input_file = fstream{filename, ios::in};
			fin >> filename;
			flux_output_file = fstream{filename, ios::out};
			flux_box = read_flux(flux_input_file, sim);
			
			assigned = true;
			
			//cout << field << " Complete \n";
		}
		
		if (field == "Stop_Step_Output") {
			sim.output_step_number = false;
			
			assigned = true;
			
			//cout << field << "\n";
		}
		
		
		if (!assigned)
			break;
		fin >> field;
	}
	
	sim.init();
	
	
	for(int i = 0; i < time_step; ++i) {
		sim.advance(fo);
	}
	
	sim.udf_output();
	
	if (n2f_box) {
		n2f_box->prepare();
		for(int i = 0; i < farfield_pos.size(); ++i) {
			auto pos = farfield_pos[i];
			n2f_box->output(farfield_output_file, pos.x, pos.y, pos.z, farfield_freq[i]);
		}
	}
	
	for (auto item : nearfield_probes) {
		item->output(nearfield_output_file);
	}
	
	if (flux_box) {
		auto res = flux_box->get();
		for(auto item : res)
			flux_output_file << item << "\n";
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
	//cout << "Simulation Finished" << endl;
	return 0;
}
