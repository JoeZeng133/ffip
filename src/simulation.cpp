#include <simulation.hpp>
#include <iostream>

namespace ffip {
	Simulation::Simulation(const real _dx, const real _dt, const iVec3 _dim): dt(_dt), dx(_dx), sim_dim(_dim) {}
	
	void Simulation::setup(const real _dx, const real _dt, const iVec3 _dim) {
		dx = _dx;
		dt = _dt;
		sim_dim = _dim;
	}

	void Simulation::probe_init() {
		for(int i = 0; i < NF_freq.size(); ++i)
			probes.push_back(new Probe_Frequency(NF_freq[i], NF_pos[i], chunk));
	}
	

	void Simulation::N2F_init() {
		if (N2F_pos.empty()) return;
		n2f_box = new N2F_Box(fVec3(N2F_p1), fVec3(N2F_p2), N2F_omega, chunk, bg_medium->get_c());
	}

	void Simulation::medium_init() {
		if (bg_medium == nullptr)
			bg_medium = make_medium(1, 0, 1, 0);
		
		prepare_medium(dt);
		const int N = 3;
		real delta = dx / (2 * N);
		std::vector<fVec3> sampled_points;
		for(int i = 1; i < 2 * N; i += 2)
			for(int j = 1; j < 2 * N; j += 2)
				for(int k = 1; k < 2 * N; k += 2) {
					sampled_points.push_back({i * delta, j * delta, k * delta});
				}
	
		auto func = [&, this](my_iterator itr) {
			for(; !itr.is_end(); itr.advance()) {
				/* assign material only when it is E or H*/
				auto p = itr.get_vec();
				auto ctype = p.get_type();
				
				if (!is_H_point(ctype) && !is_E_point(ctype))			//exclude non-material points
					continue;
				
				/* spatial averaging to get materials
				 naive sampling integration over a cube
				 can be improved using adaptive quadrature
				 and ray tracing
				 */
				fVec3 box_p1{(p.x - 1) * dx / 2, (p.y - 1) * dx / 2, (p.z - 1) * dx / 2};
				auto weights = get_zero_weights();
				
				for(auto& item : sampled_points) {
					auto sampled_point = box_p1 + item;
					bool assigned = 0;
					
					for(auto item : solids) {
						if (item->update_weights(sampled_point, weights)) {	//it is true when it is inside the solid
							assigned = 1;
							break;
						}
					}
					
					if (!assigned)			//assign background medium;
						weights[bg_medium->index] += 1;
				}
				
				if (is_E_point(ctype))
					chunk->set_medium_point(p, get_medium_ref(1, weights));
				else
					chunk->set_medium_point(p, get_medium_ref(0, weights));
			}
		};
		
		std::vector<std::thread> threads;
		for(int i = 1; i < num_proc; ++i)
			threads.push_back(std::thread{func, my_iterator(ch_p1, ch_p2, Null, i, num_proc)});
		
		func(my_iterator(ch_p1, ch_p2, Null, 0, num_proc));
		
		for(auto& item : threads)
			item.join();
	}
	
	void Simulation::source_init() {
		for(auto item : current_sources) {
			item->init(dx, ch_p1, ch_p2, chunk->get_dim(), chunk->get_origin());
			chunk->add_source_internal(item->get_source_internal());
		}
		
		for(auto item : eigen_sources) {
			item->init(tf_p1, tf_p2, chunk->get_dim(), chunk->get_origin(), ch_p1, ch_p2, dx);
			chunk->add_source_internal(item->get_source_internal());
		}
	}
	
	void Simulation::PML_init() {
		PML_init_helper(PMLs[0][0], PMLs[0][1], kx, bx, cx, sim_p1.x, sim_p2.x);
		PML_init_helper(PMLs[1][0], PMLs[1][1], ky, by, cy, sim_p1.y, sim_p2.y);
		PML_init_helper(PMLs[2][0], PMLs[2][1], kz, bz, cz, sim_p1.z, sim_p2.z);
		
		chunk->PML_init(kx, ky, kz, bx, by, bz, cx, cy, cz);
	}
	
	void Simulation::PML_init_helper(const PML& neg, const PML& pos, real_arr& k, real_arr& b, real_arr& c, const int p1, const int p2) {
		size_t dim = p2 - p1;		//[0, dim]
		
		k.resize(dim + 1, 1);
		b.resize(k.size(), 1);
		c.resize(k.size(), 0);

		for(int i = 0; i < 2 * neg.get_d(); ++i) {
			real x = neg.get_d() - i / 2.0;
			
			k[i] = neg.get_k(x);
			b[i] = neg.get_b(x, dt);
			c[i] = neg.get_c(x, dt);
		}
		
		for(int i = 0; i < 2 * pos.get_d(); ++i) {
			real x = pos.get_d() - i / 2.0;
			
			k[dim - i] = pos.get_k(x);
			b[dim - i] = pos.get_b(x, dt);
			c[dim - i] = pos.get_c(x, dt);
		}
	}
	
	void Simulation::chunk_init() {
		//dimension is in computational units, so they are two times the original values
		sim_dim = sim_dim * 2;
		
		sim_p1 = {0, 0, 0};
		sim_p2 = sim_dim;
		
		//add TFSF faces
		if (!eigen_sources.empty()) {
			tf_p1 = sim_p1;
			tf_p2 = sim_p2;
			sim_p1 = sim_p1 - iVec3{ 1, 1, 1 };
			sim_p2 = sim_p2 + iVec3{ 1, 1, 1 };
		}

		//implementations of N2F faces
		if (!N2F_pos.empty()) {
			N2F_p1 = sim_p1 - iVec3{ 1, 1, 1 };
			N2F_p2 = sim_p2 + iVec3{ 1, 1, 1 };
			sim_p1 = sim_p1 - iVec3{ 2, 2, 2 };
			sim_p2 = sim_p2 + iVec3{ 2, 2, 2 };
		}

		//add PML layers
		sim_p1.x -= int(2 * PMLs[0][0].get_d());
		sim_p1.y -= int(2 * PMLs[1][0].get_d());
		sim_p1.z -= int(2 * PMLs[2][0].get_d());
		
		sim_p2.x += int(2 * PMLs[0][1].get_d());
		sim_p2.y += int(2 * PMLs[1][1].get_d());
		sim_p2.z += int(2 * PMLs[2][1].get_d());

		//implementaions of MPI, for now 1 chunk covers the whole region
		chunk = new Chunk{sim_p1, sim_p2, sim_p1, sim_p2, dx, dt};
		ch_p1 = chunk->get_p1();
		ch_p2 = chunk->get_p2();

		chunk->set_num_proc(num_proc);
	}
	
	void Simulation::set_num_proc(const int _num_proc) {
		num_proc = _num_proc;
	}
	
	void Simulation::init() {
		chunk_init();
		PML_init();
		source_init();
		medium_init();
		N2F_init();
		probe_init();
		udf_unit();
		step = 0;
		std::cout << "Initialization Complete\n";
	}
	
	void Simulation::add_source(Source *source) {
		auto current_cast = dynamic_cast<Current_Source*>(source);
		if (current_cast) current_sources.push_back(current_cast);
		
		auto eigen_cast = dynamic_cast<Eigen_Source*>(source);
		if (eigen_cast) eigen_sources.push_back(eigen_cast);
	}
	
	void Simulation::add_solid(Solid const* solid) {
		solids.push_back(solid);
	}
	
	void Simulation::add_PML_layer(PML *PML) {
		if (PML->get_side() > 0)
			PMLs[PML->get_dir()][1] = *PML;
		else
			PMLs[PML->get_dir()][0] = *PML;
	}
	
	void Simulation::set_background_medium(Medium const*  m) {
		if (m)
			bg_medium = m;
	}

	void Simulation::add_nearfield_probe(const real freq, const fVec3& v) {
		NF_pos.push_back(v);
		NF_freq.push_back(freq);
	}

	void Simulation::add_farfield_probe(const real freq, const fVec3& pos) {
		N2F_omega.push_back(2 * pi * freq);
		N2F_pos.push_back(pos);
	}
	
	void Simulation::advance(std::ostream& os) {
		std::cout << "\n" << std::setfill('0') << std::setw(4) << step;
		real time = (step ++ ) * dt;
		
		auto func = [&, this](const int rank, const int num_proc) {
			chunk->update_Md(time, rank);
			glob_barrier->Sync();
			chunk->update_m_PML(rank);
			glob_barrier->Sync();
			chunk->update_m_source(time, rank);
			glob_barrier->Sync();
			chunk->update_B2H(time, rank);
			glob_barrier->Sync();

			chunk->update_Jd(time + 0.5 * dt, rank);
			glob_barrier->Sync();
			chunk->update_e_PML(rank);
			glob_barrier->Sync();
			chunk->update_e_source(time + 0.5 * dt, rank);
			glob_barrier->Sync();
			chunk->update_D2E(time + 0.5 * dt, rank);
			glob_barrier->Sync();

			chunk->update_padded_E(time);
			glob_barrier->Sync();
			chunk->update_padded_H(time);
			glob_barrier->Sync();

			size_t idx1, idx2;
			vector_divider(probes, rank, num_proc, idx1, idx2);
			std::for_each(probes.begin() + idx1, probes.begin() + idx2, [this](Probe* item) {item->update(step); });
			
			if (n2f_box)
				n2f_box->update(step, rank, num_proc);
			glob_barrier->Sync();
		};
		
		std::vector<std::future<void>> task_list;
		for(int i = 0; i < num_proc; ++i) {
			task_list.push_back(std::async(std::launch::async, func, i, num_proc));
		}
		for(auto &item : task_list)
			item.get();
		
		os << chunk->measure() << "\n";
		udf_advance();
	}

	real Simulation::at(const fVec3& p, const Coord_Type ctype) const{
		return chunk->at(p, ctype);
	}

	real Simulation::operator()(const fVec3& p, const Coord_Type ctype) const{
		return chunk->operator()(p, ctype);
	}

	real Simulation::operator()(const iVec3& p, const Coord_Type ctype) const {
		return chunk->operator()(p, ctype);
	}

	real Simulation::operator()(const iVec3& p) const {
		return chunk->operator()(p);
	}

	int Simulation::get_step() const{
		return step;
	}
	
	real Simulation::get_dt() const {
		return dt;
	}
	
	real Simulation::get_dx() const {
		return dx;
	}
	
	iVec3 Simulation::get_dim() const {
		return sim_dim / 2;
	}
	
	Medium const* Simulation::get_bg_medium() const {
		return bg_medium;
	}
	
	void Simulation::output_nearfield(std::ostream &os) {
		for(auto item : probes) {
			item->output(os);
		}
	}
	
	void Simulation::output_farfield(std::ostream &os) {
		if (n2f_box == nullptr) return;
		n2f_box->prepare();

		real imped = bg_medium->get_z();
		os << std::scientific;

		for (int i = 0; i < N2F_pos.size(); ++i) {
			real omega = N2F_omega[i];
			real k = omega / bg_medium->get_c();		//get the wavenumber
			real th = N2F_pos[i].x;
			real phi = N2F_pos[i].y;
			real rho = N2F_pos[i].z;
			fVec3 proj_th = {cos(th) * cos(phi), cos(th) * sin(phi), -sin(th)};
			fVec3 proj_phi = {-sin(phi), cos(phi), 0};
			Vec3<complex_num> N{0, 0, 0}, L{0, 0, 0};

			auto tmp = n2f_box->get_NL(th, phi, omega);

			complex_num Nth = inner_prod(proj_th, tmp.first);
			complex_num Nphi = inner_prod(proj_phi, tmp.first);
			complex_num Lth = inner_prod(proj_th, tmp.second);
			complex_num Lphi = inner_prod(proj_phi, tmp.second);

			complex_num arg = k / (4 * pi * rho) * complex_num{sin(k * rho), cos(k * rho)};
			complex_num Eth = -(Lphi + imped * Nth) * arg;
			complex_num Ephi = (Lth - imped * Nphi) * arg;
			complex_num Hth = (Nphi - Lth / imped) * arg;
			complex_num Hphi = -(Nth + Lphi / imped) * arg;

			os << Eth << " " << Ephi << " " << Hth << " " << Hphi << std::endl;
		}
	}
	
	void Simulation::udf_unit() {
	}
	
	void Simulation::udf_advance() {
	}

	void Simulation::udf_output() {
	}
}

