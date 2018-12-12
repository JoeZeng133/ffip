#include <simulation.hpp>
#include <iostream>

namespace ffip {
	Simulation::Simulation(const real _dx, const real _dt, const iVec3 _dim): dt(_dt), dx(_dx), sim_dim(_dim) {}
	
	void Simulation::setup(const real _dx, const real _dt, const iVec3 _dim) {
		dx = _dx;
		dt = _dt;
		sim_dim = _dim;
	}

	void Simulation::probe_init() {}

	void Simulation::N2F_init() {
		if (N2F_pos.empty()) return;

		//make frequency unique
		N2F_omega_unique = N2F_omega;
		std::sort(N2F_omega_unique.begin(), N2F_omega_unique.end());
		auto it = std::unique(N2F_omega_unique.begin(), N2F_omega_unique.end());
		N2F_omega_unique.resize(std::distance(N2F_omega_unique.begin(), it));
		
		//generate six N2F_surfaces
		auto center = (N2F_p1 + N2F_p2) * (dx / 4);
		std::pair<iVec3, iVec3> faces[6];
		faces[0] = get_face<dir_x_tag, side_low_tag>(N2F_p1, N2F_p2);
		faces[1] = get_face<dir_x_tag, side_high_tag>(N2F_p1, N2F_p2);
		faces[2] = get_face<dir_y_tag, side_low_tag>(N2F_p1, N2F_p2);
		faces[3] = get_face<dir_y_tag, side_high_tag>(N2F_p1, N2F_p2);
		faces[4] = get_face<dir_z_tag, side_low_tag>(N2F_p1, N2F_p2);
		faces[5] = get_face<dir_z_tag, side_high_tag>(N2F_p1, N2F_p2);
		
		for(int i = 0; i < 6; ++i) {
			faces[i] = get_intersection(faces[i].first, faces[i].second, ch_p1, ch_p2);
			faces[i] = get_component_interior(faces[i].first, faces[i].second, N2F_p1.get_type());
		}
		
		N2F_faces = {
			new N2F_Face<dir_x_tag>{ N2F_p1, N2F_p2, faces[0].first, faces[0].second, Side::Low, N2F_omega_unique, dx, dt, center, bg_medium},
			new N2F_Face<dir_x_tag>{ N2F_p1, N2F_p2, faces[1].first, faces[1].second, Side::High, N2F_omega_unique, dx, dt, center, bg_medium},
			new N2F_Face<dir_y_tag>{ N2F_p1, N2F_p2, faces[2].first, faces[2].second, Side::Low, N2F_omega_unique, dx, dt, center, bg_medium},
			new N2F_Face<dir_y_tag>{ N2F_p1, N2F_p2, faces[3].first, faces[3].second, Side::High, N2F_omega_unique, dx, dt, center, bg_medium},
			new N2F_Face<dir_z_tag>{ N2F_p1, N2F_p2, faces[4].first, faces[4].second, Side::Low, N2F_omega_unique, dx, dt, center, bg_medium},
			new N2F_Face<dir_z_tag>{ N2F_p1, N2F_p2, faces[5].first, faces[5].second, Side::High, N2F_omega_unique, dx, dt, center, bg_medium}
		};
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
	
		for(auto itr = my_iterator(ch_p1, ch_p2, Null); !itr.is_end(); itr.advance()) {
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
		sim_p1.x -= 2 * PMLs[0][0].get_d();
		sim_p1.y -= 2 * PMLs[1][0].get_d();
		sim_p1.z -= 2 * PMLs[2][0].get_d();
		
		sim_p2.x += 2 * PMLs[0][1].get_d();
		sim_p2.y += 2 * PMLs[1][1].get_d();
		sim_p2.z += 2 * PMLs[2][1].get_d();

		//implementaions of MPI, for now 1 chunk covers the whole region
		chunk = new Chunk{sim_p1, sim_p2, sim_p1, sim_p2, dx, dt};
		ch_p1 = chunk->get_p1();
		ch_p2 = chunk->get_p2();
	}
	
	void Simulation::init() {
		chunk_init();
		PML_init();
		source_init();
		medium_init();
		N2F_init();
		udf_unit();
		step = 0;
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

	void Simulation::add_probe(Probe* probe) {
		probes.push_back(probe);
	}

	void Simulation::add_farfield_probe(const real freq, const fVec3& pos) {
		N2F_omega.push_back(2 * pi * freq);
		N2F_pos.push_back(pos);
	}
	
	void Simulation::advance(std::ostream& os, const int num_proc) {
		std::cout << "Stepping from" << step << " to " << step + 1 << std::endl;
		real time = (step ++ ) * dt;
		chunk->update_Md(time, num_proc);
		chunk->update_B2H(time, num_proc);
		
		chunk->update_Jd(time + 0.5 * dt, num_proc);
		chunk->update_D2E(time + 0.5 * dt, num_proc);
		
		chunk->update_padded_E(time);
		chunk->update_padded_H(time);
		
		auto probes_update = [&, this](Probe* item) {
			item->update(*this);
		};
		
		task_divider(probes, probes_update, num_proc);
		
		auto N2F_update = [&, this](N2F_Face_Base* item) {
			item->update(chunk, step);
		};
		
		task_divider(N2F_faces, N2F_update, num_proc);
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
	
	void Simulation::output(std::ostream &os) {
		for(auto item : probes) {
			item->output(os);
		}
	}
	
	void Simulation::output_farfield(std::ostream &os) {
		//get impedance of the non-absorbing background medium
		real imped = bg_medium->get_z();
		os << std::scientific;
		
		for(int i = 0; i < N2F_pos.size(); ++i) {
			real omega = N2F_omega[i];
			real k = omega / bg_medium->get_c();		//get the wavenumber
			real th = N2F_pos[i].x;
			real phi = N2F_pos[i].y;
			real rho = N2F_pos[i].z;
			fVec3 proj_th = {cos(th) * cos(phi), cos(th) * sin(phi), -sin(th)};
			fVec3 proj_phi = {-sin(phi), cos(phi), 0};
			Vec3<complex_num> N{0, 0, 0}, L{0, 0, 0};
			
			for(auto item : N2F_faces) {
				auto tmp = item->get_NL(th, phi, omega);
				N = N + tmp.first;
				L = L + tmp.second;
			}
			
			complex_num Nth = inner_prod(proj_th, N);
			complex_num Nphi = inner_prod(proj_phi, N);
			complex_num Lth = inner_prod(proj_th, L);
			complex_num Lphi = inner_prod(proj_phi, L);
			
			complex_num arg = k / (4 * pi * rho) * complex_num{sin(k * rho), cos(k * rho)};
			complex_num Eth = -(Lphi + imped * Nth) * arg;
			complex_num Ephi = (Lth - imped * Nphi) * arg;
			complex_num Hth = (Nphi - Lth / imped) * arg;
			complex_num Hphi = -(Nth + Lphi / imped) * arg;
			
			os << Eth << " " << Ephi << " " << Hth << " " << Hphi << std::endl;
		}
	}
	
	std::fstream os[6];
	
	void Simulation::udf_unit() {
		/*for(int i = 0; i < 6; ++i) {
			std::string idx = std::to_string(i);
			os[i].open("face" + idx + ".out", std::fstream::out);
		}
		
		for(int i = 0; i < 6; ++i) {
			os[i] << N2F_faces[i]->get_p1() << "\n" << N2F_faces[i]->get_p2() << "\n";
			os[i] << N2F_faces[i]->get_norm_vec() << "\n";
			os[i] << std::scientific;
		}*/
	}
	
	void Simulation::udf_advance() {

		/* output fields at N2F face each (position, time step)
		for(int i = 0; i < 6; ++i) {
			auto face = N2F_faces[i];
			
			for(auto itr = my_iterator(face->get_p1(), face->get_p2(), face->get_p1().get_type()); !itr.is_end(); itr.advance()) {
				auto pos = itr.get_vec();
				os[i] << (*chunk)(pos, Ex) << " "
				<< (*chunk)(pos, Ey) << " "
				<< (*chunk)(pos, Ez) << " "
				<< (*chunk)(pos, Hx) << " "
				<< (*chunk)(pos, Hy) << " "
				<< (*chunk)(pos ,Hz) << "\n";
				
			}
		}*/
	}

	void Simulation::udf_output() {
	/*	output currents at each (position, frequency)
		for (int i = 0; i < 6; ++i)
			N2F_faces[i]->output_JM(os[i]);*/
	}
}

