#include <simulation.hpp>
#include <iostream>

namespace ffip {
	N2F_Face_Base::N2F_Face_Base(const iVec3& _p1, const iVec3& _p2, const Side _side, const real _omega, const real _dx, const real _dt, const fVec3& _center): p1(_p1), p2(_p2), side(_side), omega(_omega), dx(_dx), dt(_dt), center(_center) {
		if ((p1.x == p2.x) + (p1.y == p2.y) + (p1.z == p2.z) == 2 && ElementWise_Less_Eq(p1, p2))
			is_face = true;
	}
	
	Simulation::Simulation(const real _dx, const real _dt, const iVec3 _dim): dt(_dt), dx(_dx), sim_dim(_dim) {}

	void Simulation::probe_init() {}

	void Simulation::N2F_init() {
		auto center = sim_dim * (dx / 2 / 2);
		iVec3 N2F_p1s[6] = { 
			iVec3{ N2F_p1.x, N2F_p1.y, N2F_p1.z } ,		//x-
			iVec3{ N2F_p2.x, N2F_p1.y, N2F_p1.z } ,		//x+
			iVec3{ N2F_p1.x, N2F_p1.y, N2F_p1.z } ,		//y-
			iVec3{ N2F_p1.x, N2F_p2.y, N2F_p1.z } ,		//y+
			iVec3{ N2F_p1.x, N2F_p1.y, N2F_p1.z } ,		//z-
			iVec3{ N2F_p1.x, N2F_p1.y, N2F_p2.z } };	//z+

		iVec3 N2F_p2s[6] = {
			iVec3{ N2F_p1.x, N2F_p2.y, N2F_p2.z } ,
			iVec3{ N2F_p2.x, N2F_p2.y, N2F_p2.z } ,
			iVec3{ N2F_p2.x, N2F_p1.y, N2F_p2.z } ,
			iVec3{ N2F_p2.x, N2F_p2.y, N2F_p2.z } ,
			iVec3{ N2F_p2.x, N2F_p2.y, N2F_p1.z } ,
			iVec3{ N2F_p2.x, N2F_p2.y, N2F_p2.z } };

		std::pair<iVec3, iVec3> faces[6];

		for (int i = 0; i < 6; ++i) {
			faces[i] = get_intersection(N2F_p1s[i], N2F_p2s[i], ch_p1, ch_p2);
		}

		for (auto omega : N2F_freq) {
			if (N2F_faces.find(omega) == N2F_faces.end()) {
				
				N2F_faces[omega] = {
					new N2F_Face<dir_x_tag>{ faces[0].first, faces[0].second, Side::Low, omega, dx, dt, center },
					new N2F_Face<dir_x_tag>{ faces[1].first, faces[1].second, Side::High, omega, dx, dt, center },
					new N2F_Face<dir_y_tag>{ faces[2].first, faces[2].second, Side::Low, omega, dx, dt, center },
					new N2F_Face<dir_y_tag>{ faces[3].first, faces[3].second, Side::High, omega, dx, dt, center },
					new N2F_Face<dir_z_tag>{ faces[4].first, faces[4].second, Side::Low, omega, dx, dt, center },
					new N2F_Face<dir_z_tag>{ faces[5].first, faces[5].second, Side::High, omega, dx, dt, center } };
				}
			}
	}

	void Simulation::medium_init() {
		if (background_medium_id == -1)
			background_medium_id = make_medium(1, 0, 1, 0);
		
		prepare_medium_internal(dt);
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
			
			if (ctype == Coord_Type::Null || ctype == Coord_Type::Center || ctype == Coord_Type::Corner)			//exclude non-material points
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
					weights[background_medium_id] += 1;
			}
			
			if (ctype == Ex || ctype == Ey || ctype == Ez)
				chunk->set_medium_point(p, get_medium_internal(weights, 1));
			else
				chunk->set_medium_point(p, get_medium_internal(weights, 0));
		}
	}
	
	void Simulation::source_init() {
		for(auto item : current_sources) {
			item->init(dx, ch_p1, ch_p2, chunk->get_dim(), chunk->get_origin());
			chunk->add_source_internal(item->get_source_internal());
		}
		
		for(auto item : eigen_sources) {
			item->init({0, 0, 0}, sim_dim, chunk->get_dim(), chunk->get_origin(), ch_p1, ch_p2);
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
		sim_dim = sim_dim * 2;
		
		// TMz simulation
//		sim_p1 = {0, 0, 1};
//		sim_p2 = {sim_dim.x, sim_dim.y, 1};

		sim_p1 = {0, 0, 0};
		sim_p2 = sim_dim;
		//dimension is in computational units, so they are two times the original values
		
		//add TFSF faces
		if (!eigen_sources.empty()) {
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
		step = 0;
	}
	
	void Simulation::add_source(Source *source) {
		auto current_cast = dynamic_cast<Current_Source*>(source);
		if (current_cast) current_sources.push_back(current_cast);
		
		auto eigen_cast = dynamic_cast<Eigen_Source*>(source);
		if (eigen_cast) eigen_sources.push_back(eigen_cast);
	}
	
	void Simulation::add_solid(Solid *solid) {
		solids.push_back(solid);
	}
	
	void Simulation::add_PML_layer(PML *PML) {
		if (PML->get_side() > 0)
			PMLs[PML->get_dir()][1] = *PML;
		else
			PMLs[PML->get_dir()][0] = *PML;
	}
	
	void Simulation::set_background_medium(const int id) {
		if (id >= 0)
			background_medium_id = id;
	}

	void Simulation::add_probe(Probe* probe) {
		probes.push_back(probe);
	}

	void Simulation::add_farfield_probe(const real freq, const fVec3& p) {
		N2F_freq.push_back(freq);
		N2F_pos.push_back(p);
	}
	
	void Simulation::advance(std::ostream& os) {
		std::cout << "Stepping from" << step << " to " << step + 1 << std::endl;
		real time = (step ++ ) * dt;
		//std::cout << "testing inside simulation" << dt << std::endl;
		chunk->update_Md(time);
		chunk->update_B2H(time);
		
		chunk->update_Jd(time + 0.5 * dt);
		chunk->update_D2E(time + 0.5 * dt);
		
		chunk->update_padded_E(time);
		chunk->update_padded_H(time);

		for (auto item : probes)
			item->update(*this);
		
		
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
	
	void Simulation::output(std::ostream &o) {
		for(auto item : probes) {
			item->output(o);
		}
	}
}

