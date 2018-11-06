#include <simulation.hpp>
#include <iostream>

namespace ffip {
	N2F_Face_Base::N2F_Face_Base(const iVec3& _p1, const iVec3& _p2, const Side _side, const real _omega, const real _dx, const real _dt, const fVec3& _center): p1(_p1), p2(_p2), side(_side), omega(_omega), dx(_dx), dt(_dt), center(_center) {
		if ((p1.x == p2.x) + (p1.y == p2.y) + (p1.z == p2.z) == 2 && ElementWise_Less_Eq(p1, p2))
			is_face = true;
	}
	
	Simulation::Simulation(const real _dx, const real _dt, const iVec3 _dim): dx(_dx), dt(_dt), sim_dim(_dim) {}

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
		const int N = 3;
		if (background_medium == nullptr)
			background_medium = make_medium(1, 0);
											
		for(int x = ch_p1.x; x <= ch_p2.x; ++x)
			for(int y = ch_p1.y; y <= ch_p2.y; ++y)
				for(int z = ch_p1.z; z <= ch_p2.z; ++z) {
					/* assign material only when it is E or H*/
					iVec3 p = {x, y, z};
					if (p.get_type() == Coord_Type::Null || p.get_type() == Coord_Type::Center)
						continue;
					
					/* spatial averaging to get materials*/
					fVec3 box_p1{(x - 1) * dx / 2, (y - 1) * dx / 2, (z - 1) * dx / 2};
					real delta = dx / (2 * N);
					auto weights = get_zero_weights();
					
					for(int i = 1; i < 2 * N; i += 2)
						for(int j = 1; j < 2 * N; j += 2)
							for(int k = 1; k < 2 * N; k += 2) {
								auto sampled_point = box_p1 + fVec3{i * delta, j * delta, k * delta};
								bool assigned = 0;
								
								for(auto item : solids) {
									if (item->get_weights(sampled_point, weights)) {
										assigned = 1;
										break;
									}
								}
								
								if (!assigned)			//assign background medium;
									weights[background_medium->id] += 1;
							}
					
					chunk->set_medium_point(p, get_medium_internal(weights));
				}
	}
	
	void Simulation::source_init() {
		for(auto item : current_sources) {
			item->init(dx, chunk->get_p1(), chunk->get_p2(), chunk->get_dim(), chunk->get_origin());
			chunk->add_source_internal(item->get_source_internal());
		}
		
		for(auto item : eigen_sources) {
			item->init({0, 0, 0}, sim_dim * 2, chunk->get_dim(), chunk->get_origin(), chunk->get_p1(), chunk->get_p2());
			chunk->add_source_internal(item->get_source_internal());
		}
	}
	
	void Simulation::PML_init() {
		PML_init_helper(PMLs[0][0], PMLs[0][1], kx, bx, cx, sim_p1.x, sim_p2.x);
		PML_init_helper(PMLs[1][0], PMLs[1][1], ky, by, cy, sim_p1.y, sim_p2.y);
		PML_init_helper(PMLs[2][0], PMLs[2][1], kz, bz, cz, sim_p1.z, sim_p2.z);
		
		chunk->PML_init(kx, ky, kz, bx, by, bz, cx, cy, cz);
	}
	
	void PML_init_helper(const PML& neg, const PML& pos, real_arr& k, real_arr& b, real_arr& c, const int p1, const int p2) {
		k.resize(p2 - p1 + 1);
		b.resize(k.size());
		c.resize(k.size());
		
		for(int i = 0; i <= neg.get_d(); ++i) {
			real d = neg.get_d() - i / 2.0;
			
			k[i] = neg.get_chi(d);
			b[i] = neg.get_b(d);
			c[i] = neg.get_c(d);
		}
		
		for(int i = 0; i <= pos.get_d(); ++i) {
			real d = pos.get_d() - i / 2.0;
			
			k[p2 - i] = pos.get_chi(d);
			b[p2 - i] = pos.get_b(d);
			c[p2 - i] = pos.get_c(d);
		}
	}
	
	void Simulation::chunk_init() {
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
		chunk = new Chunk{sim_p1, sim_p2, ch_p1, ch_p2, dx, dt};
		ch_p1 = chunk->get_p1();
		ch_p2 = chunk->get_p2();
	}
	
	void Simulation::init() {
		chunk_init();
		PML_init();
		source_init();
		medium_init();
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
	
	void Simulation::set_background_medium(Medium_Type * medium) {
		background_medium = medium;
	}

	void Simulation::add_probe(Probe* probe) {
		probes.push_back(probe);
	}

	void Simulation::add_farfield_probe(const real freq, const fVec3& p) {
		N2F_freq.push_back(freq);
		N2F_pos.push_back(p);
	}
	
	void Simulation::advance() {
		real time = (step ++ ) * dt;
		chunk->update_Md(time);
		chunk->update_H2B(time);
		
		chunk->update_Jd(time);
		chunk->update_D2E(time);
		
		chunk->update_padded_E(time);
		chunk->update_padded_H(time);

		for (auto item : probes)
			item->update(*this);
	}

	real Simulation::at(const fVec3& p, const Coord_Type ctype) {
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
}
