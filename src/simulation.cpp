#include <simulation.hpp>
#include <iostream>

namespace ffip {
	void Simulation::setup(const real _dx, const real _dt, const iVec3 _dim) {
		dx = _dx;
		dt = _dt;
		sim_dim = _dim;
	}
	
	void Simulation::add_sf_layer(const int d) {
		sf_depth = d;
	}

	void Simulation::add_tf_layer(const int d) {
		tf_padded_depth = d;
	}
	
	void Simulation::add_PML_layer(const PML& pml, const Direction dir, const Side side) {
		if (side > 0)
			PMLs[dir][1] = pml;
		else
			PMLs[dir][0] = pml;
	}
	
	void Simulation::add_inc_source(Inc_Source* source) {
		Inc_Sources.push_back(source);
	}
	
	void Simulation::chunk_init() {
		//dimension is in computational units, so they are two times the original values
		sim_dim = sim_dim * 2;

		sim_p1 = { 0, 0, 0 };
		sim_p2 = sim_dim;

		//total field padded layer
		sim_p1 = sim_p1 - iVec3{ tf_padded_depth, tf_padded_depth, tf_padded_depth };
		sim_p2 = sim_p2 + iVec3{ tf_padded_depth, tf_padded_depth, tf_padded_depth };

		//scattered field layer
		tf_p1 = sim_p1;
		tf_p2 = sim_p2;
		sim_p1 = sim_p1 - iVec3{ sf_depth, sf_depth, sf_depth };
		sim_p2 = sim_p2 + iVec3{ sf_depth, sf_depth, sf_depth };

		phys_p1 = sim_p1;
		phys_p2 = sim_p2;
		
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
		chunk->set_num_proc(num_proc);
	}
	
	void Simulation::add_dipole(const real amp, const fVec3 pos, const Coord_Type ctype, const std::function<real (const real)> &profile) {
		chunk->add_dipoles(new Dipole(pos * (2 / dx), amp / (dx * dx * dx), profile, ctype, chunk));
	}
	
	void Simulation::set_background_medium(Medium const*  m) {
		if (m)
			bg_medium = m;
	}

	void Simulation::medium_init() {
		if (bg_medium == nullptr)
			bg_medium = make_medium(1, 0, 1, 0);
		
		prepare_medium(dt);
		const int N = 4;
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
				
				if (!is_M_point(ctype) && !is_E_point(ctype))			//exclude non-material points
					continue;
				
				/* staire case implementation
				 */
				//fVec3 sampled_point = p * dx / 2;
				//auto weights = get_zero_weights();
				//bool assigned = 0;
				//for(auto item : solids) {
				//	if (item->update_weights(sampled_point, weights)) {	//it is true when it is inside the solid
				//		assigned = 1;
				//		break;
				//	}
				//}

				//if (!assigned)			//assign background medium;
				//	weights[bg_medium->index] += 1;
				
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
				
				chunk->set_medium_point(p, get_medium_ref(ctype, weights));
			}
		};
		
		std::vector<std::future<void>> task_list;
		for(int i = 1; i < num_proc; ++i)
			task_list.push_back(std::async(std::launch::async, func, my_iterator(ch_p1, ch_p2, Null, i, num_proc)));
		func(my_iterator(ch_p1, ch_p2, Null, 0, num_proc));
		for(auto& item : task_list)
			item.get();
	}
	
	void Simulation::source_init() {
		for(auto item : Inc_Sources) {
			item->init(tf_p1, tf_p2, chunk->get_dim(), chunk->get_origin(), ch_p1, ch_p2, dx);
			chunk->add_inc_internal(item->get_source_internal());
		}
		
		chunk->organize_current_updates();
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
	
	void Simulation::set_num_proc(const int _num_proc) {
		num_proc = _num_proc;
		delete barrier;
		barrier = new Barrier(num_proc);
	}
	
	void Simulation::init() {
		PML_init();
		source_init();
		medium_init();
		udf_unit();
		chunk->categorize_points();
		step = 0;
		std::cout << "Initialization Complete\n";
	}
	
	void Simulation::add_solid(Solid const* solid) {
		solids.push_back(solid);
	}

	Nearfield_Probe const* Simulation::add_nearfield_probe(const real freq, fVec3 pos) {
		pos = pos * (2 / dx);
		if (!Is_Inside_Box(ch_p1, ch_p2, pos))
			throw Out_of_the_Domain{};
		
		auto res = new Nearfield_Probe{freq, pos, chunk};
		nearfield_probes.push_back(std::unique_ptr<Nearfield_Probe>{res});
		return res;
	}
	
	Flux_Box const* Simulation::add_flux_box(fVec3 p1, fVec3 p2, const real_arr& freq) {
		p1 = p1 * (2 / dx);
		p2 = p2 * (2 / dx);
		if (!Is_Inside_Box(ch_p1, ch_p2, p1) || !Is_Inside_Box(ch_p1, ch_p2, p2))
			throw Out_of_the_Domain{};
		
		auto res = new Flux_Box(p1, p2, freq, chunk);
		flux_boxes.push_back(std::unique_ptr<Flux_Box>{res});
		return res;
	}
	
	N2F_Box const* Simulation::add_n2f_box(fVec3 p1, fVec3 p2, const real_arr& freq) {
		p1 = p1 * (2 / dx);
		p2 = p2 * (2 / dx);
		if (!Is_Inside_Box(ch_p1, ch_p2, p1) || !Is_Inside_Box(ch_p1, ch_p2, p2))
			throw Out_of_the_Domain{};
		
		auto res = new N2F_Box(p1, p2, freq, chunk, bg_medium->get_c(), bg_medium->get_z());
		n2f_boxes.push_back(std::unique_ptr<N2F_Box>{res});
		return res;
	}
	
	PML Simulation::make_pml(const int d) {
		return PML(d, 0.8 * 4 / (dx * bg_medium->get_z()), 0.1, 1, 3, 1);
	}
	
	void Simulation::advance(std::ostream& os) {
		if (output_step_number)
			std::cout << "\n" << std::setfill('0') << std::setw(4) << step;
		
		real time = (step ++ ) * dt;
		
		auto func = [&, this](const int rank, const int num_proc) {
			chunk->update_Md(time, rank);
			barrier->Sync();
			chunk->update_B2H_v2(time, rank);
			barrier->Sync();

			chunk->update_Jd(time + 0.5 * dt, rank);
			barrier->Sync();
			chunk->update_D2E_v2(time + 0.5 * dt, rank);
			barrier->Sync();

			size_t idx1, idx2;
			vector_divider(nearfield_probes, rank, num_proc, idx1, idx2);
			for(auto itr = nearfield_probes.begin() + idx1; itr != nearfield_probes.begin() + idx2; ++itr)
				(*itr)->update(step);
			
			for(auto& item : n2f_boxes)
				item->update(step, rank, num_proc);
			
			for(auto& item : flux_boxes)
				item->update(step, rank, num_proc);
		};
		
		std::vector<std::future<void>> task_list;
		for(int i = 1; i < num_proc; ++i) {
			task_list.push_back(std::async(std::launch::async, func, i, num_proc));
		}
		func(0, num_proc);
		for(auto &item : task_list)
			item.get();
		
		//os << chunk->measure() << "\n";
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
	
	std::fstream os_tmp;
	iVec3 tmp_p1, tmp_p2;
	
	void Simulation::udf_unit() {
		/*os_tmp = std::fstream{"snapshot.out", std::ios::out};
		tmp_p1 = sim_p1;
		tmp_p2 = sim_p2;
		tmp_p1.x = (sim_p1.x + sim_p2.x) / 2;
		tmp_p2.x = (sim_p1.x + sim_p2.x) / 2 + 1;
		
		auto itr = my_iterator(tmp_p1, tmp_p2, Ex);
		auto dim = itr.get_dim();
		os_tmp << dim << "\n";*/

	/*	size_t num_probes = nearfield_probes.size();
		if (num_probes < 100)
			os_tmp << num_probes;
		else
			os_tmp << 100;
		os_tmp << "\n";*/
	}
	
	void Simulation::udf_advance() {
		/*if (step % 20 != 0) return;
		os_tmp << std::scientific << std::setprecision(5);
		for(auto itr = my_iterator(tmp_p1, tmp_p2, Ex); !itr.is_end(); itr.advance()) {
			os_tmp << (*chunk)(itr.get_vec()) << "\n";
		}*/

		/*size_t num_probes = nearfield_probes.size();
		if (num_probes < 100)
			for (auto& item : nearfield_probes)
				os_tmp << item->ex << " " << item->ey << " " << item->ez << " ";
		else
			for (int i = 0; i < 100; ++i) {
				auto& item = nearfield_probes[i * num_probes / 100];
				os_tmp << item->ex << " " << item->ey << " " << item->ez << " ";
			}

		os_tmp << "\n";*/
	}

	void Simulation::udf_output() {
	}
	
	
	void Simulation::prepare_medium(const real _dt) {
		for(auto& item : medium) {
			item->set_dt(_dt);
			e_medium_ref.push_back(std::make_unique<Medium_Ref>(item->get_e_medium_ref()));
			m_medium_ref.push_back(std::make_unique<Medium_Ref>(item->get_m_medium_ref()));
		}
	}
	
	std::vector<real> Simulation::get_zero_weights() {
		return std::vector<real>(medium.size(), 0);
	}
	
	Medium_Ref const* Simulation::get_medium_ref(const Coord_Type ctype, const std::vector<real>& weights) {
		
		bool is_electric_point = is_E_point(ctype);
		constexpr real tol = 1e-4;
		int nonzeros = 0;
		real total = 0;
		Medium_Ref* res = nullptr;
		auto& medium_ref = is_electric_point? e_medium_ref : m_medium_ref;
		
		if(weights.size() != medium.size())
			throw std::runtime_error("Mixing numbers have wrong length");
		
		/* rounding down parameters*/
		for(auto& x : weights) {
			if(x > tol) {
				nonzeros ++;
				total += x;
			}
		}
		
		if(nonzeros > 1) {								//for the case of mixing more than 1 material
			res = new Medium_Ref();
			for(int i = 0; i < weights.size(); ++i)
				if (weights[i] > tol) *res += *medium_ref[i] * (weights[i] / total);
			
			std::lock_guard<std::mutex> guard(medium_mutex);
			syn_medium_ref.push_back(std::unique_ptr<Medium_Ref>(res));			//garbage collecting
		}
		else if (nonzeros == 1){						//for the case of only 1 material
			for(int i = 0; i < weights.size(); ++i) {
				if (weights[i] > tol) {
					res = medium_ref[i].get();
					break;
				}
			}
		}
		
		if (res == nullptr)
			throw std::runtime_error("Illegal weights");
		return res;
	}
	
	Solid const* Simulation::make_solid(Medium const* m1, Medium const*  m2, const std::string& filename) {
		inhom_box_holder.push_back(std::make_unique<Inhomogeneous_Box>(m1, m2, filename));
		return inhom_box_holder.back().get();
	}
	
	Solid const* Simulation::make_solid(Medium const* m, const Geometry_Node& geom) {
		hom_box_holder.push_back(std::make_unique<Homogeneous_Object>(m, geom));
		return hom_box_holder.back().get();
	}
}

