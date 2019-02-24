#include <simulation.hpp>
#include <iostream>

//#include <cvmarkersobj.h>

//using namespace Concurrency;

namespace ffip {
	void Simulation::setup(const real _dx, const real _dt, const iVec3 _dim) {
		dx = _dx;
		dt = _dt;
		sim_dim = _dim;
	}

	void Simulation::add_sf_layer(const int d, const Direction dir, const Side side) {
		if (side > 0)
			bds[dir][1].sf_layer = d;
		else
			bds[dir][0].sf_layer = d;
	}

	void Simulation::add_tf_layer(const int d, const Direction dir, const Side side) {
		if (side > 0)
			bds[dir][1].tf_layer = d;
		else
			bds[dir][0].tf_layer = d;
	}
	
	void Simulation::add_sf_layer(const int d) {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 2; ++j)
				bds[i][j].sf_layer = d;
	}

	void Simulation::add_tf_layer(const int d) {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 2; ++j)
				bds[i][j].tf_layer = d;
	}
	
	void Simulation::add_PML_layer(const PML& pml, const Direction dir, const Side side) {
		if (side > 0)
			bds[dir][1].pml = pml;
		else
			bds[dir][0].pml = pml;
	}
	
	void Simulation::set_symmetry(const bool symx, const bool symy, const bool symz) {
		enable_symmetry_x = symx;
		enable_symmetry_y = symy;
		enable_symmetry_z = symz;
	}
	
	void Simulation::add_plane_wave(const Func type, const real fp, const real delay, const real amp, const real pos, const Coord_Type polarization) {
		pw_configs.push_back({ type, fp, delay, amp, pos, polarization });

		// enable symmetry in x, y direction
		enable_symmetry_x = 1;
		enable_symmetry_y = 1;
	}
	
	void Simulation::add_inc_source(Inc_Source* source) {
		if (source)
			excitation = source;
	}
	
	void Simulation::chunk_init() {
		//symmetry configurations, disable other boundary conditions
		enable_symmetry = enable_symmetry_x | enable_symmetry_y | enable_symmetry_z;
		if (enable_symmetry_x)
			bds[0][0] = bds[0][1] = Boundary{};
		
		if (enable_symmetry_y)
			bds[1][0] = bds[1][1] = Boundary{};
		
		if (enable_symmetry_z)
			bds[2][0] = bds[2][1] = Boundary{};
		
		//dimension is in computational units, so they are two times the original values
		roi_p1 = { 0, 0, 0 };
		roi_p2 = sim_dim * 2;

		//total field padded layer
		tf_p1 = roi_p1 - iVec3{ bds[0][0].tf_layer, bds[1][0].tf_layer, bds[2][0].tf_layer };
		tf_p2 = roi_p2 + iVec3{ bds[0][1].tf_layer, bds[1][1].tf_layer, bds[2][1].tf_layer };

		//scattered field layer
		phys_p1 = tf_p1 - iVec3{ bds[0][0].sf_layer, bds[1][0].sf_layer, bds[2][0].sf_layer };
		phys_p2 = tf_p2 + iVec3{ bds[0][1].sf_layer, bds[1][1].sf_layer, bds[2][1].sf_layer };
		
		//add PML layers
		sim_p1 = phys_p1 - 2 * iVec3{ bds[0][0].pml.get_d() , bds[1][0].pml.get_d(), bds[2][0].pml.get_d() };
		sim_p2 = phys_p2 + 2 * iVec3{ bds[0][1].pml.get_d() , bds[1][1].pml.get_d(), bds[2][1].pml.get_d() };
		
		//chunk region, reserved for MPI implementation
		ch_p1 = sim_p1;
		ch_p2 = sim_p2;

		config = Config{dt, dx, roi_p1, roi_p2, sim_p1, sim_p2, ch_p1, ch_p2, tf_p1, tf_p2, phys_p1, phys_p2};
		//implementaions of MPI, for now 1 chunk covers the whole region
		chunk = new Chunk{config};
	}
	
	void Simulation::add_dipole(real amp, fVec3 pos, const Coord_Type ctype, const Func func_type, const real fp, const real delay) {
		dipole_configs.push_back(Dipole_Config{ amp, pos, ctype, func_type, fp, delay });
	}
	
	void Simulation::set_background_medium(Medium const*  m) {
		if (m)
			bg_medium = m;
	}

	void Simulation::medium_init() {
		if (bg_medium == nullptr)
			bg_medium = make_medium(1, 0, 1, 0);
		
		prepare_medium(dt);
		
		//spatial average method
		std::vector<fVec3> sampled_points;
		{
			const int N = 4;
			real delta = dx / (2 * N);
			for (int i = 1; i < 2 * N; i += 2)
				for (int j = 1; j < 2 * N; j += 2)
					for (int k = 1; k < 2 * N; k += 2) {
						sampled_points.push_back({ i * delta, j * delta, k * delta });
					}
		}
		

		//R. Mittra Method (line average)
		std::vector<fVec3> sampled_points_rm[8];
		{
			const int N = 20;
			real delta = dx / (N - 1);
			for (int i = 0; i < N; ++i) {
				real x = -dx / 2 + delta * i;
				sampled_points_rm[Ex].push_back({ x, 0, 0 });
				sampled_points_rm[Ey].push_back({ 0, x, 0 });
				sampled_points_rm[Ez].push_back({ 0, 0, x });
				sampled_points_rm[Hx].push_back({ x, 0, 0 });
				sampled_points_rm[Hy].push_back({ 0, x, 0 });
				sampled_points_rm[Hz].push_back({ 0, 0, x });
			}
		}

		std::vector<std::future<void>> task_list;
		std::vector<Medium_Voxel> medium_voxels(my_iterator(ch_p1, ch_p2, All).get_size());
		auto strides = dim2stride(ch_p2 - ch_p1 + 1);

		auto get_index = [=, ch_p1 = this->ch_p1](const iVec3& p) {
			return inner_prod(p - ch_p1, strides);
		};

		//average to get property for each point
		{
			auto solid_assign = [&, this](const size_t rank) {
				for (auto itr = my_iterator(ch_p1, ch_p2, All, rank, num_proc); !itr.is_end(); itr.advance()) {
					/* assign material only when it is E or H*/
					auto p = itr.get_vec();
					auto ctype = p.get_type();

					if (!is_M_Point(ctype) && !is_E_Point(ctype))			//exclude non-material points
						continue;

					auto weights = get_zero_weights();
					/* spatial averaging */
					//{
					//	fVec3 box_p1{ (p.x - 1) * dx / 2, (p.y - 1) * dx / 2, (p.z - 1) * dx / 2 };
					//	for (auto& item : sampled_points) {
					//		auto sampled_point = box_p1 + item;
					//		bool assigned = 0;

					//		for (auto item : solids) {
					//			if (item->update_weights(sampled_point, weights)) {	//it is true when it is inside the solid
					//				assigned = 1;
					//				break;
					//			}
					//		}

					//		if (!assigned)			//assign background medium;
					//			weights[bg_medium->index] += 1;
					//	}
					//}

					/* R. Mittra averaging (line averaging) */
					{
						fVec3 fp = p * (dx / 2);
						for (auto& item : sampled_points_rm[ctype]) {
							auto sampled_point = fp + item;
							bool assigned = 0;

							for (auto item : solids) {
								if (item->update_weights(sampled_point, weights)) {	//it is true when it is inside the solid
									assigned = 1;
									break;
								}
							}

							if (!assigned)			//assign background medium;
								weights[bg_medium->index] += 1;
						}
					}
					

					medium_voxels[itr.index] = weights / (weights.sum());
				}
			};

			for (int i = 1; i < num_proc; ++i)
				task_list.push_back(std::async(std::launch::async, solid_assign, i));
			solid_assign(0);
			for (auto& item : task_list)
				item.get();
			task_list.clear();
		}
		

		// blobs material implementation 1 (source like interpolation) not recommended
		/*{
			std::vector<real> weight(8);
			interpn<3> interp(2, 2, 2);

			for (auto& item : medium_blobs) {
				for (int ctype_int = 1; ctype_int < 7; ctype_int++) {
					Coord_Type ctype = static_cast<Coord_Type>(ctype_int);
					iVec3 base = get_nearest_point<side_low_tag>(item.pos, ctype);
					fVec3 dp = (item.pos - base) / 2;

					std::fill(weight.begin(), weight.end(), 0);
					interp.put(weight, item.amp, dp.z, dp.y, dp.x);
					for (auto itr = my_iterator(base, base + iVec3{ 2, 2, 2 }, ctype); !itr.is_end(); itr.advance()) {
						medium_voxels[get_index(itr.get_vec())][item.medium->index] += weight[itr.index];
					}
				}
			}
		}*/
		
		// blobs material implementation 2 (small box implementation)
		{
			for (auto& item : medium_blobs) {
				if (!Is_Inside_Box(config.roi_p1, config.roi_p2, item.pos))
					throw Out_of_the_Domain{};

				iVec3 low(floor(item.pos - 1));
				iVec3 high(ceil(item.pos + 1));

				for (auto itr = my_iterator(low, high, All); !itr.is_end(); itr.advance()) {
					auto p = itr.get_vec();
					fVec3 dist = abs(item.pos - p);
					real weight = prod(2 - dist) / 8.0;

					size_t cur_index = get_index(p);
					if (medium_voxels[cur_index].size())
						medium_voxels[cur_index][item.medium->index] += weight * item.amp;
				}
			}
		}
		
		// inhomogeneous region (force material distribution)
		{
			for (auto inhom : inhoms) {
				auto p1 = inhom->get_p1() * (2 / dx);
				auto p2 = inhom->get_p2() * (2 / dx);
				auto region = get_component_interior(p1, p2, All);
				auto inter = get_intersection(region.first, region.second, ch_p1, ch_p2);

				for (auto itr = my_iterator(inter.first, inter.second, All); !itr.is_end(); itr.advance()) {
					auto p = itr.get_vec();
					auto ctype = p.get_type();

					if (!is_M_Point(ctype) && !is_E_Point(ctype))			//exclude non-material points
						continue;

					auto weights = get_zero_weights();

					if (inhom->update_weights(p * (dx / 2), weights))
						medium_voxels[inner_prod(p - ch_p1, strides)] = weights;
				}
			}
		}

		//{
		//	auto fo = std::fstream{ "medium.out", std::ios::out };
		//	//auto voxel = medium_voxels[inner_prod(iVec3{ 41, 40, 41 } - ch_p1, strides)];
		//	for (auto itr = my_iterator(ch_p1, ch_p2, All); !itr.is_end(); itr.advance()) {
		//		if (medium_voxels[itr.index].size() && medium_voxels[itr.index][1] > 0.2) {
		//			fo << itr.get_vec() << " " << medium_voxels[itr.index][1] << "\n";
		//		}
		//	}
		//}

		//push medium to chunk
		{
			auto push_func = [&, this](const size_t rank) {
				for (auto itr = my_iterator(ch_p1, ch_p2, All, rank, num_proc); !itr.is_end(); itr.advance()) {
					/* assign material only when it is E or H*/
					auto p = itr.get_vec();
					auto ctype = p.get_type();

					if (!is_M_Point(ctype) && !is_E_Point(ctype))			//exclude non-material points
						continue;
					chunk->set_medium_point(p, get_medium_ref(ctype, medium_voxels[itr.index]));
				}
			};

			for (int i = 1; i < num_proc; ++i)
				task_list.push_back(std::async(std::launch::async, push_func, i));
			push_func(0);
			for (auto& item : task_list)
				item.get();
		}
	}
	
	void Simulation::source_init() {
		//projector initialization
		if (excitation) {
			excitation->init(config);
			excitation->push(chunk);
		}
		
		//dipole sources initialization
		std::vector<real> weight(8);
		interpn<3> interp{ 2, 2, 2 };
		real vol = dx * dx * dx;

		for (auto& dipole : dipole_configs) {
			//reversion of amplitude for electric dipoles
			if (is_E_Point(dipole.ctype))
				dipole.amp = -dipole.amp;

			dipole.pos = dipole.pos * (2 / dx);
			if (!Is_Inside_Box(config.roi_p1 - 1e-3, config.roi_p2 + 1e-3, dipole.pos))
				throw std::runtime_error("Dipole must be inside region of interests");

			iVec3 base = get_nearest_point<side_low_tag>(dipole.pos, dipole.ctype);
			fVec3 dp = (dipole.pos - base) / 2;

			interp.put(weight, 1.0, dp.z, dp.y, dp.x);
			for (auto x : weight)
				if (x < 0 || x > 1.0)
					throw std::runtime_error("Interpolation place failed");

			for (auto itr = my_iterator(base, base + iVec3{ 2, 2, 2 }, dipole.ctype); !itr.is_end(); itr.advance()) {
				chunk->add_dipole(dipole.type, itr.get_vec(), dipole.amp * weight[itr.index] / vol, dipole.delay, dipole.fp);
			}
		}

		//symmetric plane wave initialization
		for (auto pw : pw_configs) {
			auto p1 = config.phys_p1;
			auto p2 = config.phys_p2;

			p1.z = (int)round(pw.pos / dx) * 2;
			p2.z = p1.z;

			if (!In_ClosedBounds(p1.z, config.phys_p1.z, config.phys_p2.z))
				throw Out_of_the_Domain{};

			for (auto itr = my_iterator(p1, p2, pw.polarization); !itr.is_end(); itr.advance()) {
				chunk->add_dipole(pw.type, itr.get_vec(), pw.amp, pw.delay, pw.fp);
			}
		}
	}
	
	void Simulation::PML_init() {
		std::array<real_arr, 3> k, b, c;

		PML_init_helper(bds[0][0].pml, bds[0][1].pml, k[0], b[0], c[0], sim_p1.x, sim_p2.x);
		PML_init_helper(bds[1][0].pml, bds[1][1].pml, k[1], b[1], c[1], sim_p1.y, sim_p2.y);
		PML_init_helper(bds[2][0].pml, bds[2][1].pml, k[2], b[2], c[2], sim_p1.z, sim_p2.z);
		
		chunk->PML_init(k, b, c);
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
		chunk_init();
		PML_init();
		source_init();
		medium_init();
		chunk->init();
		udf_unit();
		step = 0;
		std::cout << "Initialization Complete\n";
	}
	
	void Simulation::add_blob(fVec3 pos, real amp, Medium const* medium) {
		pos = pos * (2 / dx);
		medium_blobs.push_back(Medium_Blob{ medium, pos, amp });
	}

	void Simulation::add_solid(Solid const* solid) {
		if (solid) {
			auto ptr = dynamic_cast<Inhomogeneous_Box const*>(solid);
			if (ptr)
				inhoms.push_back(ptr);
			else
				solids.push_back(solid);
		}			
	}

	Nearfield_Probe const* Simulation::add_nearfield_probe(const real freq, fVec3 pos) {
		if (chunk == nullptr)
			throw std::runtime_error("Initialize Simulation First");

		pos = pos * (2 / dx);
		if (!Is_Inside_Box(config.sim_p1, config.sim_p2, pos))
			throw Out_of_the_Domain{};
		
		auto res = new Nearfield_Probe{freq, pos, chunk};
		nearfield_probes.push_back(std::unique_ptr<Nearfield_Probe>{res});
		return res;
	}
	
	Flux_Box const* Simulation::add_flux_box(fVec3 p1, fVec3 p2, const real_arr& freq) {
		if (chunk == nullptr)
			throw std::runtime_error("Initialize Simulation First");

		p1 = p1 * (2 / dx);
		p2 = p2 * (2 / dx);
		if (!Is_Inside_Box(config.phys_p1, config.phys_p2, p1) || !Is_Inside_Box(config.phys_p1, config.phys_p2, p2))
			throw Out_of_the_Domain{};
		
		auto res = new Flux_Box(p1, p2, freq, chunk);
		flux_boxes.push_back(std::unique_ptr<Flux_Box>{res});
		return res;
	}
	
	N2F_Box const* Simulation::add_n2f_box(fVec3 p1, fVec3 p2, const real_arr& freq) {
		if (chunk == nullptr)
			throw std::runtime_error("Initialize Simulation First");

		p1 = p1 * (2 / dx);
		p2 = p2 * (2 / dx);
		if (!Is_Inside_Box(config.phys_p1, config.phys_p2, p1) || !Is_Inside_Box(config.phys_p1, config.phys_p2, p2))
			throw Out_of_the_Domain{};
		
		auto res = new N2F_Box(p1, p2, freq, chunk, bg_medium->get_c(), bg_medium->get_z());
		n2f_boxes.push_back(std::unique_ptr<N2F_Box>{res});
		return res;
	}
	
	PML Simulation::make_pml(const int d) {
		real sigma_max = 0.8 * 4 / (dx * bg_medium->get_z());

		return PML(d, sigma_max, 1, 0, 3, 1);
	}
	
	void Simulation::advance(std::ostream& os) {
		if (output_step_number)
			std::cout << "\n" << std::setfill('0') << std::setw(4) << step;
		
		real time = (step ++ ) * dt;

//		diagnostic::marker_series pp("Simulation");
		
		auto func = [&, this](const int rank, const int num_proc) {
			{
//				diagnostic::span span(pp, "Md");
				chunk->update_Md(time, rank, barrier);
				barrier->Sync();
			}
			
			{
//				diagnostic::span span(pp, "B2H");
				chunk->update_B2H_v2(time, rank, barrier);
				barrier->Sync();
			}
			
			if (enable_symmetry) {
				if (enable_symmetry_x)
					chunk->update_ghost_helper<dir_x_tag>(rank, barrier);
				
				if (enable_symmetry_y)
					chunk->update_ghost_helper<dir_y_tag>(rank, barrier);
				
				if (enable_symmetry_z)
					chunk->update_ghost_helper<dir_z_tag>(rank, barrier);
				
				barrier->Sync();
			}

			if (excitation) {
				if (rank == 0)
					excitation->advance();
				barrier->Sync();
			}

			{
//				diagnostic::span span(pp, "Jd");
				chunk->update_Jd(time + 0.5 * dt, rank, barrier);
				barrier->Sync();
			}
			
			{
//				diagnostic::span span(pp, "D2E");
				chunk->update_D2E_v2(time + 0.5 * dt, rank, barrier);
				barrier->Sync();
			}
			
			if (enable_symmetry) {
				if (enable_symmetry_x)
					chunk->update_ghost_helper<dir_x_tag>(rank, barrier);
				
				if (enable_symmetry_y)
					chunk->update_ghost_helper<dir_y_tag>(rank, barrier);
				
				if (enable_symmetry_z)
					chunk->update_ghost_helper<dir_z_tag>(rank, barrier);
				
				barrier->Sync();
			}

			// fields calculation, does not require synchronization
			{
//				diagnostic::span span(pp, "probes");
				size_t idx1, idx2;
				vector_divider(nearfield_probes, rank, num_proc, idx1, idx2);
				for (auto itr = nearfield_probes.begin() + idx1; itr != nearfield_probes.begin() + idx2; ++itr)
					(*itr)->update(step);
			}
			
			{
//				diagnostic::span span(pp, "n2f");
				for (auto& item : n2f_boxes)
					item->update(step, rank, num_proc);
			}
			
			{
//				diagnostic::span span(pp, "flux");
				for (auto& item : flux_boxes)
					item->update(step, rank, num_proc);
			}
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
		return sim_dim;
	}
	
	Medium const* Simulation::get_bg_medium() const {
		return bg_medium; 
	}
	
	std::fstream os_tmp;
	iVec3 tmp_p1, tmp_p2;
	
	void Simulation::udf_unit() {
		os_tmp = std::fstream{"snapshot.out", std::ios::out};
		tmp_p1 = sim_p1;
		tmp_p2 = sim_p2;
		//tmp_p1.x = (sim_p1.x + sim_p2.x) / 2;
		//tmp_p2.x = (sim_p1.x + sim_p2.x) / 2 + 1;
		tmp_p1.x = 0;
		tmp_p2.x = 1;

		auto itr = my_iterator(tmp_p1, tmp_p2, Hy);
		auto dim = itr.get_dim();
		os_tmp << dim << "\n";

	/*	size_t num_probes = nearfield_probes.size();
		if (num_probes < 100)
			os_tmp << num_probes;
		else
			os_tmp << 100;
		os_tmp << "\n";*/
	}
	
	void Simulation::udf_advance() {
		if (step % 1 != 0) return;
		os_tmp << std::scientific << std::setprecision(5);
		for(auto itr = my_iterator(tmp_p1, tmp_p2, Hy); !itr.is_end(); itr.advance()) {
			os_tmp << (*chunk)(itr.get_vec()) << "\n";
		}

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
	
	Medium_Voxel Simulation::get_zero_weights() {
		return Medium_Voxel(0.0, medium.size());
	}
	
	Medium_Ref const* Simulation::get_medium_ref(const Coord_Type ctype, const Medium_Voxel& weights) {
		
		bool is_electric_point = is_E_Point(ctype);
		constexpr real tol = 1e-3;
		int nonzeros = 0;
		Medium_Ref* res = nullptr;
		auto& medium_ref = is_electric_point? e_medium_ref : m_medium_ref;
		
		if(weights.size() != medium.size())
			throw std::runtime_error("Mixing numbers have wrong length");
		
		/* one material : valid for most points */
		double val;
		int id;
		for (int i = 0; i < weights.size(); ++i) {
			if (std::abs(weights[i]) > tol) {
				nonzeros++;
				val = weights[i];
				id = i;
			}
		}

		if (nonzeros == 1 && std::abs(val - 1) < tol) {
			res = medium_ref[id].get();
		}
		else {
			res = new Medium_Ref();
			for (int i = 0; i < weights.size(); ++i)
				if (std::abs(weights[i]) > tol) *res += *medium_ref[i] * weights[i];

			std::lock_guard<std::mutex> guard(medium_mutex);
			syn_medium_ref.push_back(std::unique_ptr<Medium_Ref>(res));			//garbage collecting
		}

		if (res == nullptr)
			throw std::runtime_error("Illegal weights");

		/*if (!res->is_valid())
			throw std::runtime_error("Numerical Instability");*/

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

