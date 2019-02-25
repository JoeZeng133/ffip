#include <chunk.hpp>
#include <source.hpp>

//#include <cvmarkersobj.h>

//using namespace Concurrency;

namespace ffip {
	CU_PML::CU_PML(std::vector<real>& jmd, std::vector<real>& eh, std::array<size_t, 3> strides) : jmd(jmd), eh(eh), strides(strides) {}

	void CU_PML::update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) {}

	void CU_PML::organize() {}

	void CU_PML::update_static(const size_t rank, Barrier* barrier) {

		for (int i = 0; i < 3; ++i) {
			size_t stride1 = strides[(i + 1) % 3];
			size_t stride2 = strides[(i + 2) % 3];
			auto& rpoints = points[i];
			size_t idx1, idx2;

			vector_divider(rpoints, rank, barrier->get_num_proc(), idx1, idx2);
			for (; idx1 < idx2; ++idx1) {
				size_t index = indexes[i][idx1];
				auto& point = rpoints[idx1];

				real curl1 = (eh[index + stride1] - eh[index - stride1]);
				real curl2 = (eh[index + stride2] - eh[index - stride2]);

				point.phi1 = point.b1 * point.phi1 + point.c1 * curl1;
				point.phi2 = point.b2 * point.phi2 + point.c2 * curl2;

				jmd[index] = point.phi1 - point.phi2 + curl1 * point.k1 - curl2 * point.k2;
			}
		}
	}

	CU_Dipole::CU_Dipole(std::vector<real>& jmd, const func& f, real time, real dt) : jmd(jmd), f(f), time(time), dt(dt) {}

	void CU_Dipole::organize() {
		std::sort(points.begin(), points.end(), [](const Update_Point& a, const Update_Point& b) {return a.index < b.index; });
	}

	void CU_Dipole::update_static(const size_t rank, Barrier* barrier) {
		size_t idx1, idx2;
		vector_divider(points, rank, barrier->get_num_proc(), idx1, idx2);
		while (idx1 > 0 && idx1 < points.size() && points[idx1].index == points[idx1 - 1].index) --idx1;
		while (idx2 > 0 && idx2 < points.size() && points[idx2].index == points[idx2 - 1].index) --idx2;
		real cur_time = time + step * dt;

		for (int i = idx1; i < idx2; ++i) {
			auto& point = points[i];
			jmd[point.index] += point.c * f(cur_time, point.d, point.fp);
		}

		barrier->Sync();
		if (rank == 0)
			step++;
	}

	void CU_Dipole::add_update_point(const size_t index, const real amp, const real delay, const real fp) {
		points.push_back({ index, amp, delay, fp });
	}

	void CU_Dipole::update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) {}

	void CU_Dipole::output(std::ostream& os) {
		for (auto p : points)
			os << p.c << " " << p.d << " " << p.fp << " " << p.index << "\n";
	}

	Chunk::Chunk(const Config& config) :Chunk(config.sim_p1, config.sim_p2, config.ch_p1, config.ch_p2, config.dx, config.dt) {
		this->config = config;
	}


	Chunk::Chunk(const iVec3& sim_p1, const iVec3& sim_p2, const iVec3& ch_p1, const iVec3& ch_p2, real dx, real dt):  ch_p1(ch_p1), ch_p2(ch_p2), sim_p1(sim_p1), sim_p2(sim_p2), dx(dx), dt(dt) {
		
		ch_origin = ch_p1 - 1;
		ch_dim = ch_p2 - ch_p1 + 3;

		ch_stride_x = 1;
		ch_stride_y = ch_dim.x;
		ch_stride_z = (size_t)ch_dim.x * ch_dim.y;
		
		eh.resize((size_t)ch_dim.x * ch_dim.y * ch_dim.z, 0);
		eh1.resize(eh.size(), 0);
		jmd.resize(eh.size(), 0);
		medium_chunk.resize(eh.size(), nullptr);
		dispersive_field_chunk.resize(eh.size(), nullptr);
		
		for (int i = 0; i < 8; ++i)
			strides[i] = (i & 1? ch_stride_x * 2 : 0) + (i & 2? ch_stride_y * 2: 0) + (i & 4? ch_stride_z * 2 : 0);
		
		for(int i = 0; i < 8; ++i)
			num_ones[i ^ 0b111] = 1 << (!(i & 1) + !(i & 2) + !(i & 4));
	}
	
	real Chunk::get_dt() const {
		return dt;
	}
	
	real Chunk::get_dx() const {
		return dx;
	}
	
	iVec3 Chunk::get_dim() const{
		return ch_dim;
	}
	
	iVec3 Chunk::get_origin() const {
		return ch_origin;
	}
	
	iVec3 Chunk::get_p1() const {
		return ch_p1;
	}
	
	iVec3 Chunk::get_p2() const {
		return ch_p2;
	}
	
	size_t Chunk::get_index_ch(const iVec3& p) const {
		return (p.x - ch_origin.x) * ch_stride_x + (p.y - ch_origin.y) * ch_stride_y + (p.z - ch_origin.z) * ch_stride_z;
	}

	iVec3 Chunk::get_pos(size_t index) const {
		int x, y, z;
		
		x = ch_origin.x + (index % ch_dim.x);
		index /= ch_dim.x;
		y = ch_origin.y + (index % ch_dim.y);
		z = ch_origin.z + (index / ch_dim.y);

		return { x, y, z };
	}
	
	void Chunk::set_medium_point(const iVec3 &point, Medium_Ref const* medium_ref) {
		int index = get_index_ch(point);
		medium_chunk[index] = medium_ref;
		
		size_t num_poles = medium_ref->get_size_poles();
		if (num_poles > 0) {
			dispersive_field_chunk[index] = new Dispersive_Field{num_poles};
		}
	}
	
	void Chunk::PML_init(const std::array<real_arr, 3>& k, const std::array<real_arr, 3>& b, const std::array<real_arr, 3>& c) {
		auto p1_ch = ch_p1 - ch_origin;
		auto p2_ch = ch_p2 - ch_origin;
		auto p1_sim = ch_p1 - sim_p1;
		auto p2_sim = ch_p2 - sim_p1;

		//chop up the segment needed and add ghost points
		this->kx.insert(this->kx.end(), k[0].begin() + p1_sim.x, k[0].begin() + p2_sim.x + 1);
		this->kx.insert(this->kx.begin(), 1);
		this->kx.insert(this->kx.end(), 1);

		this->ky.insert(this->ky.end(), k[1].begin() + p1_sim.y, k[1].begin() + p2_sim.y + 1);
		this->ky.insert(this->ky.begin(), 1);
		this->ky.insert(this->ky.end(), 1);

		this->kz.insert(this->kz.end(), k[2].begin() + p1_sim.z, k[2].begin() + p2_sim.z + 1);
		this->kz.insert(this->kz.begin(), 1);
		this->kz.insert(this->kz.end(), 1);
		
		e_PML = new CU_PML(jmd, eh, { ch_stride_x, ch_stride_y, ch_stride_z });
		m_PML = new CU_PML(jmd, eh, { ch_stride_x, ch_stride_y, ch_stride_z });

		PML_init_helper<ex_tag>(k, b, c);
		PML_init_helper<ey_tag>(k, b, c);
		PML_init_helper<ez_tag>(k, b, c);
		PML_init_helper<hx_tag>(k, b, c);
		PML_init_helper<hy_tag>(k, b, c);
		PML_init_helper<hz_tag>(k, b, c);
	}
	
//	diagnostic::marker_series marker_jmd("Current Updates");
	
	void Chunk::update_Jd(const real time, const size_t rank, Barrier* barrier) {
		if (rank >= barrier->get_num_proc())
			throw std::runtime_error("Rank cannot be bigger than number of processes");

//		diagnostic::span* s;
//		s = new diagnostic::span(marker_jmd, "jx");
		update_curl<ex_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "jy");
		update_curl<ey_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "jz");
		update_curl<ez_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "e_PML");
		if (e_PML) e_PML->update_static(rank, barrier);
//		delete s;

		barrier->Sync();

//		s = new diagnostic::span(marker_jmd, "e_cu");
		if (e_tfsf) e_tfsf->update_static(rank, barrier);
//		delete s;
		
//		s = new diagnostic::span(marker_jmd, "e_dipole");
		for (auto item : e_dipoles) {
			item.second->update_static(rank, barrier);
			barrier->Sync();
		}
//		delete s;
	}
	
	void Chunk::update_Md(const real time, const size_t rank, Barrier* barrier) {
		if (rank >= barrier->get_num_proc())
			throw std::runtime_error("Rank cannot be bigger than number of processes");
		//marker_series mySeries;
//		diagnostic::span* s;
		
//		s = new diagnostic::span(marker_jmd, "mx");
		update_curl<hx_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "my");
		update_curl<hy_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "mz");
		update_curl<hz_tag>(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "m_PML");
		if (m_PML) m_PML->update_static(rank, barrier);
//		delete s;

		barrier->Sync();

//		s = new diagnostic::span(marker_jmd, "m_cu");
		if (m_tfsf) m_tfsf->update_static(rank, barrier);
//		delete s;

//		s = new diagnostic::span(marker_jmd, "m_dipole");
		for (auto item : m_dipoles) {
			item.second->update_static(rank, barrier);
			barrier->Sync();
		}
			
//		delete s;
	}
	
	void Chunk::update_jwD2E(const real time, const size_t rank, Barrier* barrier) {
		real modifier = 1 / e0 * dt;
		size_t num_proc = barrier->get_num_proc();

		for (int i = 0; i < MAX_NUM_POLES; ++i) {
			size_t idx1, idx2;
			vector_divider(e_points[i], rank, num_proc, idx1, idx2);
			for(int j = idx1; j < idx2; ++j) {
				size_t index = e_points[i][j];
				real tmp = eh[index];
				
				eh[index] = medium_chunk[index]->update_field(eh[index], eh1[index], modifier * jmd[index], dispersive_field_chunk[index]);
				eh1[index] = tmp;
			}
		}
	}
	
	void Chunk::update_jwB2H(const real time, const size_t rank, Barrier* barrier) {
		real modifier = -1 / u0 * dt;
		size_t num_proc = barrier->get_num_proc();

		for (int i = 0; i < MAX_NUM_POLES; ++i) {
			size_t idx1, idx2;
			vector_divider(m_points[i], rank, num_proc, idx1, idx2);
			for(int j = idx1; j < idx2; ++j) {
				size_t index = m_points[i][j];
				real tmp = eh[index];
				
				eh[index] = medium_chunk[index]->update_field(eh[index], eh1[index], modifier * jmd[index], dispersive_field_chunk[index]);
				eh1[index] = tmp;
			}
		}
	}
	
	void Chunk::update_ud(const real time, const size_t rank, Barrier *barrier) {
	}
	
	void Chunk::update_ghost_E(const real time, const size_t rank, Barrier *barrier) {}
	
	void Chunk::update_ghost_H(const real time, const size_t rank, Barrier *barrier) {}
	
	real Chunk::at(const fVec3& p, const Coord_Type ctype) const{
		return operator()(p / (dx / 2), ctype);
	}

	real Chunk::operator()(const iVec3 &p) const {
		if (!Is_Inside_Box(ch_p1, ch_p2, p))
			throw Out_of_the_Domain{};

		return eh[get_index_ch(p)];
	}

	real Chunk::operator()(const iVec3& p, const Coord_Type ctype) const {
		if (!Is_Inside_Box(ch_p1, ch_p2, p))
			throw Out_of_the_Domain{};

		return ave(ctype ^ p.get_type(), get_index_ch(p));
	}

	real Chunk::operator()(const fVec3& p, const Coord_Type ctype) const {
		static interp3 interp{2, 2, 2};
		if (!Is_Inside_Box(ch_p1, ch_p2, p))
			throw Out_of_the_Domain{};
 
		real_arr f(8);
		
		iVec3 base_p = get_nearest_point<side_low_tag>(p, ctype);
		size_t base_index_ch = get_index_ch(base_p);

		for (int i = 0; i < 8; ++i)
			f[i] = eh[base_index_ch + strides[i]];
		
		return interp.get(f, (p.z - base_p.z) / 2, (p.y - base_p.y) / 2, (p.x - base_p.x) / 2);
	}
	
	real Chunk::operator[](const size_t index) const {
		return eh[index];
	}

	template<>
	real Chunk::ave_helper(const int bit, const int index, const int stride) const{
		if(bit & 1) {
			return eh[index - stride] + eh[index + stride];
		}else{
			return eh[index];
		}
	}
	
	real Chunk::ave(const int bit, const int index) const{
		return ave_helper(bit, index, ch_stride_z, ch_stride_y, ch_stride_x) / num_ones[bit];
	}

	void Chunk::init() {
		//material updates
		categorize_points();
		//PML updates
		if (e_PML) e_PML->organize();
		if (m_PML) m_PML->organize();
		//TFSF updates
		if (e_tfsf) e_tfsf->organize();
		if (m_tfsf) m_tfsf->organize();
		//Dipole updates
		for (auto& item : e_dipoles)
			item.second->organize();
		for (auto& item : m_dipoles)
			item.second->organize();

		udf_init();
	}

	std::fstream tmp;
	void Chunk::udf_init() {
		tmp = std::fstream("debug.out", std::ios::out);
		for (auto& item : e_dipoles)
			item.second->output(tmp);
	}

	void Chunk::udf_update() {

	}
	
	void Chunk::categorize_points() {
		for(auto itr = my_iterator(ch_p1, ch_p2, All); !itr.is_end(); itr.advance()) {
			iVec3 pos = itr.get_vec();
			Coord_Type ctype = pos.get_type();
			size_t index = get_index_ch(pos);
			size_t num_poles = 0;
			
			if (is_E_Point(ctype)) {
				if (dispersive_field_chunk[index])
					num_poles = dispersive_field_chunk[index]->get_num_poles();
				if (num_poles > MAX_NUM_POLES) num_poles = MAX_NUM_POLES;
				
				e_points[num_poles].push_back(index);
			}
			
			if (is_M_Point(ctype)) {
				if (dispersive_field_chunk[index])
					num_poles = dispersive_field_chunk[index]->get_num_poles();
				if (num_poles > MAX_NUM_POLES) num_poles = MAX_NUM_POLES;
				
				m_points[num_poles].push_back(index);
			}
		}
	}

	real Chunk::measure() const {
		real res = 0;
		for (auto item : eh)
			res += std::abs(item);
		return res;
	}

	void Chunk::add_dipole(const Func type, const iVec3& pos, const real amp, const real delay, const real fp) {
		bool is_E = is_E_Point(pos.get_type());
		auto& dipoles = (is_E) ? e_dipoles : m_dipoles;

		switch (type)
		{
		case ricker:
			if (dipoles.find(type) == dipoles.end()) {
				dipoles[type] = new CU_Dipole(jmd, [](const real time, const real delay, const real fp) -> real {
					real arg = (pi * fp * (time - delay));
					arg *= arg;
					return (1 - 2 * arg) * exp(-arg);
				}, is_E ? 0.5 * dt : 0, dt);
			}
			dipoles[type]->add_update_point(get_index_ch(pos), amp, delay, fp);
			break;

		case sine:
			if (dipoles.find(type) == dipoles.end()) {
				dipoles[type] = new CU_Dipole(jmd, [](const real time, const real delay, const real fp) -> real {
					return sin(2 * pi * fp * (time - delay));
				}, is_E ? 0.5 * dt : 0, dt);
			}
			dipoles[type]->add_update_point(get_index_ch(pos), amp, delay, fp);
			break;

		default:
			throw std::runtime_error("Function Type not Supported");
			break;
		}
	}

	void Chunk::add_projector(Plane_Wave const* projector) {
		if (e_tfsf == nullptr) {
			e_tfsf = new CU_TFSF<Plane_Wave>(jmd);
			m_tfsf = new CU_TFSF<Plane_Wave>(jmd);
		}

		e_tfsf->add_projector(projector);
		m_tfsf->add_projector(projector);
	}

	void Chunk::add_projector_update(const iVec3& jmd_pos, const real c1, const iVec3& pr_pos) {
		bool is_E = is_E_Point(jmd_pos.get_type());
		auto tfsf = is_E ? e_tfsf : m_tfsf;

		if (tfsf) 
			tfsf->add_update_point(get_index_ch(jmd_pos), c1, pr_pos);
		else
			throw std::runtime_error("No projector has been added");
	}
}
