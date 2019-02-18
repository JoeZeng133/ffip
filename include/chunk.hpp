#pragma once

#include <medium.hpp>
//#include <source.hpp>
#include <array>
#include <map>

namespace ffip {
	class Chunk;
	class Plane_Wave;
	class Inc_Source;

	class CU {
	public:
		virtual void update_static(const size_t rank, Barrier* barrier) = 0;
		virtual void update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) = 0;
		virtual void organize() = 0;
	};

	/* for sources acting on TF/SF surface */
	template<typename F>
	class CU_Source : public CU{
	private:
		std::vector<real>& jmd;
		F& projector;
	
		struct Update_Point {
			size_t index;
			real c1;
			iVec3 c2;
		};
		std::vector<Update_Point> points;
		std::vector<size_t> index_unique;	//duplicate indexes can appear, used in update_dynamic

	public:
		CU_Source() = delete;
		CU_Source(std::vector<real>& jmd, F& projector) : jmd(jmd), projector(projector) {}

		void update_static(const size_t rank, Barrier* barrier) override;
		void update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) override;
		void organize() override;
		void add_update_point(const size_t index, const real c1, const iVec3 c2);
	};

	template<typename F>
	void CU_Source<F>::update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) {
		
	}

	template<typename F>
	void CU_Source<F>::organize() {
		std::sort(points.begin(), points.end(), [](const Update_Point& a, const Update_Point& b) {
			return a.index < b.index;
		});
	}

	template<typename F>
	void CU_Source<F>::add_update_point(const size_t index, const real c1, const iVec3 c2) {
		points.push_back(Update_Point{ index, c1, c2 });
	}

	template<typename F>
	void CU_Source<F>::update_static(const size_t rank, Barrier* barrier) {
		size_t idx1, idx2;
		vector_divider(points, rank, barrier->get_num_proc(), idx1, idx2);
		//adjust to prevent slicing the updates of a single index
		while (idx1 > 0 && idx1 < points.size() && points[idx1].index == points[idx1 - 1].index) --idx1;
		while (idx2 > 0 && idx2 < points.size() && points[idx2].index == points[idx2 - 1].index) --idx2;

		for (; idx1 < idx2; ++idx1) {
			auto& point = points[idx1];
			jmd[point.index] += point.c1 * projector(point.c2);
		}
	}

	/* for sources on PML points*/
	class CU_PML : public CU{
	private:
		std::vector<real>& jmd;
		std::vector<real>& eh;

		struct Update_Point {
			real phi1;
			real b1;
			real c1;
			real k1;

			real phi2;
			real b2;
			real c2;
			real k2;
		};

		std::array<size_t, 3> strides;
		std::vector<Update_Point> points[3];
		std::vector<size_t> indexes[3];
	public:
		CU_PML() = delete;
		CU_PML(std::vector<real>& jmd, std::vector<real>& eh, std::array<size_t, 3> strides);

		void organize() override;
		void update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) override;
		void update_static(const size_t rank, Barrier* barrier) override;

		template<typename x3>
		void add_update_point(const size_t index, real b1, real c1, real k1, real b2, real c2, real k2) {
			using x1 = typename x3::x1;
			using x2 = typename x3::x2;

			indexes[x3::val].push_back(index);
			points[x3::val].push_back(Update_Point{0, b1, c1, k1, 0, b2, c2, k2});
		}
	};

	/* for dipole currents, sin, ricker */
	class CU_Dipole : public CU{
	private:
		using func = std::function<real(const real, const real, const real)>;

		std::vector<real>& jmd;
		func f;
		real time, dt;
		int step{ 0 };

		struct Update_Point {
			size_t index;
			real c;
			real d;
			real fp;
		};
		std::vector<size_t> index_unique;	//duplicate indexes might appear
		std::vector<Update_Point> points;
	public:
		CU_Dipole() = delete;
		CU_Dipole(std::vector<real>& jmd, const func& f, const real time, const real dt);

		void organize() override;
		void update_dynamic(std::atomic<size_t>& sync_index, Barrier* barrier) override;
		void update_static(const size_t rank, Barrier* barrier) override;
		void add_update_point(const size_t index, const real amp, const real delay, const real fp);
	};
	
	class Chunk {
		/* constructor, indexing and field accessing members*/
	private:
		size_t ch_stride_x, ch_stride_y, ch_stride_z;
		iVec3 ch_dim, ch_p1, ch_p2, ch_origin, sim_p1, sim_p2;
		size_t strides[8];
		int num_ones[8];
		real dx, dt;
		Config config;
	
	public:
		Chunk(const Config& config);
		Chunk(const iVec3& _sim_p1, const iVec3& _sim_p2, const iVec3& _ch_p1, const iVec3& _ch_p2, real _dx, real _dt);
		//thread safe access at physical coordinates
		real at(const fVec3& p, const Coord_Type ctype) const;
		//thread safe access at float comp coord with ctype
		real operator()(const fVec3& p, const Coord_Type ctype) const;
		//thread safe access at integer comp coord with ctype
		real operator()(const iVec3& p, const Coord_Type ctype) const;
		//thread safe raw access
		real operator()(const iVec3& p) const;
		//thread safe indexing
		real operator[](const size_t index) const;
		/* average field according to bit patterns
		 111 = average over 8 points
		 110, 011, 101 = average over 4 points
		 100, 001, 010 = average over 2 points
		 000 = no average
		 */
		real ave(const int bit, const int index) const;
		// helper function for averaging according to bit patterns
		template<typename... Args>
		real ave_helper(const int bit, const int index, const int stride, Args... args) const;
		 
		
		/* data storage and medium*/
	private:
		std::vector<real> eh;			//E(t), H(t - 0.5dt)
		std::vector<real> eh1;			//E(t - dt), H(t - 1.5dt)
		std::vector<real> jmd;			//Jd(t), M(t - 0.5dt)
		std::vector<Dispersive_Field*> dispersive_field_chunk;	//for dispersive material only
		std::vector<Medium_Ref const*> medium_chunk;			//place holder for medium_internal
		
		static constexpr int MAX_NUM_POLES = 10;
		std::vector<size_t> e_points[MAX_NUM_POLES];	//categorize points by number of poles
		std::vector<size_t> m_points[MAX_NUM_POLES];	//used for better work schedule
		
		/* PML members*/
	private:
		std::vector<real> kx, ky, kz;
	public:
		// PML layer initializations
		void PML_init(const std::array<real_arr, 3>& k, const std::array<real_arr, 3>& b, const std::array<real_arr, 3>& c);
		template<typename F>
		void PML_init_helper(const std::array<real_arr, 3>& k, const std::array<real_arr, 3>& b, const std::array<real_arr, 3>& c);
	
	public:
		void init();
		void categorize_points();
		
		/* called by Simulation */
	public:
		// set medium to a specific point
		void set_medium_point(const iVec3& point, Medium_Ref const* medium_internal);

		/* getters */
	public:
		// get index stride in T direction
		template<typename T>
		const size_t get_ch_stride() const {
			if constexpr(std::is_same_v<T, dir_x_tag>)
				return ch_stride_x;
			else if constexpr(std::is_same_v<T, dir_y_tag>)
				return ch_stride_y;
			else if constexpr(std::is_same_v<T, dir_z_tag>)
				return ch_stride_z;
		}
		// get k_T for updating curl
		template<typename T>
		const real get_k(const int id) const {
			if constexpr(std::is_same_v<T, dir_x_tag>)
				return kx[id];
			else if constexpr(std::is_same_v<T, dir_y_tag>)
				return ky[id];
			else if constexpr(std::is_same_v<T, dir_z_tag>)
				return kz[id];
		}

		real get_dt() const;
		real get_dx() const;
		iVec3 get_dim() const;
		iVec3 get_origin() const;
		iVec3 get_p1() const;
		iVec3 get_p2() const;
		iVec3 get_pos(size_t index) const;
		size_t get_index_ch(const iVec3& p) const;			//get index relative to chunk origin
		
		/* update members*/
	public:
		/* Jd, Md (currenst, curl) updates
		  D(n + 1) - D(n) = Jd(n + 0.5) = curl(H(n + 0.5)) - Ji(n + 0.5)
		  B(n + 1) - B(n) = Md(n + 0.5) = -(curl(E(n + 0.5)) + Mi(n + 0.5))
		*/
		void update_Jd(const real time, const size_t rank, Barrier* barrier);		//curl H
		void update_Md(const real time, const size_t rank, Barrier* barrier);		//curl E
		template<typename T>
		void update_JMd_helper(const size_t rank, Barrier* barrier);	//helper function for calculating curl
		/* Material updates */
		void update_D2E_v2(const real time, const size_t rank, Barrier* barrier);
		void update_B2H_v2(const real time, const size_t rank, Barrier* barrier);
		
		void update_ud(const real time, const size_t rank, Barrier* barrier);
		
		/* ghost points updates, periodic boundaries */
		void update_ghost_H(const real time, const size_t rank, Barrier* barrier);
		void update_ghost_E(const real time, const size_t rank, Barrier* barrier);
		
		//peridoic boundary conditions, assuming normal incident angle
		template<typename Dir>
		void update_ghost_helper(const size_t rank, Barrier* barrier);
		
		/* miscellaneous*/
		real measure() const;
		
		/* Generic Currents Updates */
	private:
		CU_PML * e_PML{ nullptr }, *m_PML{ nullptr };			//PML updates
		CU_Source<Plane_Wave> * e_source{ nullptr }, *m_source{ nullptr };	//plane_wave projector
		std::map<Func, CU_Dipole*> e_dipoles;						//electric dipole currents
		std::map<Func, CU_Dipole*> m_dipoles;						//magnetic dipole currents

	public:
		//add dipole sources, currently support only ricker, sin (taks 3 arguments)
		void add_dipole(const Func type, const iVec3& pos, const real amp, const real delay, const real fp);
		void set_projector(Plane_Wave& projector);
		void add_projector_update(const iVec3& pos, const real amp, const iVec3& inc_pos);
	};
	
	template<typename Dir>
	void Chunk::update_ghost_helper(const size_t rank, Barrier *barrier) {
		auto p1 = ch_p1 - 1;
		auto p2 = ch_p2 + 1;
		auto dp = ch_p2 - ch_p1;
		iVec3 jump{0, 0, 0};
		
		choose<Dir>::get(jump) = choose<Dir>::get(dp);
		
		//low from high
		choose<Dir>::get(p1) = choose<Dir>::get(p2) = choose<Dir>::get(ch_p1 - 1);
		for(auto itr = my_iterator(p1, p2, All, rank, barrier->get_num_proc()); !itr.is_end(); itr.advance()) {
			auto vec = itr.get_vec();
			eh[get_index_ch(vec)] = eh[get_index_ch(vec + jump)];
		}
		
		//high from low
		choose<Dir>::get(p1) = choose<Dir>::get(p2) = choose<Dir>::get(ch_p2 + 1);
		for(auto itr = my_iterator(p1, p2, All, rank, barrier->get_num_proc()); !itr.is_end(); itr.advance()) {
			auto vec = itr.get_vec();
			eh[get_index_ch(vec)] = eh[get_index_ch(vec - jump)];
		}
	}

	template<typename F>
	void Chunk::PML_init_helper(const std::array<real_arr, 3>& k, const std::array<real_arr, 3>& b, const std::array<real_arr, 3>& c) {
		using x3 = typename F::dir_base;
		using x1 = typename F::dir_base::x1;
		using x2 = typename F::dir_base::x2;

		CU_PML * w_PML = (is_E_Point(F::ctype)) ? e_PML : m_PML;

		//std::cout << "PML: " << F::ctype << " " << x1::val << " " << x2::val << "\n";

		for (auto itr = my_iterator(ch_p1, ch_p2, F::ctype); !itr.is_end(); itr.advance()) {
			iVec3 pos = itr.get_vec();
			if (Is_Inside_Box(config.phys_p1, config.phys_p2, pos))
				continue;

			size_t index = get_index_ch(pos);
			int x1_index = choose<x1>::get(pos) - choose<x1>::get(ch_p1);
			int x2_index = choose<x2>::get(pos) - choose<x2>::get(ch_p1);
			
			w_PML->add_update_point<x3>(index, 
				b[x1::val][x1_index], c[x1::val][x1_index] / dx, 1 / dx / k[x1::val][x1_index],
				b[x2::val][x2_index], c[x2::val][x2_index] / dx, 1 / dx / k[x2::val][x2_index]);
		}
	}
	
	template<typename... Args>
	real Chunk::ave_helper(const int bit, const int index, const int stride, Args... args) const{
		constexpr int bit_N =1 << sizeof...(Args);
		if(bit & bit_N) {
			return ave_helper(bit, index + stride, args...) + ave_helper(bit, index - stride, args...);
		}else {
			return ave_helper(bit, index, args...);
		}
	}
	
	template<typename F>
	void Chunk::update_JMd_helper(const size_t rank, Barrier* barrier) {
		/* Curl updates without PML */
		using x3 = typename F::dir_base;
		using x1 = typename x3::x1;
		using x2 = typename x3::x2;
		
		auto tmp = get_component_interior(config.phys_p1, config.phys_p2, F::ctype);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		int ch_stride_x1 = get_ch_stride<x1>();
		int ch_stride_x2 = get_ch_stride<x2>();
		size_t num_proc = barrier->get_num_proc();
		
		for(auto itr = my_iterator(p1_ch, p2_ch, p1_ch.get_type(), rank, num_proc); !itr.is_end(); itr.advance()) {
			int index = itr.x * ch_stride_x + itr.y * ch_stride_y + itr.z * ch_stride_z;
			jmd[index] = (eh[index + ch_stride_x1] - eh[index - ch_stride_x1]) / dx / get_k<x1>(choose<x1>::get(itr)) - (eh[index + ch_stride_x2] - eh[index - ch_stride_x2]) / dx / get_k<x2>(choose<x2>::get(itr));
		}
	}
}
