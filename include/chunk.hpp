#pragma once

#include <medium.hpp>
#include <source.hpp>


namespace ffip {
	class Inc_Internal;
	class Dipole;
	
	/* Current update element */
	struct CU {
		const size_t index;
		CU(const size_t _index): index(_index) {}
		virtual void update(std::vector<real>& jmd, const real time) = 0;
		virtual ~CU() {};
	};
	
	/* Incident current update */
	template<typename F>
	struct CU_Inc : CU {
		const real c1;
		F& ft;
		const iVec3 c2;
		
		CU_Inc(const size_t index, const real c1, F& ft, const iVec3 c2)
		: CU{index}, c1(c1), ft(ft), c2(c2) {}
		
		void update(std::vector<real>& jmd, const real time) override final {
			jmd[index] += c1 * ft(c2);
		}
	};

	struct CU_PML : CU {
		const real c1;
		const real c2;
		const real* p1;
		const real* p2;
		real ct{ 0 };

		CU_PML(const size_t index, const real c1, const real c2, const real* p1, const real* p2);
		void update(std::vector<real>& jmd, const real time) override final;
	};
	
	/* Dipole current update*/
	struct CU_Dipole : CU {
		using func = std::function<real(const real)>;
		const real c;
		func f;
		
		CU_Dipole(const size_t index, const real c, const func& f): CU(index), c(c), f(f) {}
		void update(std::vector<real>& jmd, const real time) override final {
			jmd[index] += c * f(time);
		}
	};
	
	class Chunk {
		/* constructor, indexing and field accessing members*/
	private:
		size_t ch_jump_x, ch_jump_y, ch_jump_z;
		iVec3 ch_dim, ch_p1, ch_p2, ch_origin, sim_p1, sim_p2;
		size_t jumps[8];
		int num_ones[8];
		real dx, dt;
		
		my_iterator ch_itr;
	public:
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
		real ave_helper(const int bit, const int index, const int jump, Args... args) const;
		
		
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
		std::vector<PML_Point> e_PML, m_PML;
		std::vector<real> kx, ky, kz;
	public:
		// concurrent PML points update
		void e_PML_push(const size_t rank);
		void m_PML_push(const size_t rank);
		void PML_update_helper(std::vector<PML_Point>& PML, const size_t rank);
		// PML layer initializations
		void PML_init(const real_arr& kx, const real_arr& ky, const real_arr& kz,
					  const real_arr& bx, const real_arr& by, const real_arr& bz,
					  const real_arr& cx, const real_arr&cy, const real_arr&cz);
		
		/* dipole sources*/
	private:
		std::vector<std::unique_ptr<Dipole>> e_dipoles_list;
		std::vector<std::unique_ptr<Dipole>> m_dipoles_list;
	
		/* incident waves */
	private:
		std::vector<std::unique_ptr<Inc_Internal>> inc_list;
		
		/* concurrency members*/
	private:
		Barrier* barrier{ new Barrier{ 1 } };
		size_t num_proc{ 1 };
	public:
		void set_num_proc(size_t _num_proc);
		void categorize_points();
		
		
		/* called by Simulation */
	public:
		// set medium to a specific point
		void set_medium_point(const iVec3& point, Medium_Ref const* medium_internal);
		// add an incident wave source
		void add_inc_internal(Inc_Internal* source);
		// add dipoles
		void add_dipoles(Dipole* dipole);
		
		/* getters */
	public:
		// get index jump in T direction
		template<typename T>
		const size_t get_ch_jump() const;
		// get k_T for updating curl
		template<typename T>
		const real get_k(const int id) const;
		real get_dt() const;
		real get_dx() const;
		iVec3 get_dim() const;
		iVec3 get_origin() const;
		iVec3 get_p1() const;
		iVec3 get_p2() const;
		iVec3 get_pos(const size_t index) const;
		size_t get_index_ch(const iVec3& p) const;			//get index relative to chunk origin
		
		/* update members*/
	public:
		/* Jd, Md (currenst, curl) updates
		  D(n + 1) - D(n) = Jd(n + 0.5) = curl(H(n + 0.5)) - Ji(n + 0.5)
		  B(n + 1) - B(n) = Md(n + 0.5) = -(curl(E(n + 0.5)) + Mi(n + 0.5))
		*/
		void update_Jd(const real time, const size_t rank);		//curl H
		void update_Md(const real time, const size_t rank);		//curl E
		template<typename T>
		void update_JMd_helper(const iVec3 p1, const iVec3 p2, const size_t rank);	//helper function for calculating curl
		/* Material updates */
		void update_D2E(const real time, const size_t rank);
		void update_B2H(const real time, const size_t rank);
		void update_DEHB_helper(const Coord_Type F, const iVec3 p1, const iVec3 p2, const size_t rank);
		
		void update_D2E_v2(const real time, const size_t rank);
		void update_B2H_v2(const real time, const size_t rank);
		
		/* MPI updates of the boundary */
		void update_ghost_H(const real time);
		void update_ghost_E(const real time);
		
		/* miscellaneous*/
	public:
		real measure() const;
		
		/* Generic Currents Updates */
	private:
		std::vector<size_t> e_cu_indexes;
		std::vector<size_t> m_cu_indexes;
		std::vector<CU*> e_current_updates;
		std::vector<CU*> m_current_updates;
		
		std::atomic<size_t> top;
		std::mutex e_currents;
		std::mutex m_currents;
	public:
		void add_current_update(CU* cu, const iVec3& pos);
		void add_e_current_update(CU* cu);
		void add_m_current_update(CU* cu);
		void organize_current_updates();
		void dynamic_e_current_update(const real time, const size_t chunk_size);
		void dynamic_m_current_update(const real time, const size_t chunk_size);
		void reset_scheduler();
	};
	
	template<typename... Args>
	real Chunk::ave_helper(const int bit, const int index, const int jump, Args... args) const{
		constexpr int bit_N =1 << sizeof...(Args);
		if(bit & bit_N) {
			return ave_helper(bit, index + jump, args...) + ave_helper(bit, index - jump, args...);
		}else {
			return ave_helper(bit, index, args...);
		}
	}
	
	template<typename T>
	void Chunk::update_JMd_helper(const iVec3 p1, const iVec3 p2, const size_t rank) {
		/* Curl updates without PML */
		using dir_base = typename T::dir_base;
		using x1 = typename dir_base::x1;
		using x2 = typename dir_base::x2;
		
		auto tmp = get_component_interior(p1, p2, T::ctype);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		int ch_jump_x1 = get_ch_jump<x1>();
		int ch_jump_x2 = get_ch_jump<x2>();
		
		for(auto itr = my_iterator(p1_ch, p2_ch, p1_ch.get_type(), rank, num_proc); !itr.is_end(); itr.advance()) {
			int index = itr.x * ch_jump_x + itr.y * ch_jump_y + itr.z * ch_jump_z;
			jmd[index] = (eh[index + ch_jump_x1] - eh[index - ch_jump_x1]) / dx / get_k<x1>(choose<x1>::get(itr)) - (eh[index + ch_jump_x2] - eh[index - ch_jump_x2]) / dx / get_k<x2>(choose<x2>::get(itr));
		}
	}
}
