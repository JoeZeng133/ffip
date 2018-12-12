#pragma once

#include <medium.hpp>
#include <source.hpp>


namespace ffip {
	class Chunk {
	private:
		std::vector<real> eh;			//E(t), H(t - 0.5dt)
		std::vector<real> eh1;			//E(t - dt), H(t - 1.5dt)
		std::vector<real> jmd;			//Jd(t), M(t - 0.5dt)
		
		std::vector<Dispersive_Field*> dispersive_field_chunk;	//array of dispersive_field
		std::vector<Medium_Ref const*> medium_chunk;				//place holder for medium_internal
		
		/* PML members*/
		std::vector<PML_Point> e_PML, h_PML;
		std::vector<real> kx, ky, kz;
		
		/* Coordinate Members*/
		int ch_jump_x, ch_jump_y, ch_jump_z;
		iVec3 ch_dim, ch_p1, ch_p2, ch_origin, sim_p1, sim_p2;
		
		int jump[8];
		
		/* source */
		std::vector<Source_Internal*> source_list;				//
		real dx, dt;
	public:
		//given start, end domain_comp coord and initialize an empty Chunk
		Chunk(const iVec3& _sim_p1, const iVec3& _sim_p2, const iVec3& _ch_p1, const iVec3& _ch_p2, real _dx, real _dt);
		void set_medium_point(const iVec3& point, Medium_Ref const* medium_internal);
		void add_source_internal(Source_Internal* const source);
		void PML_init(const real_arr& kx, const real_arr& ky, const real_arr& kz,
					  const real_arr& bx, const real_arr& by, const real_arr& bz,
					  const real_arr& cx, const real_arr&cy, const real_arr&cz);		//provide with PML parameters in the entire simulation
		//generic getter
		template<typename T>
		const int get_ch_jump() const;
		
		template<typename T>
		const real get_k(const int id) const;
		
		//getter for dimension
		iVec3 get_dim() const;
		iVec3 get_origin() const;
		iVec3 get_p1() const;
		iVec3 get_p2() const;
		int get_index_ch(const iVec3& p) const;			//get index relative to chunk origin
		
		/* Jd, Md (currenst, curl) updates*/
		void update_Jd(const real time, const int num_proc);		//D(n + 1) - D(n) = Jd(n + 0.5) = curl(H(n + 0.5)) - Ji(n + 0.5)
		void update_Md(const real time, const int num_proc);		//B(n + 1) - B(n) = Md(n + 0.5) = -(curl(H(n + 0.5)) + Mi(n + 0.5))
		
		template<typename T>
		void update_JMd_helper(const iVec3 p1, const iVec3 p2);
		
		/* Material updates */
		void update_D2E(const real time, const int num_proc);
		void update_B2H(const real time, const int num_proc);
		void update_DEHB_helper(const Coord_Type F, const iVec3 p1, const iVec3 p2);
		
		/* MPI updates of the boundary */
		void update_padded_H(const real time);
		void update_padded_E(const real time);
		
		/* field access functions*/
		real at(const fVec3& p, const Coord_Type ctype) const;					//access at float physical coordinates
		
		/* field access at computation coordinates*/
		real operator()(const fVec3& p, const Coord_Type ctype) const;		//access at float computation coordinates
		real operator()(const iVec3& p, const Coord_Type ctype) const;		//access at integer computation coordinates
		real operator()(const iVec3& p) const;								//raw access			
		
		/* average field according to bit patterns
		 111 = average over 8 points
		 110, 011, 101 = average over 4 points
		 100, 001, 010 = average over 2 points
		 000 = no average
		 take advantage of generic programming to efficiently compute the average without too much overhead
		 */
		template<typename... Args>
		real ave_helper(const int bit, const int index, const int jump, Args... args) const{
			constexpr int bit_N =1 << sizeof...(Args);
			if(bit & bit_N) {
				return ave_helper(bit, index + jump, args...) + ave_helper(bit, index - jump, args...);
			}else {
				return ave_helper(bit, index, args...);
			}
		}
		
		real ave(const int bit, const int index) const;

		/* interpolation helper functions*/
		real interp_helper(const real* data, const real w) const;

		template<typename... Args>
		real interp_helper(const real* data, const real w, Args... args) const{
			constexpr int N = sizeof...(Args);
			constexpr int shift = 1 << N;

			return interp_helper(data, args...) * (1 - w) + interp_helper(data + shift, args...) * w;
		}
	};
	
	
	template<typename T>
	void Chunk::update_JMd_helper(const iVec3 p1, const iVec3 p2) {
		/* Curl updates without PML */
		using dir_base = typename T::dir_base;
		using x1 = typename dir_base::x1;
		using x2 = typename dir_base::x2;
		
		auto tmp = get_component_interior(p1, p2, T::ctype);
		auto p1_ch = tmp.first - ch_origin;
		auto p2_ch = tmp.second - ch_origin;
		int ch_jump_x1 = get_ch_jump<x1>();
		int ch_jump_x2 = get_ch_jump<x2>();
		
		for(auto itr = my_iterator(p1_ch, p2_ch, p1_ch.get_type()); !itr.is_end(); itr.advance()) {
			int index = itr.x * ch_jump_x + itr.y * ch_jump_y + itr.z * ch_jump_z;
			
			jmd[index] = (eh[index + ch_jump_x1] - eh[index - ch_jump_x1]) / dx / get_k<x1>(choose<x1>::get(itr)) - (eh[index + ch_jump_x2] - eh[index - ch_jump_x2]) / dx / get_k<x2>(choose<x2>::get(itr));
		}
	}
}
