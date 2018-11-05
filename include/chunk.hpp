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
		std::vector<Medium_Internal*> medium_chunk;				//place holder for medium_internal
		
		/* PML members*/
		std::vector<PML_Point> e_PML, h_PML;
		std::vector<real> kx, ky, kz;
		
		/* Coordinate Members*/
		int ch_jump_x, ch_jump_y, ch_jump_z;
		iVec3 ch_dim, ch_p1, ch_p2, ch_origin, sim_p1, sim_p2;
		
		std::vector<int> jump[8];
		
		/* source */
		std::vector<Source_Internal*> source_list;				//
		real dx, dt;
	public:
		//given start, end domain_comp coord and initialize an empty Chunk
		Chunk(const iVec3& _sim_p1, const iVec3& _sim_p2, const iVec3& _ch_p1, const iVec3& _ch_p2, real _dx, real _dt);
		void set_medium_point(const iVec3& point, Medium_Internal* const  medium_internal);
		void add_source_internal(Source_Internal* const source);
		void PML_init(const real_arr& kx, const real_arr& ky, const real_arr& kz,
					  const real_arr& bx, const real_arr& by, const real_arr& bz,
					  const real_arr& cx, const real_arr&cy, const real_arr&cz);		//provide with PML parameters in the entire simulation
		
		//getter for dimension
		iVec3 get_dim() const;
		iVec3 get_origin() const;
		iVec3 get_p1() const;
		iVec3 get_p2() const;
		int get_index_ch(const iVec3& p) const;
		
		/* Jd, Md (currenst, curl) updates*/
		void update_Jd(real time);
		void update_Md(real time);
		void update_JMd_helper(const Coord_Type Fx, const Coord_Type fy, const Coord_Type fz, std::vector<PML_Point>& PML);
		
		/* Material updates */
		void update_D2E(real time);
		void update_H2B(real time);
		void update_DEHB_helper(const Coord_Type F);
		
		/* MPI members */
		void update_padded_H(real time);
		void update_padded_E(real time);
		
		/* */
		real operator[](const iVec3& p) const;		//access field at a domain computation
		real get_field(const iVec3& p, const Coord_Type ctype) const;	//get field of ctype at a domain comp coordinate
		real get_prev_field(const iVec3& p, const Coord_Type ctype) const;
		
		
		/* average field according to bit patterns
		 111 = average over 8 points
		 110, 011, 101 = average over 4 points
		 100, 001, 010 = average over 2 points
		 000 = no average
		 take advantage of generic programming to efficiently compute the average without too much overhead
		 */
		template<typename... Args>
		int ave_helper(const int bit, const int index, const int jump, Args... args) const{
			if(bit && (1 << sizeof...(Args))) {
				return ave_helper(bit, index + jump, args...) + ave_helper(bit, index - jump, args...);
			}else {
				return ave_helper(bit, index, args...);
			}
		}
		
		int ave(const int bit, const int index) const;
	};
}
