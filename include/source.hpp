#pragma once

#include <utility.hpp>
#include <chunk.hpp>

namespace ffip {
	/* Plane_Wave projector with different excitation functions
	 it is obligation of script to make sure it is compatible with the main simulation so that it can be used as a projector
	 it can used as a stand-alone 1D simulation
	 exciation at -dim_neg * dx, the overall effective region is [-dim_neg, dim_pos] * dx (excluding PML region)
	 */
	class Chunk;
	class Plane_Wave;
	
	class Plane_Wave {
	private:
		// constructor parameters
		real dx, dt;
		int dim_neg, dim_pos, n;
		PML pml;
		real amp;
		Direction polorization;
		
		std::function<real(const real)> f; //excitation functor
		real courant;
		int dim, origin;
		
		int time_step{0};
		
		std::vector<real> eh;
		std::vector<real> psi_pos;
		std::vector<real> psi_neg;
		
		std::vector<real> b_z;
		std::vector<real> c_z;
		std::vector<real> k_z;

		real er{1}, ur{1};
	public:
		Plane_Wave() = delete;										//no default constructor
		Plane_Wave(const real _dx, const real _dt, const int _dim_neg, const int _dim_pos);
		Plane_Wave(const Plane_Wave&) = default;					//copy
		Plane_Wave& operator=(Plane_Wave&) = default;
		Plane_Wave(Plane_Wave&&) = default;							//move
		Plane_Wave& operator=(Plane_Wave&&) = default;
		
		void set_excitation(const std::function<real(const real)>& _f, const real _amp = 1, const Direction _polorization = Direction::X);	//set up excitation position, function choice, amplitude, excitation position

		void set_PML(const PML& _pml);							//set up PML layer in x+
		void init();														//after set-up, it has to be initialized
		void set_medium(const real er, const real ur);
		
		//access at float physical coordinates
		real at(const fVec3& p, const Coord_Type ctype) const;					
		//access at float computation coordinates, with type specified
		real operator()(const fVec3& p, const Coord_Type ctype) const;
		//raw access, no type guaranteed
		real operator()(const iVec3& p) const;
		
		void hard_E(real time);		//hard E source to excite wave
		void advance();		//H(t-0.5dt), E(t) -> H(t+0.5dt), E(t+dt)
		void update_H();			//update Hy field
		void update_E();			//update Ex field
		int get_time_step();
		
		void udf_init();
		void udf_advance();
	};
	
	/* TFSF surfaces */
	struct TFSF_Surface {
		iVec3 d1, d2;			//surface corner points in domain computation coordinates
		Direction dir;			//surface direction x, y, z
		Side side;				//surface side +1, -1
		int sign_correction;	//=1, when TF surface, =-1 when SF
		
		iVec3 get_pos_offset() const;	//position offset for a particular TF/SF surface
		void TF2SF();					//transit from TF surface to a SF surface
		
		TFSF_Surface(const iVec3& _d1, const iVec3& _d2, const Direction _dir, const Side side, int _sign_correction);
		TFSF_Surface(const std::pair<iVec3, iVec3>& d, const Direction dir, const Side side, int sign_correction);
	};
	
	/* Projector_Source provided with different types of projectors, plane wave so far */
	class Projector_Source {
	public:
		Projector_Source() = default;
		void init(const Config& config);
		void add_projector(const Plane_Wave& projector);
		void push(Chunk* chunk);

		template<typename x1>
		void push_helper(Chunk* chunk, const TFSF_Surface& face);
		void advance();
		
	private:
		std::vector<Plane_Wave*> projectors;
		Config config;
		std::vector<TFSF_Surface> tfsf_list;
	};

	template<typename x1>
	void Projector_Source::push_helper(Chunk* chunk, const TFSF_Surface& face) {
		using x2 = typename x1::x1;
		using x3 = typename x1::x2;

		// direction should be paralell to the face
		if (x1::val == face.dir)
			return;

		real amp = face.side * ((face.dir == x2::val) ? 1 : -1) / config.dx;
		iVec3 pos_offset = face.get_pos_offset();

		for (auto itr = my_iterator(face.d1, face.d2, x1::E); !itr.is_end(); itr.advance()) {
			chunk->add_projector_update(itr.get_vec(), amp, itr.get_vec() + pos_offset);
		}

		for (auto itr = my_iterator(face.d1, face.d2, x1::H); !itr.is_end(); itr.advance()) {
			chunk->add_projector_update(itr.get_vec(), amp, itr.get_vec() + pos_offset);
		}
	}
}


