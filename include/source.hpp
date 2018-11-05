#pragma once

#include <utility.hpp>


namespace ffip {
	/* Plane_Wave projector with different excitation functions (currently using Sinuosuidal
	 it is obligation of script to make sure it is compatible with the main simulation so that it can be used as a projector
	 it can used as a stand-alone 1D simulation
	 */
	class Plane_Wave {
	private:
		// constructor parameters
		real dx, dt;
		int n;
		
		real excite_p{0};
		PML PML_pos, PML_neg;
		
		real amp;
		Direction polorization;
		Sinuosuidal_Func f{0};										//excitation functor
		
		real courant;
		int n_neg, n_pos, n_tot;
		int excite_p1, excite_p2, c1, c2;
		
		std::vector<real> eh;
		std::vector<real> psi_pos;
		std::vector<real> psi_neg;
		
		std::vector<real> b_z;
		std::vector<real> c_z;
		std::vector<real> k_z;
		
	public:
		Plane_Wave() = delete;										//no default constructor
		Plane_Wave(const real _dx, const real _dt, const int _n);	//constructor
		Plane_Wave(const Plane_Wave&) = default;					//copy
		Plane_Wave& operator=(Plane_Wave&) = default;
		Plane_Wave(Plane_Wave&&) = default;							//move
		Plane_Wave& operator=(Plane_Wave&&) = default;
		
		void set_excitation(const Sinuosuidal_Func& _f,
							const real phys_excite_p = 0,
							const real _amp = 1,
							const Direction _polorization = Direction::X);	//set up excitation position, function choice, amplitude, excitation position
		void set_PML_neg(const PML& _PML_neg);							//set up PML layer in x-
		void set_PML_pos(const PML& _PML_pos);							//set up PML layer in x+
		void init();														//after set-up, it has to be initialized
		
		real operator[](const iVec3& p) const;								//allow fields access through domain computation coordinates
		
		void hard_E(real time);		//hard E source to excite wave
		void update(real time);		//H(t-0.5dt), E(t) -> H(t+0.5dt), E(t+dt)
		void update_H();			//update Hy field
		void update_E();			//update Ex field
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
		
		static iVec3 vec3_base[3];	//bases vector in 3 dimension
	};
	
	/* source interfaces for use in chunk*/
	class Source_Internal {
	public:
		virtual ~Source_Internal() {}
		virtual void update_Jd(std::vector<real>& jmd) = 0;
		virtual void update_Md(std::vector<real>& jmd) = 0;
		virtual void get_Jd(real time) = 0;
		virtual void get_Md(real time) = 0;
	};
	
	/*allow homogeneous storage of sources
	  source types include eigen source and current source
	  they are taken care of separately since they are different in nature
	 prefix ch = chunk's = property of chunk
	 postfix ch = with respect to chunk
	 */
	
	class Source {
	public:
		virtual ~Source() {}
		virtual Source_Internal* get_source_internal() = 0;
	};
	
	/* current source J(x) = A(x) * f(t)*/
	class Current_Source : public Source {
	public:
		/* constructor access by scripts*/
		Current_Source(const GriddedInterp& _interp, const Sinuosuidal_Func& _phase, const Coord_Type _ctype);		//constructor from interpolation class
		Current_Source(const std::string& filename, const Sinuosuidal_Func& _phase, const Coord_Type _ctype);		//constructor from interpolation file
		
		/* override function */
		Source_Internal* get_source_internal() override;
		
		/* initialization access by simulation class*/
		void init(const real dx,
				  const iVec3 _ch_p1, const iVec3 _ch_p2,
				  const iVec3 _ch_dim, const iVec3 _ch_origin);		//pass chunk information and simulation information
		
	private:
		GriddedInterp interp;			//A(x) for the current
		Sinuosuidal_Func phase;		//f(t) for the current
		Coord_Type ctype;			//Jx, Jy, Jz, Mx, My, Mz
		
		iVec3 p1, p2;					//specify region to be updated
		iVec3 ch_p1, ch_p2;				//corner coordinates of chunk
		iVec3 ch_dim, ch_origin;	//
		real interp_ratio;			//interpolation ratio
		std::vector<real> amp;
	};
	
	/* current source used in chunk*/
	class Current_Internal : public Source_Internal{
	public:
		Current_Internal(const Sinuosuidal_Func _phase, const std::vector<real>& _amp,
						 const iVec3 _ch_dim, const iVec3 _ch_origin,
						 const iVec3 _p1, const iVec3 _p2);
		
		/* override functions*/
		void update_Jd(std::vector<real>& jmd) override;
		void update_Md(std::vector<real>& jmd) override;
		void get_Jd(real time) override;
		void get_Md(real time) override;
		
		void update_helper(std::vector<real> &jmd);	//update function
		
	private:
		Sinuosuidal_Func phase;		//f(t) for the current
		std::vector<real> amp;		//modified amplitude
		iVec3 ch_dim, ch_origin;			//chunk dimension
		iVec3 p1, p2;				//p1, p2 specify the region in chunk to be updated
		
		Coord_Type ctype;						//Coord type
		iVec3 dp;								//p2 - p1
		int ch_jump_x, ch_jump_y, ch_jump_z;	//jumps in chunk
		int amp_jump_x, amp_jump_y, amp_jump_z;	//jumps in amp
		int base_index_ch;						//chunk p1 index
		real cur_phase;							//f(t) updated by get_Jd or get_Md
	};
	
	/* eigen_source provided with different types of projectors, plane wave so far */
	class Eigen_Source : public Source {
	public:
		Eigen_Source(const Plane_Wave& _projector);
		Source_Internal* get_source_internal() override;
		void init(const iVec3 _tf1, const iVec3 _tf2,
				  const iVec3& _ch_dim, const iVec3& _ch_origin,
				  const iVec3 _ch_p1, const iVec3 _ch_p2);
		
	private:
		Plane_Wave projector;
		iVec3 ch_p1, ch_p2, tf1, tf2, ch_dim, ch_origin;
	};
	
	/* Eigen source used in chunk*/
	class Eigen_Internal : public Source_Internal {
		static int TFSF_Mat[3][3];
		
	public:
		Eigen_Internal(const std::vector<TFSF_Surface> _tsfs_list, const Plane_Wave& _projector,
					   const iVec3& _ch_dim, const iVec3& _ch_origin);	//copy every thing inside
		Eigen_Internal(const Eigen_Internal&) = delete; //no copy
		Eigen_Internal& operator=(const Eigen_Internal&) = delete;
		Eigen_Internal(Eigen_Internal&&) = delete; //no move
		Eigen_Internal& operator=(Eigen_Internal&&) = delete;
		
		/* overriding functions*/
		void update_Jd(std::vector<real>& jmd) override;
		void update_Md(std::vector<real>& jmd) override;
		void get_Jd(double time) override;
		void get_Md(double time) override;
		
		void update_helper(std::vector<real>& jmd, const TFSF_Surface face, Coord_Type type);	//helper functions for update_Jd, update_Md
	private:
		std::vector<TFSF_Surface> tfsf_list;
		Plane_Wave projector;
		iVec3 ch_dim, ch_origin;	//information about chunk
		int jump_x, jump_y, jump_z;
	};
}


