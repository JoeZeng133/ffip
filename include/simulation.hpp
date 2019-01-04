#pragma once

#include <utility.hpp>
#include <geometry.hpp>
#include <medium.hpp>
#include <source.hpp>
#include <chunk.hpp>
#include <analysis.hpp>
#include <map>

#include <memory>
#include <vector>

namespace ffip {
	class Simulation;
	class Nearfield_Probe;
	class N2F_Box;
	class Flux_Box;
	
	class Simulation {
	private:
		real dt, dx;
		iVec3 sim_dim;
		
		int step;
		iVec3 sim_p1, sim_p2;			//domain coordinates of lower, upper points
		iVec3 ch_p1, ch_p2;		//domain coordinates of lower, upper points
		
		int num_proc{1};
		Chunk* chunk;
		Medium const* bg_medium{nullptr};
		
		std::vector<const Solid*> solids;
		std::vector<Current_Source*> current_sources;
		std::vector<Eigen_Source*> eigen_sources;
		
		std::vector<real> NF_freq;
		std::vector<fVec3> NF_pos;
		std::vector<Nearfield_Probe*> probes;
		
		//far fields calculation
		fVec3 N2F_p1, N2F_p2;	//coordinates of N2F box
		std::vector<real> N2F_omega;
		std::vector<fVec3> N2F_pos;
		N2F_Box* n2f_box{nullptr};

		//scattering coefficient calculation
		real_arr c_scat_omega;
		Flux_Box* scat_flux_box{nullptr};
		
		iVec3 tf_p1, tf_p2;

		PML PMLs[3][2];
		std::vector<real> kx, ky, kz, bx, by, bz, cx, cy, cz;
		
		void probe_init();
		void N2F_init();
		void medium_init();
		void source_init();
		void PML_init();
		void PML_init_helper(const PML& neg, const PML& pos, real_arr& k, real_arr& b, real_arr& c, const int p1, const int p2);
		void chunk_init();
		
	public:
		Simulation() = default;
		Simulation(const real _dx, const real _dt, const iVec3 _dim);
		
		void setup(const real _dx, const real _dt, const iVec3 _dim);
		void add_solid(Solid const* solid);
		void add_source(Source* source);
		void add_PML_layer(PML* PML);
		void add_nearfield_probe(const real freq, const fVec3& pos);
		void add_farfield_probe(const real freq, const fVec3& pos);
		void add_c_scat_freq(const real freq);
		
		void set_background_medium(Medium const* m);
		void set_num_proc(const int _num_proc);
		void init();
		
		void udf_unit();
		void udf_advance();
		void udf_output();
		void advance(std::ostream& os);
		/* field access functions*/
		real at(const fVec3& p, const Coord_Type ctype) const;					//access at float physical coordinates

		/* field access at computation coordinates*/
		real operator()(const fVec3& p, const Coord_Type ctype) const;		//access at float computation coordinates
		real operator()(const iVec3& p, const Coord_Type ctype) const;		//access at integer computation coordinates
		real operator()(const iVec3& p) const;								//raw access	

		/* getter */
		int get_step() const;
		real get_dt() const;
		real get_dx() const;
		iVec3 get_dim() const;
		Medium const* get_bg_medium() const;
		
		void output_nearfield(std::ostream& os);
		void output_farfield(std::ostream& os);
		void output_c_scat(std::ostream& os);
	};
}
