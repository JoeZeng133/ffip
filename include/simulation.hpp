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
		iVec3 sim_p1, sim_p2;	//simulation region corner points
		iVec3 ch_p1, ch_p2;		//chunk region corner points
		iVec3 tf_p1, tf_p2;		//total field region corner points
		iVec3 phys_p1, phys_p2; //physical region corner points
		int sf_depth{0};		//Scattered field depth in computation unit
		int tf_padded_depth{ 0 };//total field padded depth in computation unit
		
		Barrier* barrier{new Barrier{1}};		//local barrier
		int num_proc{1};		//number of processes to use
		
		Chunk* chunk{nullptr};
		Medium const* bg_medium{nullptr};
		
		
		
		std::vector<const Solid*> solids;
		std::vector<Inc_Source*> Inc_Sources;
		
		//nearfield probes, near to far field conversion, flux box
		std::vector<std::unique_ptr<Nearfield_Probe>> nearfield_probes;
		std::vector<std::unique_ptr<N2F_Box>> n2f_boxes;
		std::vector<std::unique_ptr<Flux_Box>> flux_boxes;

		//CPML members
		PML PMLs[3][2];
		std::vector<real> kx, ky, kz, bx, by, bz, cx, cy, cz;
		
		//medium factory
		std::vector<std::unique_ptr<Medium>> medium;
		std::vector<std::unique_ptr<Medium_Ref>> e_medium_ref;
		std::vector<std::unique_ptr<Medium_Ref>> m_medium_ref;
		std::vector<std::unique_ptr<Medium_Ref>> syn_medium_ref;
		std::mutex medium_mutex;
		
		// geometry factory
		std::vector<std::unique_ptr<Inhomogeneous_Box>> inhom_box_holder;
		std::vector<std::unique_ptr<Homogeneous_Object>> hom_box_holder;
		std::vector<std::unique_ptr<Primitive>> primitive_holder;
		
		//initialization
		void medium_init();
		void source_init();
		void PML_init();
		void PML_init_helper(const PML& neg, const PML& pos, real_arr& k, real_arr& b, real_arr& c, const int p1, const int p2);
		
	public:
		bool output_step_number{1};
		
		Simulation() = default;		//allow delayed initializations
		
		//function members to be used before chunk is initialized
		void add_inc_source(Inc_Source* source);
		void add_sf_layer(const int d);
		void add_tf_layer(const int d);
		void add_PML_layer(const PML& pml, const Direction dir, const Side side);
		void setup(const real _dx, const real _dt, const iVec3 _dim);
		void chunk_init();
		
		//functin members to be used that doesn't care whether chunk is initialized
		void add_solid(Solid const* solid);
		
		//function members to be used after chunk is initialized
		Nearfield_Probe const* add_nearfield_probe(const real freq, fVec3 pos);
		PML make_pml(const int d);
		Flux_Box const* add_flux_box(fVec3 p1, fVec3 p2, const real_arr& freq);
		N2F_Box const* add_n2f_box(fVec3 p1, fVec3 p2, const real_arr& freq);
		
		void add_dipole(const real amp, const fVec3 pos, const Coord_Type ctype, const std::function<real(const real)>& profile);
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
		
		void output_farfield(std::ostream& os);
		
		//medium factory
		template<typename... Args>
		Medium* make_medium(Args&&... args);
		
		Medium_Ref const* get_medium_ref(const Coord_Type ctype, const std::vector<real>& weights);
		std::vector<real> get_zero_weights();
		void prepare_medium(const real _dt);
		
		//geometry factory
		template<typename... Args>
		Geometry_Node make_sphere(Args&&... args);
		
		template<typename... Args>
		Geometry_Node make_box(Args&&... args);

		template<typename... Args>
		Geometry_Node make_disk(Args&&... args);
		
		//make inhomogeneous regions
		Solid const* make_solid(Medium const* m1, Medium const* m2, const std::string& filename);
		//make homogeneous regions
		Solid const* make_solid(Medium const* m, const Geometry_Node& geom);
	};
	
	template<typename... Args>
	Geometry_Node Simulation::make_sphere(Args&&... args) {
		primitive_holder.push_back(std::make_unique<Sphere>(std::forward<Args>(args)...));
		return Geometry_Node{primitive_holder.back().get()};
	}
	
	template<typename... Args>
	Geometry_Node Simulation::make_box(Args&&... args) {
		primitive_holder.push_back(std::make_unique<Box>(std::forward<Args>(args)...));
		return Geometry_Node{primitive_holder.back().get()};
	}

	template<typename... Args>
	Geometry_Node Simulation::make_disk(Args&&... args) {
		primitive_holder.push_back(std::make_unique<Disk>(std::forward<Args>(args)...));
		return Geometry_Node{ primitive_holder.back().get() };
	}
	
	template<typename... Args>
	Medium* Simulation::make_medium(Args&&... args) {
		size_t top = medium.size();
		
		medium.push_back(std::make_unique<Medium>(std::forward<Args>(args)...));
		medium.back()->set_index(top);
		return medium.back().get();
	}
}
