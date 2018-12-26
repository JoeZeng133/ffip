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
	class Probe;
	/* near to far field struct, resides in Simulation class*/
	class N2F_Face_Base {
	public:
		N2F_Face_Base() = default;
		virtual void update(Chunk* const chunk, const int n) = 0;
		virtual std::pair<Vec3<complex_num>, Vec3<complex_num>> get_NL(const real theta, const real phi, const real w) const = 0;
		virtual iVec3 get_norm_vec() const = 0;
		virtual iVec3 get_p1() const = 0;
		virtual iVec3 get_p2() const = 0;
		virtual void output_JM(std::ostream& os) const = 0;
	};
	
	/* D in {dir_x_tag, dir_y_tag, dir_z_tag}
	 local reference frame = {x1, x2, x3} where x3 = D and {x1, x2, x3} conforms right hand rule
	 */
	template<typename D>
	class N2F_Face : public N2F_Face_Base {
	private:
		iVec3 p1, p2;			//face lower and upper points in domain comp coordinates
		Side side;				//face side +-
		real_arr& omega;		//list of frequency (strictly ascending order)
		real dx, dt;			//dx, dt from Simulation
		fVec3 center;			//center in domain phys coordinates used in defining observation points
		Medium const* bg_medium;		//background medium
		
		mutable std::vector<complex_arr> j1_list, j2_list, m1_list, m2_list;	//{j1, j2, j3}, {m1, m2, m3} are surface currents where j3, m3 are always zeros
		std::vector<int> x1, x2;	//{x1, x2} are gridded points coordinates that cover the whole N2F face
		real_arr fx1, fx2;			//{fx1, fx2} are {x1, x2} converetd to domain phys coordinates and subtrated by center coordinates
		real_arr int_weights;		//integration weights for outpuing N,L in the end
		mutable bool is_phase_corrected{0}; //phase needs to be corrected for H at the first time get_NL is called
		void validity_check();
		void correct_phase() const;
		
	public:
		N2F_Face(iVec3 N2F_p1, iVec3 N2F_p2, const iVec3& _p1, const iVec3& _p2, const Side _side, real_arr& _omega, const real _dx, const real _dt, const fVec3& _center, Medium const* _bg_medium);
		
		void update(Chunk* const chunk, const int time_step) override;
		std::pair<Vec3<complex_num>, Vec3<complex_num>> get_NL(const real theta, const real phi, const real w) const override;
		iVec3 get_norm_vec() const override;
		
		iVec3 get_p1() const override {return rotate_frame(p1, typename D::z{});}
		iVec3 get_p2() const override {return rotate_frame(p2, typename D::z{});}

		//debugging function
		void output_JM(std::ostream& os) const override{
			os << std::scientific;
			for (int i = 0;  i < omega.size(); ++i) {
				auto& j1 = j1_list[i];
				auto& j2 = j2_list[i];
				auto& m1 = m1_list[i];
				auto& m2 = m2_list[i];

				for (int k = 0; k < j1.size(); ++k) {
					auto J = cVec3{ j1[k], j2[k], 0 };
					auto M = cVec3{ m1[k], m2[k], 0 };
					os << rotate_frame(J, typename D::z{}) << " " << rotate_frame(M, typename D::z{}) << "\n";
				}
			}
		}
	};


	
	
	template<typename D>
	iVec3 N2F_Face<D>::get_norm_vec() const {
		return vec3_base[D::val] * (int)side;
	}
	
	template<typename D>
	void N2F_Face<D>::validity_check() {
		if(!ElementWise_Less_Eq(p1, p2))
			throw std::runtime_error("Invalid lower and upper points for N2F face");
		
		if(ElementWise_Less(p1, p2))
			throw std::runtime_error("p1 p2 specifies a box not a face");
		
		if(p1.z != p2.z)
			throw std::runtime_error("Face direction is wrong");
		
		if(p1.get_type() != p2.get_type())
			throw std::runtime_error("Two corners should have the same ctype");
		
		for(int i = 1; i < omega.size(); ++i)
			if(omega[i] <= omega[i - 1])
				throw std::runtime_error("Frequency lists should be strictly monotone");
	}
	
	template<typename D>
	N2F_Face<D>::N2F_Face(iVec3 N2F_p1, iVec3 N2F_p2, const iVec3& _p1, const iVec3& _p2, const Side _side, real_arr& _omega, const real _dx, const real _dt, const fVec3& _center, Medium const* _bg_medium): p1(_p1), p2(_p2), side(_side), omega(_omega), dx(_dx), dt(_dt), center(_center), bg_medium(_bg_medium) {
		//change frame to local coordinates
		N2F_p1 = rotate_frame(N2F_p1, D{});
		N2F_p2 = rotate_frame(N2F_p2, D{});
		p1 = rotate_frame(p1, D{});
		p2 = rotate_frame(p2, D{});
		center = rotate_frame(center, D{});
		
		validity_check();
		
		for(int i = p1.x; i <= p2.x; i += 2) {
			x1.push_back(i);
			fx1.push_back(i * dx / 2 - center.x);
		}

		for(int j = p1.y; j <= p2.y; j += 2) {
			x2.push_back(j);
			fx2.push_back(j * dx / 2 - center.y);
		}

		for(int j = p1.y; j <= p2.y; j += 2)
			for(int i = p1.x; i <= p2.x; i += 2) {
				real side1 = (i == N2F_p1.x || i == N2F_p2.x)? 0.5 : 1;
				real side2 = (j == N2F_p1.y || j == N2F_p2.y)? 0.5 : 1;
				int_weights.push_back(side1 * side2 * dx * dx);
		}
		
		j1_list.resize(omega.size());
		j2_list.resize(omega.size());
		m1_list.resize(omega.size());
		m2_list.resize(omega.size());
		
		size_t N = x1.size() * x2.size();
		
		auto f_resize = [=](complex_arr& arr){
			arr.resize(N);
		};
		
		for_each(j1_list.begin(), j1_list.end(), f_resize);
		for_each(j2_list.begin(), j2_list.end(), f_resize);
		for_each(m1_list.begin(), m1_list.end(), f_resize);
		for_each(m2_list.begin(), m2_list.end(), f_resize);
	}
	
	template<typename D>
	void N2F_Face<D>::update(Chunk* const chunk, const int n) {
		int index, k = p1.z;
		/*int side_int = static_cast<int>(side);*/

		for(int f = 0; f < omega.size(); ++f) {
			complex_num exp_omega_n(cos(omega[f] * n * dt), -sin(omega[f] * n * dt));
			index = 0;
			auto& j1 = j1_list[f];
			auto& j2 = j2_list[f];
			auto& m1 = m1_list[f];
			auto& m2 = m2_list[f];
			
			for(auto j : x2)
				for(auto i : x1){
					iVec3 pos(i, j, k);
					pos = rotate_frame(pos, typename D::z{});	//get fields using the coordinates in the original ref frame
					
					/* J = (0, 0, side) cross H, M = (0, 0, -side) cross E */
					j1[index] += -side * (*chunk)(pos, dir_traits<D>::x2::H) * exp_omega_n;
					j2[index] += side * (*chunk)(pos, dir_traits<D>::x1::H) * exp_omega_n;
					
					m1[index] += side * (*chunk)(pos, dir_traits<D>::x2::E) * exp_omega_n;
					m2[index] += -side * (*chunk)(pos, dir_traits<D>::x1::E) * exp_omega_n;
					
					++index;
				}
		}
	}
	
	template<typename D>
	void N2F_Face<D>::correct_phase() const{
		for(int f = 0; f < omega.size(); ++f) {
			complex_num phase_correction{cos(omega[f] * dt * 0.5), sin(omega[f] * dt * 0.5)};
			for(auto& x : j1_list[f]) {
				x *= phase_correction;
			}
			
			for(auto& x : j2_list[f]) {
				x *= phase_correction;
			}
		}
		
		is_phase_corrected = 1;
	}
	
	template<typename D>
	std::pair<Vec3<complex_num>, Vec3<complex_num>> N2F_Face<D>::get_NL(const real theta, const real phi, const real w)  const{
		if(!is_phase_corrected) correct_phase();
		
		auto low = std::lower_bound(omega.begin(), omega.end(), w);
		if (low == omega.end() || *low != w)
			throw std::runtime_error("Invalid requested frequency");
		
		size_t f = low - omega.begin();
		auto& j1 = j1_list[f];
		auto& j2 = j2_list[f];
		auto& m1 = m1_list[f];
		auto& m2 = m2_list[f];
		
		fVec3 op{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
		op = rotate_frame(op, D{});		//observation direction
		
		Vec3<complex_num> N{0, 0, 0};
		Vec3<complex_num> L{0, 0, 0};
		
		real z = p1.z * dx / 2 - center.z;	//z
		real k = omega[f] / bg_medium->get_c();		//wave number assumed
		int index = 0;
		
		for(auto y : fx2)
			for(auto x : fx1) {
				real tmp = (op.x * x + op.y * y + op.z * z) * k;
				auto phase = complex_num{cos(tmp), sin(tmp)};
				
				N.x += j1[index] * phase * int_weights[index];
				N.y += j2[index] * phase * int_weights[index];
				L.x += m1[index] * phase * int_weights[index];
				L.y += m2[index] * phase * int_weights[index];
				++index;
			}
		
		return {rotate_frame(N, typename D::z{}), rotate_frame(L, typename D::z{})};
	}
	
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
		std::vector<Probe*> probes;
		
		std::vector<real> N2F_omega;
		std::vector<real> N2F_omega_unique;
		std::vector<fVec3> N2F_pos;
		std::vector<N2F_Face_Base*> N2F_faces;
		
		iVec3 tf_p1, tf_p2;
		iVec3 N2F_p1, N2F_p2;	//coordinates of N2F box

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
		void advance(std::ostream& os);
		void add_solid(Solid const* solid);
		void add_source(Source* source);
		void add_PML_layer(PML* PML);
		void add_probe(Probe* probe);
		void add_farfield_probe(const real freq, const fVec3& pos);
		void set_background_medium(Medium const* m);
		void set_num_proc(const int _num_proc);
		void init();
		
		void udf_unit();
		void udf_advance();
		void udf_output();

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
		
		void output(std::ostream& os);
		void output_farfield(std::ostream& os);
	};
}
