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
	protected:
		
		iVec3 p1, p2;		//face lower and upper points in domain comp coordinates
		Side side;			//face side +-
		real omega;			//frequency in rad/s
		real dx, dt;		//dx, dt from Simulation
		fVec3 center;		//center in domain phys coordinates used in defining observation points
		
		mutable complex_arr j1, j2, m1, m2;	//{j1, j2, j3}, {m1, m2, m3} are surface currents where j3, m3 are always zeros
		std::vector<int> x1, x2;	//{x1, x2} are gridded points coordinates that cover the whole N2F face
		real_arr fx1, fx2;			//{fx1, fx2} are {x1, x2} converetd to domain phys coordinates and subtrated by center coordinates
		mutable bool is_phase_corrected{0}; //phase needs to be corrected for H at the first time get_NL is called
		bool is_face{ 0 };			//if it is not a face, ignore it and always output 0 to get_NL

	public:
		N2F_Face_Base(const iVec3& _p1, const iVec3& _p2, const Side _side, const real _omega, const real _dx, const real _dt, const fVec3& _center);
		virtual void update(Chunk* const chunk, const int n) = 0;
		virtual std::pair<Vec3<complex_num>, Vec3<complex_num>> get_NL(const real theta, const real phi) const = 0;
	};
	
	/* D in {dir_x_tag, dir_y_tag, dir_z_tag}
	 local reference frame = {x1, x2, x3} where x3 = D and {x1, x2, x3} conforms right hand rule
	 */
	template<typename D>
	class N2F_Face : public N2F_Face_Base {
	private:
	public:
		N2F_Face(const iVec3& _p1, const iVec3& _p2, const Side _side, const real _omega, const real _dx, const real _dt, const fVec3& _center);
		void update(Chunk* const chunk, const int n) override;
		std::pair<Vec3<complex_num>, Vec3<complex_num>> get_NL(const real theta, const real phi) const override;
	};
	
	template<typename D>
	N2F_Face<D>::N2F_Face(const iVec3& _p1, const iVec3& _p2, const Side _side, const real _omega, const real _dx, const real _dt, const fVec3& _center): N2F_Face_Base(_p1, _p2, _side, _omega, _dx, _dt, _center) {
		//change frame to local coordinates
		p1 = rotate_frame(p1, D{});
		p2 = rotate_frame(p2, D{});
		center = rotate_frame(center, D{});
		
		if(!ElementWise_Less_Eq(p1, p2))
			throw std::runtime_error("Invalid lower and upper points for N2F face");
		
		if(ElementWise_Less(p1, p2))
			throw std::runtime_error("p1 p2 specifies a box not a face");
		
		if(p1.z != p2.z)
			throw std::runtime_error("Face direction is wrong");
		
		//add extra point to cover the entire region [p1, p2]
		for(int i = p1.x; i <= p2.x; i += 2)
			x1.push_back(i);
		if((p2.x - p1.x) & 1)
			x1.push_back(p2.x);
		
		for(int i = p1.y; i <= p2.y; i += 2)
			x2.push_back(i);
		if((p2.y - p2.y) & 1)
			x2.push_back(p2.y);
		
		fx1.resize(x1.size());
		fx2.resize(x2.size());
		
		for(auto i : x1)
			fx1.push_back(i * dx / 2 - center.x);
		for(auto j : x2)
			fx2.push_back(j * dx / 2 - center.y);
		
		
		j1.resize(x1.size() * x2.size(), 0);
		j2.resize(j1.size(), 0);
		m1.resize(j1.size(), 0);
		m2.resize(j1.size(), 0);
	}
	
	template<typename D>
	void N2F_Face<D>::update(Chunk* const chunk, const int n) {
		if (!is_face)
			return;

		int index, k = p1.z;
		complex_num exp_omega_n;
	
		exp_omega_n = complex_num(cos(omega * n * dt), -sin(omega * n * dt));
		index = 0;
		
		for(auto i : x1)
			for(auto j : x2) {
				iVec3 pos{i, j, k};
				rotate_frame(pos, typename D::z{});	//get fields using the coordinates in the original ref frame
				
				/* J = (0, 0, side) cross H, M = (0, 0, -side) cross E */
				j1[index] += -side * (*chunk)(pos, dir_traits<D>::x2::H) * exp_omega_n;
				j2[index] += side * (*chunk)(pos, dir_traits<D>::x1::H) * exp_omega_n;
				m1[index] += side * (*chunk)(pos, dir_traits<D>::x2::E) * exp_omega_n;
				m2[index] += -side * (*chunk)(pos, dir_traits<D>::x1::E) * exp_omega_n;
				
				++index;
			}
	}
	
	template<typename D>
	std::pair<Vec3<complex_num>, Vec3<complex_num>> N2F_Face<D>::get_NL(const real theta, const real phi)  const{
		if (!is_face)
			return { {0, 0, 0}, {0, 0, 0} };

		if(!is_phase_corrected) {
			auto phase_correction = exp(complex_num(0, omega * dt / 2));
			
			for(auto& x : m1) {
				x = x * phase_correction;
			}
			
			for(auto& x : m2) {
				x = x * phase_correction;
			}
			is_phase_corrected = 1;
		}
		
		Vec3 op{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};							//observation point
		/* integrand of j1, j2, m1, m2 at each point
		 j_int = j * exp(j * k * (r' dot r0));
		 */
		complex_arr j_int_1(j1.size()), j_int_2(j1.size()), m_int_1(j1.size()), m_int_2(j1.size());
		Vec3<complex_num> N, L;
		
		real z = p1.z * dx / 2 - center.z;	//z
		real beta = omega / ffip::c0;		//wave number
		int index = 0;
		
		op = rotate_frame(op, D{});
		for(auto x : fx1)
			for(auto y : fx2) {
				auto phase = exp(complex_num(0, (op.x * x + op.y * y + op.z * z) * beta));
				
				j_int_1[index] = j1[index] * phase;
				j_int_2[index] = j2[index] * phase;
				m_int_1[index] = m1[index] * phase;
				m_int_2[index] = m2[index] * phase;
				++index;
			}
		
		N.x = integral_ndim(j_int_1.data(), fx1, fx2);
		N.y = integral_ndim(j_int_2.data(), fx1, fx2);
		N.z = 0;
		
		L.x = integral_ndim(m_int_1.data(), fx1, fx2);
		L.y = integral_ndim(m_int_2.data(), fx1, fx2);
		L.z = 0;
		
		return {N, L};
	}
	
	class Simulation {
	private:

		real dt, dx;
		iVec3 sim_dim;
		
		int step;
		iVec3 sim_p1, sim_p2;			//domain coordinates of lower, upper points
		iVec3 ch_p1, ch_p2;		//domain coordinates of lower, upper points
		
		Chunk* chunk;
		int background_medium_id{-1};
		
		std::vector<Solid*> solids;
		std::vector<Current_Source*> current_sources;
		std::vector<Eigen_Source*> eigen_sources;
		std::vector<Probe*> probes;
		
		std::vector<real> N2F_freq;
		std::vector<fVec3> N2F_pos;
		std::map<real, std::vector<N2F_Face_Base*>> N2F_faces;

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
		Simulation(const real _dx, const real _dt, const iVec3 _dim);
		void advance(std::ostream& os);
		void add_solid(Solid* solid);
		void add_source(Source* source);
		void add_PML_layer(PML* PML);
		void add_probe(Probe* probe);
		void add_farfield_probe(const real freq, const fVec3& p);
		void set_background_medium(const int id);
		void init();

		/* field access functions*/
		real at(const fVec3& p, const Coord_Type ctype) const;					//access at float physical coordinates

		/* field access at computation coordinates*/
		real operator()(const fVec3& p, const Coord_Type ctype) const;		//access at float computation coordinates
		real operator()(const iVec3& p, const Coord_Type ctype) const;		//access at integer computation coordinates
		real operator()(const iVec3& p) const;								//raw access	

		/* */
		int get_step() const;
		real get_dt() const;
		void output(std::ostream&);
	};
}
