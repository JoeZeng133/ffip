#pragma once

#include <utility.hpp>
#include <simulation.hpp>
#include <iostream>

namespace ffip {
	class Probe {
	public:
		virtual ~Probe() {}
		virtual void update(const int time_step) = 0;
		virtual void output(std::ostream&) = 0;
	};
	
	class Probe_Frequency : public Probe{
	private:
		real omega;
		fVec3 pos;
		Chunk const* chunk;
		std::atomic<bool> phase_corrected{0};
	public:
		complex_num ex{0}, ey{0}, ez{0}, hx{0}, hy{0}, hz{0};
		
		Probe_Frequency() = delete;
		Probe_Frequency(const real freq, const fVec3& _pos, Chunk const* chunk);
		void update(const int time_step) override;
		void output(std::ostream&) override;
	};
	
	/* Used internally inside simulation
	 	fourier transforms of a particualr field for a box region
	 	field inside the box region can be accessed once simulation is done
	 */
	class Box_Freq_Req {
	private:
		//constructor needed
		fVec3 p1, p2;								//corner points for the requested region
		Coord_Type F;								//field coord type
		const real_arr& omega_list;					//list of omegas (frequencies)
		Chunk const* chunk{nullptr};
		
		//other private data memebrs
		real dt, dx;
		iVec3 near_p1, near_p2;
		sVec3 dim;									//dimension of the array
		interpn<3> interp;
		
		std::vector<complex_arr> raw_fields;		//fourier transform of raw fields from chunk
		std::atomic<bool> phase_corrected{0};		//phase correction for H field
	public:
		Box_Freq_Req() = default;
		Box_Freq_Req(const fVec3& p1, const fVec3& p2, const Coord_Type F, const real_arr& omega_list, Chunk const* chunk);
		
		Box_Freq_Req(const Box_Freq_Req&) = default;
		Box_Freq_Req& operator=(const Box_Freq_Req&) = default;
		Box_Freq_Req(Box_Freq_Req&&) = default;
		Box_Freq_Req& operator=(Box_Freq_Req&&) = default;
		
		void validity_check() const;
		void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1);
		void correct_phase();						//thread safe function that performs phase correction once
		
		complex_arr operator()(const fVec3& p) const;
//		complex_arr operator()(const iVec3& p) const;
	};
	
	/* N2F region class*/
	class N2F_Face_Base {
	public:
		virtual void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1) = 0;
		virtual std::pair<cVec3, cVec3> get_NL(const real theta, const real phi, const int f_index) const = 0;
		virtual ~N2F_Face_Base() {};
		virtual void prepare(const size_t rank = 0, const size_t num_proc = 1) = 0;
	};
	
	template<typename Dir>
	class N2F_Face : public N2F_Face_Base {
		using x1 = typename Dir::x1;
		using x2 = typename Dir::x2;
		using x1_a = typename Dir::x1_a;
		using x2_a = typename Dir::x2_a;
		
	private:
		fVec3 p1, p2, ref_p;			//corner points specified in computational coordinates
		Side side;
		const real_arr& omega_list;
		Chunk const* chunk{nullptr};
		real c;
		
		Box_Freq_Req* e1_req, *e2_req, *h1_req, *h2_req;
		std::vector<complex_arr> j1_list, j2_list, m1_list, m2_list;
		std::atomic<bool> prepared{0};
		
		real lx1, lx2;
		size_t nx1, nx2;
		real dx1, dx2;
		
		void validity_check() const;
	public:
		
		
		N2F_Face(const std::pair<fVec3, fVec3>& corners, const fVec3& ref_p, const Side side, const real_arr& omega_list, Chunk const* chunk, const real c);
		N2F_Face(const fVec3& p1, const fVec3& p2, const fVec3& ref_p, const Side side, const real_arr& freq, Chunk const* chunk, const real c);
		~N2F_Face();
		
		void prepare(const size_t rank = 0, const size_t num_proc = 1) override;
		void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1) override;
		std::pair<cVec3, cVec3> get_NL(const real theta, const real phi, const int f_index) const override;
	};
	
	template<typename Dir>
	N2F_Face<Dir>::~N2F_Face() {
		delete e1_req;
		delete e2_req;
		delete h1_req;
		delete h2_req;
	}
	
	template<typename Dir>
	N2F_Face<Dir>::N2F_Face(const std::pair<fVec3, fVec3>& corners, const fVec3& ref_p, const Side side, const real_arr& omega_list, Chunk const* chunk, const real c): N2F_Face(corners.first, corners.second, ref_p, side, omega_list, chunk, c) {}
	
	template<typename Dir>
	N2F_Face<Dir>::N2F_Face(const fVec3& _p1, const fVec3& _p2, const fVec3& _ref_p, const Side _side, const real_arr& _omega_list, Chunk const* _chunk, const real _c): p1(_p1), p2(_p2), ref_p(_ref_p),  side(_side), omega_list(_omega_list), chunk(_chunk), c(_c) {
		
		validity_check();
		fVec3 dp = p2 - p1;
		lx1 = choose<x1>::get(dp);
		lx2 = choose<x2>::get(dp);
		nx1 = (size_t)round(lx1 / 2);
		nx2 = (size_t)round(lx2 / 2);
		dx1 = lx1 / nx1;
		dx2 = lx2 / nx2;
		
		auto init_func = [&, this](decltype(j1_list)& list) {
			list.resize(omega_list.size());
			for(auto& item : list)
				item.resize((nx1 + 1) * (nx2 + 1));
		};
		
		init_func(j1_list);
		init_func(j2_list);
		init_func(m1_list);
		init_func(m2_list);
		
		e1_req = new Box_Freq_Req(p1, p2, x1::E, omega_list, chunk);
		e2_req = new Box_Freq_Req(p1, p2, x2::E, omega_list, chunk);
		h1_req = new Box_Freq_Req(p1, p2, x1::H, omega_list, chunk);
		h2_req = new Box_Freq_Req(p1, p2, x2::H, omega_list, chunk);
	}
	
	template<typename Dir>
	void N2F_Face<Dir>::validity_check() const {
		if (choose<Dir>::get(p1) != choose<Dir>::get(p2))
			throw Invalid_Direction{};
		
		if (!ElementWise_Less_Eq(p1, p2))
			throw Invalid_Corner_Points{};
	}
	
	template<typename Dir>
	void N2F_Face<Dir>::update(const int time_step, const size_t rank, const size_t num_proc) {
		e1_req->update(time_step, rank, num_proc);
		e2_req->update(time_step, rank, num_proc);
		h1_req->update(time_step, rank, num_proc);
		h2_req->update(time_step, rank, num_proc);
	}
	
	template<typename Dir>
	void N2F_Face<Dir>::prepare(const size_t rank, const size_t num_proc){
		real side_float = side;
		h1_req->correct_phase();
		h2_req->correct_phase();
		
		for(auto itr = my_iterator(iVec3(0, 0, 0), iVec3(nx1, nx2, 0), Null, rank, num_proc); !itr.is_end(); itr.advance()) {
			fVec3 sample_point = p1;
			choose<x1>::get(sample_point) += itr.x * dx1;
			choose<x2>::get(sample_point) += itr.y * dx2;
			
			auto e1 = (*e1_req)(sample_point);
			auto e2 = (*e2_req)(sample_point);
			auto h1 = (*h1_req)(sample_point);
			auto h2 = (*h2_req)(sample_point);
			
			for(int k = 0; k < omega_list.size(); ++k) {
				j1_list[k][itr.index] = -side_float * h2[k];
				j2_list[k][itr.index] = side_float * h1[k];
				m1_list[k][itr.index] = side_float * e2[k];
				m2_list[k][itr.index] = -side_float * e1[k];
			}
		}
	}

	template<typename Dir>
	std::pair<cVec3, cVec3> N2F_Face<Dir>::get_NL(const real theta, const real phi, const int f_index) const {
		
		real dx = chunk->get_dx();
		fVec3 op{sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};
		fVec3 tmp = p1 - ref_p;
		cVec3 N{0, 0, 0};
		cVec3 L{0, 0, 0};
		real k = omega_list[f_index] / c;
		auto& j1 = j1_list[f_index];
		auto& j2 = j2_list[f_index];
		auto& m1 = m1_list[f_index];
		auto& m2 = m2_list[f_index];
		
		for(auto itr = my_iterator(iVec3(0, 0, 0), iVec3(nx1, nx2, 0), Null); !itr.is_end(); itr.advance()) {
			
			fVec3 sample_point = tmp;
			choose<x1>::get(sample_point) += itr.x * dx1;
			choose<x2>::get(sample_point) += itr.y * dx2;
			
			real phase = (op.x * sample_point.x + op.y * sample_point.y + op.z * sample_point.z) * k * dx / 2;
			auto exp_phase = complex_num(cos(phase), sin(phase));
			real trapz_weight = (itr.x == 0 || itr.x == nx1 ? 0.5 : 1) * (itr.y == 0 || itr.y == nx2 ? 0.5 : 1);
			
			//trapezoidal weighted sum
			choose<x1>::get(N) += j1[itr.index] * exp_phase * trapz_weight;
			choose<x2>::get(N) += j2[itr.index] * exp_phase * trapz_weight;
			choose<x1>::get(L) += m1[itr.index] * exp_phase * trapz_weight;
			choose<x2>::get(L) += m2[itr.index] * exp_phase * trapz_weight;
		}
		
		//scale to real integration
		real dS = dx * dx * (dx1 / 2) * (dx2 / 2);
		N = N * dS;
		L = L * dS;
		return {N, L};
	}
	
	/* N2F Box class*/
	class N2F_Box {
	private:
		fVec3 p1, p2;
		real_arr omega_list;
		Chunk* const chunk {nullptr};
		real c;
		std::vector<N2F_Face_Base*> n2f_faces;
		
	public:
		N2F_Box(const fVec3& p1, const fVec3& p2, const real_arr& omega, Chunk* const chunk, const real c);
		~N2F_Box();
		
		N2F_Box(const N2F_Box&) = delete;
		N2F_Box& operator=(const N2F_Box&) = delete;
		N2F_Box(N2F_Box&&) = default;
		N2F_Box& operator=(N2F_Box&&) = default;
		
		void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1);
		void prepare(const size_t rank = 0, const size_t num_proc = 1);
		std::pair<cVec3, cVec3> get_NL(const real theta, const real phi, const real omega);
	};
	
	/* Flux region class*/
//	class Flux_Face_Base {
//	public:
//		virtual void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1) = 0;
//		virtual complex_arr get() const = 0;
//		virtual ~Flux_Face_Base() {};
//	};
//
//	template<typename Dir>
//	class Flux_Face : public Flux_Face_Base {
//		using x1 = typename Dir::x1;
//		using x2 = typename Dir::x2;
//		using x1_a = typename Dir::x1_a;
//		using x2_a = typename Dir::x2_a;
//
//	private:
//		fVec3 p1, p2;
//		Side side;
//		const real_arr& freq;
//		Chunk const* chunk{nullptr};
//		real dx, dt;
//
//		Box_Freq_Req* e1_req, *e2_req, *h1_req, *h2_req;
//
//		void validity_check() const;
//	public:
//		Flux_Face(const std::pair<fVec3, fVec3>& corners, const Side side, const real_arr& freq, Chunk const* chunk);
//		Flux_Face(const fVec3& p1, const fVec3& p2, const Side side, const real_arr& freq, Chunk const* chunk);
//		void update(const int time_step, const size_t rank = 0, const size_t num_proc = 1) override;
//		complex_arr get() const override;
//	};
//
//	template<typename Dir>
//	void Flux_Face<Dir>::validity_check() const {
//		if (choose<Dir>::get(p1) != choose<Dir>::get(p2))
//			throw Invalid_Direction{};
//
//		if (!ElementWise_Less_Eq(p1, p2))
//			throw Invalid_Corner_Points{};
//	}
//
//	template<typename Dir>
//	complex_arr Flux_Face<Dir>::get() const {
//		h1_req->correct_phase();
//		h2_req->correct_phase();
//
//		auto dp = p2 - p1;
//		real lx1 = choose<x1>::get(dp);
//		real lx2 = choose<x2>::get(dp);
//		size_t nx1 = (int)ceil(lx1 / dx) + 1;
//		size_t nx2 = (int)ceil(lx2 / dx) + 1;
//		real dx1 = lx1 / nx1;
//		real dx2 = lx2 / nx2;
//
//		complex_arr res(freq.size(), 0);
//		for(int j = 0; j <= nx2; ++j)
//			for(int i = 0; i <= nx1; ++i) {
//				auto sample_point = dp;
//				choose<x1>::get(sample_point) += i * dx1;
//				choose<x2>::get(sample_point) += j * dx2;
//
//				auto e1 = (*e1_req)(sample_point);
//				auto e2 = (*e2_req)(sample_point);
//				auto h1 = (*h1_req)(sample_point);
//				auto h2 = (*h2_req)(sample_point);
//
//				real modifier = ((i == 0) ^ (i == nx1) ? 0.5 : 1) * ((j == 0) ^ (j == nx2) ? 0.5 : 1);
//				for(int f = 0; f < freq.size(); ++f) {
//					res[f] += int(side) * modifier * (e1[f] * std::conj(h2[f]) - e2[f] * std::conj(h1[f])) ;
//				}
//			}
//		for(auto& item : res)
//			item *= dx1 * dx2;
//
//		return res;
//	}
//
//	template<typename Dir>
//	Flux_Face<Dir>::Flux_Face(const std::pair<fVec3, fVec3>& corners, const Side side, const real_arr& freq, Chunk const* chunk): Flux_Face(corners.first, corners.second, side, freq, chunk) {}
//
//	template<typename Dir>
//	Flux_Face<Dir>::Flux_Face(const fVec3& _p1, const fVec3& _p2, const Side _side, const real_arr& _freq, Chunk const* _chunk):p1(_p1), p2(_p2), side(_side), freq(_freq), chunk(_chunk) {
//		validity_check();
//		dx = chunk->get_dx();
//		dt = chunk->get_dt();
//		e1_req = new Box_Freq_Req(p1, p2, Dir::x1::E, freq, chunk);
//		e2_req = new Box_Freq_Req(p1, p2, Dir::x2::E, freq, chunk);
//		h1_req = new Box_Freq_Req(p1, p2, Dir::x1::H, freq, chunk);
//		h2_req = new Box_Freq_Req(p1, p2, Dir::x2::H, freq, chunk);
//	}
//
//	template<typename Dir>
//	void Flux_Face<Dir>::update(const int time_step, const size_t rank, const size_t num_proc) {
//		e1_req->update(time_step, rank, num_proc);
//		e2_req->update(time_step, rank, num_proc);
//		h1_req->update(time_step, rank, num_proc);
//		h2_req->update(time_step, rank, num_proc);
//	}
//
//	class Flux_Box {
//	private:
//		fVec3 p1, p2;
//		real_arr omega;
//		Chunk* const chunk {nullptr};
//		std::vector<Flux_Face_Base*> flux_faces;
//
//	public:
//		Flux_Box(const fVec3& p1, const fVec3& p2, const real_arr& omega, Chunk* const chunk);
//		~Flux_Box();
//
//		Flux_Box(const Flux_Box&) = delete;
//		Flux_Box& operator=(const Flux_Box&) = delete;
//		Flux_Box(Flux_Box&&) = default;
//		Flux_Box& operator=(Flux_Box&&) = default;
//		void update(const int time_step, const size_t rank, const size_t num_proc);
//		complex_arr get() const;
//	};
}
