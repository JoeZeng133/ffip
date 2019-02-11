#include <utility.hpp>

using namespace std;

namespace ffip {
	const real e0 = 8.854187817e-12;
	const real u0 = 1.2566370614e-6;
	const real pi = 3.141592653589793e+00;
	const real z0 = 376.73031346;
	const real c0 = 3e8;
	
	const int side_low_tag::val = -1;
	
	int side_low_tag::get_next_int(int x) {
		return x - 1;
	}
	
	int side_low_tag::round(real x) {
		return int(std::floor(x));
	}
	
	const int side_high_tag::val = 1;
	
	int side_high_tag::get_next_int(int x) {
		return x + 1;
	}
	
	int side_high_tag::round(real x) {
		return int(std::ceil(x));
	}
	
	const int odd_tag::val = 1;
	const int even_tag::val = 0;
	
	const int dir_x_tag::val = 0;
	const Coord_Type dir_x_tag::E = Ex;
	const Coord_Type dir_x_tag::H = Hx;
	
	const int dir_y_tag::val = 1;
	const Coord_Type dir_y_tag::E = Ey;
	const Coord_Type dir_y_tag::H = Hy;
	
	const int dir_z_tag::val = 2;
	const Coord_Type dir_z_tag::E = Ez;
	const Coord_Type dir_z_tag::H = Hz;
	
	const Coord_Type ex_tag::ctype = Ex;
	const Coord_Type ey_tag::ctype = Ey;
	const Coord_Type ez_tag::ctype = Ez;
	const Coord_Type hx_tag::ctype = Hx;
	const Coord_Type hy_tag::ctype = Hy;
	const Coord_Type hz_tag::ctype = Hz;
	
	
	 
	auto make_gaussian_func(real sigma_t, real d) -> std::function<real(const real)> {
		real a = 0.5 / (sigma_t * sigma_t);
		
		return [=](const real time) -> real {
			return exp(-a * (time - d) * (time - d));
		};
	}
	
	auto make_sin_func(real freq, real d) -> std::function<real(const real)> {
		real a = 2 * pi * freq;
		
		return [=](const real time) -> real {
			return sin(a* (time - d));
		};
	}
	
	auto make_ricker_func(real fp, real d) -> std::function<real(const real)> {
		real a = (pi * fp) * (pi * fp);
		
		return [=](const real time) -> real {
			real arg = a * (time - d) * (time - d);
			return (1 - 2 * arg) * exp(-arg);
		};
	}
	
	/* Gridded Interpolation class */
	GriddedInterp::GriddedInterp(string _file_name):file_name(_file_name) {
		fstream fin{file_name, ios::in};
		
		if (!fin.is_open()) {
			throw std::runtime_error("Cannot open the interpolation file");
		}
		
		fin >> dimx >> dimy >> dimz;
		
		fin >> x0 >> dx;
		fin >> y0 >> dy;
		fin >> z0 >> dz;
		
		//pad interpolation with zeros
		{
			v.resize(size_t(dimx + 2) * (dimy + 2) * (dimz + 2));
			auto strides = dim2stride(iVec3(dimx + 2, dimy + 2, dimz + 2));

			for (auto itr = my_iterator(iVec3{ 1, 1, 1 }, iVec3{ dimx, dimy, dimz }, All); !itr.is_end(); itr.advance()) {
				size_t index = inner_prod(itr.get_vec(), strides);
				fin >> v[index];
			}

			dimx += 2;	x0 -= dx;
			dimy += 2;	y0 -= dy;
			dimz += 2;	z0 -= dz;
		}
		
		//no padding
		/*{
			v.resize(size_t(dimx) * dimy * dimz);
			for (auto itr = my_iterator(iVec3{ 0, 0, 0 }, iVec3{ dimx - 1, dimy - 1, dimz - 1 }, All); !itr.is_end(); itr.advance()) {
				fin >> v[itr.index];
			}
		}*/

		interp = interp3(dimx, dimy, dimz);
	}
	
	/*GriddedInterp::GriddedInterp(const iVec3& dim, const fVec3& w0, const fVec3& dw, const real_arr& _v):
	v(_v), dx(dw.x), dy(dw.y), dz(dw.z) {
		
		for(int i = 0; i < dim.x; ++i)
			x.push_back(w0.x + dw.x * i);
		
		for(int j = 0; j < dim.y; ++j)
			y.push_back(w0.y + dw.y * j);
		
		for(int k = 0; k < dim.z; ++k)
			z.push_back(w0.z + dw.z * k);
		
		interp = interpn<3>(z.size(), y.size(), x.size());
	}*/

	real GriddedInterp::get(const fVec3& p) const{
		return interp.get(v, (p.z - z0) / dz, (p.y - y0) / dy, (p.x - x0) / dx);
	}
	
	fVec3 GriddedInterp::get_p1() const{
		return {x0, y0, z0};
	}
	
	fVec3 GriddedInterp::get_p2() const{
		return {x0 + (dimx - 1) * dx, y0 + (dimy - 1) * dy, z0 + (dimz - 1) * dz};
	}
	
	/* PML */
	PML::PML(int d, real sigma_max, real k_max, real a_max, int m, int m_a):d(d), sigma_max(sigma_max), k_max(k_max), a_max(a_max), m(m), m_a(m_a) {};
	
	int PML::get_d() const {
		return d;
	}
	
	real PML::get_sigma(const real x) const {
		return pow(x / d, m) * sigma_max;
	}
	
	real PML::get_k(const real x) const {
		return 1 + (k_max - 1) * pow(x / d, m);
	}
	
	real PML::get_a(const real x) const {
		return a_max * pow(1 - x / d, m_a);
	}
	
	real PML::get_b(const real x, const real dt) const {
		return exp(-(get_sigma(x) / get_k(x) + get_a(x)) * dt / e0);
	}
	
	real PML::get_c(const real x, const real dt) const {
		return get_sigma(x) * (get_b(x, dt) - 1) / (get_sigma(x) * get_k(x) + get_k(x) * get_k(x) * get_a(x));
	}
	
	real PML::optimal_sigma_max(real m_k, real dx, real er, real ur) {
		return 0.8 * (m_k + 1) / (z0 * dx * sqrt(er * ur));
	}
	
	Dispersive_Field::Dispersive_Field(const size_t num_poles) {
		eh2 = 0;
		p.resize(num_poles, 0);
		p1.resize(num_poles, 0);
	}
	
	size_t Dispersive_Field::get_num_poles() const {
		return p1.size();
	}
		
	/* ############################################### */
	my_iterator::my_iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype): my_iterator{corners.first, corners.second, ctype} {}
	
	my_iterator::my_iterator(const fVec3& p1, const fVec3& p2, const Coord_Type ctype): my_iterator{get_component_interior(p1, p2, ctype), ctype} {}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype) {
		std::pair<iVec3, iVec3> tmp;
		
		// loop through points inside [p1, p2]
		if (ctype == All) {
			tmp = get_component_interior(p1, p2, ctype);
			stride = 1;
		}
		else if (ctype == None) {
			throw std::runtime_error("Invalid Ctype");
		}
		else {
			tmp = get_component_interior(p1, p2, ctype);
			stride = 2;
		}
		
		x = x0 = tmp.first.x;
		y = y0 = tmp.first.y;
		z = z0 = tmp.first.z;
		x1 = tmp.second.x;
		y1 = tmp.second.y;
		z1 = tmp.second.z;
		
		index = 0;
		end = size = get_size();
	}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype, const size_t rank, const size_t num): my_iterator(p1, p2, ctype) {
		if (rank >= num)
			throw std::runtime_error("Rank is larger than number of threads");
		
		if (size <= num)
			throw std::runtime_error("Region is too small to divide, Better use non-parallel version");
		
		size_t idx1 = (rank * size) / num;
		size_t idx2 = ((rank + 1) * size) / num;
		
		index = idx1;
		end = idx2;
		
		sVec3 dim = get_dim();
		x = x0 + (idx1 % dim.x) * stride;
		idx1 /= dim.x;
		y = y0 + (idx1 % dim.y) * stride;
		idx1 /= dim.y;
		z = z0 + (idx1 % dim.z) * stride;
	}
	
	iVec3 my_iterator::get_vec(size_t index) const {
		int x, y, z;
		x = x0 + (index % ((x1 - x0) / stride + 1)) * stride;
		index /= (x1 - x0) / stride + 1;
		y = y0 + (index % ((y1 - y0) / stride + 1)) * stride;
		index /= ((y1 - y0) / stride + 1);
		z = z0 + (index % ((z1 - z0) / stride + 1)) * stride;
		return {x, y, z};
	}
	
	void my_iterator::advance() {
		index++;
		if((x += stride) > x1) {
			x = x0;
			if((y += stride) > y1) {
				y = y0;
				z += stride;
			}
		}
	}
	
	bool my_iterator::is_end() const{
		return index >= end;
	}
	
	bool my_iterator::is_empty() const{
		return x0 > x1 || y0 > y1 || z0 > z1;
	}
	
	size_t my_iterator::get_size() const{
		return is_empty()? 0 : (size_t)((x1 - x0) / stride + 1) * ((y1 - y0) / stride + 1) * ((z1 - z0) / stride + 1);
	}
	
	sVec3 my_iterator::get_dim() const {
		return is_empty()? Vec3<size_t>(0, 0, 0) : Vec3<size_t>((x1 - x0) / stride + 1, (y1 - y0) / stride + 1, (z1 - z0) / stride + 1);
	}
	
	iVec3 my_iterator::get_vec() const {
		return {x, y, z};
	}
	
	iVec3 vec3_base[3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	
	std::ostream& operator<<(std::ostream& os, const complex_num& c) {
		os << c.real() << " " << c.imag();
		return os;
	}

	Barrier::Barrier(std::size_t count) : m_count{ count }, m_initial{ count }, m_state{ State::Down } { }

	void Barrier::Sync()
	{
		std::unique_lock<std::mutex> lock{ m_mutex };

		if (m_state == State::Down)
		{
			// Counting down the number of syncing threads
			if (--m_count == 0) {
				m_state = State::Up;
				m_cv.notify_all();
			}
			else {
				m_cv.wait(lock, [this] { return m_state == State::Up; });
			}
		}

		else // (m_state == State::Up)
		{
			// Counting back up for Auto reset
			if (++m_count == m_initial) {
				m_state = State::Down;
				m_cv.notify_all();
			}
			else {
				m_cv.wait(lock, [this] { return m_state == State::Down; });
			}
		}
	}

	size_t Barrier::get_num_proc() const {
		return m_initial;
	}

	Barrier* glob_barrier{ new Barrier{ 1 } };
	void set_num_proc(const size_t num_proc) {
		delete glob_barrier;
		glob_barrier = new Barrier(num_proc);
	}

	
}
