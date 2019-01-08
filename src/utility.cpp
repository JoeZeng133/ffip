#include <utility.hpp>

using namespace std;

namespace ffip {
	const real e0 = 8.854187817e-12;
	const real u0 = 1.2566370614e-6;
	const real pi = 3.141592653589793e+00;
	const real z0 = 376.73031346;
	const real c0 = 3e8;
	
	int Ctype2DirInt(const Coord_Type ctype) {
		static int look_up_table[8] = {-1, 0, 1, 2, 2, 1, 0, -1}; //000, 001, 010, 011, 100, 101, 110, 111
		if (ctype == Null)
			throw runtime_error("Ctype cannot be Null");
		
		return look_up_table[static_cast<int>(ctype)];
	}
	
	const int side_low_tag::val = -1;
	
	int side_low_tag::get_next_int(int x) {
		return x - 1;
	}
	
	int side_low_tag::round(real x) {
		return int(floor(x));
	}
	
	const int side_high_tag::val = 1;
	
	int side_high_tag::get_next_int(int x) {
		return x + 1;
	}
	
	int side_high_tag::round(real x) {
		return int(ceil(x));
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
	
	/* specializations of Vec3*/
	template<>
	Coord_Type iVec3::get_type() const {
		return static_cast<Coord_Type>((x & 1) | ((y & 1) << 1) | ((z & 1) << 2));
	}
	
	template<>
	Coord_Type iVec3::get_type(const Coord_Type other) const {
		return static_cast<Coord_Type>(other ^ get_type());
	}
	
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
	
	/* Gridded Interpolation class*/
	GriddedInterp::GriddedInterp(string _file_name):file_name(_file_name) {
		fstream fin{file_name, ios::in};
		
		if (!fin.is_open()) {
			throw std::runtime_error("Cannot open the interpolation file");
		}
		
		int n, m, p;
		real x0, y0, z0;
		
		fin >> n >> m >> p;
		
		fin >> x0 >> dx;
		fin >> y0 >> dy;
		fin >> z0 >> dz;
		
		for(int i = 0; i < n; ++i)
			x.push_back(x0 + dx * i);
		
		for(int j = 0; j < m; ++j)
			y.push_back(y0 + dy * j);
		
		for(int k = 0; k < p; ++k)
			z.push_back(z0 + dz * k);
		
		for(int i = 0; i < n * m * p; ++i) {
			double tmp;
			fin >> tmp;
			v.push_back(tmp);
		}
		
		interp = interpn<3>(p, m, n);
	}
	
	GriddedInterp::GriddedInterp(const iVec3& dim, const fVec3& w0, const fVec3& dw, const real_arr& _v):
	v(_v), dx(dw.x), dy(dw.y), dz(dw.z) {
		
		for(int i = 0; i < dim.x; ++i)
			x.push_back(w0.x + dw.x * i);
		
		for(int j = 0; j < dim.y; ++j)
			y.push_back(w0.y + dw.y * j);
		
		for(int k = 0; k < dim.z; ++k)
			z.push_back(w0.z + dw.z * k);
		
		interp = interpn<3>(z.size(), y.size(), x.size());
	}

	real GriddedInterp::request_value(const fVec3& p) const{
		return interp.get(v, (p.z - z.front()) / dz, (p.y - y.front()) / dy, (p.x - x.front()) / dx);
	}
	
	void GriddedInterp::expand_dim_helper(std::vector<real> &w, const real lo, const real hi, const Direction dir) {
		if(w.size() > 1 || lo > w[0] || hi < w[0] || lo >= hi) return;		//no need to expand if it is already non-zero dimension or w[0] is not inside [lo, hi]
		
		real ratio_lo, ratio_hi;
		vector<real> new_v(v.size() * 2);
		
		ratio_lo = 2 * (hi - w[0]) / (hi - lo) / (hi - lo);			//the lower elements need to be multiplied by a linear interpolation ratio and divided by the length spanned by [lo, hi] to maintain the same integral
		ratio_hi = 2 * (w[0] - lo) / (hi - lo) / (hi - lo);			//the same goes for high elements
		
		Vec3<size_t> new_dim = {x.size() + (dir == X), y.size() + (dir == Y), z.size() + (dir == Z)};
		Vec3<size_t> jump = {1, x.size(), x.size() * y.size()};
		Vec3<size_t> new_jump = {1, new_dim.x, new_dim.x * new_dim.y};
		
		auto get_index = [](const size_t x, const size_t y, const size_t z, const Vec3<size_t>& jump) -> size_t {
			return x * jump.x + y * jump.y + z * jump.z;
		};//used for calculation of index
		
		//get new value arrays based on interpolation
		for(int i = 0; i < x.size(); ++i)
			for(int j = 0; j < y.size(); ++j)
				for(int k = 0; k < z.size(); ++k) {
					size_t index = get_index(i, j, k, jump);
					size_t new_index_lo = get_index(i, j, k, new_jump);
					size_t new_index_hi = get_index(i + (dir == X), j + (dir == Y), k + (dir == Z), new_jump);
					
					new_v[new_index_lo] = ratio_lo * v[index];
					new_v[new_index_hi] = ratio_hi * v[index];
				}
		
		//final udpates
		w[0] = lo;
		w.push_back(hi);
		v = std::move(new_v);
	}
	
	void GriddedInterp::expand_dim(const real lo, const real hi, const Direction dir) {
		switch (dir) {
			case Direction::X:
				expand_dim_helper(x, lo, hi, dir);
				break;
				
			case Direction::Y:
				expand_dim_helper(y, lo, hi, dir);
				break;
				
			case Direction::Z:
				expand_dim_helper(z, lo, hi, dir);
				break;
				
			default:
				break;
		}
	}
	
	void GriddedInterp::scale_xyz(const real sx, const real sy, const real sz) {
		if (sx == 0 || sy == 0 || sz == 0)
			throw runtime_error("The scale factor cannot be zero");
		
		for_each(x.begin(), x.end(), [=](real &x){x *= sx;});
		for_each(y.begin(), y.end(), [=](real &y){y *= sy;});
		for_each(z.begin(), z.end(), [=](real &z){z *= sz;});
	}
	
	real GriddedInterp::request_integral() const {
		return integral_ndim(v.data(), z, y, x);
	}
	
	fVec3 GriddedInterp::get_p1() const{
		return {x.front(), y.front(), z.front()};
	}
	
	fVec3 GriddedInterp::get_p2() const{
		return {x.back(), y.back(), z.back()};
	}
	
	inline real interp_helper(const real* data, const real w) {
		return data[0] * w + data[1] * w;
	}
	
	/* PML */
	PML::PML(int _d, real _sigma_max): d(_d), sigma_max(_sigma_max) {}
	PML::PML(int _d, real _sigma_max, real _k_max, real _a_max, real _m, real _m_a):d(_d), sigma_max(_sigma_max), k_max(_k_max), a_max(_a_max), m(_m), m_a(_m_a) {};
	
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
	
	/* PML point constructor*/
	PML_Point::PML_Point(const int _index, const int _jump_pos, const int _jump_neg, const real _b_pos, const real _b_neg, const real _c_pos, const real _c_neg):
	index(_index), jump_pos(_jump_pos), jump_neg(_jump_neg), b_pos(_b_pos), b_neg(_b_neg), c_pos(_c_pos), c_neg(_c_neg) {}
	
	/* ############################################### */
	my_iterator::my_iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype): my_iterator{corners.first, corners.second, ctype} {}
	
	my_iterator::my_iterator(const fVec3& p1, const fVec3& p2, const Coord_Type ctype): my_iterator{get_component_interior(p1, p2, ctype), ctype} {}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype) {
		std::pair<iVec3, iVec3> tmp;
		
		// loop through points inside [p1, p2]. Null means loop through every types of points
		if (ctype != Coord_Type::Null)  {
			tmp = get_component_interior(p1, p2, ctype);
			jump = 2;
		} else {
			tmp = {p1, p2};
			jump = 1;
		}
		
		x = x0 = tmp.first.x;
		y = y0 = tmp.first.y;
		z = z0 = tmp.first.z;
		x1 = tmp.second.x;
		y1 = tmp.second.y;
		z1 = tmp.second.z;
		
		index = 0;
		end = size = get_size();
		if (size == 0)
			z = z1 + 1;
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
		
		x = x0 + (idx1 % ((x1 - x0) / jump + 1)) * jump;
		idx1 /= (x1 - x0) / jump + 1;
		y = y0 + (idx1 % ((y1 - y0) / jump + 1)) * jump;
		idx1 /= ((y1 - y0) / jump + 1);
		z = z0 + (idx1 % ((z1 - z0) / jump + 1)) * jump;
	}
	
	iVec3 my_iterator::get_vec(size_t index) const {
		int x, y, z;
		x = x0 + (index % ((x1 - x0) / jump + 1)) * jump;
		index /= (x1 - x0) / jump + 1;
		y = y0 + (index % ((y1 - y0) / jump + 1)) * jump;
		index /= ((y1 - y0) / jump + 1);
		z = z0 + (index % ((z1 - z0) / jump + 1)) * jump;
		return {x, y, z};
	}
	
	void my_iterator::advance() {
		index++;
		if((x += jump) > x1) {
			x = x0;
			if((y += jump) > y1) {
				y = y0;
				z += jump;
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
		return is_empty()? 0 : (size_t)((x1 - x0) / jump + 1) * ((y1 - y0) / jump + 1) * ((z1 - z0) / jump + 1);
	}
	
	Vec3<size_t> my_iterator::get_dim() const {
		return is_empty()? Vec3<size_t>(0, 0, 0) : Vec3<size_t>((x1 - x0) / jump + 1, (y1 - y0) / jump + 1, (z1 - z0) / jump + 1);
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

	Barrier* glob_barrier{ new Barrier{ 1 } };
	void set_num_proc(const size_t num_proc) {
		delete glob_barrier;
		glob_barrier = new Barrier(num_proc);
	}

	size_t Barrier::get_num_proc() const {
		return m_initial;
	}
}
