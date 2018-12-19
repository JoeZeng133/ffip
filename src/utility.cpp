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
		return floor(x);
	}
	
	const int side_high_tag::val = 1;
	
	int side_high_tag::get_next_int(int x) {
		return x + 1;
	}
	
	int side_high_tag::round(real x) {
		return ceil(x);
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
	
	/* Gaussian functions */
	Gaussian_Func::Gaussian_Func(real sigma_t, real _mu):mu(_mu) {
		a = 0.5 / (sigma_t * sigma_t);
	}
	
	Gaussian_Func::Gaussian_Func(real sigma_f) {
		real sigma_t = 1 / (2 * pi) / sigma_f;
		a = 0.5 / (sigma_t * sigma_t);
	}
	
	void Gaussian_Func::set_dt(real _dt) {
		dt = _dt;
	}
	
	real Gaussian_Func::operator()(real time) const{
		return exp(-a * time * time);
	}
	
	real Gaussian_Func::operator[](int time) const{
		return operator()(time * dt);
	}
	
	auto Gaussian_Func::get_functor() -> std::function<real(const real)> const {
		
		return [a=this->a](const real time) -> real {
			return exp(-a * time * time);
		};
	}
	
	/* Sinuosuidal Functions*/
	Sinuosuidal_Func::Sinuosuidal_Func(real freq, real d): a(2 * pi * freq), phase(d / (2 * pi * freq)) {}
	
	real Sinuosuidal_Func::operator()(real time) const{
		return sin(a * time - phase);
	}
	
	real Sinuosuidal_Func::operator[](int time) const{
		return operator()(time * dt);
	}
	
	void Sinuosuidal_Func::set_dt(real _dt) {
		dt = _dt;
	}
	
	auto Sinuosuidal_Func::get_functor() -> std::function<real(const real)> const {
		return [a = this->a, phase = this->phase](const real time) -> real {
			return sin(a * time - phase);
		};
	}
	
	/* ricker wavelet function*/
	Rickerwavelet_Func::Rickerwavelet_Func(real fp, real _d) {
		a = (pi * fp) * (pi * fp);
		d = _d;
	}
	
	real Rickerwavelet_Func::operator()(real time) const{
		real arg = a * (time - d) * (time - d);
		//std::cout << "testing inside rickerwavelet function " << a  << " " << (time - d) << std::endl;
		return (1 - 2 * arg) * exp(-arg);
	}
	
	real Rickerwavelet_Func::operator[](int time) const{
		return operator()(time * dt);
	}
	
	void Rickerwavelet_Func::set_dt(real _dt) {
		dt = _dt;
	}
	
	auto Rickerwavelet_Func::get_functor() -> std::function<real(const real)> const {
		return [a = this->a, d = this->d](const real time) -> real {
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
		real dx, dy, dz;
		
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
	}
	
	GriddedInterp::GriddedInterp(const iVec3& dim, const fVec3& w0, const fVec3& dw, const real_arr& _v):
	v(_v) {
		
		for(int i = 0; i < dim.x; ++i)
			x.push_back(w0.x + dw.x * i);
		
		for(int j = 0; j < dim.y; ++j)
			y.push_back(w0.y + dw.y * j);
		
		for(int k = 0; k < dim.z; ++k)
			z.push_back(w0.z + dw.z * k);
	}

	
	real GriddedInterp::request_value(const fVec3& p) const{
		return interp_ndim(v.data(), p.z, z, p.y, y, p.x, x);
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
	PML::PML(Direction _dir, Side _side): dir(_dir), side(_side) {}
	
	PML::PML(Direction _dir, Side _side, real _d): d(_d), dir(_dir), side(_side) {}
	
	PML::PML(Direction _dir, Side _side, real _d, real _sigma_max): d(_d), sigma_max(_sigma_max), dir(_dir), side(_side) {}
	
	Direction PML::get_dir() const {
		return dir;
	}
	
	Side PML::get_side() const {
		return side;
	}
	
	real PML::get_d() const {
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
	
	/* PML point constructor*/
	PML_Point::PML_Point(const int _index, const int _jump_pos, const int _jump_neg, const real _b_pos, const real _b_neg, const real _c_pos, const real _c_neg):
	index(_index), jump_pos(_jump_pos), jump_neg(_jump_neg), b_pos(_b_pos), b_neg(_b_neg), c_pos(_c_pos), c_neg(_c_neg) {}
	
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, Coord_Type ctype) {
		decltype(get_component_interior(p1, p2, ctype)) tmp;
		
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
		
		size = get_size();
		if (size == 0)
			z = z1 + 1;
	}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype, const int rank, const int num): my_iterator(p1, p2, ctype) {
		if (rank >= num)
			throw std::runtime_error("Rank is larger than number of threads");
		
		if (size <= num)
			throw std::runtime_error("Region is too small to divide");
		
		int idx1 = (rank * size) / num;
		int idx2 = ((rank + 1) * size) / num;
		size = idx2 - idx1;
		
		x = x0 + (idx1 % ((x1 - x0) / jump + 1)) * jump;
		idx1 /= (x1 - x0) / jump + 1;
		y = y0 + (idx1 % ((y1 - y0) / jump + 1)) * jump;
		idx1 /= ((y1 - y0) / jump + 1);
		z = z0 + (idx1 % ((z1 - z0) / jump + 1)) * jump;
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
		return index >= size;
	}
	
	bool my_iterator::is_empty() const{
		return x0 > x1 || y0 > y1 || z0 > z1;
	}
	
	size_t my_iterator::get_size() const{
		return is_empty()? 0 : (size_t)((x1 - x0) / jump + 1) * ((y1 - y0) / jump + 1) * ((z1 - z0) / jump + 1);
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

	std::pair<iVec3, iVec3> divide_region(iVec3 p1, iVec3 p2, const int r, const int n) {
		if (!ElementWise_Less_Eq(p1, p2))
			throw std::runtime_error("The region cannot be divided");

		iVec3 dp = (p2 - p1) + iVec3{ 1, 1, 1 };
		if (dp.z >= dp.y && dp.z >= dp.x) {
			if (n > dp.z)
				throw std::runtime_error("The region cannot be divided");

			p2.z = ((r + 1) * dp.z) / n - 1 + p1.z;
			p1.z = (r * dp.z) / n + p1.z;
			
		}
		else if (dp.y >= dp.x && dp.y >= dp.z) {
			if (n > dp.y)
				throw std::runtime_error("The region cannot be divided");

			p2.y = ((r + 1) * dp.y) / n - 1 + p1.y;
			p1.y = (r * dp.y) / n + p1.y;
			
		}
		else {
			if (n > dp.x)
				throw std::runtime_error("The region cannot be divided");

			p2.x = ((r + 1) * dp.x) / n - 1 + p1.x;
			p1.x = (r * dp.x) / n + p1.x;
		}

		return { p1, p2 };
	}
}
