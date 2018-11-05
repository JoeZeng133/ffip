#include <utility.hpp>

using namespace std;

namespace ffip {
	const real e0 = 8.854187817e-12;
	const real u0 = 1.2566370614e-6;
	const real pi = boost::math::constants::pi<double>();
	const real z0 = 376.73031346;
	const real c = 3e8;
	
	const int side_low_tag::val = -1;
	
	int side_low_tag::get_next_int(int x) {
		return x - 1;
	}
	
	int side_low_tag::round(real x) {
		return floor(x);
	}
	
	int side_low_tag::round_wtol(real x) {
		int res = floor(x);
		if(x - res < 0.99)
			return res;
		else
			return res + 1;
	}
	
	const int side_high_tag::val = 1;
	
	int side_high_tag::get_next_int(int x) {
		return x + 1;
	}
	
	int side_high_tag::round(real x) {
		return ceil(x);
	}
	
	int side_high_tag::round_wtol(real x) {
		int res = ceil(x);
		if(res - x < 0.99)
			return res;
		else
			return res - 1;
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
	
	
	/* specializations of Vec3*/
	template<>
	Coord_Type iVec3::get_type() const {
		return static_cast<Coord_Type>((x & 1) | ((y & 1) << 1) | ((z & 1) << 2));
	}
	
	template<>
	Coord_Type iVec3::get_type(const Coord_Type other) const {
		return static_cast<Coord_Type>(other ^ get_type());
	}
	
	fVec3 operator/(const fVec3& a, const real b) {
		return fVec3{a.x / b, a.y / b, a.z / b};
	}
	
	fVec3 operator*(const iVec3& a, const real b) {
		return {a.x * b, a.y * b, a.z * b};
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
	
	real Gaussian_Func::operator()(real time) {
		return exp(-a * time * time);
	}
	
	real Gaussian_Func::operator[](int time) {
		return operator()(time * dt);
	}
	
	/* Sinuosuidal Functions*/
	Sinuosuidal_Func::Sinuosuidal_Func(real freq) {
		a = 2 * pi * freq;
	}
	
	real Sinuosuidal_Func::operator()(real time) {
		return sin(a * time);
	}
	
	real Sinuosuidal_Func::operator[](int time) {
		return operator()(time * dt);
	}
	
	void Sinuosuidal_Func::set_dt(real _dt) {
		dt = _dt;
	}
	
	/* ricker wavelet function*/
	Rickerwavelet_Func::Rickerwavelet_Func(real fp, real d) {
		a = (pi * fp) * (pi * fp);
	}
	
	real Rickerwavelet_Func::operator()(real time) {
		real arg = a * (time - d) * (time - d);
		return (1 - 2 * arg) * exp(-arg);
	}
	
	real Rickerwavelet_Func::operator[](int time) {
		return operator()(time * dt);
	}
	
	void Rickerwavelet_Func::set_dt(real _dt) {
		dt = _dt;
	}
	
	/* Gridded Interpolation class*/
	GriddedInterp::GriddedInterp(string _file_name):file_name(_file_name) {
		fstream f1(_file_name, ios::in);
		if(!f1.is_open())
			throw runtime_error("File name invalid");
		
		int n, m, p;
		f1 >> n >> m >> p;
		for(int i = 0; i < n; ++i) {
			real tmp;
			f1 >> tmp;
			x.push_back(tmp);
		}
		
		for(int i = 0;i < m; ++i) {
			real tmp;
			f1 >> tmp;
			y.push_back(tmp);
		}
		
		for(int i = 0; i < p; ++i) {
			real tmp;
			f1 >> tmp;
			y.push_back(tmp);
		}
		
		for(int i = 0; i < n; ++i)
			for(int j = 0; j < m; ++j)
				for(int k = 0; k < p; ++k) {
					real tmp;
					f1 >> tmp;
					v.push_back(tmp);
				}
		
		init();
	}
	
	GriddedInterp::GriddedInterp(const std::vector<real>& _x, const std::vector<real>& _y, const std::vector<real>& _z, const std::vector<real>& _v):x(_x), y(_y), z(_z), v(_v) {
		init();
	}
	
	void GriddedInterp::set_tol(real _tol) {
		tol = _tol;
	}
	
	void GriddedInterp::init() {
		for(int i = 1; i < x.size(); ++i)
			if(x[i] < x[i - 1])
				throw runtime_error("Monotone series x is requried");
		
		for(int i = 1; i < y.size(); ++i)
			if(y[i] < y[i - 1])
				throw runtime_error("Monotone series y is requried");
		
		for(int i = 1; i < z.size(); ++i)
			if(z[i] < z[i - 1])
				throw runtime_error("Monotone series z is requried");
		
		if (x.size() == 1)
			active_dim[0] = 0;
		else
			active_dim[0] = 1;
		
		if (y.size() == 1)
			active_dim[1] = 0;
		else
			active_dim[1] = 1;
		
		if (z.size() == 1)
			active_dim[2] = 0;
		else
			active_dim[2] = 1;
		
		dim = active_dim[0] + active_dim[1] + active_dim[2];
		
		for(int i = 0; i < active_dim[0] * 2; ++i)
			for(int j = 0; j < active_dim[1] * 2; ++j)
				for(int k = 0; k < active_dim[2] * 2; ++k) {
					int tmp = i * y.size() * z.size() + j * z.size() + k;
					jump_arr.push_back(tmp);
				}
		
		size = jump_arr.size();
		val_arr.resize(size);
		loc_jump_x = y.size() * z.size();
		loc_jump_y = z.size();
		loc_jump_z = 1;
	}
	
	int GriddedInterp::add_val_arr(const std::vector<real> &v, real x, int &index) const{
		/* igore 0 dimension*/
		if(v.size() == 1) return 0;
		
		if(x < v.front() - tol || x > v.back() + tol)
			return -1;
		
		x = max(x, v.front());				//geometry tolorance is important
		x = min(x, v.back());
		
		auto itr = lower_bound(v.begin(), v.end(), x);
		/* -1 if out of bound*/
		
		if(*itr == x)
			if(itr == v.begin())
				itr++;
		
		real v1 = *itr - x;
		real v0 = x - *(itr - 1);
		
		if (index) {
			for(int i = index; i < (index << 1); ++i)
				val_arr[i] = val_arr[i - index] * v1;
			
			for(int i = 0; i < index; ++i)
				val_arr[i] *= v0;
			
			index = (index << 1);
		} else {
			index = 2;
			val_arr[0] = v0;
			val_arr[1] = v1;
		}
		
		return distance(itr - 1, v.begin());
	}
	
	real GriddedInterp::request_value(const fVec3& p) const{
		real res = 0;
		int tail = 0;
		
		int idx_z = add_val_arr(z, p.z, tail);
		if(idx_z < 0) return 0;
		
		int idx_y = add_val_arr(y, p.y, tail);
		if(idx_y < 0) return 0;
		
		int idx_x = add_val_arr(x, p.x, tail);
		if(idx_x < 0) return 0;
		
		int base_index = idx_x * loc_jump_x + idx_y * loc_jump_y + idx_z * loc_jump_z;
		
		for(int i = 0; i < size; ++i) {
			res += val_arr[i] * v[jump_arr[i ^ (size - 1)] + base_index];
		}
		
		return res;
	}
	
	void GriddedInterp::expand_dim(real lo, real hi, Direction dir) {
		real ratio_lo, ratio_hi;
		vector<real> new_v;
		
		switch (dir) {
			case Direction::X:
				if(x.size() > 1 || lo > x[0] || hi < x[0]) break;		//no need to expand if it is already non-zero dimension or x[0] is not inside [lo, hi]
				
				ratio_lo = (hi - x[0]) / (hi - lo) / (hi - lo);			//the lower elements need to be multiplied by a linear interpolation ratio and divided by the length spanned by [lo, hi] to maintain the same integral
				ratio_hi = (x[0] - lo) / (hi - lo) / (hi - lo);			//the same goes for high elements
				new_v.resize(v.size() * 2);
				
				for(int j = 0; j < y.size(); ++j)
					for(int k = 0; k < z.size(); ++k) {
						int index = j * loc_jump_y + k * loc_jump_z;
						new_v[index] = ratio_lo * v[index];
						new_v[index + loc_jump_x] = ratio_hi * v[index];
					}
				v = std::move(new_v);
				break;
				
			case Direction::Y:
				if(y.size() > 1 || lo > y[0] || hi < y[0]) break;
				
				ratio_lo = (hi - y[0]) / (hi - lo) / (hi - lo);
				ratio_hi = (y[0] - lo) / (hi - lo) / (hi - lo);
				new_v.resize(v.size() * 2);
				
				for(int i = 0; i < x.size(); ++i)
					for(int k = 0; k < z.size(); ++k) {
						int index = i * loc_jump_x + k * loc_jump_z;
						new_v[index] = ratio_lo * v[index];
						new_v[index + loc_jump_y] = ratio_hi * v[index];
					}
				v = std::move(new_v);
				break;
				
			case Direction::Z:
				if(z.size() > 1 || lo > z[0] || hi < z[0]) break;
				
				ratio_lo = (hi - z[0]) / (hi - lo) / (hi - lo);
				ratio_hi = (z[0] - lo) / (hi - lo) / (hi - lo);
				new_v.resize(v.size() * 2);
				
				for(int i = 0; i < x.size(); ++i)
					for(int j = 0; j < y.size(); ++j) {
						int index = i * loc_jump_x + j * loc_jump_y;
						new_v[index] = ratio_lo * v[index];
						new_v[index + loc_jump_z] = ratio_hi * v[index];
					}
				v = std::move(new_v);
				break;
				
			default:
				break;
		}
	}
	
	real GriddedInterp::request_integral() const {
		return integral3(x, y, z, v.data());
	}
	
	fVec3 GriddedInterp::get_p1() const{
		return {x.front(), y.front(), z.front()};
	}
	
	fVec3 GriddedInterp::get_p2() const{
		return {x.back(), y.back(), z.back()};
	}
	
	
	/* PML */
	PML::PML(Direction _dir, Side _side): dir(_dir), side(_side) {}
	
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
	
	real PML::get_chi(const real x) const {
		return 1 + (chi_max - 1) * pow(x / d, m);
	}
	
	real PML::get_a(const real x) const {
		return a_max * pow(1 - x / d, m_a);
	}
	
	real PML::get_b(const real x) const {
		return exp(-dt * (get_sigma(x) / e0 / get_chi(x) + get_a(x) / e0));
	}
	
	real PML::get_c(const real x) const {
		return get_sigma(x) * (get_b(x) - 1) / (get_sigma(x) * get_chi(x) + get_chi(x) * get_chi(x) * get_a(x));
	}
	
	real PML::optimal_sigma_max(real m_k, real dt, real er, real ur) {
		return 0.8 * (m_k + 1) / (z0 * dt * sqrt(er * ur));
	}
	
	/* PML point constructor*/
	PML_Point::PML_Point(const int _index, const int _jump_pos, const int _jump_neg, const real _b_pos, const real _b_neg, const real _c_pos, const real _c_neg):
	index(_index), jump_pos(_jump_pos), jump_neg(_jump_neg), b_pos(_b_pos), b_neg(_b_neg), c_pos(_c_pos), c_neg(_c_neg) {}
	
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, Coord_Type ctype) {
		auto tmp = get_component_interior(p1, p2, ctype);
		if(ElementWise_Less_Eq(p1, p2)) {
			x = x0 = tmp.first.x;
			y = y0 = tmp.first.y;
			z = z0 = tmp.first.z;
			x1 = tmp.second.x;
			y1 = tmp.second.y;
			z1 = tmp.second.z;
		} else {
			x = x1 + 1;
		}
	}
	
	void my_iterator::advance() {
		index++;
		if((z += 2) > z1) {
			z = z0;
			if((y += 2) > y1) {
				y = 0;
				x += 2;
			}
		}
	}
	
	bool my_iterator::is_end() {
		return x > x1;
	}
	
	bool my_iterator::is_empty() {
		return x0 > x1 || y0 > y1 || z0 > z1;
	}
	
	size_t my_iterator::size() {
		return is_empty()? 0 : (size_t)((x1 - x0 + 2) >> 1) * ((y1 - y0 + 2) >> 1) * ((z1 - z0 + 2) >> 1);
	}

}
