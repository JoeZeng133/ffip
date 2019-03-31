#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>
#include <thread>
#include <complex>
#include <iomanip>
#include <numeric>
#include <future>
#include <valarray>
#include <type_traits>
#include <mpi.h>
#include <medium.hpp>

namespace ffip {
	//exceptions
	struct Invalid_Corner_Points {};
	struct Invalid_Direction {};
	struct Chunk_is_Not_Initialized{};
	struct Invalide_Coord_Type{};
	struct Out_of_the_Domain{};
	
	//type aliases
	using double_arr = std::vector<double>;
	using complex_pair = std::pair<std::complex<double>, std::complex<double> >;
	using complex_num = std::complex<double>;
	using complex_arr = std::vector<complex_num>;

	//constants
	extern const double pi;
	extern unsigned int hex_bit_count[16];
	extern unsigned int dir_of_ctypes[16];

	//enums with bit pattern indication
	enum Direction {X = 0, Y = 1, Z = 2};
	enum Side {Negative = -1, Positive = 1};
	enum Coord_Type {
		Ex = 0b001, Ey = 0b010, Ez = 0b100, Hx = 0b110, Hy = 0b101, Hz = 0b011, 
		Dx = 0b1001, Dy = 0b1010, Dz = 0b1100, Bx = 0b1110, By = 0b1101, Bz = 0b1011,
		Corner = 0b000, Center = 0b111, All = 0b1000};

	//initializations of utility (calculating constants)
	void init();

	//count the number of 1s in the first 4 bits
	constexpr unsigned int count_hex_bits(int x) {

		return (x & 1) + ((x & 2) >> 1) + ((x & 4) >> 2) + ((x & 8) >> 3);
	}

	//get direction integer from coord_type
	constexpr int get_dir_from_ctype(Coord_Type ctype) {
		return dir_of_ctypes[ctype & 0b111];
	}

	//get normal vector in certain direction with sign
	constexpr iVec3 get_norm_vec(Direction dir, Side side) {
		iVec3 res;
		res[dir] = side;
		return res;
	}

	//return whether it is electrical point
	constexpr bool is_e_point(const Coord_Type ctype) {

		return count_hex_bits(ctype & 0b111) == 1;
	}
	
	//return wehether it is magnetic point
	constexpr bool is_m_point(const Coord_Type ctype) {

		return count_hex_bits(ctype & 0b111) == 2;
	}
	
	// light weight 3 element vector
	template<typename T = double>
	struct Vec3{
		using value_type = T;

		T x{}, y{}, z{};

		Vec3() = default;
		Vec3(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}

		// copy construction
		template<typename T2>
		explicit Vec3(const Vec3<T2>& other): x(other.x), y(other.y), z(other.z) {}

		// copy assignment
		template<typename T2>
		Vec3& operator=(const Vec3<T2>& other) {
			return Vec3{other};
		};

		//Unary -
		Vec3 operator-() const {
			return {-x, -y, -z};
		}

		//Unary+
		Vec3 operator+() const {
			return {x, y, z};
		}

		//dot product
		template<typename T2>
		std::common_type_t<T, T2> dot(const Vec3<T2>& other) const {
			return x * other.x + y * other.y + z * other.z;
		}

		//scalar division
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator/(const T2 mult) const {
			return {x / mult, y / mult, z / mult};
		}

		//scalar division assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator/=(const T2 mult) const {
			x /= mult; y /= mult; z /= mult;
		}

		//scalar multiplication
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const T2 mult) const {
			return {x * mult, y * mult, z * mult};
		}

		//scalar multiplication assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator*=(const T2 mult) const {
			x *= mult; y *= mult; z *= mult;
		}

		//pointwise multiplication
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const Vec3<T2>& other) const {
			return {x * other.x, y * other.y, z * other.z};
		}

		//pointwise multiplication assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator*=(const Vec3<T2>& other) const {
			x *= other.x; y *= other.y; z *= other.z;
		}

		//scalar addition
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const T2 add) const {
			return {x + add, y + add, z + add};
		}

		//pointwise scalar addition
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const Vec3<T2>& other) const {
			return {x + other.x, y + other.y, z + other.z};
		}

		//scalar addition assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator+=(const Vec3<T2>& other) const {
			x += other.x; y += other.y; z += other.z;
		}

		//scalar subtraction
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const T2 minus) const {
			return {x - minus, y - minus, z - minus};
		}

		//scalar subtraction assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator-=(const T2 minus) const {
			x -= minus; y -= minus; z -= minus;
		}

		//pointwise scalar subtraction
		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const Vec3<T2>& other) const {
			return {x - other.x, y - other.y, z - other.z};
		}

		//pointwise scalar subtraction assignment
		template<typename T2>
		Vec3<std::common_type_t<T, T2>>& operator-=(const Vec3<T2>& other) const {
			x -= other.x; y -= other.y; z -= other.z;
		}
		
		//ceil
		Vec3 ceil() const {
			return Vec3( std::ceil(x), std::ceil(y), std::ceil(z) );
		}

		//floor
		Vec3 floor() const {
			return Vec3( std::floor(x), std::floor(y), std::floor(z) );
		}

		//absolute value
		Vec3 abs() const {
			return Vec3( std::abs(x),  std::abs(y), std::abs(z) );
		}

		//return x * y * z
		T prod() const {
			return x * y * z;	
		}

		//get<N> style access
		template<unsigned int N>
		T& get() {
			if constexpr(N == 0)
				return x;
			else if constexpr(N == 1)
				return y;
			else if constexpr(N == 2)
				return z;
		}

		//get<N> style access, const
		template<unsigned int N>
		const T& get() const {
			if constexpr(N == 0)
				return x;
			else if constexpr(N == 1)
				return y;
			else if constexpr(N == 2)
				return z;
		}

		//[] style access
		T& operator[](const size_t n) {
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		//[] style access
		const T& operator[](const size_t n) const {
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		//get Coord_Type for the given grid coordinates
		template<typename Dummy = T, typename = std::enable_if_t<std::is_integral<Dummy>::value>>
		Coord_Type get_type() const {
			return static_cast<Coord_Type>((x & 1) | ((y & 1) << 1) | ((z & 1) << 2));
		}
	};

	using iVec3 = Vec3<int>;			
	using fVec3 = Vec3<double>;	
	using cVec3 = Vec3<complex_num>;
	using sVec3 = Vec3<long long>;

	//strictly pointwise less
	template<typename T1, typename T2>
	constexpr bool le_vec3(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x < b.x && a.y < b.y && a.z < b.z;
	}

	//pointwise less equal
	template<typename T1, typename T2>
	constexpr bool leq_vec3(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x <= b.x && a.y <= b.y && a.z <= b.z;
	}

	//array dimension to stride
	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	constexpr Vec3<T> dim2stride(const Vec3<T>& dim) {
		return { 1, dim.x, dim.x * dim.y }; 
	}
	
	// template<int N, typename T>
	// constexpr Vec3<T> rotate_frame(const Vec3<T>& p) {
	// 	if constexpr(N == 0)
	// 		return {p.y, p.z, p.x};
	// 	else if constexpr(N == 1)
	// 		return {p.z, p.x, p.y};
	// 	else if constexpr(N == 2)
	// 		return p;
	// }
	
	//return nearest odd(even) integer of a integer
	template<int S, unsigned int T>
	constexpr int get_nearest_int(const int x) {
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");
		
		return (x & 1) ^ T? x + S : x;
	}
	
	//return nearest odd(even) integer of a double
	template<int S, unsigned int T>
	constexpr int get_nearest_int(const double x) {
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");
		
		if constexpr(S == -1)
			return get_nearest_int<S, T>(std::floor(x));
		else
			return get_nearest_int<S, T>(std::ceil(x));
	}

	//get nearest point of a coord type
	template<int S, typename T>
	constexpr iVec3 get_nearest_point(const Vec3<T>& p, const Coord_Type ctype) {
		iVec3 res;

		if (ctype & 1) {
			res.x = get_nearest_int<S, 1>(p.x);
		} else {
			res.x = get_nearest_int<S, 0>(p.x);
		}

		if (ctype & 2) {
			res.y = get_nearest_int<S, 1>(p.y);
		} else {
			res.y = get_nearest_int<S, 0>(p.y);
		}

		if (ctype & 4) {
			res.z = get_nearest_int<S, 1>(p.z);
		} else {
			res.z = get_nearest_int<S, 0>(p.z);
		}

		return res;
	}

	//return grid points inside a box
	template<typename T1, typename T2>
    constexpr std::pair<iVec3, iVec3> 
	get_component_interior(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type ctype) {

		if (ctype == Coord_Type::All)
			return {iVec3(p1.ceil()), iVec3(p2.floor())};

		return std::make_pair(	get_nearest_point<Positive>(p1, ctype), get_nearest_point<Negative>(p2, ctype));
	}

	//return grid points enclosing the box
	template<typename T1, typename T2>
	constexpr std::pair<iVec3, iVec3> 
	get_component_closure(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {

		if (type == Coord_Type::All)
			return {iVec3(p1.floor()), iVec3(p2.ceil())};

		return std::make_pair(	get_nearest_point<Negative>(p1, type), get_nearest_point<Positive>(p2, type));
	}

	// constexpr std::pair<iVec3, iVec3>
	// 	get_component_closure_exact(const fVec3& p1, const fVec3& p2, const Coord_Type ctype, double dx) {
	// 		auto res = get_component_closure(p1/dx, p2/dx, ctype);
	// 		int stride = (ctype == Coord_Type::All)? 1 : 2;
			
	// 		res.first -= stride * 2;
	// 		res.second += stride * 2;
	// 		while( res.first.x * dx <= p1.x ) res.first.x++;
	// 		while( res.first.y * dx <= p1.y ) res.first.y++;
	// 		while( res.first.z * dx <= p1.z ) res.first.z++;

	// 		while( res.second.x * dx >= p2.x ) res.first.x++;
	// 		while( res.second.y * dx >= p2.y ) res.first.y++;
	// 		while( res.second.z * dx >= p2.z ) res.first.z++;

	// 		res.first -= 1;
	// 		res.second += 1;

	// 		return res;
	// 	}

	//return intersection of two intervals
	template<typename T>
	constexpr std::pair<T, T> 
	get_intersection(const T p1, const T p2, const T q1, const T q2) {

		return {std::max(p1, q1), std::min(p2, q2)};
	}

	//return intersection of two closed cubes
	template<typename T>
	constexpr std::pair<Vec3<T>, Vec3<T>> 
	get_intersection(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {

		return	{{std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
				{std::min(p2.x, q2.x),std::min(p2.y, q2.y), std::min(p2.z, q2.z)}};
	}

	//determine if two closed interval intersects, [p1, p2] ^ [q1, q2]
	template<typename T>
	constexpr bool is_intersect(const T p1, const T p2, const T q1, const T q2) {
		return p1 <= q2 && q1 <= p2;
	}

	//determine if two closed cubes intersects
	template<typename T>
	constexpr bool 
	is_intersect(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {

		return	is_intersect(p1.x, p2.x, q1.x, q2.x) && 
				is_intersect(p1.y, p2.y, q1.y, q2.y) && 
				is_intersect(p1.z, p2.z, q1.z, q2.z);
	}

	//return one face of a box
	template<unsigned int D, int S, typename T>
	std::pair<Vec3<T>, Vec3<T>> get_face(Vec3<T> p1, Vec3<T> p2) {
		
		static_assert(S == -1 || S == 1, "S must be -1 or 1");
		
		if constexpr(S == -1) {
			p2.template get<D>() = p1.template get<D>();
		}
		else {
			p1.template get<D>() = p2.template get<D>();
		}
		
		return {p1, p2};
	}

	//return one face of a box, non-template function
	template<typename T>
	std::pair<Vec3<T>, Vec3<T>> get_face(Vec3<T> p1, Vec3<T> p2, Direction dir, Side side) {

		switch(dir) {
			case X:
				if(side == Positive)
					p1.template get<X>() = p2.template get<X>();
				else
					p2.template get<X>() = p1.template get<X>();
				break;

			case Y:
				if(side == Positive)
					p1.template get<Y>() = p2.template get<Y>();
				else
					p2.template get<Y>() = p1.template get<Y>();
				break;

			case Z:
				if(side == Positive)
					p1.template get<Z>() = p2.template get<Z>();
				else
					p2.template get<Z>() = p1.template get<Z>();
				break;
		}

		return {p1, p2};
	}

	//return whether it is in a closed interval
	template <typename T>
	bool is_in_interval(const T& value, const T& Negative, const T& Positive) {

		return !(value < Negative) && !(Positive < value);
	}
	
	//interpn, linear interpolation with nearest extrapolation
	template<int N>
	class interpn : public interpn<N - 1> {
	public:
		using base_class = interpn<N - 1>;
		size_t dimn{1}, stride{1};
		
		interpn() = default;
		
		template<typename... Args>
		interpn(const size_t _dimn, Args... args): base_class(args...), dimn(_dimn) {

			if (dimn < 1) throw std::runtime_error("Invalid Dimension");

			stride = base_class::get_size();
		}
		
		//get the array size
		size_t get_size() const {
			return dimn * base_class::get_size();
		}
		
		//interpolate xn, xn-1, ..., x1
		template<typename T, typename... Args>
		T get(T const* data, double xq, Args&&... args) const {
			//ignore this dimension if it is 1
			if (dimn == 1)
				return base_class::get(data, args...);
			
			//nearest extrapolation
			xq = std::clamp<int>(0, dimn - 1, xq);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (index + 1 < dimn)
				return	tx * base_class::get(data + stride * (index + 1), args...) +
						(1 - tx) * base_class::get(data + stride * index, args...);
			else
				return	base_class::get(data + stride * index, args...);
		}
		
		//interpolate xn, xn-1, ..., x1
		template<typename T, typename... Args>
		T operator()(const std::vector<T>& vec, Args&&... args) const {

			static_assert(sizeof...(args) == N, "Invalid Request");
			return get(vec.data(), args...);
		}
		
		//transpose_helper 
		template<typename T, typename... Args>
		void transpose_helper(T* data, const T& val, double xq, Args&&... args) const {

			//ignore 1 dimension
			if ( dimn == 1) {
				base_class::get_helper(data, args...);
				return;
			}
			
			//nearest extrapolation
			xq = std::clamp<int>(0, dimn - 1, xq);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (index + 1 < dimn)
				base_class::transpose_helper(data + stride * (index + 1), val * tx, args...);

			base_class::transpose_helper(data + stride * index, val * (1 - tx), args...);
		}
		
		// the put function might skip placing 0s in some elements of the array
		template<typename T, typename... Args>
		void transpose(std::vector<T>& vec, const T& val, Args&&... args) const {

			static_assert(sizeof...(args) == N, "Invalid Request");
			transpose_helper(vec.data(), val, args...);
		}
	};
	
	template<>
	class interpn<1> {
	public:
		size_t dimn{1}, stride{1};
		
		interpn() = default;
		
		interpn(const size_t _dimn): dimn(_dimn) {

			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");
			
			stride = 1;
		}
		
		size_t get_size() const {
			return dimn;
		}
		
		template<typename T>
		T get(const T* data, double xq) const {

			//ignore this dimension if it is 1
			if (dimn == 1)		
				return data[0];
			
			//nearest extrapolation
			xq = std::clamp<int>(0, dimn - 1, xq);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (index + 1 < dimn)
				return tx * data[index + 1] + (1 - tx) * data[index];
			else
				return data[index];
			
		}
		
		template<typename T>
		T operator()(const std::vector<T>& vec, double xq) const {
			return operator()(vec.data(), xq);
		}
		
		template<typename T>
		void transpose(std::vector<T>& vec, const T& val, double xq) const{
			return transpose_helper(vec.data(), val, xq);
		}
		
		template<typename T>
		void transpose_helper(T* data, const T& val, double xq) const{

			//ignore this dimension if it is 1
			if ( dimn == 1) {		
				data[0] = val;
				return;
			}
			
			xq = std::clamp<int>(0, dimn - 1, xq);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (index + 1 < dimn) {
				data[index + 1] = val * tx;
				data[index] = val * (1 - tx);
			}
			else {
				data[index] = val;
			}
		}
	};

	//Iterator through points in Yee Cells
	struct Yee_Iterator {

		int x0, y0, z0;
		int lenx, leny, lenz;
		int i, j, k;
		int stride;
		
		Yee_Iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype);
		Yee_Iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype);

		//move to next coordinate
		void next();

		//does it reach the end
		bool is_end() const;

		//is it empty
		bool is_empty() const;

		//return the current coordinate
		iVec3 get_coord() const;

		//return total num of coordinates
		size_t get_size() const;

		//return dimension of the region
		iVec3 get_dim() const;
	};

	//Yee3 grid
	//Ghost points aware interpolation
    struct Yee3 {
		
		//grid points coordinates
        iVec3 grid_p1, grid_p2;
		//ghost point coordinates
		iVec3 ghost_p1, ghost_p2;
		//dimensions of the array
		iVec3 dim;
		//strides of the array
		sVec3 stride;

		Yee3(iVec3 ghost_p1, iVec3 ghost_p2);

		//get grid lower corner, excluding ghost point
		iVec3 get_grid_p1() const;

		//get grid uppper corner, excluding ghost point
		iVec3 get_grid_p2() const;

		//return size of array compatible with the grid
		size_t get_size() const;

		//return strides in array
		sVec3 get_stride() const;

		//an alias for get_nearest_point<Negative>(pos, ctype)
		iVec3 get_base_point(const fVec3& pos, Coord_Type ctype) const;

		//get coordinate from a index, no range checking
		iVec3 get_coord_from_index(size_t n) const;

		//get index for a given point, no range checking
		size_t get_index_from_coord(const iVec3& pos) const;
		size_t get_index_from_coord(int i, int j, int k) const;

		//inside grids, excluding ghost points
		bool is_inside(const fVec3& pos) const;
		bool is_inside(int i, int j, int k) const;
		bool is_inside(const iVec3& pos) const;

		//get index offset, same as get_index(p + offset) - get_index(p)
		long long get_index_offset(const iVec3& offset) const;
		long long get_index_offset(int i, int j, int k) const;

		//get corresponding dimension
		template<int N>
		size_t get_dim() const {
			return dim.get<N>();
		}

		

		//access the raw value from a grid coordinate, no range checking
		template<typename T>
		T get_raw_val(const std::vector<T>& data, const iVec3& pos) const {
			return data[get_index_from_coord(pos)];
		}

		//return 3 dim interpolation weights based on sx, sy, sz
		std::vector<double> get_interp_weights(double sx, double sy, double sz) const;

		//s.x, s.y, s.z
		std::vector<double> get_interp_weights(const fVec3& s) const;

		//return interpolation weights
        std::vector<double> get_interp_weights(fVec3 pos, Coord_Type ctype) const;

		//trilinear interpolate a value from a 3d gridded data
		//in MPI reduce, it makes sure interpolation added together is correct
		template<typename T>
		T interp(const std::vector<T>& data, fVec3 pos, const Coord_Type ctype) const {
			iVec3 base = get_base_point(pos, ctype);
			T res{};

			//if all points fall inside
			if (leq_vec3(grid_p1, base) && leq_vec3(base + 2, grid_p2)) {
				
				auto weights = get_interp_weights((pos - base) / 2);
				size_t base_index = get_index_from_coord(base);
				size_t weight_index = 0;

				for(int k = 0; k < 4; k += 2)
				for(int j = 0; j < 4; j += 2)
				for(int i = 0; i < 4; i += 2) {
					res += data[base_index + get_index_offset(i, j, k)] * weights[weight_index];
					++weight_index;
				}
			}
			//if some points intersect, add range checking
			else if (is_intersect(base, base + 2, grid_p1, grid_p2)) {
				auto weights = get_interp_weights((pos - base) / 2);
				size_t weight_index = 0;

				for(int k = 0; k < 4; k += 2)
				for(int j = 0; j < 4; j += 2)
				for(int i = 0; i < 4; i += 2) {
					auto p = base + iVec3{i, j, k};
					if (is_inside(p))
						res += data[get_index_from_coord(p)] * weights[weight_index];
					++weight_index;
				}
			}
		
			return res;
		}
    };

	//output overloading
	std::ostream& operator<<(std::ostream& os, const std::complex<double>& c);
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec3<T>& c) {
		os << c.x << " " << c.y << " " << c.z;
		return os;
	}

	//scalar functions depending on positions
	template<typename T>
	using Scalar_Function = std::function<T(const fVec3&)>;

	//Material functions for use in user-defined geometry
	using Material_Function = Scalar_Function<Medium>;

	//Gaussian first derivative
	double Gaussian1(double t, double width);

	//Gaussian second derivative
	double Gaussian2(double t, double width);
}
