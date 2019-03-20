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
	/* Exception*/
	struct Invalid_Corner_Points {};
	struct Invalid_Direction {};
	struct Chunk_is_Not_Initialized{};
	struct Invalide_Coord_Type{};
	struct Out_of_the_Domain{};
	
	/* type aliases */
	using double_arr = std::vector<double>;
	using complex_pair = std::pair<std::complex<double>, std::complex<double>>;
	using complex_num = std::complex<double>;
	using complex_arr = std::vector<complex_num>;

	extern const double pi;

	/* enums with implied bit pattern for calculations*/
	enum Direction {X = 0, Y = 1, Z = 2};
	enum Side {Low = -1, High = 1};
	enum Coord_Type {Ex = 0b001, Ey = 0b010, Ez = 0b100, Hx = 0b110, Hy = 0b101, Hz = 0b011, 
		Dx = 0b1001, Dy = 0b1010, Dz = 0b1100, Bx = 0b1110, By = 0b1101, Bz = 0b1011,
		Corner = 0b000, Center = 0b111, All = 0b1000};

	//count the number of 1s in the first 4 bits
	constexpr unsigned int count_hex_bits(int x) {

		return (x & 1) + ((x & 2) >> 1) + ((x & 4) >> 2) + ((x & 8) >> 3);
	}

	//get direction integer from coord_type
	constexpr int get_dir_from_ctype(Coord_Type ctype) {

		ctype = Coord_Type(ctype & 0b111);
		if (ctype == Ex || ctype == Hx) return 0;
		if (ctype == Ey || ctype == Hy) return 1;
		if (ctype == Ez || ctype == Hz) return 2;

		return -1;
	}

	//return whether it is electrical point
	constexpr bool is_e_point(const Coord_Type ctype) {

		return (ctype == Ex || ctype == Ey || ctype == Ez);
	}
	
	//return wehether it is magnetic point
	constexpr bool is_m_point(const Coord_Type ctype) {

		return (ctype == Hx || ctype == Hy || ctype == Hz);
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

		/* arithmetic operations */
		template<typename T2>
		std::comon_type_t<T, T2> dot(const Vec3<T2>& other) const {
			return x * other.x + y * other.y + z * other.z;
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const T2 mult) const {
			return {x * mult, y * mult, z * mult};
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const Vec3<T2>& other) const {
			return {x * other.x, y * other.y, z * other.z};
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const T2 add) const {
			return {x + add, y + add, z + add};
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const Vec3<T2>& other) const {
			return {x + other.x, y + other.y, z + other.z};
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const T2 minus) const {
			return {x - minus, y - minus, z - minus};
		}

		template<typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const Vec3<T2>& other) const {
			return {x - other.x, y - other.y, z - other.z};
		}
		
		Vec3 ceil() const {
			return { std::ceil(x), std::ceil(y), std::ceil(z) };
		}

		Vec3 floor() const {
			return { std::floor(x), std::floor(y), std::floor(z) };
		}

		Vec3 abs() const {
			return { std::abs(x),  std::abs(y), std::abs(z) };
		}

		T prod() const {
			return x * y * z;	
		}

		/* generic opeartions */
		template<int N>
		T& get() {
			if constexpr(N == 0)
				return x;
			else if constexpr(N == 1)
				return y;
			else if constexpr(N == 2)
				return z;
		}

		template<int N>
		const T& get() const {
			if constexpr(N == 0)
				return x;
			else if constexpr(N == 1)
				return y;
			else if constexpr(N == 2)
				return z;
		}

		T& operator[](const size_t n) {
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		const T& operator[](const size_t n) {
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		//get Coord_Type for the given comp coordinates
		template<typename Dummy = T, typename = std::enable_if_t<std::is_integral<Dummy>::value>>
		Coord_Type get_type() const {
			return static_cast<Coord_Type>((x & 1) | ((y & 1) << 1) | ((z & 1) << 2));
		}
	};

	template<typename T1, typename T2>
	bool le_vec3(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x < b.x && a.y < b.y && a.z < b.z;
	}

	template<typename T1, typename T2>
	bool leq_vec3(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x <= b.x && a.y <= b.y && a.z <= b.z;
	}

	using iVec3 = Vec3<int>;			
	using fVec3 = Vec3<double>;	
	using cVec3 = Vec3<complex_num>;
	using sVec3 = Vec3<size_t>;

	template<typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
	constexpr Vec3<T> dim2stride(const Vec3<T>& dim) {
		return { 1, dim.x, dim.x * dim.y }; 
	}

	/* testing whether a point is inside of a box*/
	template<typename T1, typename T2>
	constexpr bool Is_Inside_Box(const Vec3<T1>& p1, const Vec3<T1>& p2, const Vec3<T2>& tp, const T2 tol = T2{}) {
		return (leq_vec3(p1, tp + tol) && leq_vec3(tp - tol, p2));
	}
	
	/* make xN the x3 axis */
	template<int N, typename T>
	constexpr Vec3<T> rotate_frame(const Vec3<T>& p) {
		if constexpr(N == 0)
			return {p.y, p.z, p.x};
		else if constexpr(N == 1)
			return {p.z, p.x, p.y};
		else if constexpr(N == 2)
			return p;
	}
	
	/* return the nearest integer (self included) that is even or odd*/
	template<int S, unsigned int T>
	constexpr int get_nearest_int(const int x) {
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");
		
		return (x & 1) ^ T? x + S : x;
	}
	
	template<int S, unsigned int T>
	constexpr int get_nearest_int(const double x) {
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");
		
		if constexpr(S == -1)
			return get_nearest_int<S, T>(std::floor(x));
		else
			return get_nearest_int<S, T>(std::ceil(x));
	}

	/* get the nearest point of coord_type to a global computation point on side S (low or high)*/
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

	/* p1, p2 are computational coordinates, T1 , T2 in {floating types  and integer types}
	 return the corner points of fields of coord_type type inside [p1, p2]
	 if coord_type is null, return fields of all types
	 */
	template<typename T1, typename T2>
    constexpr std::pair<iVec3, iVec3> 
		get_component_interior(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {

		if (type == Coord_Type::All)
			return {iVec3(p1.ceil()), iVec3(p2.floor())};

		return std::make_pair(	get_nearest_point<1>(p1, type), get_nearest_point<-1>(p2, type));
	}

	/* return the smallest box*/
	template<typename T1, typename T2>
	constexpr std::pair<iVec3, iVec3> 
	get_component_closure(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {

		if (type == Coord_Type::All)
			return {iVec3(p1.floor()), iVec3(p2.ceil())};

		return std::make_pair(	get_nearest_point<-1>(p1, type), get_nearest_point<1>(p2, type));
	}

	/* return the intersection of two boxes specified by (p1, p2) and (q1, q2)*/
	template<typename T>
	constexpr std::pair<Vec3<T>, Vec3<T>> get_intersection(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {
		return	{{std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
				{std::min(p2.x, q2.x),std::min(p2.y, q2.y), std::min(p2.z, q2.z)}};
	}
	
	/* return a particular face of a box region specified by (p1, p2)*/
	template<unsigned int D, int S, typename T>
	constexpr std::pair<Vec3<T>, Vec3<T>> 
		get_face(Vec3<T> p1, Vec3<T> p2) {
		
		static_assert(S == -1 || S == 1, "S must be -1 or 1");
		if constexpr(S == 1) {
			p1.get<D>() = p2.get<D>();
		}
		else {
			p2.get<D>() = p1.get<D>();
		}
		
		return {p1, p2};
	}

	/* check whether a value is in closed interval [low, high]*/
	template <typename T>
	bool In_ClosedBounds(const T& value, const T& low, const T& high) {
		return !(value < low) && !(high < value);
	}

#define tol_interp 1e-5
	
	template<int N>
	class interpn : public interpn<N - 1> {
	public:
		using base_class = interpn<N - 1>;
		size_t dimn{1}, stride{1};
		
		interpn() = default;
		
		template<typename... Args>
		interpn(const size_t _dimn, Args... args): base_class(args...), dimn(_dimn) {
			static_assert(sizeof...(args) == N - 1, "Invalid Number of Size");
			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");
			
			stride = base_class::get_size();
		}
		
		size_t get_size() const {
			return dimn * base_class::get_size();
		}
		
		template<typename T, typename... Args>
		T get_helper(T const* data, double xq, Args&&... args) const{
			if ( dimn == 1)		//ignore this dimension if it is 1
				return base_class::get_helper(data, args...);
			
			//force points to be inside of region
			if (xq > dimn - 1) return base_class::get_helper(data + stride * (dimn - 1), args...);
			if (xq < 0) return base_class::get_helper(data, args...);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (tx < tol_interp)
				return base_class::get_helper(data + stride * index, args...);
			else
				return	tx * base_class::get_helper(data + stride * (index + 1), args...) +
				(1 - tx) * base_class::get_helper(data + stride * index, args...);
		}
		
		template<typename T, typename... Args>
		T get(const std::vector<T>& vec, Args&&... args) const{
			static_assert(sizeof...(args) == N, "Invalid Request");
			return get_helper(vec.data(), args...);
		}
		
		template<typename T, typename... Args>
		void put_helper(T* data, const T& val, double xq, Args&&... args) const{
			if ( dimn == 1) {	//ignore 1 dimension
				base_class::get_helper(data, args...);
				return;
			}
			
			//clamping coordinates
			if (xq > dimn - 1) return base_class::put_helper(data + stride * (dimn - 1), val, args...);
			if (xq < 0) return base_class::put_helper(data, val, args...);
			
			size_t index = xq;
			double tx = xq - index;
			
			if (tx < tol_interp)
				base_class::put_helper(data + stride * index, val, args...);
			else {
				base_class::put_helper(data + stride * (index + 1), val * tx, args...);
				base_class::put_helper(data + stride * index, val * (1 - tx), args...);
			}
		}
		
		// the put function might skip placing 0s in some elements of the array
		template<typename T, typename... Args>
		void put(std::vector<T>& vec, const T& val, Args&&... args) const {
			static_assert(sizeof...(args) == N, "Invalid Request");
			return put_helper(vec.data(), val, args...);
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
		T get(const std::vector<T>& vec, double xq) const{
			return get_helper(vec.data(), xq);
		}
		
		template<typename T>
		T get_helper(T const* data, double xq) const{
			if (dimn == 1)		//ignore this dimension if it is 1
				return data[0];
			
			if (xq > dimn - 1) return data[dimn - 1];
			if (xq < 0) return data[0];
			
			size_t index = xq;
			double tx = xq - index;
			
			if (tx < tol_interp)
				return data[index];
			else
				return tx * data[index + 1] + (1 - tx) * data[index];
		}
		
		template<typename T>
		void put(std::vector<T>& vec, const T& val, double xq) const{
			return put_helper(vec.data(), val, xq);
		}
		
		template<typename T>
		void put_helper(T* data, const T& val, double xq) const{
			if ( dimn == 1) {		//ignore this dimension if it is 1
				data[0] = val;
				return;
			}
			//clamping coordinates
			if (xq > dimn - 1) {
				data[dimn - 1] = val;
				return;
			}
			
			if (xq < 0) {
				data[0] = val;
				return;
			}
			
			size_t index = xq;
			double tx = xq - index;
			
			if (tx < tol_interp)
				data[index] = val;
			else {
				data[index + 1] = val * tx;
				data[index] = val * (1 - tx);
			}
		}
	};

	/*
		Iterator through grid coordinates
	*/
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

	/*
		3D Yee-Cells data
		for dipole source creation and fields interpolation
	*/
    struct Yee3 {
		
		//physical coordinates
        fVec3 phys_p1, phys_p2;
		//grid coordinates
        iVec3 grid_p1, grid_p2;
		//grid distance times 2
        double dx;
		//dimensions of the array
		iVec3 dim;
		//strides of the array
		sVec3 stride;
		//[grid_p1 - tol,grid_p2 + tol] will be deemed as inside the cell
		double tol = 1e-3;

		Yee3(iVec3 grid_p1, iVec3 grid_p2, fVec3 phys_p1, double dx);

		//get grid corner points
		iVec3 get_grid_p1() const;
		iVec3 get_grid_p2() const;

		//get coordinate from a index, throws error if outside
		iVec3 get_coord_from_index(size_t n) const;

		//get index for a given point, undefined if pos is outside
		size_t get_index_from_coord(const iVec3& pos) const;
		size_t get_index_from_coord(int i, int j, int k) const;

		//get index offset at (i, j, k) relative to lower corner points
		size_t get_index_offset(const iVec3& pos) const;
		size_t get_index_offset(int i, int j, int k) const;

		//get corresponding dimension
		template<int N>
		size_t get_dim() const {
			return dim.get<N>();
		}

		//return interpolation details
		std::vector<std::pair<size_t, double>> get_interp_info(bool is_grid_coordinate, fVec3 pos, Coord_Type ctype);

		//access the raw value from a grid coordinate, undefined if not insdie the grid
		template<typename T>
		T get_raw_val(const std::vector<T>& data, const iVec3& pos) const {
			return data[(pos - grid_p1).dot(stride)];
		}

		//determine whether a point is inside a grid
		bool is_in_grid(bool is_grid_coordinate, const fVec3& pos) const;
		
		//convert a physical coordinate to a grid coordinate
		fVec3 phys_coord2grid_coord(const fVec3& pos) const;

		//transpose value into a 3d gridded data, return whether successful
        template<typename T>
        bool transpose_interp
		(std::vector<T>& data, bool is_grid_coordinate, fVec3 pos, T amp, Coord_Type ctype) const {

			if (!is_grid_coordinate)
				pos = (pos - phys_p1) * (2 / dx);

			if(!Is_Inside_Box(grid_p1, grid_p2, pos))
				return false;

			iVec3 base = get_nearest_point<-1>(pos, ctype);

			double sx = (pos.x - base.x) / 2;
			double sy = (pos.y - base.y) / 2;
			double sz = (pos.z - base.z) / 2;
			T* base_ptr = (base - grid_p1).dot(stride) + data.data();

			base_ptr[get_index_offset(0, 0, 0)] += (1 - sx) * (1 - sy) * (1 - sz)	* amp;
			base_ptr[get_index_offset(2, 0, 0)] += sx * (1 - sy) * (1 - sz)			* amp;
			base_ptr[get_index_offset(0, 2, 0)] += (1 - sx) * sy * (1 - sz)			* amp;
			base_ptr[get_index_offset(2, 2, 0)] += sx * sy * (1 - sz)				* amp;
			base_ptr[get_index_offset(0, 0, 2)] += (1 - sx) * (1 - sy) * sz			* amp;
			base_ptr[get_index_offset(2, 0, 2)] += sx * (1 - sy) * sz				* amp;
			base_ptr[get_index_offset(0, 2, 2)] += (1 - sx) * sy * sz 				* amp;
			base_ptr[get_index_offset(2, 2, 2)] += sx * sy * sz						* amp;

			return true;
		}

		//trilinear interpolate a value from a 3d gridded data, undefined if not inside the grid
		template<typename T>
		T interp
		(const std::vector<T>& data, bool is_grid_coordinate, const fVec3& pos, const Coord_Type ctype) const {
			if (!is_grid_coordinate)
				pos = (pos - phys_p1) * (2 / dx);

			if(!Is_Inside_Box(grid_p1, grid_p2, pos))
				return T{};

			iVec3 base = get_nearest_point<-1>(pos, ctype);

			double sx = (pos.x - base.x) / 2;
			double sy = (pos.y - base.y) / 2;
			double sz = (pos.z - base.z) / 2;
			T* base_ptr = (base - grid_p1).dot(stride) + data.data();

			return (1 - sx) * (1 - sy) * (1 - sz) * base_ptr[get_index_offset(0, 0, 0)] + 
					sx * (1 - sy) * (1 - sz) * 		base_ptr[get_index_offset(2, 0, 0)] +
					(1 - sx) * sy * (1 - sz) * 		base_ptr[get_index_offset(0, 2, 0)] +
					sx * sy * (1 - sz) * 			base_ptr[get_index_offset(2, 2, 0)] +
					(1 - sx) * (1 - sy) * sz * 		base_ptr[get_index_offset(0, 0, 2)] +
					sx * (1 - sy) * sz * 			base_ptr[get_index_offset(2, 0, 2)] +
					(1 - sx) * sy * sz * 			base_ptr[get_index_offset(2, 2, 2)] +
					sx * sy * sz * 					base_ptr[get_index_offset(2, 2, 2)];
		}
    };

	//output overloading
	std::ostream& operator<<(std::ostream& os, const std::complex<double>& c);
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec3<T>& c) {
		os << c.x << " " << c.y << " " << c.z;
		return os;
	}

	// to pass around classes and stuff
	struct Config {
		double dt, dx;			//discretized length and time
		iVec3 roi_p1, roi_p2;	//region of interest, 0->dim * 2
		iVec3 sim_p1, sim_p2;	//region of simulation, everything included
		iVec3 ch_p1, ch_p2;		//region of chunk, collection of chunks is a disjoint cover of region of simulation
		iVec3 tf_p1, tf_p2;		//region of total field
		iVec3 phys_p1, phys_p2;	//region of physical field, region of simulation minus PML layer
	};

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
