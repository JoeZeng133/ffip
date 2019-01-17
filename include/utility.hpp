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
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <numeric>
#include <future>

namespace ffip {
	/* Exception*/
	struct Invalid_Corner_Points {};
	struct Invalid_Direction {};
	struct Chunk_is_Not_Initialized{};
	struct Invalide_Coord_Type{};
	struct Out_of_the_Domain{};
	
	/* type aliases */
    using real = double;
	using real_arr = std::vector<real>;
	using complex_pair = std::pair<std::complex<real>, std::complex<real>>;
	using complex_num = std::complex<real>;
	using complex_arr = std::vector<complex_num>;

    /* constants */
	extern const real e0, u0, z0, c0;
	extern const real pi;

	/* enums with implied bit pattern for calculations*/
	enum Direction {X = 0, Y = 1, Z = 2};
	enum Side {Low = -1, High = 1};
	enum Coord_Type {Ex = 0b001, Ey = 0b010, Ez = 0b100, Hx = 0b110, Hy = 0b101, Hz = 0b011, Corner = 0b000, Center = 0b111, All = -2, Null = -1};
	
	int Ctype2DirInt(const Coord_Type ctype);
	
	constexpr bool is_E_point(const Coord_Type ctype) {
		return (ctype == Ex || ctype == Ey || ctype == Ez);
	}
	
	constexpr bool is_M_point(const Coord_Type ctype) {
		return (ctype == Hx || ctype == Hy || ctype == Hz);
	}
	/* forward declrations */
	template<typename T> struct Vec3;

    /* tags */
	struct side_low_tag;
	struct side_high_tag;
	
	struct dir_x_tag;
	struct dir_y_tag;
	struct dir_z_tag;
	
	struct e_tag;
	struct h_tag;
	
	struct hx_tag;
	struct hy_tag;
	struct hz_tag;
	struct ex_tag;
	struct ey_tag;
	struct ez_tag;
	
	/* lower than something */
    struct side_low_tag{
		using opposite_side = side_high_tag;
		static const int val;
		static int get_next_int(int x);
		static int round(real x);
		
		template<typename T>
		static T& choose_between(T& a, T& b) {
			return a;
		}
	};

	/* higher than sth */
    struct side_high_tag {
		using opposite_side = side_low_tag;
		static const int val;
		static int get_next_int(int x);
		static int round(real x);
		
		template<typename T>
		static T& choose_between(T& a, T& b) {
			return b;
		}
	};

	/* odd number, val = 1 */
    struct odd_tag{
		static const int val;
	};
	
	/* even number, val = 0 */
    struct even_tag{
		static const int val;
	};
	
	/* x direction, x1(y), x2(z), x3(x) */
	struct dir_x_tag{
		using x1 = dir_y_tag;
		using x2 = dir_z_tag;
		
		using x1_a = dir_y_tag;
		using x2_a = dir_z_tag;
		
		using z = dir_y_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
		
	};
	
	/* y direction, x1(z), x2(x), x3(y)*/
	struct dir_y_tag{
		using x1 = dir_z_tag;
		using x2 = dir_x_tag;
		
		using x1_a = dir_x_tag;
		using x2_a = dir_z_tag;
		
		using z = dir_x_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	/* z direction, x1(x), x2(y), x3(z)*/
	struct dir_z_tag{
		using x1 = dir_x_tag;
		using x2 = dir_y_tag;
		
		using x1_a = dir_x_tag;
		using x2_a = dir_y_tag;
		
		using z = dir_z_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	struct e_tag {};
	struct h_tag {};
	
	struct hx_tag : public dir_x_tag, public h_tag {
		using dir_base = dir_x_tag;
		using f_base = h_tag;
		static const Coord_Type ctype;
	};
	
	struct hy_tag : public dir_y_tag, public h_tag {
		using dir_base = dir_y_tag;
		using f_base = h_tag;
		static const Coord_Type ctype;
	};
	
	struct hz_tag : public dir_z_tag, public h_tag {
		using dir_base = dir_z_tag;
		using f_base = h_tag;
		static const Coord_Type ctype;
	};
	
	struct ex_tag : public dir_x_tag, public e_tag{
		using dir_base = dir_x_tag;
		using f_base = e_tag;
		static const Coord_Type ctype;
	};
	
	struct ey_tag : public dir_y_tag, public e_tag {
		using dir_base = dir_y_tag;
		using f_base = e_tag;
		static const Coord_Type ctype;
	};
	
	struct ez_tag : public dir_z_tag, public e_tag {
		using dir_base = dir_z_tag;
		using f_base = e_tag;
		static const Coord_Type ctype;
	};
	
	template<typename D>
	struct dir_traits {
		//orientation preserving
		using x1 = typename D::x1;
		using x2 = typename D::x2;
		using x3 = D;
		using z = typename D::z;		//rotate_frame<D::z> (rotate_frame<D>() ) = Identity map
		
		//array arrangement preserving
		using x1_a = typename D::x1_a;
		using x2_a = typename D::x2_a;
	};
	
	template<typename T = real>
	struct Vec3{
		using value_type = T;
		/* data members */
		T x{}, y{}, z{};

		/* Semiregular members*/
		Vec3() = default;
		Vec3(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}
		
		template<typename T2>
		explicit Vec3(const Vec3<T2>&);			//copy construction
		
		template<typename T2>
		Vec3& operator=(const Vec3<T2>&) ;		//copy assignment

		/* function members */
		Coord_Type get_type() const;						//get Coord_Type for the given comp coordinates
		Coord_Type get_type(const Coord_Type other) const;	//get Coord_Type for the given relative comp coordinates
	};
	
	using iVec3 = Vec3<int>;
	using fVec3 = Vec3<real>;
	using cVec3 = Vec3<complex_num>;
	using sVec3 = Vec3<size_t>;
	
	template<typename T>
		template<typename T2>
	Vec3<T>::Vec3(const Vec3<T2>& other): x(other.x), y(other.y), z(other.z) {}
	
	template<typename T>
		template<typename T2>
	Vec3<T>& Vec3<T>::operator=(const Vec3<T2>& other) {
		return Vec3<T>{other};
	}
	
	/* point wise less */
	template<typename T1, typename T2>
	constexpr bool ElementWise_Less(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x < b.x && a.y < b.y && a.z < b.z;
	}
	
	/* point wise less than*/
	template<typename T1, typename T2>
	constexpr bool ElementWise_Less_Eq(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x <= b.x && a.y <= b.y && a.z <= b.z;
	}
	
	/* testing whether a point is inside of a box*/
	template<typename T1, typename T2>
	constexpr bool Is_Inside_Box(const Vec3<T1>& p1, const Vec3<T1>& p2, const Vec3<T2>& tp) {
		return (ElementWise_Less_Eq(p1, tp) && ElementWise_Less_Eq(tp, p2));
	}
	
	/* dot product */
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} * T2{})> operator*(const Vec3<T1>& a, const Vec3<T2>& b) {
		return {a.x * b.x, a.y * b.y, a.z * b.z};
	}
	
	/* inner product*/
	template<typename T1, typename T2>
	constexpr decltype(T1{} * T2{}) inner_prod(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x * b.x + a.y * b.y +a.z * b.z;
	}
	
	/* right scalar product */
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} * T2{})> operator*(const Vec3<T1>& a, const T2 b) {
		return {a.x * b, a.y * b, a.z * b};
	}
	/* left scalar product*/
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} * T2{})> operator*(const T1 b, const Vec3<T2>& a) {
		return {a.x * b, a.y * b, a.z * b};
	}
	
	/* point wise addtion */
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} + T2{})> operator+(const Vec3<T1>& a, const Vec3<T2>& b) {
		return {a.x + b.x, a.y + b.y, a.z + b.z};
	}
	
	/* point wise substraction */
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} - T2{})> operator-(const Vec3<T1>& a, const Vec3<T2>& b) {
		return {a.x - b.x, a.y - b.y, a.z - b.z};
	}
	
	/* right scaler division*/
	template<typename T1, typename T2>
	constexpr Vec3<decltype(T1{} / T2{})> operator/(const Vec3<T1>& a, const T2 b) {
		return {a.x / b, a.y / b, a.z / b};
	}
	
	/* make dir_x_tag the x3 axis */
	template<typename T>
	constexpr Vec3<T> rotate_frame(const Vec3<T>& p, dir_x_tag) {
		return {p.y, p.z, p.x};
	}
	
	/* make dir_y_tag the x3 axis */
	template<typename T>
	constexpr Vec3<T> rotate_frame(const Vec3<T>& p, dir_y_tag) {
		return {p.z, p.x, p.y};
	}
	
	/* make dir_z_tag the x3 axis */
	template<typename T>
	constexpr Vec3<T> rotate_frame(const Vec3<T>& p, dir_z_tag) {
		return p;
	}
	
	/* return the nearest integer (self included) that is even or odd*/
	template<typename S, typename T>
	constexpr int get_nearest_int(const int x) {
		return ((x & 1) ^ T::val)? S::get_next_int(x) : x;
	}
	
	template<typename S, typename T>
	constexpr int get_nearest_int(const real x) {
		return get_nearest_int<S, T>(S::round(x));
	}
	
	/* helper empty struct for accessing x,y,z elements of sth generically*/
	template<typename D>
	struct choose {
		
		template<typename T>
		static typename T::value_type& get(T&);
		template<typename T>
		static typename T::value_type const& get(const T&);
	};
	
	template<>
	template<typename T>
	typename T::value_type& choose<dir_x_tag>::get(T& v) {
		return v.x;
	}
	
	template<>
	template<typename T>
	typename T::value_type& choose<dir_y_tag>::get(T& v) {
		return v.y;
	}
	
	template<>
	template<typename T>
	typename T::value_type& choose<dir_z_tag>::get(T& v) {
		return v.z;
	}
	
	template<>
	template<typename T>
	typename T::value_type const& choose<dir_x_tag>::get(const T& v) {
		return v.x;
	}
	
	template<>
	template<typename T>
	typename T::value_type const& choose<dir_y_tag>::get(const T& v) {
		return v.y;
	}
	
	template<>
	template<typename T>
	typename T::value_type const& choose<dir_z_tag>::get(const T& v) {
		return v.z;
	}
	
	/* get the nearest point of coord_type to a global computation point on side S (low or high)*/
	template<typename S, typename T>
	constexpr iVec3 get_nearest_point(const Vec3<T>& p, const Coord_Type type) {
		iVec3 res;
		if (type & 1) {
			res.x = get_nearest_int<S, odd_tag>(p.x);
		} else {
			res.x = get_nearest_int<S, even_tag>(p.x);
		}

		if (type & 2) {
			res.y = get_nearest_int<S, odd_tag>(p.y);
		} else {
			res.y = get_nearest_int<S, even_tag>(p.y);
		}

		if (type & 4) {
			res.z = get_nearest_int<S, odd_tag>(p.z);
		} else {
			res.z = get_nearest_int<S, even_tag>(p.z);
		}

		return res;
	}

	/* p1, p2 are computational coordinates, T1 , T2 in {floating types  and integer types}
	 return the corner points of fields of coord_type type inside [p1, p2]
	 if coord_type is null, return fields of all types
	 */
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_interior(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			return {iVec3(ceil(p1.x), ceil(p1.y), ceil(p1.z)), iVec3(floor(p2.x), floor(p2.y), floor(p2.z))};

		return std::make_pair(	get_nearest_point<side_high_tag>(p1, type),
								get_nearest_point<side_low_tag>(p2, type));
	}

	/* return the smallest box*/
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_closure(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			return {iVec3(floor(p1.x), floor(p1.y), floor(p1.z)), iVec3(ceil(p2.x), ceil(p2.y), ceil(p2.z))};

		return std::make_pair(	get_nearest_point<side_low_tag>(p1, type),
								get_nearest_point<side_high_tag>(p2, type));
	}

	/* return the intersection of two boxes specified by (p1, p2) and (q1, q2)*/
	template<typename T>
	std::pair<Vec3<T>, Vec3<T>> get_intersection(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {
		return std::make_pair<Vec3<T>, Vec3<T>>({std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
									{std::min(p2.x, q2.x),std::min(p2.y, q2.y), std::min(p2.z, q2.z)});
	}
	
	/* return a particular face of a box region specified by (p1, p2)*/
	template<typename D, typename S, typename T>
	std::pair<Vec3<T>, Vec3<T>> get_face(Vec3<T> p1, Vec3<T> p2) {
		choose<D>::get(S::choose_between(p2, p1)) = choose<D>::get(S::choose_between(p1, p2));
		return {p1, p2};
	}

    /* source functions */
	auto make_gaussian_func(real sigma_t, real d = 0) -> std::function<real(const real)>;
	auto make_sin_func(real freq, real d = 0) -> std::function<real(const real)>;
	auto make_ricker_func(real fp, real d = 0) -> std::function<real(const real)>;
	
	/* fast 1 dimensional linear interpolation */
	class interp1 {
	public:
		size_t dimx;
		interp1(const size_t _dimx): dimx(_dimx) {
			if (dimx <= 1)
				throw std::runtime_error("Invalid Side");
		}
		
		template<typename T>
		T get(const std::vector<T>& data, const real x) {
			if (data.size()!= dimx)
				throw std::runtime_error("Invalid Data Size");
			
			if (x < 0 || x >= dimx - 1)
				throw Out_of_the_Domain{};
			
			size_t x_b = floor(x);
			return data[x_b] * (x_b + 1 - x) + data[x_b + 1] * (x - x_b);
		}
	};
	
	
	/* fast 2 dimensional linear interpolation */
	class interp2 {
	public:
		size_t dimx, dimy;
		size_t size, jump_y;
		
		interp2(const size_t _dimx, const size_t _dimy): dimx(_dimx), dimy(_dimy) {
			if (dimx <= 1 || dimy <= 1)
				throw std::runtime_error("Invalid Size");
			
			size = dimx * dimy;
			jump_y = dimx;
		}
		
		template<typename T>
		T get(const std::vector<T>& data, const real y, const real x) {
			if (data.size()!= size)
				throw  std::runtime_error("Invalid Data Size");
			
			if (x < 0 || x > dimx - 1 || y < 0 || y > dimy - 1)
				throw Out_of_the_Domain{};
			
			size_t xb = floor(x);
			size_t yb = floor(y);
			real tx = x - xb;
			real ty = y - yb;
			T const* ptr = &data[xb + yb * jump_y];
			
			return (1 - tx) * (1 - ty) * ptr[0] +
			tx * (1 - ty) * ptr[1] +
			(1 - tx) * ty * ptr[jump_y] +
			tx * ty * ptr[1 + jump_y];
		}
	};
	
	/* fast 3 dimensional linear interpolation
	 */
	class interp3 {
	public:
		size_t dimx, dimy, dimz;
		size_t size, jump_y, jump_z;
		
		interp3(const size_t _dimx, const size_t _dimy, const size_t _dimz): dimx(_dimx), dimy(_dimy), dimz(_dimz) {
			if (dimx <= 1 || dimy <= 1 || dimz <= 1)
				throw std::runtime_error("Invalid Size");
			
			size = dimx * dimy * dimz;
			jump_y = dimx;
			jump_z = dimx * dimy;
		}
		
		template<typename T>
		T get(const std::vector<T>& data, const real z, const real y, const real x) {
			if (data.size()!= size)
				throw std::runtime_error("Invalid Data Size");
			
			if (x < 0 || x > dimx - 1 || y < 0 || y > dimy - 1|| z < 0 || z > dimz - 1)
				throw Out_of_the_Domain{};
			
			size_t xb = floor(x);
			size_t yb = floor(y);
			size_t zb = floor(z);
			real tx = x - xb;
			real ty = y - yb;
			real tz = z - zb;
			T const* ptr = &data[xb + yb * jump_y + zb * jump_z];
			
			const T& c000 = ptr[0];
			const T& c001 = ptr[1];
			const T& c010 = ptr[jump_y];
			const T& c011 = ptr[1 + jump_y];
			const T& c100 = ptr[jump_z];
			const T& c101 = ptr[jump_z + 1];
			const T& c110 = ptr[jump_z + jump_y];
			const T& c111 = ptr[jump_z + jump_y + 1];
			
			return (1 - tx) * (1 - ty) * (1 - tz) * c000 +
			tx * (1 - ty) * (1 - tz) * c001 +
			(1 - tx) * ty * (1 - tz) * c010 +
			tx * ty * (1 - tz) * c011 +
			(1 - tx) * (1 - ty) * tz * c100 +
			tx * (1 - ty) * tz * c101 +
			(1 - tx) * ty * tz * c110 +
			tx * ty * tz * c111;
		}
	};
	
	/* Generic Interpolation class
	 	two times slower than the fast one
	 */
#define tol_interp 1e-5
	
	template<int N>
	class interpn : public interpn<N - 1> {
	public:
		using base_class = interpn<N - 1>;
		size_t dimn{1}, jump{1};
		
		interpn() = default;
		
		template<typename... Args>
		interpn(const size_t _dimn, Args... args): base_class(args...), dimn(_dimn) {
			static_assert(sizeof...(args) == N - 1, "Invalid Number of Size");
			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");
			
			jump = base_class::get_size();
		}
		
		size_t get_size() const {
			return dimn * base_class::get_size();
		}
		
		template<typename T, typename... Args>
		T get_helper(T const* data, real xq, Args&&... args) const{
			if ( dimn == 1)		//ignore this dimension if it is 1
				return base_class::get_helper(data, args...);
			
			//force points to be inside of region
			if (xq > dimn - 1) return base_class::get_helper(data + jump * (dimn - 1), args...);
			if (xq < 0) return base_class::get_helper(data, args...);
			
			size_t index = xq;
			real tx = xq - index;
			
			if (tx < tol_interp)
				return base_class::get_helper(data + jump * index, args...);
			else
				return	tx * base_class::get_helper(data + jump * (index + 1), args...) +
				(1 - tx) * base_class::get_helper(data + jump * index, args...);
		}
		
		template<typename T, typename... Args>
		T get(const std::vector<T>& vec, Args&&... args) const{
			static_assert(sizeof...(args) == N, "Invalid Request");
			return get_helper(vec.data(), args...);
		}
		
		template<typename T, typename... Args>
		void put_helper(T* data, const T& val, real xq, Args&&... args) const{
			if ( dimn == 1) {	//ignore 1 dimension
				base_class::get_helper(data, args...);
				return;
			}
			
			//force points to be inside of region
			if (xq > dimn - 1) return base_class::put_helper(data + jump * (dimn - 1), val, args...);
			if (xq < 0) return base_class::put_helper(data, val, args...);
			
			size_t index = xq;
			real tx = xq - index;
			
			if (tx < tol_interp)
				base_class::put_helper(data + jump * index, val, args...);
			else {
				base_class::put_helper(data + jump * (index + 1), val * tx, args...);
				base_class::put_helper(data + jump * index, val * (1 - tx), args...);
			}
		}
		
		template<typename T, typename... Args>
		void put(std::vector<T>& vec, const T& val, Args&&... args) const {
			static_assert(sizeof...(args) == N, "Invalid Request");
			return put_helper(vec.data(), val, args...);
		}
	};
	
	template<>
	class interpn<1> {
	public:
		size_t dimn{1}, jump{1};
		
		interpn() = default;
		
		interpn(const size_t _dimn): dimn(_dimn) {
			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");
			
			jump = 1;
		}
		
		size_t get_size() const {
			return dimn;
		}
		
		template<typename T>
		T get(const std::vector<T>& vec, real xq) const{
			return get_helper(vec.data(), xq);
		}
		
		template<typename T>
		T get_helper(T const* data, real xq) const{
			if (dimn == 1)		//ignore this dimension if it is 1
				return data[0];
			
			if (xq > dimn - 1) return data[dimn - 1];
			if (xq < 0) return data[0];
			
			size_t index = xq;
			real tx = xq - index;
			
			if (tx < tol_interp)
				return data[index];
			else
				return tx * data[index + 1] + (1 - tx) * data[index];
		}
		
		template<typename T>
		void put(std::vector<T>& vec, const T& val, real xq) const{
			return put_helper(vec.data(), val, xq);
		}
		
		template<typename T>
		void put_helper(T* data, const T& val, real xq) const{
			if ( dimn == 1) {		//ignore this dimension if it is 1
				data[0] = val;
				return;
			}
			//force points to be inside of region
			if (xq > dimn - 1) {
				data[dimn - 1] = val;
				return;
			}
			
			if (xq < 0) {
				data[0] = val;
				return;
			}
			
			size_t index = xq;
			real tx = xq - index;
			
			if (tx < tol_interp)
				data[index] = val;
			else {
				data[index + 1] = val * tx;
				data[index] = val * (1 - tx);
			}
		}
	};

    /* linear interpolation on fixed interval grids*/
    class GriddedInterp {
    private:
		real_arr x, y, z, v;
		real dx, dy, dz;
        std::string file_name;
		interpn<3> interp{1, 1, 1};
		
    public:
        /* Semiregular members*/
        GriddedInterp(std::string file_name);
		GriddedInterp(const iVec3& dim, const fVec3& w0, const fVec3& dw, const real_arr& _v);
		
		GriddedInterp(const GriddedInterp&) = default;				//copy
		GriddedInterp& operator=(const GriddedInterp&) = default;
		GriddedInterp(GriddedInterp&&) = default;					//move
		GriddedInterp& operator=(GriddedInterp&&) = default;

        /* Function Members */
        real request_value(const fVec3& p) const;			//return value at a point, 0 dimension is ignored
        real request_integral() const;						//return the integral over the entire volume
		
		void expand_dim(const real lo, const real hi, const Direction dir);	//expand zero dimension without changing the integral, this is for use in current source
		void scale_xyz(const real sx, const real sy, const real sz);	//x *= sx, y *= sy, z *= sz, to scale all dimensions, caution: it changes the integral by sx * sy * sz
		
		void expand_dim_helper(std::vector<real>& w, const real lo, const real hi, const Direction dir);
        fVec3 get_p1() const;								//return the lower corner point
        fVec3 get_p2() const;								//return the upper corner point

		//public functions for calculating generic line, face and volume integral
    };
	
	
	
	/* n dimensional linear integration, 1 dimension,
	 xn does not need to have fixed interval*/
	template<typename T>
	T integral_ndim(const T* data, const std::vector<real> &xn) {
		if(xn.size() == 0)
			throw std::runtime_error("Invalid array size");
		
		// ignore zero dimension
		if(xn.size() == 1)
			return data[0];
		
		int size = xn.size();
		T res = 0.5 * (xn[1] - xn[0]) * data[0] + 0.5 * (xn[size - 1] - xn[size - 2]) * data[size-1];
		
		for(int i = 1; i < size - 1; ++i) {
			res += 0.5 * data[i] * (xn[i + 1] - xn[i - 1]);
		}
		
		return res;
	}
	
	/* empty struct to get size of an object*/
	template<typename T>
	struct size_trait {
		static size_t size(const T&) {return 1;};
	};
	
	/* return size if it is a vector */
	template<typename T>
	struct size_trait<std::vector<T>> {
		static size_t size(const std::vector<T>& v) {return v.size();}
	};
	
	/* get size of series of vectors, one vector */
	template<typename T>
	inline size_t get_vec_size(const T& vn) {
		return size_trait<T>::size(vn);
	}
	
	/* get size of series of vectors, template */
	template<typename T, typename... Args>
	inline size_t get_vec_size(const T& vn, Args... args) {
		return size_trait<T>::size(vn) * get_vec_size(args...);
	}
	
	/* n dimensional linear integration, template,
	 xn does not need to have fixed interval*/
	template<typename T, typename... Args>
	T integral_ndim(const T* data, const std::vector<real> &xn, Args&&... args) {
		if (xn.size() == 0)
			throw std::runtime_error("Invalid array size");
		
		// ignore zero dimension
		if (xn.size() == 1)
			return integral_ndim(data, args...);
		
		
		std::vector<T> tmp(xn.size());
		size_t shift = get_vec_size(args...);
		
		for(int i = 0; i < xn.size(); ++i) {
			tmp[i] = integral_ndim(data + shift * i, args...);
		}
		
		return integral_ndim(tmp.data(), xn);
	}
	
	/* PML class for calculating k, b, c arrays used in updating perfectly matched layer using CPML*/
	class PML {
	private:
		int d;
		real sigma_max;
		real k_max{1};
		real a_max{0};
		real m{3};
		real m_a{1};

	public:
		/* constructors*/
		PML() = default;
		PML(int d, real sigma_max);
		PML(int d, real sigma_max, real k_max, real a_max, real m, real m_a);

		PML(const PML&) = default;				//copy
		PML& operator=(const PML&) = default;
		PML(PML&&) = default;						//move
		PML& operator=(PML&&) = default;

		int get_d() const;
		real get_sigma(const real x) const;
		real get_a(const real x) const;
		
		real get_k(const real x) const;
		real get_b(const real x, const real dt) const;
		real get_c(const real x, const real dt) const;

		static real optimal_sigma_max(real m_k, real dx, real er = 1, real ur = 1);
	};
	
	/* extra field information for use in dispersive field updates
	stores e(n-2) and polarization p(n) p(n - 1)
	 */
	struct Dispersive_Field {
		union {real ex2, ey2, ez2, e2, hx2, hy2, hz2, h2, eh2;};
		std::vector<real> p1, p;
		Dispersive_Field(const size_t num_poles);
		size_t get_num_poles() const;
	};
	
	/* field information for points inside PML layers*/
	struct PML_Point {
		int index, jump_pos, jump_neg;
		real b_pos, b_neg, c_pos, c_neg;
		real psi_pos{0}, psi_neg{0};
		
		PML_Point(const int _index, const int _jump_pos, const int _jump_neg, const real _b_pos, const real _b_neg, const real _c_pos, const real _c_neg);
	};
	
	template<typename P, int N, int M>
	inline typename P::value_type get(const P& p);
	
	/* an iterator to iterate through a box region specified by two corner regions [p1, p2]
	 it allows iteration of particular coord_type points: ex, ..., ez, hx, ..., hz
	 for null, it will loop through all points
	 */
	struct my_iterator {
		using value_type = int;
		
		size_t size, index, end;
		int x0, y0, z0, x1, y1, z1;
		int x, y, z;
		int jump;
		
		/* non-parallel version constructor*/
		my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype);
		my_iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype);
		my_iterator(const fVec3& p1, const fVec3& p2, const Coord_Type ctype);
		/* used for paralellization so that the itr only loops through rank/num of the total points*/
		my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype, const size_t rank, const size_t num);
		
		void advance();					//increase index by 1
		bool is_end() const;			//does it reach the end?
		bool is_empty() const;			//is the region illegal?
		iVec3 get_vec() const;			//return the current x, y, z in iVec3
		iVec3 get_vec(size_t index) const;	//return x, y, z specified by index
		size_t get_index(const int x, const int y, const int z) const;	//return the index specified by x, y, z
		size_t get_size() const;			//return size of the region
		Vec3<size_t> get_dim() const;		//return the dimension of the region
	};
	
	/* Barrier implementation with reset at the end*/
	class Barrier
	{
	private:
		std::mutex m_mutex;
		std::condition_variable m_cv;

		size_t m_count;
		const size_t m_initial;

		enum State : unsigned char {
			Up, Down
		};
		State m_state;
		   
	public:
		explicit Barrier(std::size_t count);
		Barrier& operator=(Barrier&&) = delete;
		Barrier& operator=(const Barrier&) = delete;
		/// Blocks until all N threads reach here
		void Sync();
		size_t get_num_proc() const;
	};

	/* Global Barrier for use in anywhere*/
	extern Barrier* glob_barrier;
	void set_num_proc(const size_t num_proc);

	// divide a list evenly
	template<typename T>
	inline void vector_divider(const std::vector<T>& list, const size_t rank, const size_t num_proc, size_t& idx1, size_t& idx2) {
		if (num_proc > list.size()) {
			if (rank >= list.size()) {
				idx1 = idx2 = 0;
			}
			else {
				idx1 = rank;
				idx2 = rank + 1;
			}
		}
		else {
			idx1 = rank * list.size() / num_proc;
			idx2 = (rank + 1) * list.size() / num_proc;
		}
	}
	extern iVec3 vec3_base[3];
	
	//output overloading
	std::ostream& operator<<(std::ostream& os, const complex_num& c);
	
	template<typename T>
	std::ostream& operator<<(std::ostream& os, const Vec3<T>& c) {
		os << c.x << " " << c.y << " " << c.z;
		return os;
	}

	// to pass around classes and stuff
	struct Config {
		real dt, dx;
		iVec3 sim_p1, sim_p2;
		iVec3 ch_p1, ch_p2;
		iVec3 tf_p1, tf_p2;
		iVec3 phys_p1, phys_p2;
	};
}
