#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <boost/math/constants/constants.hpp>

namespace ffip {
    using real = double;
	using real_arr = std::vector<double>;
	using complex_pair = std::pair<std::complex<real>, std::complex<real>>;
	using complex_num = std::complex<real>;
	using complex_arr = std::vector<complex_num>;

    /* constants */
	extern const real e0, u0, z0, c;
	extern const real pi;

	/* some enums */
	enum Direction {X = 0, Y = 1, Z = 2};
	enum Side {Low = -1, High = 1};
	enum Coord_Type {Ex = 0b001, Ey = 0b010, Ez = 0b100, Hx = 0b110, Hy = 0b101, Hz = 0b011, Corner = 0b000, Center = 0b111, Null = -1};

    /* tags */
	struct dir_x_tag;
	struct dir_y_tag;
	struct dir_z_tag;
	
    struct side_low_tag{
		static const int val;
		static int get_next_int(int x);
		static int round(real x);
		static int round_wtol(real x);
	};

    struct side_high_tag{
		static const int val;
		static int get_next_int(int x);
		static int round(real x);
		static int round_wtol(real x);
	};

    struct odd_tag{
		static const int val;
	};
    struct even_tag{
		static const int val;
	};
	
	struct dir_x_tag{
		using x1 = dir_y_tag;
		using x2 = dir_z_tag;
		using z = dir_y_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	struct dir_y_tag{
		using x1 = dir_z_tag;
		using x2 = dir_x_tag;
		using z = dir_x_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	struct dir_z_tag{
		using x1 = dir_x_tag;
		using x2 = dir_y_tag;
		using z = dir_z_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	
	/*rotational symmetry of cartesian coordinates (x, y, z)
	(y, z, x) or (z, x, y) or (x, y, z)
	*/
	template<typename D>
	struct dir_traits {
		using x1 = typename D::x1;
		using x2 = typename D::x2;
		using x3 = D;
		using z = typename D::z;		//rotate_frame<D::z> (rotate_frame<D>() ) = Identity map
		
	};
	
	template<typename T = real>
	struct Vec3{
		using value_type = T;
		/* data members */
		T x, y, z;

		/* Semiregular members*/
		Vec3() = default;
		Vec3(T _x, T _y, T _z): x(_x), y(_y), z(_z) {}
		Vec3(const Vec3&) = default;			//copy
		Vec3& operator=(const Vec3&) = default;

		/* function members */
		Coord_Type get_type() const;						//get Coord_Type for the given comp coordinates
		Coord_Type get_type(const Coord_Type other) const;	//get Coord_Type for the given relative comp coordinates
	};
	
	using iVec3 = Vec3<int>;
	using fVec3 = Vec3<real>;
	using cVec3 = Vec3<complex_num>;
	
	/* point wise less */
	template<typename T1, typename T2>
	inline bool ElementWise_Less(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x < b.x && a.y < b.y && a.z < b.z;
	}
	
	/* point wise less than*/
	template<typename T1, typename T2>
	inline bool ElementWise_Less_Eq(const Vec3<T1>& a, const Vec3<T2>& b) {
		return a.x <= b.x && a.y <= b.y && a.z <= b.z;
	}
	
	/* dot product */
	template<typename T>
	inline Vec3<T> operator*(const Vec3<T>& a, const Vec3<T>& b) {
		return Vec3<T>{a.x * b.x, a.y * b.y, a.z * b.z};
	}
	
	/* right scalar product */
	template<typename T>
	inline Vec3<T> operator*(const Vec3<T>& a, const T b) {
		return Vec3<T>{a.x * b, a.y * b, a.z * b};
	}
	/* left scalar product*/
	template<typename T>
	inline Vec3<T> operator*(const T b, const Vec3<T>& a) {
		return Vec3<T>{a.x * b, a.y * b, a.z * b};
	}
	
	/* point wise addtion */
	template<typename T>
	inline Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b) {
		return Vec3<T>{a.x + b.x, a.y + b.y, a.z + b.z};
	}
	
	/* point wise substraction */
	template<typename T>
	inline Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b) {
		return Vec3<T>{a.x - b.x, a.y - b.y, a.z - b.z};
	}
	
	/* make dir_x_tag the x3 axis */
	template<typename T>
	inline Vec3<T> rotate_frame(const Vec3<T>& p, dir_x_tag) {
		return Vec3<T>{p.y, p.z, p.x};
	}
	
	/* make dir_y_tag the x3 axis */
	template<typename T>
	inline Vec3<T> rotate_frame(const Vec3<T>& p, dir_y_tag) {
		return Vec3<T>{p.z, p.x, p.y};
	}
	
	/* make dir_z_tag the x3 axis */
	template<typename T>
	inline Vec3<T> rotate_frame(const Vec3<T>& p, dir_z_tag) {
		return p;
	}
	
	/* overloading functions of Vec3*/
	fVec3 operator/(const fVec3& a, const real b);
	fVec3 operator*(const iVec3& a, const real b);

	/* return the nearest integer (self included) that is even or odd*/
	template<typename S, typename T>
	inline int get_nearest_int(const int x) {
		return ((x & 1) ^ T::val)? S::get_next_int(x) : x;
	}
	
	template<typename S, typename T>
	inline int get_nearest_int(const real x) {
		return get_nearest_int<S, T>(S::round(x));
	}
	
	/* get the nearest point of coord_type to a global computation point on side S (low or high)*/
	template<typename S, typename T>
	inline iVec3 get_nearest_point(const Vec3<T>& p, const Coord_Type type) {
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

	/* p1, p2 are domain comp coordinates, return the largest box(surface, line) containing all the points of type component inside specified by the start and end point*/
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_interior(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			throw std::runtime_error("Coord_Type cannot be Null");

		return std::make_pair(	get_nearest_point<side_high_tag, Vec3<T1>>(p1, type),
								get_nearest_point<side_low_tag, Vec3<T2>>(p2, type));
	}

	/* return the smallest box*/
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_closure(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			throw std::runtime_error("Coord_Type cannot be Null");

		return std::make_pair(	get_nearest_point<side_low_tag, Vec3<T1>>(p1, type),
								get_nearest_point<side_high_tag, Vec3<T2>>(p2, type));
	}

	/* return the intersection of two boxes specified by (p1, p2) and (q1, q2)*/
	template<typename T>
	std::pair<Vec3<T>, Vec3<T>> get_intersection(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {
		return std::make_pair<Vec3<T>, Vec3<T>>({std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
									{std::min(p2.x, q2.x),std::min(p2.y, q2.y), std::min(p2.z, q2.z)});
	}

    /* source functions */
    struct Gaussian_Func {
        real a, mu;									// exp(-a * t^2)s
        real dt;

		Gaussian_Func() = delete;					//delete default constructors
        Gaussian_Func(real sigma_t, real mu = 0);	//exp(-0.5 (x - mu)^2 / sigma_t^2)
        explicit Gaussian_Func(real sigma_f);		//fourier = exp(-0.5 f^2 / sigma_f^2)

		Gaussian_Func(const Gaussian_Func&) = default;	//copy
		Gaussian_Func& operator=(const Gaussian_Func&) = default;


        real operator() (real time);
        real operator[] (int time);					//(time * dt)

        void set_dt(real _dt);
    };

    struct Sinuosuidal_Func {
        real a; //sin(a * t)
        real dt;

		Sinuosuidal_Func() = delete;							//no default constructors
        explicit Sinuosuidal_Func(real freq); 					//sin(2 * pi * freq * t)
		Sinuosuidal_Func(const Sinuosuidal_Func&) = default;	//copy
		Sinuosuidal_Func& operator=(const Sinuosuidal_Func&) = default;

        real operator() (real time);
        real operator[] (int time);

        void set_dt(real _dt);
    };

    struct Rickerwavelet_Func {
        real a, d;											//(1 - 2 * a * (t - d)) * exp(-a * (t - d)^2))
        real dt;

        Rickerwavelet_Func(real fp, real d);				//(1 - 2(pi * fp * (t - d))^2) * exp(-(pi * fp * (t - d))^2)
		Rickerwavelet_Func(const Rickerwavelet_Func&) = default;
		Rickerwavelet_Func& operator=(const Rickerwavelet_Func&) = default;

        real operator() (real time);
        real operator[] (int time);

        void set_dt(real dt);
    };

    /* linear interpolation on fixed interval grids*/
    class GriddedInterp {
    private:
        std::vector<real> x, y, z, v;
		real x0, y0, z0;
		real dx, dy, dz;
        std::string file_name;

		size_t size, loc_jump_x, loc_jump_y, loc_jump_z;

    public:
        /* Semiregular members*/
        GriddedInterp(std::string file_name);
        GriddedInterp(const std::vector<real>& _x, const std::vector<real>& _y, const std::vector<real>& _z, const std::vector<real>& _v);
		GriddedInterp(const GriddedInterp&) = default;				//copy
		GriddedInterp& operator=(const GriddedInterp&) = default;
		GriddedInterp(GriddedInterp&&) = default;					//move
		GriddedInterp& operator=(GriddedInterp&&) = default;

        /* Function Members */
		void init();
        real request_value(const fVec3& p) const;			//return value at a point, 0 dimension is ignored
        real request_integral() const;						//return the integral over the entire volume
		void expand_dim(real lo, real hi, Direction dir);	//expand zero dimension without changing the integral, this is for use in current source
		
		
        fVec3 get_p1() const;								//return the lower corner point
        fVec3 get_p2() const;								//return the upper corner point

		//public functions for calculating generic line, face and volume integral
    };



	/* get size of series of vectors, one vector */
	template<typename T>
	inline size_t get_vec_size(const std::vector<T>& vn) {
		return vn.size();
	}

	/* get size of series of vectors, template */
	template<typename T, typename... Args>
	inline size_t get_vec_size(const std::vector<T>& vn, Args... args) {
		return vn.size() * get_vec_size(args...);
	}

	/* n dimensional linear interpolation, 1 dimension
	   Xn must have fixed interval, otherwise the result does not make sense*/
	template<typename T>
	inline T interp_ndim(const T* data, const real xn, const std::vector<real>& Xn) {
		if (Xn.size() == 0)
			throw std::runtime_error("number of coordinates are zero");

		if (Xn.size() == 1)		//ignore this dimension if it is 1
			return data[0];

		if (xn > Xn.back() || xn < Xn.front())
			return 0;			//zero if it is outside of the region

		real index = (xn - Xn[0]) / (Xn.back() - Xn.front()) * (Xn.size() - 1);

		return (index - floor(index)) * data[(int)index + 1] + (ceil(index) - index) * data[(int)index];
	}

	/* n dimensional linear interpolation, template 
	   Xn must have fixed interval, otherwise the result does not make sense*/
	template<typename T, typename... Args>
	inline T interp_ndim(const T* data, const real xn, const std::vector<real>& Xn, Args... args) {
		if (Xn.size() == 0)
			throw std::runtime_error("number of coordinates are zero");

		if (Xn.size() == 1)		//ignore this dimension if it is 1
			return interp_ndim(data, args...);

		if (xn > Xn.back() || xn < Xn.front())
			return 0;			//zero if it is outside of the region

		real index = (xn - Xn[0]) / (Xn.back() - Xn.front()) * (Xn.size() - 1);
		size_t shift = get_vec_size(args...);

		return (index - floor(index)) * interp_dim(data + shift, args...) + (ceil(index) - index) *  interp_dim(data, args...);
	}
	
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
	
	/* n dimensional linear integration, template,
	xn does not need to have fixed interval*/
	template<typename T, typename... Args>
	T integral_ndim(const T* data, const std::vector<real> &xn, Args... args) {
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
		int d{0};
		real dt;
		real sigma_max{0}, chi_max{1}, a_max{0.15};
		real m_a{1}, m{3};
		Direction dir;
		Side side;

	public:
		/* constructors*/
		PML() = default;
		PML(Direction _dir, Side _side);
		PML(Direction _dir, Side _side, real _d, real _sigma_max);

		PML(const PML&) = default;				//copy
		PML& operator=(const PML&) = default;
		PML(PML&&) = default;						//move
		PML& operator=(PML&&);

		void set_dt(const real _dt);

		real get_d() const;
		real get_sigma(const real x) const;
		real get_chi(const real x) const;
		real get_a(const real x) const;
		real get_b(const real x) const;
		real get_c(const real x) const;
		Direction get_dir() const;
		Side get_side() const;

		static real optimal_sigma_max(real m_k, real dt, real er = 1, real ur = 1);
	};
	
	/* extra field information for use in dispersive field updates*/
	struct Dispersive_Field {
		union {real ex2, ey2, ez2, e2, hx2, hy2, hz2, h2, eh2;};
		std::vector<real> jp1, jp;
	};
	
	/* all points that are in PML layers*/
	struct PML_Point {
		int index, jump_pos, jump_neg;
		real b_pos, b_neg, c_pos, c_neg;
		real psi_pos{0}, psi_neg{0};
		
		PML_Point(const int _index, const int _jump_pos, const int _jump_neg, const real _b_pos, const real _b_neg, const real _c_pos, const real _c_neg);
	};
	
	template<typename P, int N, int M>
	inline typename P::value_type get(const P& p);
	
	
	struct my_iterator {
		int x0, y0, z0, x1, y1, z1;
		int x, y, z, index{0};
		
		my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype);
		void advance();
		bool is_end();
		bool is_empty();
		size_t size();
	};
}
