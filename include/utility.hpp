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

namespace ffip {
	/* type aliases */
    using real = double;
	using real_arr = std::vector<real>;
	using complex_pair = std::pair<std::complex<real>, std::complex<real>>;
	using complex_num = std::complex<real>;
	using complex_arr = std::vector<complex_num>;

    /* constants */
	extern const real e0, u0, z0, c0;
	extern const real pi;

	/* enums with underlying bit pattern for calculations*/
	enum Direction {X = 0, Y = 1, Z = 2};
	enum Side {Low = -1, High = 1};
	enum Coord_Type {Ex = 0b001, Ey = 0b010, Ez = 0b100, Hx = 0b110, Hy = 0b101, Hz = 0b011, Corner = 0b000, Center = 0b111, Null = -1};
	
	int Ctype2DirInt(const Coord_Type ctype);
	
	constexpr bool is_E_point(const Coord_Type ctype) {
		return (ctype == Ex || ctype == Ey || ctype == Ez);
	}
	
	constexpr bool is_H_point(const Coord_Type ctype) {
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
    struct side_high_tag{
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
		using z = dir_y_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
		
	};
	
	/* y direction, x1(z), x2(x), x3(y)*/
	struct dir_y_tag{
		using x1 = dir_z_tag;
		using x2 = dir_x_tag;
		using z = dir_x_tag;
		
		static const int val;
		static const Coord_Type E;
		static const Coord_Type H;
	};
	
	/* z direction, x1(x), x2(y), x3(z)*/
	struct dir_z_tag{
		using x1 = dir_x_tag;
		using x2 = dir_y_tag;
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
	
	using iVec3 = Vec3<int>;
	using fVec3 = Vec3<real>;
	using cVec3 = Vec3<complex_num>;
	
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

	/* p1, p2 are domain comp coordinates, return the largest box(surface, line) containing all the points of type component inside specified by the start and end point*/
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_interior(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			throw std::runtime_error("Coord_Type cannot be Null");

		return std::make_pair(	get_nearest_point<side_high_tag>(p1, type),
								get_nearest_point<side_low_tag>(p2, type));
	}

	/* return the smallest box*/
	template<typename T1, typename T2>
	inline std::pair<iVec3, iVec3> get_component_closure(const Vec3<T1>& p1, const Vec3<T2>& p2, const Coord_Type type) {
		if (type == Coord_Type::Null)
			throw std::runtime_error("Coord_Type cannot be Null");

		return std::make_pair(	get_nearest_point<side_low_tag>(p1, type),
								get_nearest_point<side_high_tag>(p2, type));
	}

	/* return the intersection of two boxes specified by (p1, p2) and (q1, q2)*/
	template<typename T>
	std::pair<Vec3<T>, Vec3<T>> get_intersection(const Vec3<T>& p1, const Vec3<T>& p2, const Vec3<T>& q1, const Vec3<T>& q2) {
		return std::make_pair<Vec3<T>, Vec3<T>>({std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
									{std::min(p2.x, q2.x),std::min(p2.y, q2.y), std::min(p2.z, q2.z)});
	}
	
	/* return a particular face of a box region specified by (p1, p2), very nasty*/
	template<typename D, typename S, typename T>
	std::pair<Vec3<T>, Vec3<T>> get_face(Vec3<T> p1, Vec3<T> p2) {
		choose<D>::get(S::choose_between(p2, p1)) = choose<D>::get(S::choose_between(p1, p2));
		return {p1, p2};
	}

    /* source functions */
    struct Gaussian_Func {
        real a, mu;									// exp(-a * t^2)s
        real dt;

		Gaussian_Func() = delete;					//delete default constructors
        Gaussian_Func(real sigma_t, real mu);	//exp(-0.5 (x - mu)^2 / sigma_t^2)
        explicit Gaussian_Func(real sigma_f);		//fourier = exp(-0.5 f^2 / sigma_f^2)

		Gaussian_Func(const Gaussian_Func&) = default;	//copy
		Gaussian_Func& operator=(const Gaussian_Func&) = default;


        real operator() (real time) const;
        real operator[] (int time) const;					//(time * dt)
		
		auto get_functor() -> std::function<real(const real)> const;
        void set_dt(real _dt);
    };

    struct Sinuosuidal_Func {
        real a; //sin(a * t)
		real phase{0};
        real dt;

		Sinuosuidal_Func() = delete;							//no default constructors
        Sinuosuidal_Func(real freq, real d = 0); 					//sin(2 * pi * freq * (t - d))
		Sinuosuidal_Func(const Sinuosuidal_Func&) = default;	//copy
		Sinuosuidal_Func& operator=(const Sinuosuidal_Func&) = default;

        real operator() (real time) const;
        real operator[] (int time) const;

		auto get_functor() -> std::function<real(const real)> const;
        void set_dt(real _dt);
    };

    struct Rickerwavelet_Func {
        real a, d;											//(1 - 2 * a * (t - d)) * exp(-a * (t - d)^2))
        real dt;

        Rickerwavelet_Func(real fp, real _d);				//(1 - 2(pi * fp * (t - d))^2) * exp(-(pi * fp * (t - d))^2)
		Rickerwavelet_Func(const Rickerwavelet_Func&) = default;
		Rickerwavelet_Func& operator=(const Rickerwavelet_Func&) = default;

        real operator() (real time) const;
        real operator[] (int time) const;

		auto get_functor() -> std::function<real(const real)> const;
        void set_dt(real dt);
    };

    /* linear interpolation on fixed interval grids*/
    class GriddedInterp {
    private:
		std::vector<real> x, y, z, v;
        std::string file_name;

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
		
		if (index == int(index))
			return data[int(index)];
		else
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
		//std::cout << index << " ";
		size_t shift = get_vec_size(args...);

		if (index == int(index))
			return interp_ndim(data + shift * int(index), args...);
		else
			return (index - floor(index)) * interp_ndim(data + shift * int(index + 1), args...) + (ceil(index) - index) *  interp_ndim(data + shift * int(index), args...);
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
		real sigma_max{0}, k_max{1}, a_max{0.1};	//it was found out a_max is critical in absorbing waves in 1D simulation
		real m_a{1}, m{3};
		Direction dir;
		Side side;

	public:
		/* constructors*/
		PML() = default;
		PML(Direction _dir, Side _side);
		PML(Direction _dir, Side _side, int _d);
		PML(Direction _dir, Side _side, int _d, real _sigma_max);

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
		Direction get_dir() const;
		Side get_side() const;

		static real optimal_sigma_max(real m_k, real dx, real er = 1, real ur = 1);
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
	
	
	/* an iterator to iterate through a box region specified by two corner regions [p1, p2]
	 it allows iteration of particular coord_type points: ex, ..., ez, hx, ..., hz
	 for null, it will loop through all points
	 it allows looping through a particular part of the region
	 */
	struct my_iterator {
		using value_type = int;
		
		size_t size{0}, index{0};
		int x0, y0, z0, x1, y1, z1;
		int x, y, z;
		int jump;
		
		my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype);
		my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype, const size_t rank, const size_t num);
		
		void advance();
		bool is_end() const;
		bool is_empty() const;
		iVec3 get_vec() const;
		size_t get_size() const;
	};

	template<typename T, typename F>
	void task_divider(std::vector<T>& list, F& func, const int num_proc) {
		if (list.size() < num_proc)
			for(auto& item : list) {
				func(item);
			}
		else {
			std::vector<std::thread> threads;

			for(int i = 1; i < num_proc; ++i) {
				auto itr1 = list.begin() + (i * list.size()) / num_proc;
				auto itr2 = list.begin() + ((i + 1) * list.size()) / num_proc;
				threads.push_back(std::thread(std::for_each<decltype(itr1), F>, itr1, itr2, std::ref(func)));
			}
			auto itr1 = list.begin() + (0 * list.size()) / num_proc;
			auto itr2 = list.begin() + (1 * list.size()) / num_proc;
			std::for_each(itr1, itr2, func);
			
			for(auto& item : threads)
				item.join();
		}
	}

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
}
