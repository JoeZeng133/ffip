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
#include <nlohmann/json.hpp>

namespace ffip
{
	using json = nlohmann::json;

	//type aliases
	using double_arr = std::vector<double>;
	using complex_pair = std::pair<std::complex<double>, std::complex<double>>;
	using complex_num = std::complex<double>;
	using complex_arr = std::vector<complex_num>;

	//constants
	extern const double pi;
	extern unsigned int hex_bit_count[16];
	extern unsigned int dir_of_ctypes[16];

	//enums with bit pattern indication
	enum Direction
	{
		X = 0,
		Y = 1,
		Z = 2
	};
	enum Side
	{
		Negative = -1,
		Positive = 1
	};
	enum Coord_Type
	{
		Ex = 0b001,
		Ey = 0b010,
		Ez = 0b100,
		Hx = 0b110,
		Hy = 0b101,
		Hz = 0b011,
		Dx = 0b1001,
		Dy = 0b1010,
		Dz = 0b1100,
		Bx = 0b1110,
		By = 0b1101,
		Bz = 0b1011,
		Corner = 0b000,
		Center = 0b111
	};

	//initializations of utility (calculating constants)
	void init();

	//count the number of 1s in the first 4 bits
	constexpr unsigned int count_hex_bits(unsigned int x)
	{

		return (x & 1) + ((x & 2) >> 1) + ((x & 4) >> 2) + ((x & 8) >> 3);
	}

	//return whether it is electrical point
	constexpr bool is_e_point(Coord_Type ctype)
	{

		return count_hex_bits(ctype & 0b111) == 1;
	}

	//return wehether it is magnetic point
	constexpr bool is_m_point(Coord_Type ctype)
	{

		return count_hex_bits(ctype & 0b111) == 2;
	}

	//return whether it is an eh point
	constexpr bool is_eh_point(Coord_Type ctype)
	{

		return ctype < 0b111 && ctype > 0;
	}

	//return whether it is a db point
	constexpr bool is_db_point(Coord_Type ctype)
	{

		return ctype > 0b111 && ctype < 0b1111;
	}

	//return Ed
	constexpr Coord_Type get_e_ctype_from_dir(Direction d)
	{
		return Coord_Type(1 << d);
	}

	//return Hd
	constexpr Coord_Type get_m_ctype_from_dir(Direction d)
	{
		return Coord_Type(0b111 ^ (1 << d));
	}

	//get direction integer from coord_type, 3 means no direction
	constexpr Direction get_dir_from_ctype(Coord_Type ctype) {
		int ctype_int = ctype & 0b111;
		if (is_m_point(ctype))
			ctype_int = 0b111 ^ ctype_int;
			
		if(ctype_int == 0b001) return X;
		if(ctype_int == 0b010) return Y;
		if(ctype_int == 0b100) return Z;
		return Direction(3);
	}

	//return x2 in (x1, x2, x3) 
	constexpr Direction get_x2_from_x3(Direction x3) {
		return static_cast<Direction>((x3 + 2) % 3);
	}

	//return x1 in (x1, x2, x3)
	constexpr Direction get_x1_from_x3(Direction x3) {
		return static_cast<Direction>((x3 + 1) % 3);
	}

	// light weight 3 element vector
	template <typename T = double>
	struct Vec3
	{
		using value_type = T;

		T x{}, y{}, z{};

		Vec3() = default;

		//construction
		Vec3(T x, T y, T z) : x(x), y(y), z(z) {}

		//implicit copy construction
		template <typename T2>
		Vec3(const Vec3<T2> &other) : x(other.x), y(other.y), z(other.z) {}

		//copy assignment
		template <typename T2>
		Vec3 &operator=(const Vec3<T2> &other)
		{
			x = other.x;
			y = other.y;
			z = other.z;
			return *this;
		};

		//pointwise division assiginment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator/=(const Vec3<T2> &div) const
		{
			x /= div.x;
			y /= div.y;
			z /= div.z;
			return *this;
		}

		//scalar division assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator/=(const T2 div) const
		{
			x /= div;
			y /= div;
			z /= div;
			return *this;
		}

		//scalar multiplication assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator*=(const T2 mult) const
		{
			x *= mult;
			y *= mult;
			z *= mult;
			return *this;
		}

		//pointwise multiplication assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator*=(const Vec3<T2> &mult) const
		{
			x *= mult.x;
			y *= mult.y;
			z *= mult.z;
			return *this;
		}

		//pointwise addition assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator+=(const Vec3<T2> &add) const
		{
			x += add.x;
			y += add.y;
			z += add.z;
			return *this;
		}

		//scalar addition assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator+=(const T2 add) const
		{
			x += add;
			y += add;
			z += add;
			return *this;
		}

		//scalar subtraction assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator-=(const T2 minus) const
		{
			x -= minus;
			y -= minus;
			z -= minus;
			return *this;
		}

		//pointwise scalar subtraction assignment
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> &operator-=(const Vec3<T2> &minus) const
		{
			x -= minus.x;
			y -= minus.y;
			z -= minus.z;
			return *this;
		}

		//pointwise division
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator/(const Vec3<T2> &div) const
		{
			return {x / div.x, y / div.y, z / div.z};
		}

		//scalar division
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator/(const T2 div) const
		{
			return {x / div, y / div, z / div};
		}

		//scalar multiplication
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const T2 mult) const
		{
			return {x * mult, y * mult, z * mult};
		}

		//pointwise multiplication
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator*(const Vec3<T2> &mult) const
		{
			return {x * mult.x, y * mult.y, z * mult.z};
		}

		//scalar addition
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const T2 add) const
		{
			return {x + add, y + add, z + add};
		}

		//pointwise scalar addition
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator+(const Vec3<T2> &add) const
		{
			return {x + add.x, y + add.y, z + add.z};
		}

		//pointwise scalar subtraction
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const Vec3<T2> &minus) const
		{
			return {x - minus.x, y - minus.y, z - minus.z};
		}

		//scalar subtraction
		template <typename T2>
		Vec3<std::common_type_t<T, T2>> operator-(const T2 minus) const
		{
			return {x - minus, y - minus, z - minus};
		}

		//to vector container
		std::vector<T> to_vector() const
		{
			return {x, y, z};
		}

		//to vector container and reverse
		std::vector<T> to_vector_reverse() const
		{
			return {z, y, x};
		}

		//to fixed length array
		std::array<T, 3> to_array() const
		{
			return {x, y, z};
		}

		//to fixed length array reverse
		std::array<T, 3> to_array_reverse() const
		{
			return {z, y, x};
		}

		//serialize to string
		std::string to_string() const
		{
			return std::to_string(x) + std::to_string(y) + std::to_string(z);
		}

		//pointwise comparisons
		bool operator==(const Vec3 &other) const
		{
			return x == other.x && y == other.y && z == other.z;
		}

		//Unary -
		Vec3 operator-() const
		{
			return {-x, -y, -z};
		}

		//Unary+
		Vec3 operator+() const
		{
			return {x, y, z};
		}

		//dot product
		template <typename T2>
		std::common_type_t<T, T2> dot(const Vec3<T2> &other) const
		{
			return x * other.x + y * other.y + z * other.z;
		}

		//ceil
		Vec3 ceil() const
		{
			return Vec3(std::ceil(x), std::ceil(y), std::ceil(z));
		}

		//floor
		Vec3 floor() const
		{
			return Vec3(std::floor(x), std::floor(y), std::floor(z));
		}

		Vec3 round() const
		{
			return Vec3(std::round(x), std::round(y), std::round(z));
		}

		//absolute value
		Vec3 abs() const
		{
			return Vec3(std::abs(x), std::abs(y), std::abs(z));
		}

		//return x * y * z
		T prod() const
		{
			return x * y * z;
		}

		//get<N> style access
		template <unsigned int N>
		T &get()
		{
			if constexpr (N == 0)
				return x;
			else if constexpr (N == 1)
				return y;
			else if constexpr (N == 2)
				return z;
		}

		//get<N> style access, const
		template <unsigned int N>
		const T &get() const
		{
			if constexpr (N == 0)
				return x;
			else if constexpr (N == 1)
				return y;
			else if constexpr (N == 2)
				return z;
		}

		//[] style access
		T &operator[](const size_t n)
		{
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		//[] style access
		const T &operator[](const size_t n) const
		{
			if (n == 0)
				return x;
			else if (n == 1)
				return y;
			return z;
		}

		//get Coord_Type for the given grid coordinates
		template <typename Dummy = T, typename = std::enable_if_t<std::is_integral<Dummy>::value>>
		Coord_Type get_type() const
		{
			return static_cast<Coord_Type>((x & 1) | ((y & 1) << 1) | ((z & 1) << 2));
		}
	};

	using iVec3 = Vec3<int>;
	using fVec3 = Vec3<double>;
	using cVec3 = Vec3<complex_num>;
	using sVec3 = Vec3<long long>;

	template<typename T>
	using pair_Vec3 = std::pair<Vec3<T>, Vec3<T>>;

	//get normal vector in certain direction with sign
	iVec3 get_norm_vec(Direction dir, Side side);

	//strictly pointwise less
	template <typename T1, typename T2>
	constexpr bool le_vec3(const Vec3<T1> &a, const Vec3<T2> &b)
	{
		return a.x < b.x && a.y < b.y && a.z < b.z;
	}

	//pointwise less equal
	template <typename T1, typename T2>
	constexpr bool leq_vec3(const Vec3<T1> &a, const Vec3<T2> &b)
	{
		return a.x <= b.x && a.y <= b.y && a.z <= b.z;
	}

	//return nearest odd(even) integer of a integer
	template <int S, unsigned int T>
	constexpr int get_nearest_int(const int x)
	{
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");

		return (x & 1) ^ T ? x + S : x;
	}

	//return nearest odd(even) integer of a double
	template <int S, unsigned int T>
	constexpr int get_nearest_int(const double x)
	{
		static_assert(T < 2, "T should be 0(even) or 1(odd)");
		static_assert(S == 1 || S == -1, "S should be -1 or 1");

		if constexpr (S == -1)
			return get_nearest_int<S, T>((int)std::floor(x));
		else
			return get_nearest_int<S, T>((int)std::ceil(x));
	}

	//get nearest point of a coord type
	template <int S, typename T>
	constexpr iVec3 get_nearest_point(const Vec3<T> &p, const Coord_Type ctype)
	{
		iVec3 res;

		if (ctype & 1)
		{
			res.x = get_nearest_int<S, 1>(p.x);
		}
		else
		{
			res.x = get_nearest_int<S, 0>(p.x);
		}

		if (ctype & 2)
		{
			res.y = get_nearest_int<S, 1>(p.y);
		}
		else
		{
			res.y = get_nearest_int<S, 0>(p.y);
		}

		if (ctype & 4)
		{
			res.z = get_nearest_int<S, 1>(p.z);
		}
		else
		{
			res.z = get_nearest_int<S, 0>(p.z);
		}

		return res;
	}

	//return grid points inside a box
	template <typename T1, typename T2>
	constexpr pair_Vec3<int>
	get_component_interior(const Vec3<T1> &p1, const Vec3<T2> &p2, const Coord_Type ctype)
	{

		return {get_nearest_point<Positive>(p1, ctype), get_nearest_point<Negative>(p2, ctype)};
	}

	//return grid points enclosing the box
	template <typename T1, typename T2>
	constexpr pair_Vec3<int>
	get_component_closure(const Vec3<T1> &p1, const Vec3<T2> &p2, const Coord_Type type)
	{

		return {get_nearest_point<Negative>(p1, type), get_nearest_point<Positive>(p2, type)};
	}

	//return intersection of two closed intervals
	template <typename T>
	constexpr std::pair<T, T>
	get_intersection(const T p1, const T p2, const T q1, const T q2)
	{

		return {std::max(p1, q1), std::min(p2, q2)};
	}

	//return intersection of two closed cubes
	template <typename T>
	constexpr pair_Vec3<T>
	get_intersection(const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &q1, const Vec3<T> &q2)
	{

		return {{std::max(p1.x, q1.x), std::max(p1.y, q1.y), std::max(p1.z, q1.z)},
				{std::min(p2.x, q2.x), std::min(p2.y, q2.y), std::min(p2.z, q2.z)}};
	}

	//determine if two closed interval intersects
	template <typename T>
	constexpr bool is_intersect(const T p1, const T p2, const T q1, const T q2)
	{
		return p1 <= q2 && q1 <= p2;
	}

	//determine if two closed cubes intersects
	template <typename T>
	constexpr bool
	is_intersect(const Vec3<T> &p1, const Vec3<T> &p2, const Vec3<T> &q1, const Vec3<T> &q2)
	{

		return is_intersect(p1.x, p2.x, q1.x, q2.x) &&
			is_intersect(p1.y, p2.y, q1.y, q2.y) &&
			is_intersect(p1.z, p2.z, q1.z, q2.z);
	}

	//return one face of a box
	template <unsigned int D, int S, typename T>
	pair_Vec3<T> get_face(Vec3<T> p1, Vec3<T> p2)
	{

		static_assert(D < 3, "D should be 0(even) or 1(odd)");
		static_assert(S == -1 || S == 1, "S must be -1 or 1");

		if constexpr (S == -1)
		{
			p2.template get<D>() = p1.template get<D>();
		}
		else
		{
			p1.template get<D>() = p2.template get<D>();
		}

		return {p1, p2};
	}

	//return once face of a box (overload)
	template <unsigned int D, int S, typename T>
	pair_Vec3<T> get_face(const pair_Vec3<T>& p)
	{
		return get_face<D, S, T>(p.first, p.second);
	}

	//return one face of a box, non-template function
	template <typename T>
	pair_Vec3<T> get_face(Vec3<T> p1, Vec3<T> p2, Direction dir, Side side)
	{

		switch (dir)
		{
		case X:
			if (side == Positive)
				p1.template get<X>() = p2.template get<X>();
			else
				p2.template get<X>() = p1.template get<X>();
			break;

		case Y:
			if (side == Positive)
				p1.template get<Y>() = p2.template get<Y>();
			else
				p2.template get<Y>() = p1.template get<Y>();
			break;

		case Z:
			if (side == Positive)
				p1.template get<Z>() = p2.template get<Z>();
			else
				p2.template get<Z>() = p1.template get<Z>();
			break;
		}

		return {p1, p2};
	}

	//return one face of a box, non-template function (overload)
	template <typename T>
	pair_Vec3<T> get_face(const pair_Vec3<T>& p, Direction dir, Side side) {
		return get_face(p.first, p.second, dir, side);
	}

	//interpn, linear interpolation with 0 padding extrapolation
	//0 padding extrapolation is compatible with paralell interpolation
	//ignoring dimension of 1
	template <int N>
	class interpn : public interpn<N - 1>
	{
	public:
		using base_class = interpn<N - 1>;
		size_t dimn{1}, stride{1};

		interpn() = default;

		template <typename... Args>
		interpn(const size_t dimn, Args... args) : base_class(args...), dimn(dimn)
		{

			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");

			stride = base_class::get_size();
		}

		//get the array size
		size_t get_size() const
		{
			return dimn * base_class::get_size();
		}

		//interpolate xn, xn-1, ..., x1
		template <typename T, typename... Args>
		T operator()(T const *data, double xq, Args &&... args) const
		{
			//ignore this dimension if it is 1
			if (dimn == 1)
				return base_class::operator()(data, args...);

			//0 padding extrapolation
			if (xq <= 0)
			{
				if (xq <= -1)
					return T{};
				return (xq + 1) * base_class::operator()(data, args...);
			}

			if (xq >= dimn - 1)
			{
				if (xq >= dimn)
					return T{};
				return (dimn - xq) * base_class::operator()(data + stride * (dimn - 1), args...);
			}

			size_t index = xq;
			double tx = xq - index;

			return tx * base_class::operator()(data + stride * (index + 1), args...) +
				(1 - tx) * base_class::operator()(data + stride * index, args...);
		}

		//interpolate xn, xn-1, ..., x1
		template <typename T, typename... Args>
		T operator()(const std::vector<T> &vec, Args &&... args) const
		{
			return operator()(vec.data(), args...);
		}

		//transpose_helper
		template <typename T, typename... Args>
		void transpose(T *data, const T &val, double xq, Args &&... args) const
		{

			//ignore 1 dimension
			if (dimn == 1)
			{
				base_class::transpose(data, args...);
				return;
			}

			//0 padding extrapolation
			if (xq <= 0)
			{
				if (xq <= -1)
					return;
				base_class::transpose(data, val * (1 + xq), args...);
				return;
			}

			if (xq >= dimn - 1)
			{
				if (xq >= dimn)
					return;
				base_class::transpose(data + stride * (dimn - 1), val * (dimn - xq), args...);
				return;
			}

			size_t index = xq;
			double tx = xq - index;

			base_class::transpose(data + stride * (index + 1), val * tx, args...);
			base_class::transpose(data + stride * index, val * (1 - tx), args...);

			return;
		}

		//transpose val additively back to original grid
		template <typename T, typename... Args>
		void transpose(std::vector<T> &vec, const T &val, Args &&... args) const
		{
			transpose(vec.data(), val, args...);
		}
	};

	template <>
	class interpn<1>
	{
	public:
		size_t dimn{1}, stride{1};

		interpn() = default;

		interpn(const size_t _dimn) : dimn(_dimn)
		{

			if (dimn < 1)
				throw std::runtime_error("Invalid Dimension");

			stride = 1;
		}

		size_t get_size() const
		{
			return dimn;
		}

		template <typename T>
		T operator()(const T *data, double xq) const
		{

			//ignore this dimension if it is 1
			if (dimn == 1)
				return data[0];

			//0 padding extrapolation
			if (xq <= 0)
			{
				if (xq <= -1)
					return T{};
				return (xq + 1) * data[0];
			}

			if (xq >= dimn - 1)
			{
				if (xq >= dimn)
					return T{};
				return (dimn - xq) * data[dimn - 1];
			}

			size_t index = xq;
			double tx = xq - index;

			return tx * data[index + 1] + (1 - tx) * data[index];
		}

		template <typename T>
		T operator()(const std::vector<T> &vec, double xq) const
		{
			return operator()(vec.data(), xq);
		}

		template <typename T>
		void transpose(std::vector<T> &vec, const T &val, double xq) const
		{
			return transpose(vec.data(), val, xq);
		}

		template <typename T>
		void transpose(T *data, const T &val, double xq) const
		{

			//ignore this dimension if it is 1
			if (dimn == 1)
			{
				data[0] = val;
				return;
			}

			//0 padding extrapolation
			if (xq <= 0)
			{
				if (xq <= -1)
					return;
				data[0] += val * (1 + xq);
				return;
			}

			if (xq >= dimn - 1)
			{
				if (xq >= dimn)
					return;
				data[dimn - 1] += val * (dimn - xq);
				return;
			}

			size_t index = xq;
			double tx = xq - index;

			data[index + 1] += val * tx;
			data[index] += val * (1 - tx);
		}
	};

	//Iterator through points in Yee Cells
	struct Yee_Iterator
	{

		int x0, y0, z0;
		int lenx, leny, lenz;
		int i, j, k;
		int stride;

		//floating point constructor
		Yee_Iterator(const fVec3 &p1, const fVec3 &p2, const Coord_Type ctype);

		//floating point constructor
		Yee_Iterator(const fVec3 &p1, const fVec3 &p2);

		//loop through specific coord type
		Yee_Iterator(const iVec3 &p1, const iVec3 &p2, const Coord_Type ctype);

		//loop through specific coord type, delegating constructor
		Yee_Iterator(const pair_Vec3<int> &corners, const Coord_Type ctype);

		//loop through all points
		Yee_Iterator(const iVec3 &p1, const iVec3 &p2);

		//loop through all points, delegating constructor
		Yee_Iterator(const pair_Vec3<int> &corners);

		void init(const iVec3 &p1, const iVec3 &p2);

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

		//return get_size without constructing an instance
		static size_t get_size(const iVec3 &p1, const iVec3 &p2, const Coord_Type ctype);

		//return get_size without constructing an instance
		static size_t get_size(const iVec3 &p1, const iVec3 &p2);

		//return dimension of the region
		iVec3 get_dim() const;
	};

	//Yee3 grid
	//Ghost points aware interpolation
	struct Yee3
	{

		//grid points coordinates
		iVec3 grid_p1, grid_p2;
		//ghost point coordinates
		iVec3 ghost_p1, ghost_p2;
		//dimensions of the array
		iVec3 dim;
		//strides of the array
		sVec3 stride;

		Yee3(const iVec3& ghost_p1, const iVec3& ghost_p2);

		//get grid lower corner, excluding ghost point
		iVec3 get_grid_p1() const;

		//get grid uppper corner, excluding ghost point
		iVec3 get_grid_p2() const;

		//return size of array compatible with the grid
		size_t get_size() const;

		//return strides in array
		sVec3 get_stride() const;

		//intersect with
		pair_Vec3<int> intersect_with(const iVec3 &p1, const iVec3 &p2) const;
		pair_Vec3<int> intersect_with(const pair_Vec3<int> &box) const;

		//an alias for get_nearest_point<Negative>(pos, ctype)
		iVec3 get_base_point(const fVec3 &pos, Coord_Type ctype) const;

		//get coordinate from a index, no range checking
		iVec3 get_coord_from_index(size_t n) const;

		//get index for a given point, no range checking
		size_t get_index_from_coord(const iVec3 &pos) const;
		size_t get_index_from_coord(int i, int j, int k) const;

		//inside grids, excluding ghost points
		bool is_inside(const fVec3 &pos) const;
		bool is_inside(int i, int j, int k) const;
		bool is_inside(const iVec3 &pos) const;

		//get index offset, same as get_index(p + offset) - get_index(p)
		long long get_index_offset(const iVec3 &offset) const;
		long long get_index_offset(int i, int j, int k) const;

		//get corresponding dimension
		template <int N>
		size_t get_dim() const
		{
			return dim.get<N>();
		}

		//access the raw value from a grid coordinate, no range checking
		template <typename T>
		T get_raw_val(const std::vector<T> &data, const iVec3 &pos) const
		{
			return data[get_index_from_coord(pos)];
		}

		//return 3 dim interpolation weights based on sx, sy, sz
		std::array<double, 8> get_interp_weights(double sx, double sy, double sz) const;

		//s.x, s.y, s.z
		std::array<double, 8> get_interp_weights(const fVec3 &s) const;

		//return interpolation weights
		std::array<double, 8> get_interp_weights(const fVec3 &pos, Coord_Type ctype) const;

		//trilinear interpolate a value from a 3d gridded data
		//in MPI reduce, it makes sure interpolation added together is correct
		template <typename T>
		T interp(const std::vector<T> &data, const fVec3 &pos, const Coord_Type ctype) const
		{
			iVec3 base = get_base_point(pos, ctype);
			T res{};

			//if all points fall inside
			if (leq_vec3(grid_p1, base) && leq_vec3(base + 2, grid_p2))
			{

				auto weights = get_interp_weights((pos - base) / 2);
				size_t base_index = get_index_from_coord(base);
				size_t weight_index = 0;

				for (int k = 0; k < 4; k += 2)
					for (int j = 0; j < 4; j += 2)
						for (int i = 0; i < 4; i += 2)
						{
							res += data[base_index + get_index_offset(i, j, k)] * weights[weight_index++];
						}
			}
			//if some points intersect, add range checking (slower)
			else if (is_intersect(base, base + 2, grid_p1, grid_p2))
			{
				auto weights = get_interp_weights((pos - base) / 2);
				size_t weight_index = 0;

				for (int k = 0; k < 4; k += 2)
					for (int j = 0; j < 4; j += 2)
						for (int i = 0; i < 4; i += 2)
						{
							auto p = base + iVec3{i, j, k};
							if (is_inside(p))
								res += data[get_index_from_coord(p)] * weights[weight_index];
							++weight_index;
						}
			}

			return res;
		}
	};

	//for interpolation of a particular coord_type
	struct Grid_3
	{
		iVec3 p1, p2;

		size_t size;
		sVec3 dim;
		interpn<3> interpolant;

		Grid_3(const iVec3 &p1, const iVec3 &p2);

		template <typename T>
		T operator()(const std::vector<T> &data, const fVec3 &pos) const
		{
			auto rel_pos = (pos - p1) / 2.0;
			return interpolant(data, rel_pos.z, rel_pos.y, rel_pos.x);
		}

		template <typename T>
		T operator()(const T *data, const fVec3 &pos) const
		{
			auto rel_pos = (pos - p1) / 2.0;
			return interpolant(data, rel_pos.z, rel_pos.y, rel_pos.x);
		}

		size_t get_size() const;

		sVec3 get_dim() const;
	};

	//output overloading
	std::ostream &operator<<(std::ostream &os, const std::complex<double> &c);

	template <typename T>
	std::ostream &operator<<(std::ostream &os, const Vec3<T> &c)
	{
		os << c.x << " " << c.y << " " << c.z;
		return os;
	}

	//Gaussian first derivative
	double Gaussian1(double t, double width);

	//Gaussian second derivative
	double Gaussian2(double t, double width);

	//linspace integer
	std::vector<double> linspace(int s, int e, int stride = 1);

	//return maximum size of chunk among chunks in domain divided by [num1, num2, num3]
	size_t get_max_size_chunk(const iVec3 &dim, const iVec3 &num);

	//return optimal decomposition of domain given number of processes and dimension of domain
	iVec3 decompose_domain(const iVec3 &dim, int np);

	//return chunk corner points given its coordinates
	pair_Vec3<int> get_chunk_from_coords(const iVec3 &dim, const iVec3 &num, const iVec3 &coord);
} // namespace ffip
