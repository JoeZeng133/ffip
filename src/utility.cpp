#include <utility.hpp>

namespace ffip
	{
	const double pi = 3.141592653589793e+00;
	unsigned int hex_bit_count[16];
	unsigned int dir_of_ctypes[16];

	void init()
	{
		for (int i = 0; i < 16; ++i)
			hex_bit_count[i] = (i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3);

		std::fill(dir_of_ctypes, dir_of_ctypes + 16, 3);
		dir_of_ctypes[Ex] = dir_of_ctypes[Hx] = dir_of_ctypes[Dx] = dir_of_ctypes[Bx] = 0;
		dir_of_ctypes[Ey] = dir_of_ctypes[Hy] = dir_of_ctypes[Dy] = dir_of_ctypes[By] = 1;
		dir_of_ctypes[Ez] = dir_of_ctypes[Hz] = dir_of_ctypes[Dz] = dir_of_ctypes[Bz] = 2;
	}

	iVec3 get_norm_vec(Direction dir, Side side)
	{
		iVec3 res;
		res[(size_t)dir] = side;
		return res;
	}

	//Yee iterator
	Yee_Iterator::Yee_Iterator(const fVec3 &p1, const fVec3 &p2, const Coord_Type ctype)
	{

		stride = 2;
		auto tmp = get_component_interior(p1, p2, ctype);
		init(tmp.first, tmp.second);
	}

	Yee_Iterator::Yee_Iterator(const fVec3 &p1, const fVec3 &p2)
	{

		stride = 1;
		init(iVec3(p1.ceil()), iVec3(p2.floor()));
	}

	Yee_Iterator::Yee_Iterator(const pair_Vec3<int> &corners, const Coord_Type ctype) : Yee_Iterator{corners.first, corners.second, ctype} {}

	Yee_Iterator::Yee_Iterator(const iVec3 &p1, const iVec3 &p2, const Coord_Type ctype)
	{

		stride = 2;
		auto tmp = get_component_interior(p1, p2, ctype);
		init(tmp.first, tmp.second);
	}

	Yee_Iterator::Yee_Iterator(const pair_Vec3<int> &corners) : Yee_Iterator{corners.first, corners.second} {}

	Yee_Iterator::Yee_Iterator(const iVec3 &p1, const iVec3 &p2)
	{

		stride = 1;
		init(p1, p2);
	}

	size_t Yee_Iterator::get_size(const iVec3 &p1, const iVec3 &p2, const Coord_Type ctype)
	{
		auto tmp = get_component_interior(p1, p2, ctype);

		return (((p2.x - p1.x) / 2 + 1) * ((p2.y - p1.y) / 2 + 1) * ((p2.z - p1.z) / 2 + 1));
	}

	size_t Yee_Iterator::get_size(const iVec3 &p1, const iVec3 &p2)
	{

		return ((p2.x - p1.x + 1) * (p2.y - p1.y + 1) * (p2.z - p1.z + 1));
	}

	void Yee_Iterator::init(const iVec3 &p1, const iVec3 &p2)
	{

		auto len = (p2 - p1) / stride + 1;

		x0 = p1.x;
		y0 = p1.y;
		z0 = p1.z;
		lenx = len.x;
		leny = len.y;
		lenz = len.z;

		i = j = k = 0;

		if (is_empty())
			lenx = leny = lenz = 0;
	}

	void Yee_Iterator::next()
	{
		if (++i == lenx)
		{
			i = 0;
			if (++j == leny)
			{
				j = 0;
				++k;
			}
		}
	}

	bool Yee_Iterator::is_end() const
	{
		return k >= lenz;
	}

	bool Yee_Iterator::is_empty() const
	{
		return lenx <= 0 || leny <= 0 || lenz <= 0;
	}

	size_t Yee_Iterator::get_size() const
	{
		return (size_t)lenx * leny * lenz;
	}

	iVec3 Yee_Iterator::get_dim() const
	{
		return iVec3(lenx, leny, lenz);
	}

	iVec3 Yee_Iterator::get_coord() const
	{
		return {x0 + i * stride, y0 + j * stride, z0 + k * stride};
	}

	//Yee3
	Yee3::Yee3(const iVec3& ghost_p1, const iVec3& ghost_p2) : ghost_p1(ghost_p1), ghost_p2(ghost_p2)
	{

		dim = ghost_p2 - ghost_p1 + 1;
		grid_p1 = ghost_p1 + 1;
		grid_p2 = ghost_p2 - 1;

		stride.z = dim.x * dim.y;
		stride.y = dim.x;
		stride.x = 1;
	}

	iVec3 Yee3::get_grid_p1() const
	{
		return grid_p1;
	}

	iVec3 Yee3::get_grid_p2() const
	{
		return grid_p2;
	}

	size_t Yee3::get_size() const
	{
		return (size_t)dim.x * dim.y * dim.z;
	}

	sVec3 Yee3::get_stride() const
	{
		return stride;
	}

	pair_Vec3<int> Yee3::intersect_with(const iVec3 &p1, const iVec3 &p2) const
	{
		return get_intersection(grid_p1, grid_p2, p1, p2);
	}

	pair_Vec3<int> Yee3::intersect_with(const pair_Vec3<int> &box) const
	{
		return get_intersection(grid_p1, grid_p2, box.first, box.second);
	}

	iVec3 Yee3::get_base_point(const fVec3 &pos, Coord_Type ctype) const
	{
		return get_nearest_point<Negative>(pos, ctype);
	}

	iVec3 Yee3::get_coord_from_index(size_t index) const
	{

		int x, y, z;
		x = ghost_p1.x + (index % dim.x);
		index /= dim.x;
		y = ghost_p1.y + (index % dim.y);
		z = ghost_p1.z + (index / dim.y);
		return {x, y, z};
	}

	long long Yee3::get_index_offset(int i, int j, int k) const
	{

		return (long long)i * stride.x + (long long)j * stride.y + (long long)k * stride.z;
	}

	long long Yee3::get_index_offset(const iVec3 &pos) const
	{

		return get_index_offset(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(const iVec3 &pos) const
	{

		return get_index_from_coord(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(int i, int j, int k) const
	{

		return (i - ghost_p1.x) * stride.x + (j - ghost_p1.y) * stride.y + (k - ghost_p1.z) * stride.z;
	}

	bool Yee3::is_inside(const fVec3 &pos) const
	{

		return leq_vec3(grid_p1, pos) && leq_vec3(pos, grid_p2);
	}

	bool Yee3::is_inside(int i, int j, int k) const
	{

		return i >= grid_p1.x && i <= grid_p2.x &&
			j >= grid_p1.y && j <= grid_p2.y &&
			k >= grid_p1.z && k <= grid_p2.z;
	}

	bool Yee3::is_inside(const iVec3 &pos) const
	{
		return is_inside(pos.x, pos.y, pos.z);
	}

	std::array<double, 8> Yee3::get_interp_weights(double sx, double sy, double sz) const
	{
		return {(1 - sx) * (1 - sy) * (1 - sz),
				sx * (1 - sy) * (1 - sz),
				(1 - sx) * sy * (1 - sz),
				sx * sy * (1 - sz),
				(1 - sx) * (1 - sy) * sz,
				sx * (1 - sy) * sz,
				(1 - sx) * sy * sz,
				sx * sy * sz};
	}

	std::array<double, 8> Yee3::get_interp_weights(const fVec3 &s) const
	{
		return get_interp_weights(s.x, s.y, s.z);
	}

	std::array<double, 8> Yee3::get_interp_weights(const fVec3 &pos, Coord_Type ctype) const
	{

		iVec3 base = get_base_point(pos, ctype);
		double sx = (pos.x - base.x) / 2;
		double sy = (pos.y - base.y) / 2;
		double sz = (pos.z - base.z) / 2;

		return get_interp_weights(sx, sy, sz);
	}

	//Gaussian functions
	double Gaussian1(double t, double width)
	{
		return t * std::exp(-(t * t / 2 / width / width));
	}

	double Gaussian2(double t, double width)
	{
		double arg = (t / width);
		arg = arg * arg;

		return (1 - arg) * exp(-arg);
	}

	//Output overloading
	std::ostream &operator<<(std::ostream &os, const std::complex<double> &c)
	{
		os << c.real() << " " << c.imag();
		return os;
	}

	Grid_3::Grid_3(const iVec3 &p1, const iVec3 &p2) : p1(p1), p2(p2)
	{
		if (p1.get_type() != p2.get_type())
			throw std::runtime_error("Invalid Coord Types");

		if (!leq_vec3(p1, p2))
			throw std::runtime_error("Invalid Volume");

		dim = ((p2 - p1) / 2 + 1);
		interpolant = interpn<3>(dim.z, dim.y, dim.x);
		size = interpolant.get_size();
	}

	size_t Grid_3::get_size() const
	{
		return size;
	}

	sVec3 Grid_3::get_dim() const
	{
		return dim;
	}

	std::vector<double> linspace(int s, int e, int stride)
	{
		std::vector<double> res;
		for (; s <= e; s += stride)
			res.push_back(s);

		return res;
	}

	size_t get_max_size_chunk(const iVec3 &dim, const iVec3 &num)
	{
		//each dimension has to be less than dim
		if (!(leq_vec3(num, dim)))
			return std::numeric_limits<size_t>::max();

		return (ceil((double)dim.x / num.x)) *
			(ceil((double)dim.y / num.y)) *
			(ceil((double)dim.z / num.z));
	};

	iVec3 decompose_domain(const iVec3 &dim, int np)
	{

		iVec3 res = {1, 1, np};
		size_t min_size = get_max_size_chunk(dim, {1, 1, np});

		for (int i = 1; i <= std::sqrt(np); ++i)
			if (np % i == 0)
			{

				int rest = np / i;

				for (int j = 1; j <= std::sqrt(rest); ++j)
					if (rest % j == 0)
					{
						int k = rest / j;
						//(i, j, k)
						if (auto size = get_max_size_chunk(dim, {i, j, k}); size < min_size)
						{
							min_size = size;
							res = {i, j, k};
						}

						//(i, k, j)
						if (auto size = get_max_size_chunk(dim, {i, k, j}); size < min_size)
						{
							min_size = size;
							res = {i, k, j};
						}
					}

				for (int j = 1; j <= std::sqrt(i); ++j)
					if (i % j == 0)
					{
						int k = i / j;
						//(rest, j, k)
						if (auto size = get_max_size_chunk(dim, {rest, j, k}); size < min_size)
						{
							min_size = size;
							res = {rest, j, k};
						}

						//(rest, k, j)
						if (auto size = get_max_size_chunk(dim, {rest, k, j}); size < min_size)
						{
							min_size = size;
							res = {rest, k, j};
						}
					}
			}

		return res;
	}

	pair_Vec3<int> get_chunk_from_coords(const iVec3 &dim, const iVec3 &num, const iVec3 &coord)
	{
		iVec3 p1 = (dim * coord) / num;
		iVec3 p2 = (dim * (coord + 1)) / num;

		return {p1, p2};
	}
} // namespace ffip
