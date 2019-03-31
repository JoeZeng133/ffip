#include <utility.hpp>

namespace ffip {
    const double pi = 3.141592653589793e+00;
	unsigned int hex_bit_count[16];
	unsigned int dir_of_ctypes[16];

	void init() {
		for(int i = 0; i < 16; ++i)
			hex_bit_count[i] = (i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3);
			

		std::fill(dir_of_ctypes, dir_of_ctypes + 16, 3);
		dir_of_ctypes[Ex] = dir_of_ctypes[Hx] = dir_of_ctypes[Dx] = dir_of_ctypes[Bx] = 0;
		dir_of_ctypes[Ey] = dir_of_ctypes[Hy] = dir_of_ctypes[Dy] = dir_of_ctypes[By] = 1;
		dir_of_ctypes[Ez] = dir_of_ctypes[Hz] = dir_of_ctypes[Dz] = dir_of_ctypes[Bz] = 2;

	}

	//Yee iterator
	Yee_Iterator::Yee_Iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype): Yee_Iterator{corners.first, corners.second, ctype} {}
	
	Yee_Iterator::Yee_Iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype) {
		
		// loop through points inside [p1, p2]
		auto tmp = get_component_interior(p1, p2, ctype);
		if (ctype == All) {
			stride = 1;
		} else {
			stride = 2;
		}

		auto len = (tmp.second - tmp.first) / stride + 1;
		
		x0 = tmp.first.x;
		y0 = tmp.first.y;
		z0 = tmp.first.z;
		lenx = len.x;
		leny = len.y;
		lenz = len.z;

		i = j = k = 0;

		if (is_empty()) lenx = leny = lenz = 0;
	}
	
	void Yee_Iterator::next() {
		if (++i == lenx) {
			i = 0;
			if (++j == leny) {
				j = 0;
				++k;
			}
		}
	}
	
	bool Yee_Iterator::is_end() const{
		return k >= lenz;
	}
	
	bool Yee_Iterator::is_empty() const{
		return lenx <= 0 || leny <= 0 || lenz <= 0;
	}
	
	size_t Yee_Iterator::get_size() const{
		return (size_t)lenx * leny * lenz;
	}
	
	iVec3 Yee_Iterator::get_dim() const {
		return iVec3(lenx, leny, lenz);
	}

	
	iVec3 Yee_Iterator::get_coord() const {
		return {x0 + i * stride, y0 + j * stride, z0 + k * stride};
	}

	//Yee3
	Yee3::Yee3(iVec3 ghost_p1, iVec3 ghost_p2):
		ghost_p1(ghost_p1), ghost_p2(ghost_p2) {

		dim = ghost_p2 - ghost_p1 + 1;
		grid_p1 = ghost_p1 + 1;
		grid_p2 = ghost_p2 - 1;

		stride.z = dim.x * dim.y;
		stride.y = dim.x;
		stride.x = 1;
	}

	iVec3 Yee3::get_grid_p1() const {
		return grid_p1;
	}

	iVec3 Yee3::get_grid_p2() const {
		return grid_p2;
	}

	size_t Yee3::get_size() const {
		return (size_t)dim.x * dim.y * dim.z;
	}

	sVec3 Yee3::get_stride() const {
		return stride;
	}

	iVec3 Yee3::get_base_point(const fVec3& pos, Coord_Type ctype) const {
		return get_nearest_point<Negative>(pos, ctype);
	}

	iVec3 Yee3::get_coord_from_index(size_t index) const {

		int x, y, z;
		x = ghost_p1.x + (index % dim.x);
		index /= dim.x;
		y = ghost_p1.y + (index % dim.y);
		z = ghost_p1.z + (index / dim.y);
		return {x, y, z};
	}

	long long Yee3::get_index_offset(int i, int j, int k) const {

		return (long long)i * stride.x + (long long)j * stride.y + (long long)k * stride.z;
	}

	long long Yee3::get_index_offset(const iVec3& pos) const {

		return get_index_offset(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(const iVec3& pos) const {

		return get_index_from_coord(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(int i, int j, int k) const {

		return (i - ghost_p1.x) * stride.x + (j - ghost_p1.y) * stride.y + (k - ghost_p1.z) * stride.z;
	}

	bool Yee3::is_inside(const fVec3& pos) const {
		
		return leq_vec3(grid_p1, pos) && leq_vec3(pos, grid_p2);
	}

	bool Yee3::is_inside(int i, int j, int k) const {

		return	i >= grid_p1.x && i <= grid_p2.x &&
				j >= grid_p1.y && j <= grid_p2.y &&
				k >= grid_p1.z && k <= grid_p2.z;
	}

	bool Yee3::is_inside(const iVec3& pos) const {
		return	is_inside(pos.x, pos.y, pos.z);
	}


	std::vector<double> Yee3::get_interp_weights(double sx, double sy, double sz) const {
		return {(1 - sx) * (1 - sy) * (1 - sz),
				sx * (1 - sy) * (1 - sz),
				(1 - sx) * sy * (1 - sz),
				sx * sy * (1 - sz),
				(1 - sx) * (1 - sy) * sz,
				sx * (1 - sy) * sz,
				(1 - sx) * sy * sz,
				sx * sy * sz};
	}

	std::vector<double> Yee3::get_interp_weights(const fVec3& s) const {
		return get_interp_weights(s.x, s.y, s.z);
	}

	std::vector<double> Yee3::get_interp_weights
	(fVec3 pos, Coord_Type ctype) const {

		iVec3 base = get_base_point(pos, ctype);
		double sx = (pos.x - base.x) / 2;
		double sy = (pos.y - base.y) / 2;
		double sz = (pos.z - base.z) / 2;
		
		return get_interp_weights(sx, sy, sz);
	}
	
	//Gaussian functions
	double Gaussian1(double t, double width) {
		return t * std::exp(-(t * t / 2 / width / width));
	}

	double Gaussian2(double t, double width) {
		double arg = (t / width);
		arg = arg * arg;

		return (1 - arg) * exp(-arg);
	}

	//Output overloading
	std::ostream& operator<<(std::ostream& os, const std::complex<double>& c) {
		os << c.real() << " " << c.imag();
		return os;
	}
}
