#include <utility.hpp>

namespace ffip {
    const double pi = 3.141592653589793e+00;

	/*
		Yee3
	*/
	Yee3::Yee3(iVec3 grid_p1, iVec3 grid_p2, fVec3 phys_p1, double dx):

	grid_p1(grid_p1), grid_p2(grid_p2), phys_p1(phys_p1), dx(dx) {

		phys_p2 = phys_p1 + (grid_p2 - grid_p1) * dx / 2;
		dim = grid_p2 - grid_p1 + 1;

		stride.z = dim.x * dim.y;
		stride.y = dim.x;
		stride.x = 1;
	}

	iVec3 Yee3::get_coord_from_index(size_t index) const {

		int x, y, z;
		x = grid_p1.x + (index % dim.x);
		index /= dim.x;
		y = grid_p1.y + (index % dim.y);
		z = grid_p1.z + (index / dim.y);
		return {x, y, z};
	}

	size_t Yee3::get_index_offset(int i, int j, int k) const {

		return i * stride.x + j * stride.y + k * stride.z;
	}

	size_t Yee3::get_index_offset(const iVec3& pos) const {

		return get_index_offset(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(const iVec3& pos) const {

		return get_index_from_coord(pos.x, pos.y, pos.z);
	}

	size_t Yee3::get_index_from_coord(int i, int j, int k) const {

		return (i - grid_p1.x) * stride.x + (j - grid_p1.y) * stride.y + (k - grid_p1.z) * stride.z;
	}

	bool Yee3::is_in_grid(bool is_grid_coordinate, const fVec3& pos) const {
		
		if (!is_grid_coordinate)
			return Is_Inside_Box(grid_p1, grid_p2, pos);
		else
			return Is_Inside_Box(phys_p1, phys_p2, pos);
	}

	fVec3 Yee3::phys_coord2grid_coord(const fVec3& pos) const {

		return (pos - phys_p1) * (2 / dx);
	}

	std::vector<std::pair<size_t, double>> Yee3::get_interp_info
	(bool is_grid_coordinate, fVec3 pos, Coord_Type ctype) {

		if (!is_grid_coordinate)
				pos = phys_coord2grid_coord(pos);

		if(!Is_Inside_Box(grid_p1, grid_p2, pos))
			return {};

		iVec3 base = get_nearest_point<-1>(pos, ctype);

		double sx = (pos.x - base.x) / 2;
		double sy = (pos.y - base.y) / 2;
		double sz = (pos.z - base.z) / 2;
		size_t base_index = get_index_from_coord(base);

		return {{get_index_offset(0, 0, 0),	(1 - sx) * (1 - sy) * (1 - sz)},
				{get_index_offset(2, 0, 0), sx * (1 - sy) * (1 - sz)},
				{get_index_offset(0, 2, 0),	(1 - sx) * sy * (1 - sz)},
				{get_index_offset(2, 2, 0), sx * sy * (1 - sz)},
				{get_index_offset(0, 0, 2), (1 - sx) * (1 - sy) * sz},
				{get_index_offset(2, 0, 2), sx * (1 - sy) * sz},
				{get_index_offset(0, 2, 2), (1 - sx) * sy * sz },
				{get_index_offset(2, 2, 2), sx * sy * sz}};
	}
	
	/* 
		Yee_Iterator
	*/
	Yee_Iterator::Yee_Iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype): Yee_Iterator{corners.first, corners.second, ctype} {}
	
	Yee_Iterator::Yee_Iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype) {
		std::pair<iVec3, iVec3> tmp;
		
		// loop through points inside [p1, p2]
		tmp = get_component_interior(p1, p2, ctype);
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

	/*
		Gaussian functions
	*/
	double Gaussian1(double t, double width) {
		return t * std::exp(-(t * t / 2 / width / width));
	}

	double Gaussian2(double t, double width) {
		double arg = (t / width);
		arg = arg * arg;

		return (1 - arg) * exp(-arg);
	}

	/*
		operaor overloading
	*/
	std::ostream& operator<<(std::ostream& os, const std::complex<double>& c) {
		os << c.re() << " " << c.im();
	}
}
