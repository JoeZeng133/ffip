#include <utility.hpp>

namespace ffip {
    const double pi = 3.141592653589793e+00;

	/*
		Grid3
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
		my_iteartor
	*/
	my_iterator::my_iterator(const std::pair<iVec3, iVec3> corners, const Coord_Type ctype): my_iterator{corners.first, corners.second, ctype} {}
	
	my_iterator::my_iterator(const fVec3& p1, const fVec3& p2, const Coord_Type ctype): my_iterator{get_component_interior(p1, p2, ctype), ctype} {}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype) {
		std::pair<iVec3, iVec3> tmp;
		
		// loop through points inside [p1, p2]
		if (ctype == All) {
			tmp = get_component_interior(p1, p2, ctype);
			stride = 1;
		}
		else if (ctype == None) {
			throw std::runtime_error("Invalid Ctype");
		}
		else {
			tmp = get_component_interior(p1, p2, ctype);
			stride = 2;
		}
		
		x = x0 = tmp.first.x;
		y = y0 = tmp.first.y;
		z = z0 = tmp.first.z;
		x1 = tmp.second.x;
		y1 = tmp.second.y;
		z1 = tmp.second.z;
		
		index = 0;
		end = size = get_size();
	}
	
	my_iterator::my_iterator(const iVec3& p1, const iVec3& p2, const Coord_Type ctype, const size_t rank, const size_t num): my_iterator(p1, p2, ctype) {
		if (rank >= num)
			throw std::runtime_error("Rank is larger than number of threads");
		
		if (size <= num)
			throw std::runtime_error("Region is too small to divide, Better use non-parallel version");
		
		size_t idx1 = (rank * size) / num;
		size_t idx2 = ((rank + 1) * size) / num;
		
		index = idx1;
		end = idx2;
		
		sVec3 dim = get_dim();
		x = x0 + (idx1 % dim.x) * stride;
		idx1 /= dim.x;
		y = y0 + (idx1 % dim.y) * stride;
		idx1 /= dim.y;
		z = z0 + (idx1 % dim.z) * stride;
	}
	
	iVec3 my_iterator::get_vec(size_t index) const {
		int x, y, z;
		x = x0 + (index % ((x1 - x0) / stride + 1)) * stride;
		index /= (x1 - x0) / stride + 1;
		y = y0 + (index % ((y1 - y0) / stride + 1)) * stride;
		index /= ((y1 - y0) / stride + 1);
		z = z0 + (index % ((z1 - z0) / stride + 1)) * stride;
		return {x, y, z};
	}
	
	void my_iterator::advance() {
		index++;
		if((x += stride) > x1) {
			x = x0;
			if((y += stride) > y1) {
				y = y0;
				z += stride;
			}
		}
	}
	
	bool my_iterator::is_end() const{
		return index >= end;
	}
	
	bool my_iterator::is_empty() const{
		return x0 > x1 || y0 > y1 || z0 > z1;
	}
	
	size_t my_iterator::get_size() const{
		return is_empty()? 0 : (size_t)((x1 - x0) / stride + 1) * ((y1 - y0) / stride + 1) * ((z1 - z0) / stride + 1);
	}
	
	sVec3 my_iterator::get_dim() const {
		return is_empty()? Vec3<size_t>(0, 0, 0) : Vec3<size_t>((x1 - x0) / stride + 1, (y1 - y0) / stride + 1, (z1 - z0) / stride + 1);
	}
	
	iVec3 my_iterator::get_vec() const {
		return {x, y, z};
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
