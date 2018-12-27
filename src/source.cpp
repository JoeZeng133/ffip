#include <utility.hpp>
#include <source.hpp>

namespace ffip {
	/* Plane_Wave definitions*/
	Plane_Wave::Plane_Wave(const real _dx, const real _dt, const int _n): dx(_dx), dt(_dt), n(_n) {}
	
	void Plane_Wave::set_excitation(const std::function<real(const real)>& _f, const real _amp, const Direction _polorization) {
		f = _f;
		polorization = _polorization;
		amp = _amp;
		
		if (polorization == Direction::Z)
			throw std::runtime_error("Z polorization not allowed");
	}
	
	void Plane_Wave::set_PML(const PML &_pml) {
		pml = _pml;
	}
	
	void Plane_Wave::hard_E(real time)  {
		real val = f(time) * amp;
		eh[0] = val;
	}
	
	void Plane_Wave::advance(std::ostream& os) {
		update_H();
		update_E();
		hard_E((++time_step) * dt);
		
		//os << operator()(fVec3(1, 0, n), Ex) << "\n";
	}
	
	int Plane_Wave::get_time_step() {
		return time_step;
	}
	
	void Plane_Wave::init() {
		courant = c0 * dt / dx;
		if (courant > 1)
			throw std::runtime_error("Courant number has to be less than 1");
		
		dim = int((n + pml.get_d() + 2) * 2);	//[-0.5, n + 0.5] padded to be [-1, n + 1 + pml]
		origin = 2;
		
		eh.resize(dim + 1, 0);
		psi_pos.resize(dim + 1, 0);
		psi_neg.resize(dim + 1, 0);
		
		k_z.resize(dim + 1, 1);
		b_z.resize(dim + 1, 1);
		c_z.resize(dim + 1, 0);
		
		for(int i = 0; i < pml.get_d() * 2; ++i) {
			real x = pml.get_d() - i / 2.0;
			
			b_z[dim - i] = pml.get_b(x, dt);
			c_z[dim - i] = pml.get_c(x, dt);
			k_z[dim - i] = pml.get_k(x);
		}
		
		time_step = 0;
	}
	
	void Plane_Wave::set_medium(const real _er, const real _ur) {
		er = _er;
		ur = _ur;
	}
	
	void Plane_Wave::update_H(){
		for (int i = 1; i < dim; i += 2) {
			psi_pos[i] = b_z[i] * psi_pos[i] + c_z[i] * (eh[i + 1] - eh[i - 1]) / dx;
			eh[i] += -((eh[i + 1] - eh[i - 1]) / dx / k_z[i] + psi_pos[i]) / (ur * u0) * dt;
		}
	}
	
	void Plane_Wave::update_E() {
		for (int i = 2; i < dim; i += 2) {
			psi_neg[i] = b_z[i] * psi_neg[i] + c_z[i] * (eh[i + 1] - eh[i - 1]) / dx;
			eh[i] += -((eh[i + 1] - eh[i - 1]) / dx / k_z[i] + psi_neg[i]) / (e0 * er) * dt;
		}
	}
	
	real Plane_Wave::operator()(const iVec3 &p) const {
		if (p.z < -1 || p.z > (2 * n + 1))
			throw std::runtime_error("Access exceeds domain");
		
		if(polorization == Direction::X) {
			if ((p.x & 1) && !(p.y & 1))
				return eh[p.z + origin];
		} else {
			if (!(p.x & 1) && (p.y & 1))
				return eh[p.z + origin];
		}
		
		return 0;
	}
	
	real Plane_Wave::operator()(const fVec3 &p, const Coord_Type ctype) const {
		return ((*this)(get_nearest_point<side_low_tag>(p, ctype)) + (*this)(get_nearest_point<side_high_tag>(p, ctype))) / 2;
	}
	
	real Plane_Wave::at(const fVec3 &p, const Coord_Type ctype) const {
		return operator()(p / (dx / 2), ctype);
	}
	
	/* Current_Source definitions */
	Current_Source::Current_Source(const GriddedInterp& _interp, const std::function<real(const real)>& _phase, const Coord_Type _ctype): interp(_interp), phase(_phase), ctype(_ctype) {
		
		//std::cout << "test phase at the construction of current_source " << phase(0) << std::endl;
	}
	
	Current_Source::Current_Source(const std::string& filename, const std::function<real(const real)>& _phase, const Coord_Type _ctype): interp(filename), phase(_phase), ctype(_ctype) {
	}
	
	Source_Internal* Current_Source::get_source_internal() {
		
		
		return new Current_Internal{phase, amp, ch_dim, ch_origin, p1, p2};
	}
	
	void Current_Source::init(const real dx,
							  const iVec3 _ch_p1, const iVec3 _ch_p2,
							  const iVec3 _ch_dim, const iVec3 _ch_origin) {
		ch_p1 = _ch_p1;
		ch_p2 = _ch_p2;
		ch_dim = _ch_dim;
		ch_origin = _ch_origin;
		
		real original_integral = interp.request_integral();			//the original integral is defined by user (line, surface or volume current), the integral should be invariant when applied to FDTD grid points
		
		interp.scale_xyz((2 / dx), (2 / dx), (2 / dx));				//scale interpolant from physical coordinates to computation coordinates, this allows direct value access from computation coordinates
		
		/* expanding dimensions, this is neccessary to model point, line and surface source and treat all of them as volumetric current source*/
		auto fp1 = interp.get_p1();
		auto fp2 = interp.get_p2();
		auto closure = get_component_closure(fp1, fp2, ctype);	//enclose the area with ctype points
		auto cp1 = closure.first;
		auto cp2 = closure.second;
		
		/* expand region to its closure. */
		if(fp1.x == fp2.x && ch_p1.x != ch_p2.x) {
			interp.expand_dim(cp1.x, cp2.x, Direction::X);
		}
		
		if(fp1.y == fp2.y && ch_p1.y != ch_p2.y) {
			interp.expand_dim(cp1.y, cp2.y, Direction::Y);
		}
		
		if(fp1.z == fp2.z && ch_p1.z != ch_p2.z) {
			interp.expand_dim(cp1.z, cp2.z, Direction::Z);
		}
		
		/* initialization of amp, find the interp_ratio
		   interp_ratio gaurantees that integral of the effective current is equivalent to the original current integral
		 */
		fp1 = interp.get_p1();
		fp2 = interp.get_p2();
		closure = get_component_closure(fp1, fp2, ctype);
		cp1 = closure.first;
		cp2 = closure.second;
		
		real integral = 0;

		for(my_iterator itr{cp1, cp2, ctype}; !itr.is_end(); itr.advance()) {
			real temp = interp.request_value(fVec3(itr.x, itr.y, itr.z));
			integral += temp;
		}
		
		integral *= dx * dx * dx;
		real interp_ratio = original_integral / integral;
		
		/* considers only the intersection between chunk and current region*/
		auto intersection = get_intersection(cp1, cp2, ch_p1, ch_p2);
		p1 = intersection.first;
		p2 = intersection.second;
		
		for(my_iterator itr{p1, p2, ctype}; !itr.is_end(); itr.advance()) {
			real temp = interp.request_value(fVec3(itr.x, itr.y, itr.z));
			if (ctype == Coord_Type::Ex || ctype == Coord_Type::Ey || ctype == Coord_Type::Ez)		//Jd = curl - Js, Md = curl + Ms
				temp = -temp;
			
			if(isnan(temp * interp_ratio))
				throw std::runtime_error("Current amplitidu is nan");

			amp.push_back(temp * interp_ratio);
		}
	}
	
	/* Current_Internal definitions*/
	Current_Internal::Current_Internal(const std::function<real(const real)> _phase, const std::vector<real>& _amp,
									   const iVec3 _ch_dim, const iVec3 _ch_origin,
									   const iVec3 _p1, const iVec3 _p2):
	phase(_phase), amp(_amp), ch_dim(_ch_dim), ch_origin(_ch_origin), p1(_p1), p2(_p2) {
		
		ctype = p1.get_type();
		
		if(ctype != p2.get_type())
			throw std::runtime_error("Corner points has to be the same type");
		
		auto p1_ch = p1 - ch_origin;
		ch_jump_x = 1;
		ch_jump_y = ch_dim.x;
		ch_jump_z = ch_dim.x * ch_dim.y;
		base_index_ch = p1_ch.x * ch_jump_x + p1_ch.y * ch_jump_y + p1_ch.z * ch_jump_z;	//p1's chunk index
		
		dp = p2 - p1;
		amp_jump_x = 1;
		amp_jump_y = dp.x / 2 + 1;
		amp_jump_z = (dp.x / 2 + 1) * (dp.y / 2 + 1);
		
		//std::cout << "testing phase at the contruction of current_internal " << phase(0) << std::endl;
	}
	
	void Current_Internal::update_helper(std::vector<real> &jmd) {
		for(auto itr = my_iterator({0, 0, 0}, dp, Corner); !itr.is_end(); itr.advance()) {
			int index_ch = base_index_ch + itr.x * ch_jump_x + itr.y * ch_jump_y + itr.z * ch_jump_z;
			
			jmd[index_ch] += amp[itr.index] * cur_phase;
		}
		
		//std::cout << "testing inside current_internal" << amp[0] * cur_phase << " " << amp[0] * cur_phase << std::endl;
	}
	
	void Current_Internal::update_Jd(std::vector<real> &jmd, const size_t rank) {
		//std::cout << "testing phase at update_Jd" << phase(0) << std::endl;
		if(ctype != Coord_Type::Ex && ctype != Coord_Type::Ey && ctype != Coord_Type::Ez) return;
		update_helper(jmd);
	}
	
	void Current_Internal::update_Md(std::vector<real> &jmd, const size_t rank) {
		if(ctype != Coord_Type::Hx && ctype != Coord_Type::Hy && ctype != Coord_Type::Hz) return;
		update_helper(jmd);
	}
	
	void Current_Internal::get_Jd(real time, const size_t rank) {
		cur_phase = phase(time);

		//std::cout << cur_phase << std::endl;
	}
	
	void Current_Internal::get_Md(real time, const size_t rank) {
		cur_phase = phase(time);
	}
	
	/* TFSF_Surface definitions*/
	TFSF_Surface::TFSF_Surface(const iVec3& _d1, const iVec3& _d2, const Direction _dir, const Side _side, int _sign_correction): d1(_d1), d2(_d2), dir(_dir), side(_side), sign_correction(_sign_correction) {
		
		if(side != Side::High && side != Side::Low)
			throw std::runtime_error("Invalid TFSF Side");
		
		if(dir != Direction::X && dir != Direction::Y && dir != Direction::Z)
			throw std::runtime_error("Invalid TFSF direction");
		
		if(sign_correction != -1 && sign_correction != 1)
			throw std::runtime_error("Invalid TFSF sign correction");
	}
	
	TFSF_Surface::TFSF_Surface(const std::pair<iVec3, iVec3>& d, const Direction dir, const Side side, int sign_correction): TFSF_Surface(d.first, d.second, dir, side, sign_correction) {}
	
	iVec3 TFSF_Surface::get_pos_offset() const{
		int dir_int = static_cast<int>(dir);
		return vec3_base[dir_int] * (side * sign_correction);
	}
	
	void TFSF_Surface::TF2SF() {
		int dir_int = static_cast<int>(dir);
		iVec3 offset = vec3_base[dir_int] * (int)side;
		
		d1 = d1 + offset;
		d2 = d2 + offset;
		sign_correction = -1;
	}
	
	/* Eigen_Source definitions*/
	Eigen_Source::Eigen_Source(const Plane_Wave& _projector): projector(_projector) {}
	
	void Eigen_Source::init(const iVec3& _tf1, const iVec3& _tf2,
							const iVec3& _ch_dim, const iVec3& _ch_origin,
							const iVec3& _ch_p1, const iVec3& _ch_p2, const real _dx) {
		ch_p1 = _ch_p1;
		ch_p2 = _ch_p2;
		tf1 = _tf1;
		tf2 = _tf2;
		ch_dim = _ch_dim;
		ch_origin = _ch_origin;
		dx = _dx;
	}
	
	Source_Internal* Eigen_Source::get_source_internal() {
		std::vector<TFSF_Surface> tfsf_list;
		
		/* generate TF surfaces*/
		tfsf_list.push_back(TFSF_Surface{get_face<dir_x_tag, side_high_tag>(tf1, tf2), Direction::X, Side::High, 1});	//x+
		tfsf_list.push_back(TFSF_Surface{get_face<dir_x_tag, side_low_tag>(tf1, tf2), Direction::X, Side::Low, 1}); //x-
		tfsf_list.push_back(TFSF_Surface{get_face<dir_y_tag, side_high_tag>(tf1, tf2), Direction::Y, Side::High, 1}); //y+
		tfsf_list.push_back(TFSF_Surface{get_face<dir_y_tag, side_low_tag>(tf1, tf2), Direction::Y, Side::Low, 1}); //y-
		tfsf_list.push_back(TFSF_Surface{get_face<dir_z_tag, side_high_tag>(tf1, tf2), Direction::Z, Side::High, 1});	//z+
		tfsf_list.push_back(TFSF_Surface{get_face<dir_z_tag, side_low_tag>(tf1, tf2), Direction::Z, Side::Low, 1});	//z-
		
		/* generate SF surfaces */
		for(int i = 0; i < 6; ++i) {
			auto tmp = tfsf_list[i];
			tmp.TF2SF();
			tfsf_list.push_back(tmp);
		}
		
		/* intersect with chunk*/
		for(auto& i : tfsf_list) {
			auto tmp = get_intersection(i.d1, i.d2, ch_p1, ch_p2);
			i.d1 = tmp.first;
			i.d2 = tmp.second;
		}
		
		return new Eigen_Internal{tfsf_list, projector, ch_dim, ch_origin, dx};
	}
	
	/* Eigen_Internal definitions*/
	int Eigen_Internal::TFSF_Mat[3][3] = {
		{0, -1, 1},
		{1, 0, -1},
		{-1, 1, 0}};
	
	Eigen_Internal::Eigen_Internal(const std::vector<TFSF_Surface>& _tsfs_list, const Plane_Wave& _projector, const iVec3& _ch_dim, const iVec3& _ch_origin, const real _dx): tfsf_list(_tsfs_list), projector(_projector), ch_dim(_ch_dim), ch_origin(_ch_origin), dx(_dx) {
		
		jump_x = 1;
		jump_y = ch_dim.x;
		jump_z = ch_dim.x * ch_dim.y;
		projector.init();
	}
	
	void Eigen_Internal::update_helper(std::vector<real> &jmd, const TFSF_Surface face, Coord_Type ctype, const size_t rank) {
		/* return if type is parellel to face direction*/
		switch (face.dir) {
			case Direction::X:
				if(ctype == Coord_Type::Ex || ctype == Coord_Type::Hx)
					return;
				break;
				
			case Direction::Y:
				if(ctype == Coord_Type::Ey || ctype == Coord_Type::Hy)
					return;
				break;
				
			case Direction::Z:
				if(ctype == Coord_Type::Ez || ctype == Coord_Type::Hz)
					return;
				break;
				
			default:
				throw std::runtime_error("Direction of the TF/SF surface is illegal");
				break;
		}
		
		auto interior = get_component_interior(face.d1, face.d2, ctype);
		if(!ElementWise_Less_Eq(interior.first, interior.second))	//return if the interior is empty
			return;
		
		int side = static_cast<int>(face.side);
		int dir = static_cast<int>(face.dir);
		int type_dir_int = Ctype2DirInt(ctype);
		iVec3 pos_offset = face.get_pos_offset();
		iVec3 int1_ch = interior.first - ch_origin;			//chunk coordinate
		iVec3 int2_ch = interior.second - ch_origin;
		
		//loop over the coord type with chunk, domain coordinates updating on the fly
		for(auto ch_itr = my_iterator(int1_ch, int2_ch, int1_ch.get_type()), d_itr = my_iterator(interior.first, interior.second, ctype); !ch_itr.is_end(); ch_itr.advance(), d_itr.advance()) {
			int index = ch_itr.x * jump_x + ch_itr.y * jump_y + ch_itr.z * jump_z;
			
			jmd[index] += side * TFSF_Mat[dir][type_dir_int] * projector(d_itr.get_vec() + pos_offset) / dx;
		}
	}
	
	void Eigen_Internal::update_Jd(std::vector<real> &jmd, const size_t rank) {
		for(auto& face : tfsf_list) {
			update_helper(jmd, face, Coord_Type::Ex, rank);
			update_helper(jmd, face, Coord_Type::Ey, rank);
			update_helper(jmd, face, Coord_Type::Ez, rank);
		}
	}
	
	void Eigen_Internal::update_Md(std::vector<real> &jmd, const size_t rank) {
		for(auto& face : tfsf_list) {
			update_helper(jmd, face, Coord_Type::Hx, rank);
			update_helper(jmd, face, Coord_Type::Hy, rank);
			update_helper(jmd, face, Coord_Type::Hz, rank);
		}
	}
	
	void Eigen_Internal::get_Jd(real time, const size_t rank) {
		if (rank == 0)
			projector.advance(std::cout);	//this updates
	}
	
	void Eigen_Internal::get_Md(real time, const size_t rank) {}
	
}
