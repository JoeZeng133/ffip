#include <utility.hpp>
#include <source.hpp>

namespace ffip {
	/* Plane_Wave definitions*/
	Plane_Wave::Plane_Wave(const real _dx, const real _dt, const int _n): dx(_dx), dt(_dt), n(_n) {}
	
	void Plane_Wave::set_excitation(const Sinuosuidal_Func& _f, const real phys_excite_p, const real _amp, const Direction _polorization) {
		
		excite_p = phys_excite_p / dx;
		f = _f;
		polorization = _polorization;
		amp = _amp;
		
		if (polorization == Direction::Z)
			throw std::runtime_error("Z polorization not allowed");
	}
	
	void Plane_Wave::set_PML_neg(const PML &_PML_neg) {
		PML_neg = _PML_neg;
	}
	
	void Plane_Wave::set_PML_pos(const PML &_PML_pos) {
		PML_pos = _PML_pos;
	}
	
	void Plane_Wave::hard_E(real time)  {
		real val = f(time) * amp;
		eh[excite_p1] = c1 * val;
		eh[excite_p2] = c2 * val;
	}
	
	void Plane_Wave::update(real time) {
		update_H();
		update_E();
		hard_E(time + dt);
	}
	
	void Plane_Wave::init() {
		if (excite_p > n || excite_p < 0)
			throw std::runtime_error("Excitement position cannot be outside of domain");
		
		courant = c * dt / dx;
		if (courant > 1)
			throw std::runtime_error("Courant number has to be less than 1");
		n_neg = PML_neg.get_d();
		n_pos = PML_pos.get_d();
		n_tot = n + n_neg + n_pos;
		
		excite_p1 = floor(excite_p);
		excite_p2 = ceil(excite_p);
		c1 = excite_p2 - excite_p;
		c2 = excite_p - excite_p1;
		excite_p1 = (excite_p1 + n_neg) * 2;
		excite_p2 = (excite_p2 + n_neg) * 2;
		
		eh.resize(n_tot * 2 + 1, 0);
		psi_pos.resize(n_tot * 2 + 1, 0);
		psi_neg.resize(n_tot * 2 + 1, 0);
		
		k_z.resize(n_tot * 2 + 1, 1);
		b_z.resize(n_tot * 2 + 1, 1);
		c_z.resize(n_tot * 2 + 1, 0);
		
		for(int i = 0; i <= n_neg * 2; ++i) {
			real x = n_neg - i / 2.0;
			b_z[i] = PML_neg.get_b(x);
			c_z[i] = PML_neg.get_c(x);
			k_z[i] = PML_neg.get_chi(x);
		}
		
		for(int i = 0; i <= n_pos * 2; ++i) {
			real x = i / 2.0;
			b_z[i + (n + n_neg) * 2] = PML_pos.get_b(x);
			c_z[i + (n + n_neg) * 2] = PML_pos.get_c(x);
			k_z[i + (n + n_neg) * 2] = PML_pos.get_chi(x);
		}
	}
	
	void Plane_Wave::update_H(){
		for (int i = 1; i < 2 * n_tot; i += 2) {
			psi_pos[i] = b_z[i] * psi_pos[i] + c_z[i] * (eh[i + 1] - eh[i - 1]) / dx;
			eh[i] += -((eh[i + 1] - eh[i - 1]) / dx / k_z[i] + psi_pos[i]) / u0 * dt;
		}
	}
	
	void Plane_Wave::update_E() {
		for (int i = 2; i < 2 * n_tot; i += 2) {
			psi_neg[i] = b_z[i] * psi_neg[i] + c_z[i] * (eh[i + 1] - eh[i - 1]) / dx;
			eh[i] += -((eh[i + 1] - eh[i - 1]) / dx / k_z[i] + psi_neg[i]) / e0 * dt;
		}
	}
	
	real Plane_Wave::operator()(const iVec3 &p) const {
		if(p.z < 0 || p.z > (n << 1))
			throw std::runtime_error("Access exceeds domain");
		
		if(polorization == Direction::X) {
			if ((p.x & 1) == 1 && (p.y & 1) == 0)
				return eh[p.z + (n_neg << 1)];
		} else {
			if ((p.x & 1) == 0 && (p.y & 1) == 1)
				return eh[p.z + (n_neg << 1)];
		}
		
		return 0;
	}
	
	/* TFSF_Surface definitions*/
	iVec3 TFSF_Surface::vec3_base[3] = {
		{1, 0, 0},
		{0, 1, 0},
		{0, 0, 1}
	};
	
	TFSF_Surface::TFSF_Surface(const iVec3& _d1, const iVec3& _d2, const Direction _dir, const Side _side, int _sign_correction): d1(_d1), d2(_d2), dir(_dir), side(_side), sign_correction(_sign_correction) {
		
		if(side != Side::High && side != Side::Low)
			throw std::runtime_error("Invalid Side");
		
		if(dir != Direction::X && dir != Direction::Y && dir != Direction::Z)
			throw std::runtime_error("Invalid TFSF direction");
		
		if(sign_correction != -1 && sign_correction != 1)
			throw std::runtime_error("Invalid sign correction");
	}
	
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
	
	/* Current_Source definitions */
	Current_Source::Current_Source(const GriddedInterp& _interp, const Sinuosuidal_Func& _phase, const Coord_Type _ctype): interp(_interp), phase(_phase), ctype(_ctype) {}
	
	Current_Source::Current_Source(const std::string& filename, const Sinuosuidal_Func& _phase, const Coord_Type _ctype): interp(filename), phase(_phase), ctype(_ctype) {}
	
	Source_Internal* Current_Source::get_source_internal() {
		return new Current_Internal{phase, amp, ch_dim, ch_origin, p1, p2};
	}
	
	void Current_Source::init(const real dx,
							  const iVec3 _ch_p1, const iVec3 _ch_p2,
							  const iVec3 _ch_dim, const iVec3 _ch_origin) {
		ch_p1 = _ch_p1;
		ch_p1 = _ch_p2;
		ch_dim = ch_dim;
		ch_origin = _ch_origin;
		
		auto to_phys = [ratio = (dx / 2)](const auto& p) -> fVec3
		{return {p.x * ratio, p.y * ratio, p.z * ratio};};			//lambda converter from comp to phys
		
		auto to_comp = [ratio = (2 / dx)](const auto& p) -> fVec3
		{return {p.x * ratio, p.y * ratio, p.z * ratio};};			//lambda converter from phys to comp
		
		
		/* expanding dimensions, this is neccessary to model point, line and surface source and treat all of them as volumetric current source*/
		auto fp1 = interp.get_p1();
		auto fp2 = interp.get_p2();
		
		auto closure = get_component_closure(to_comp(fp1), to_comp(fp2), ctype);	//enclose the area with ctype points
		auto fclosure_p1 = to_phys(closure.first);		//convert to physical coordinates
		auto fclosure_p2 = to_phys(closure.second);
		
		/* expand region to its closure. */
		if(fp1.x == fp2.x) {
			interp.expand_dim(fclosure_p1.x, fclosure_p2.x, Direction::X);
		}
		
		if(fp1.y == fp2.y) {
			interp.expand_dim(fclosure_p1.y, fclosure_p2.y, Direction::Y);
		}
		
		if(fp1.z == fp2.z) {
			interp.expand_dim(fclosure_p1.z, fclosure_p2.z, Direction::Z);
		}
		
		/* initialization of amp, find the interp_ratio*/
		fp1 = interp.get_p1();
		fp2 = interp.get_p1();
		closure = get_component_closure(to_comp(fp1), to_comp(fp2), ctype);
		p1 = closure.first;
		p2 = closure.second;
		
		real integral = 0;
		for(int i = p1.x; i <= p2.x; i += 2)
			for(int j = p1.y; j <= p2.y; j += 2)
				for(int k = p1.z; k <= p2.z; k += 2) {
					real temp = interp.request_value(to_phys(iVec3{i, j, k}));
					integral += temp;
				}
		integral *= dx * dx * dx;
		interp_ratio = integral / interp.request_integral();
		
		/* considers only the intersection between chunk and current region*/
		auto intersection = get_intersection(p1, p2, ch_p1, ch_p2);
		p1 = intersection.first;
		p2 = intersection.second;
		for(int i = p1.x; i <= p2.x; i += 2)
			for(int j = p1.y; j <= p2.y; j += 2)
				for(int k = p1.z; k <= p2.z; k += 2) {
					real temp = interp.request_value(to_phys(iVec3{i, j, k}));
					if (ctype == Coord_Type::Ex || ctype == Coord_Type::Ey || ctype == Coord_Type::Ez)		//Jd = curl - Js, Md = curl + Ms
						temp = -temp;
					
					amp.push_back(temp * interp_ratio);
				}
	}
	
	/* Current_Internal definitions*/
	Current_Internal::Current_Internal(const Sinuosuidal_Func _phase, const std::vector<real>& _amp,
									   const iVec3 _ch_dim, const iVec3 _ch_origin,
									   const iVec3 _p1, const iVec3 _p2):
	phase(_phase), amp(_amp), ch_dim(_ch_dim), ch_origin(_ch_origin), p1(_p1), p2(_p2) {
		
		ctype = p1.get_type();
		
		if(ctype != p2.get_type())
			throw std::runtime_error("Corner points has to be the same type");
		
		auto p1_ch = p1 - ch_origin;
		ch_jump_x = ch_dim.y * ch_dim.z;
		ch_jump_y = ch_dim.z;
		ch_jump_z = 1;
		base_index_ch = p1_ch.x * ch_jump_x + p1_ch.y * ch_jump_y + p1_ch.z * ch_jump_z;	//p1's chunk index
		
		dp = p2 - p1;
		amp_jump_x = (dp.y / 2 + 1) * (dp.z / 2 + 1);
		amp_jump_y = dp.z / 2 + 1;
		amp_jump_z = 1;
	}
	
	void Current_Internal::update_helper(std::vector<real> &jmd) {
		for(int i = 0; i <= dp.x; i += 2)
			for(int j = 0; j <= dp.y; j += 2)
				for(int k = 0; k <= dp.z; k += 2) {
					int index_ch = base_index_ch + i * ch_jump_x + j * ch_jump_y + k * ch_jump_z;
					int index_amp = (i >> 1) * amp_jump_x + (j >> 1) * amp_jump_y + (k >> 1) * amp_jump_z;
					
					jmd[index_ch] += amp[index_amp] * cur_phase;
				}
	}
	
	void Current_Internal::update_Jd(std::vector<real> &jmd) {
		if(ctype != Coord_Type::Ex && ctype != Coord_Type::Ey && ctype != Coord_Type::Ez) return;
		update_helper(jmd);
	}
	
	void Current_Internal::update_Md(std::vector<real> &jmd) {
		if(ctype != Coord_Type::Hx && ctype != Coord_Type::Hy && ctype != Coord_Type::Hz) return;
		update_helper(jmd);
	}
	
	void Current_Internal::get_Jd(real time) {
		cur_phase = phase(time);
	}
	
	void Current_Internal::get_Md(real time) {
		cur_phase = phase(time);
	}
	
	/* Eigen_Source definitions*/
	Eigen_Source::Eigen_Source(const Plane_Wave& _projector): projector(_projector) {}
	
	void Eigen_Source::init(const iVec3 _tf1, const iVec3 _tf2,
							const iVec3& _ch_dim, const iVec3& _ch_origin,
							const iVec3 _ch_p1, const iVec3 _ch_p2) {
		ch_p1 = _ch_p1;
		ch_p2 = _ch_p2;
		tf1 = _tf1;
		tf2 = _tf2;
		ch_dim = _ch_dim;
		ch_origin = _ch_origin;
	}
	
	Source_Internal* Eigen_Source::get_source_internal() {
		std::vector<TFSF_Surface> tfsf_list;
		
		/* generate TF surfaces*/
		tfsf_list.push_back(TFSF_Surface{iVec3{tf2.x, tf1.y, tf1.z}, iVec3{tf2.x, tf2.y, tf2.z}, Direction::X, Side::High, 1});	//x+
		tfsf_list.push_back(TFSF_Surface{iVec3{tf1.x, tf1.y, tf1.z}, iVec3{tf1.x, tf2.y, tf2.z}, Direction::X, Side::Low, 1}); //x-
		tfsf_list.push_back(TFSF_Surface{iVec3{tf1.x, tf2.y, tf1.z}, iVec3{tf2.x, tf2.y, tf2.z}, Direction::Y, Side::High, 1}); //y+
		tfsf_list.push_back(TFSF_Surface{iVec3{tf1.x, tf1.y, tf1.z}, iVec3{tf2.x, tf1.y, tf2.z}, Direction::Y, Side::Low, 1}); //y-
		tfsf_list.push_back(TFSF_Surface{iVec3{tf1.x, tf1.y, tf2.z}, iVec3{tf2.x, tf2.y, tf2.z}, Direction::Z, Side::High, 1});	//z+
		tfsf_list.push_back(TFSF_Surface{iVec3{tf1.x, tf1.y, tf1.z}, iVec3{tf2.x, tf2.y, tf1.z}, Direction::Z, Side::Low, 1});	//z-
		
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
		
		return new Eigen_Internal{tfsf_list, projector, ch_dim, ch_origin};
	}
	
	/* Eigen_Internal definitions*/
	int Eigen_Internal::TFSF_Mat[3][3] = {{0, -1, 1}, {1, 0, -1}, {-1, 1, 0}};
	
	Eigen_Internal::Eigen_Internal(const std::vector<TFSF_Surface> _tsfs_list, const Plane_Wave& _projector,
								   const iVec3& _ch_dim, const iVec3& _ch_origin): tfsf_list(_tsfs_list), projector(_projector), ch_dim(_ch_dim), ch_origin(_ch_origin) {
		
		jump_x = ch_dim.y * ch_dim.z;
		jump_y = ch_dim.z;
		jump_z = 1;
		projector.init();
	}
	
	void Eigen_Internal::update_helper(std::vector<real> &jmd, const TFSF_Surface face, Coord_Type type) {
		/* return if type is parellel to face direction*/
		switch (face.dir) {
			case Direction::X:
				if(type == Coord_Type::Ex || type == Coord_Type::Hx)
					return;
				break;
				
			case Direction::Y:
				if(type == Coord_Type::Ey || type == Coord_Type::Hy)
					return;
				break;
				
				
			case Direction::Z:
				if(type == Coord_Type::Ez || type == Coord_Type::Hz)
					return;
				break;
				
			default:
				throw std::runtime_error("Direction of the TF/SF surface is illegal");
				break;
		}
		
		auto interior = get_component_interior(face.d1, face.d2, type);
		if(ElementWise_Less_Eq(interior.first, interior.second))	//return if the interior is empty
			return;
		
		int side = static_cast<int>(face.side);
		int dir = static_cast<int>(face.dir);
		int type_int = static_cast<int>(type);
		iVec3 pos_offset = face.get_pos_offset();
		iVec3 c1 = interior.first - ch_origin;			//chunk coordinate
		iVec3 c2 = interior.second - ch_origin;
		iVec3 d1 = interior.first;						//domain coordinate
		
		//loop over the coord type with chunk, domain coordinates updating on the fly
		for(int i = c1.x, gi = d1.x; i <= c2.x; i += 2, gi += 2)
			for(int j = c1.y, gj = d1.y; j <= c2.y; j += 2, gj += 2)
				for (int k = c1.z, gk = d1.z; k <= c2.z; k += 2, gk += 2) {
					int index = i * jump_x + j * jump_y + k * jump_z;
					
					jmd[index] += side * TFSF_Mat[dir][type_int] * projector(iVec3{gi + pos_offset.x, gj + pos_offset.y, gk + pos_offset.z});
				}
	}
	
	void Eigen_Internal::update_Jd(std::vector<real> &jmd) {
		for(auto& face : tfsf_list) {
			update_helper(jmd, face, Coord_Type::Ex);
			update_helper(jmd, face, Coord_Type::Ey);
			update_helper(jmd, face, Coord_Type::Ez);
		}
	}
	
	void Eigen_Internal::update_Md(std::vector<real> &jmd) {
		for(auto& face : tfsf_list) {
			update_helper(jmd, face, Coord_Type::Hx);
			update_helper(jmd, face, Coord_Type::Hy);
			update_helper(jmd, face, Coord_Type::Hz);
		}
	}
	
	void Eigen_Internal::get_Jd(real time) {
		projector.update(time);	//this updates
	}
	
	void Eigen_Internal::get_Md(real time) {}
	
}
