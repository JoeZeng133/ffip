#include <utility.hpp>
#include <source.hpp>

namespace ffip {
	/* Plane_Wave definitions*/
	Plane_Wave::Plane_Wave(const real _dx, const real _dt, const int _dim_neg, const int _dim_pos) : dx(_dx), dt(_dt), dim_neg(_dim_neg), dim_pos(_dim_pos) {}
	
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
	
	void Plane_Wave::advance() {
		update_H();
		update_E();
		hard_E((++time_step) * dt);
		udf_advance();
	}
	
	int Plane_Wave::get_time_step() {
		return time_step;
	}
	
	std::fstream os;
	
	void Plane_Wave::udf_init() {
		/*os = std::fstream{"projector.out", std::ios::out};*/
	}
	
	void Plane_Wave::udf_advance() {
		/*if (time_step % 10 != 0) return;
		for (int i = 2; i < dim; i += 2) {
			os << eh[i] << " ";
		}
		os << "\n";*/
	}
	
	void Plane_Wave::init() {
		courant = c0 * dt / dx;
		if (courant > 1)
			throw std::runtime_error("Courant number has to be less than 1");
		
		n = dim_neg + dim_pos;
		dim = (n + pml.get_d()) * 2;	//[-dim_neg, dim_pos + pml] * dx
		origin = dim_neg * 2;
		
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
		udf_init();
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
		int z = p.z + origin;
		Coord_Type ctype = p.get_type();
		if (z < 0 || z > n * 2)
			throw Out_of_the_Domain{};
		
		if(polorization == Direction::X) {
			if (ctype == Ex || ctype == Hy)
				return eh[z];
		} else {
			if (ctype == Ey || ctype == Hx)
				return ctype == Hx ? -eh[z] : eh[z];
		}
		
		return 0;
	}
	
	real Plane_Wave::operator()(const fVec3 &p, const Coord_Type ctype) const {
		//interpolate to get field intensity, x, y coordinates are neglected
		iVec3 p1 = get_nearest_point<side_low_tag>(p, ctype); 
		real f = (p.z - p1.z) / 2;
		real val1 = operator()(p1);
		p1.z += 2;
		real val2 = operator()(p1);

		return val1 * (1 - f) + val2 * f;
	}
	
	real Plane_Wave::at(const fVec3 &p, const Coord_Type ctype) const {
		return operator()(p / (dx / 2), ctype);
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
	
	/* Projector_Source definitions*/
	void Projector_Source::add_projector(const Plane_Wave& projector) {
		projectors.push_back(new Plane_Wave{ projector });
	}

	void Projector_Source::init(const Config& config) {
		this->config = config;

		auto tf1 = config.tf_p1;
		auto tf2 = config.tf_p2;
		auto sf_layer1 = tf1 - config.phys_p1;
		auto sf_layer2 = config.phys_p2 - tf2;

		/* generate TF surfaces*/
		if (sf_layer2.x)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_x_tag, side_high_tag>(tf1, tf2), Direction::X, Side::High, 1 });	//x+
		if (sf_layer1.x)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_x_tag, side_low_tag>(tf1, tf2), Direction::X, Side::Low, 1 }); //x-
		if (sf_layer2.y)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_y_tag, side_high_tag>(tf1, tf2), Direction::Y, Side::High, 1 }); //y+
		if (sf_layer1.y)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_y_tag, side_low_tag>(tf1, tf2), Direction::Y, Side::Low, 1 }); //y-
		if (sf_layer2.z)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_z_tag, side_high_tag>(tf1, tf2), Direction::Z, Side::High, 1 });	//z+
		if (sf_layer1.z)
			tfsf_list.push_back(TFSF_Surface{ get_face<dir_z_tag, side_low_tag>(tf1, tf2), Direction::Z, Side::Low, 1 });	//z-

		for (int i = 0, len = tfsf_list.size(); i < len; ++i) {
			auto tmp = tfsf_list[i];
			tmp.TF2SF();
			tfsf_list.push_back(tmp);
		}

		for (auto pr : projectors)
			pr->init();
	}

	void Projector_Source::advance() {
		for (auto pr : projectors)
			pr->advance();
	}
	
	void Projector_Source::push(Chunk* chunk) {
		for (auto pr : projectors)
			chunk->add_projector(pr);

		for (auto& face : tfsf_list) {
			push_helper<dir_x_tag>(chunk, face);
			push_helper<dir_y_tag>(chunk, face);
			push_helper<dir_z_tag>(chunk, face);
		}
	}
}
