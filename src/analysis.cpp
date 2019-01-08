#include <analysis.hpp>

namespace ffip {
	/* ################################ */
	Nearfield_Probe::Nearfield_Probe(const real freq, const fVec3& _pos, Chunk const* _chunk): omega(2 * pi * freq), pos(_pos), chunk(_chunk) {
		real phase = 0.5 * chunk->get_dt() * omega;
		correction = complex_num{cos(phase), sin(phase)};
	}

	void Nearfield_Probe::set_chunk(const Chunk *_chunk) {
		chunk = _chunk;
	}
	
	void Nearfield_Probe::update(const int time_step) {
		real phase = -time_step * chunk->get_dt() * omega;
		complex_num exp_omega_n = { cos(phase), sin(phase) };

		ex += (*chunk)(pos, Ex) * exp_omega_n;
		ey += (*chunk)(pos, Ey) * exp_omega_n;
		ez += (*chunk)(pos, Ez) * exp_omega_n;

		hx += (*chunk)(pos, Hx) * exp_omega_n;
		hy += (*chunk)(pos, Hy) * exp_omega_n;
		hz += (*chunk)(pos, Hz) * exp_omega_n;
	}

	void Nearfield_Probe::output(std::ostream &o) const{
		o << std::scientific
		<< ex << " " << ey << " " << ez << " "
		<< hx * correction << " " << hy * correction << " " << hz * correction << std::endl;
	}
	
	/* ################################ */
	Box_Freq_Req::Box_Freq_Req(const fVec3& _p1, const fVec3& _p2, const Coord_Type _F,  const real_arr& _omega_list, Chunk const* _chunk): p1(_p1), p2(_p2), F(_F), omega_list(_omega_list), chunk(_chunk) {
		
		validity_check();
		
		dt = chunk->get_dt();
		dx = chunk->get_dx();
		near_p1 = get_nearest_point<side_low_tag>(p1, F);
		near_p2 = get_nearest_point<side_high_tag>(p2, F);
		
		auto itr = my_iterator(near_p1, near_p2, F);
		dim = itr.get_dim();
		interp = interpn<3>(dim.z, dim.y, dim.x);
		raw_fields.resize(omega_list.size());
		for(auto& item : raw_fields)
			item.resize(itr.size, 0);
		
		if (is_M_point(F))
			step_correction = -0.5;
		else
			step_correction = 0;
	}
	
	void Box_Freq_Req::validity_check() const {
		if (F == Null)
			throw Invalide_Coord_Type{};
		
		if (!Is_Inside_Box(chunk->get_p1(), chunk->get_p2(), p1) || !Is_Inside_Box(chunk->get_p1(), chunk->get_p2(), p2))
			throw Out_of_the_Domain{};
	}

	void Box_Freq_Req::update(const int time_step, const size_t rank, const size_t num_proc) {
		complex_arr fourier_coef;
		for(auto omega : omega_list) {
			real phase = -omega * (time_step + step_correction) * dt;
			fourier_coef.emplace_back(cos(phase), sin(phase));
		}
		
		for(auto itr = my_iterator(near_p1, near_p2, F, rank, num_proc); !itr.is_end(); itr.advance()) {
			real val = (*chunk)(itr.get_vec());

			for(int i = 0; i < omega_list.size(); ++i) {
				raw_fields[i][itr.index] += val * fourier_coef[i];
			} 
		}
	}

	complex_arr Box_Freq_Req::operator()(const fVec3 &p) const{
		if (!ElementWise_Less_Eq(p1, p) || !ElementWise_Less_Eq(p, p2))
			throw Out_of_the_Domain{};
		
		real x = (p.x - near_p1.x) / 2;
		real y = (p.y - near_p1.y) / 2;
		real z = (p.z - near_p1.z) / 2;

		complex_arr res(omega_list.size(), 0);
		for (int f = 0; f < omega_list.size(); ++f)
			res[f] = interp.get(raw_fields[f], z, y, x);

		return res;
	}
	
	/* ################################ */
	N2F_Box::N2F_Box(const fVec3& _p1, const fVec3& _p2, const real_arr& _freq_list, Chunk* const _chunk, const real _c, const real _z):p1(_p1), p2(_p2), omega_list(_freq_list), chunk(_chunk), c(_c), z(_z) {
		
		for(auto& item : omega_list)
			item *= 2 * pi;
		
		//unique frequency
		std::sort(omega_list.begin(), omega_list.end());
		auto it = std::unique(omega_list.begin(), omega_list.end());
		omega_list.resize(std::distance(omega_list.begin(), it));
		
		fVec3 ref_p = (p1 + p2) / 2;			//reference point is taken as the center of the n2f box
		n2f_faces.push_back(new N2F_Face<dir_x_tag>(get_face<dir_x_tag, side_low_tag>(p1, p2), ref_p, Low, omega_list, chunk, c));
		n2f_faces.push_back(new N2F_Face<dir_y_tag>(get_face<dir_y_tag, side_low_tag>(p1, p2), ref_p, Low, omega_list, chunk, c));
		n2f_faces.push_back(new N2F_Face<dir_z_tag>(get_face<dir_z_tag, side_low_tag>(p1, p2), ref_p, Low, omega_list, chunk, c));
		n2f_faces.push_back(new N2F_Face<dir_x_tag>(get_face<dir_x_tag, side_high_tag>(p1, p2), ref_p, High, omega_list, chunk, c));
		n2f_faces.push_back(new N2F_Face<dir_y_tag>(get_face<dir_y_tag, side_high_tag>(p1, p2), ref_p, High, omega_list, chunk, c));
		n2f_faces.push_back(new N2F_Face<dir_z_tag>(get_face<dir_z_tag, side_high_tag>(p1, p2), ref_p, High, omega_list, chunk, c));
	}
	
	N2F_Box::~N2F_Box() {
		for(auto item : n2f_faces)
			delete item;
	}
	
	void N2F_Box::update(const int time_step, const size_t rank, const size_t num_proc) {
		for(auto item : n2f_faces)
			item->update(time_step, rank, num_proc);
	}
	
	std::pair<cVec3, cVec3> N2F_Box::get_NL(const real theta, const real phi, const real omega) const{
		//fetch the index of the frequency
		auto low = std::lower_bound(omega_list.begin(), omega_list.end(), omega);
		if (low == omega_list.end() || *low != omega)
			throw std::runtime_error("Invalid requested frequency");
		
		size_t omega_index = low - omega_list.begin();
		
		cVec3 N{0, 0, 0};
		cVec3 L{0, 0, 0};
		
		for(auto item : n2f_faces) {
			auto tmp = item->get_NL(theta, phi, omega_index);
			N = N + tmp.first;
			L = L + tmp.second;
		}
		
		return {N, L};
	}
	
	void N2F_Box::prepare(const size_t rank, const size_t num_proc) const{
		for(auto item : n2f_faces)
			item->prepare(rank, num_proc);
	}
	
	void N2F_Box::output(std::ostream& os, const real th, const real phi, const real rho, const real freq) const {
	
		real omega = 2 * pi * freq;
		real k = omega / c;		//get the wavenumber
		fVec3 proj_th = {cos(th) * cos(phi), cos(th) * sin(phi), -sin(th)};
		fVec3 proj_phi = {-sin(phi), cos(phi), 0};
		Vec3<complex_num> N{0, 0, 0}, L{0, 0, 0};
		
		auto tmp = get_NL(th, phi, omega);
		
		complex_num Nth = inner_prod(proj_th, tmp.first);
		complex_num Nphi = inner_prod(proj_phi, tmp.first);
		complex_num Lth = inner_prod(proj_th, tmp.second);
		complex_num Lphi = inner_prod(proj_phi, tmp.second);
		
		complex_num arg = k / (4 * pi * rho) * complex_num{sin(k * rho), cos(k * rho)};
		complex_num Eth = -(Lphi + z * Nth) * arg;
		complex_num Ephi = (Lth - z * Nphi) * arg;
		complex_num Hth = (Nphi - Lth / z) * arg;
		complex_num Hphi = -(Nth + Lphi / z) * arg;
		
		os << std::scientific << Eth << " " << Ephi << " " << Hth << " " << Hphi << std::endl;
	}
	
	/* ################################ */
	Flux_Box::Flux_Box(const fVec3& _p1, const fVec3& _p2, const real_arr& _freq_list, Chunk* const _chunk): p1(_p1), p2(_p2), omega_list(_freq_list), chunk(_chunk) {
		
		flux_faces.push_back(new Flux_Face<dir_x_tag>{get_face<dir_x_tag, side_low_tag>(p1, p2), Low, omega_list, chunk});
		flux_faces.push_back(new Flux_Face<dir_y_tag>{get_face<dir_y_tag, side_low_tag>(p1, p2), Low, omega_list, chunk});
		flux_faces.push_back(new Flux_Face<dir_z_tag>{get_face<dir_z_tag, side_low_tag>(p1, p2), Low, omega_list, chunk});
		flux_faces.push_back(new Flux_Face<dir_x_tag>{get_face<dir_x_tag, side_high_tag>(p1, p2), High, omega_list, chunk});
		flux_faces.push_back(new Flux_Face<dir_y_tag>{get_face<dir_y_tag, side_high_tag>(p1, p2), High, omega_list, chunk});
		flux_faces.push_back(new Flux_Face<dir_z_tag>{get_face<dir_z_tag, side_high_tag>(p1, p2), High, omega_list, chunk});
		
		for(auto& item : omega_list)
			item *= 2 * pi;
	}
	
	Flux_Box::~Flux_Box() {
		for(auto item : flux_faces)
			delete item;
	}
	
	void Flux_Box::update(const int time_step, const size_t rank, const size_t num_proc) {
		for(auto item : flux_faces)
			item->update(time_step, rank, num_proc);
	}
	
	real_arr Flux_Box::get() const {
		real_arr res(omega_list.size(), 0);
		
		for(auto item : flux_faces) {
			auto tmp = item->get();
			for(int i = 0; i < tmp.size(); ++i)
				res[i] += tmp[i];
		}
		
		return res;
	}
	
}
