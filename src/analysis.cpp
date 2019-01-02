#include <analysis.hpp>

namespace ffip {
	/* ################################ */
	Probe_Frequency::Probe_Frequency(const real freq, const fVec3& _pos, Chunk const* _chunk): omega(2 * pi * freq), pos(_pos), chunk(_chunk) {}

	void Probe_Frequency::update(const int time_step) {
		real phase = -time_step * chunk->get_dt() * omega;
		complex_num exp_omega_n = { cos(phase), sin(phase) };

		ex += chunk->at(pos, Ex) * exp_omega_n;
		ey += chunk->at(pos, Ey) * exp_omega_n;
		ez += chunk->at(pos, Ez) * exp_omega_n;

		hx += chunk->at(pos, Hx) * exp_omega_n;
		hy += chunk->at(pos, Hy) * exp_omega_n;
		hz += chunk->at(pos, Hz) * exp_omega_n;
	}

	void Probe_Frequency::output(std::ostream &o) {
		if (!phase_corrected.exchange(true)) {
			real phase = 0.5 * chunk->get_dt() * omega;
			complex_num correction = complex_num{cos(phase), sin(phase)};		//exp(0.5dtw)
			
			hx *= correction;
			hy *= correction;
			hz *= correction;
		}
		
		o << std::scientific;
		o << ex << " " << ey << " " << ez << " "
		<< hx << " " << hy << " " << hz << std::endl;
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
	}
	
	void Box_Freq_Req::validity_check() const {
		if (F == Null)
			throw Invalide_Coord_Type{};
		
		if (!ElementWise_Less_Eq(chunk->get_p1(), p1)
			|| !ElementWise_Less_Eq(p1, chunk->get_p2())
			|| !ElementWise_Less_Eq(chunk->get_p1(), p2)
			|| !ElementWise_Less_Eq(p2, chunk->get_p2()))
			throw Out_of_the_Domain{};
	}

	void Box_Freq_Req::update(const int time_step, const size_t rank, const size_t num_proc) {
		complex_arr fourier_coef;
		for(auto omega : omega_list)
			fourier_coef.push_back(complex_num(cos(omega * time_step * dt), -sin(omega * time_step * dt)));

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
	
	void Box_Freq_Req::correct_phase() {
		if (phase_corrected.exchange(true)) return;
		
		if (is_H_point(F)) {
			for(int f = 0; f < omega_list.size(); ++f) {
				real phase = omega_list[f] * dt * 0.5;
				complex_num exp_phase{cos(phase), sin(phase)};
				for(auto& val : raw_fields[f])
					val *= exp_phase;
			}
		}
	}
	
	/* ################################ */
	N2F_Box::N2F_Box(const fVec3& _p1, const fVec3& _p2, const real_arr& _omega_list, Chunk* const _chunk, const real _c):p1(_p1), p2(_p2), omega_list(_omega_list), chunk(_chunk), c(_c) {
		
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
	
	void N2F_Box::update(const int time_step, const size_t rank, const size_t num_proc) {
		for(auto item : n2f_faces)
			item->update(time_step, rank, num_proc);
	}
	
	std::pair<cVec3, cVec3> N2F_Box::get_NL(const real theta, const real phi, const real omega) {
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
	
	void N2F_Box::prepare(const size_t rank, const size_t num_proc) {
		for(auto item : n2f_faces)
			item->prepare(rank, num_proc);
	}
}
