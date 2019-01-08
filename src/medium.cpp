#include <medium.hpp>

namespace ffip {
	real Medium_Ref::get_d1() const {
		return d1;
	}
	
	real Medium_Ref::get_d2() const {
		return d2;
	}
	
	real Medium_Ref::get_d3() const {
		return d3;
	}
	
	size_t Medium_Ref::get_size_poles() const {
		return poles.size();
	}
	
	real Medium_Ref::update_field(real eh, real eh1, real x, Dispersive_Field *f2) const{
		
		if (poles.empty()) {		/* non dispersive material updates*/
			return (x - eh * d2) / d1;
			
		} else {					/* update current contributions*/
			auto& pn = f2->p;
			auto& pn1 = f2->p1;
			
			real sum_J = 0;
			for(int i = 0; i < pn.size(); ++i) {
				auto pole = poles[i];
				sum_J += (pole->a1 - 1) * pn[i] + pole->a2 * pn1[i];
			}
			
			real res = (x - eh * d2 - eh1 * d3 - sum_J) / d1;
			
			for(int i = 0; i < pn.size(); ++i) {
				auto pole = poles[i];
				real tmp = pn[i];
				pn[i] = pole->a1 * pn[i] + pole->a2 * pn1[i] + (pole->b0 * res + pole->b1 * eh + pole->b2 * eh1) * weights[i];
				pn1[i] = tmp;
			}
			
			return res;
		}
	}
	
	Medium_Ref& Medium_Ref::operator*=(const real f) {
		for(auto& item : weights)
			item *= f;
		
		d1 *= f;
		d2 *= f;
		d3 *= f;
		return *this;
	}
	
	Medium_Ref& Medium_Ref::operator+=(const Medium_Ref &other) {
		poles.insert(poles.end(), other.poles.begin(), other.poles.end());
		weights.insert(weights.end(), other.weights.begin(), other.weights.end());
		
		d1 += other.d1;
		d2 += other.d2;
		d3 += other.d3;
		return *this;
	}
	
	Medium_Ref operator+(const Medium_Ref& A, const Medium_Ref& B) {
		Medium_Ref res{A};
		return res+=B;
	}
	
	Medium_Ref operator*(const Medium_Ref& A, const real f) {
		Medium_Ref res{A};
		return res*=f;
	}
	
	Medium_Ref operator*(const real f, const Medium_Ref& A) {
		Medium_Ref res{A};
		return res*=f;
	}
	
	/* Pole_base and derived classes*/
	void Pole_Base::set_dt(const real _dt) {dt = _dt;}
	void Pole_Base::set_index(const int _index) {index = _index;}
	Pole_Ref Pole_Base::get_ref() const {
		return Pole_Ref{get_a1(), get_a2(), get_b0(), get_b1(), get_b2()};
	}
	
	/* Critical Points Pole*/
	CP_Pole::CP_Pole(real A, real phi, real Omega, real Gamma) {
		Omega *= 2 * pi;
		Gamma *= 2 * pi;
		
		a0 = 2 * A * Omega * (Omega * cos(phi) - Gamma * sin(phi));
		a1 = -2 * A * Omega * sin(phi);
		b0 = Omega * Omega + Gamma * Gamma;
		b1 = 2 * Gamma;
		b2 = 1;
	}
	
	real CP_Pole::get_cp() const {
		return b2 / (dt * dt) + b1 / (2 * dt) + b0 / 4;
	}
	
	real CP_Pole::get_a1() const {
		return (2 * b2 / (dt * dt) - b0 / 2) / get_cp();
	}
	
	real CP_Pole::get_a2() const {
		return (b1 / (2 * dt) - b2 / (dt * dt) - b0 / 4) / get_cp();
	}
	
	real CP_Pole::get_b0() const {
		return (a0 / 4 + a1 / (2 * dt)) / get_cp();
	}
	
	real CP_Pole::get_b1() const {
		return (a0 / 2) / get_cp();
	}
	
	real CP_Pole::get_b2() const {
		return (a0 / 4 - a1 / (2 * dt)) / get_cp();
	}
	
	/* Lorentz Pole */
	Lorentz_Pole::Lorentz_Pole(real rel_perm, real freq, real damp) {
		a0 = rel_perm * (2 * pi * freq) * (2 * pi * freq);
		a1 = 0;
		b0 = (2 * pi * freq) * (2 * pi * freq);
		b1 = 2 * pi * damp;
		b2 = 1;
	}
	
	/* Deybe pole */
	Deybe_Pole::Deybe_Pole(real rel_perm, real relaxation) {
		a0 = rel_perm;
		a1 = 0;
		b0 = 1;
		b1 = relaxation;
		b2 = 0;
	}
	
	/* Drude pole*/
	Drude_Pole::Drude_Pole(real freq, real inv_relaxation) {
		a0 = (2 * pi * freq) * (2 * pi * freq);
		a1 = 0;
		b0 = 0;
		b1 = inv_relaxation * 2 * pi;
		b2 = 1;
	}
	
	/* Medium used for scripting*/
	Medium::~Medium() {
		for(auto item : e_poles)
			delete item;
		
		for(auto item : m_poles)
			delete item;
		
		for(auto item : e_poles_ref)
			delete item;
		
		for(auto item : m_poles_ref)
			delete item;
	}
	
	Medium::Medium(const real _e_inf, const real _sigma_e): e_inf{_e_inf}, sigma_e{_sigma_e} {}
	
	Medium::Medium(const real _e_inf, const real _sigma_e, const real _u_inf, const real _sigma_u): Medium{_e_inf, _sigma_e} {
		u_inf = _u_inf;
		sigma_u = _sigma_u;
	}
	
	void Medium::add_e_poles(Pole_Base* const  pole) {
		e_poles.push_back(pole);
	}
	
	void Medium::add_m_poles(Pole_Base* const pole) {
		m_poles.push_back(pole);
	}
	
	void Medium::set_index(const int _index) {
		index = _index;
	}
	
	void Medium::set_dt(const real _dt) {
		dt = _dt;
		for(auto item : e_poles) {
			item->set_dt(dt);
			e_poles_ref.push_back(new Pole_Ref{item->get_ref()});
		}
		
		for(auto item : m_poles) {
			item->set_dt(dt);
			m_poles_ref.push_back(new Pole_Ref{item->get_ref()});
		}
	}
	
	real Medium::get_e_inf() const {
		return e_inf;
	}
	
	real Medium::get_u_inf() const {
		return u_inf;
	}
	
	real Medium::get_c() const {
		return c0 / sqrt(e_inf * u_inf);
	}
	
	real Medium::get_z() const {
		return z0 * sqrt(u_inf / e_inf);
	}
	
	Medium_Ref Medium::get_e_medium_ref() const {
		Medium_Ref res;
		
		res.weights.resize(e_poles_ref.size(), 1);
		res.d1 = e_inf + sigma_e * dt / (2 * e0);
		res.d2 = -e_inf + sigma_e * dt / (2 * e0);
		res.d3 = 0;
		
		for(auto& pole : e_poles_ref) {
			res.d1 += pole->b0;
			res.d2 += pole->b1;
			res.d3 += pole->b2;
			res.poles.push_back(pole);
		}
		return res;
	}
	
	Medium_Ref Medium::get_m_medium_ref() const {
		Medium_Ref res;
		
		res.weights.resize(m_poles_ref.size(), 1);
		res.d1 = u_inf + sigma_u * dt / (2 * u0);
		res.d2 = -u_inf + sigma_u * dt / (2 * u0);
		res.d3 = 0;
		
		for(auto& pole : m_poles_ref) {
			res.d1 += pole->b0;
			res.d2 += pole->b1;
			res.d3 += pole->b2;
			res.poles.push_back(pole);
		}
		
		return res;
	}
}
