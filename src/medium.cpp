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
	
	real Medium_Ref::update_field(real eh, real eh1, real norm_db, Dispersive_Field *f2) const{
		
		if (poles.empty()) {		/* non dispersive material updates*/
			return (2 * norm_db - eh * d2) / d1;
			
		} else {					/* update current contributions*/
			
			real sum_J = 0;
			for(int i = 0; i < poles.size(); ++i) {
				auto pole = poles[i];
				
				if (pole->a2 == 0) {				//deybe, drue update
					f2->jp1[i] = f2->jp[i] * pole->a1;
					sum_J += f2->jp[i] + f2->jp1[i];
				}
				else {										//lorentz pole
					f2->jp1[i] =f2->jp[i] * pole->a1 + f2->jp1[i] * pole->a2;
					sum_J += f2->jp[i] + f2->jp1[i];
				}
			}
			
			/* buffer E(n + 1) */
			real ehp1 = (2 * norm_db - sum_J - eh1 * d3 - eh * d2) / d1;
			
			/* update currents */
			for(int i = 0; i < poles.size(); ++i) {
				auto pole = poles[i];
				
				if(pole->a2 == 0) {			//Deybe or Drude update
					f2->jp[i] = f2->jp1[i] + (ehp1 * pole->b0 + eh * pole->b1) * weights[i];
					
				} else {						//Lorentz update
					real tmp = f2->jp[i];		//buffer J(n)
					f2->jp[i] = f2->jp1[i] + (ehp1 * pole->b0 + eh1 * pole->b2) * weights[i];
					f2->jp1[i] = tmp;
				}
			}
			return ehp1;
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
	
	/* Lorentz Pole */
	Lorentz_Pole::Lorentz_Pole(real rel_perm, real freq, real damp): epsilon(rel_perm), omega(freq * 2 * pi), delta(damp * pi) {}
	
	real Lorentz_Pole::get_a1() const {
		return (2 - omega * omega * dt * dt) / (1 + delta * dt);
	}
	
	real Lorentz_Pole::get_a2() const {
		return (delta * dt - 1) / (delta * dt + 1);
	}
	
	real Lorentz_Pole::get_b0() const {
		return epsilon * omega * omega * dt / 2 / (1 + delta * dt);
	}
	
	real Lorentz_Pole::get_b1() const {
		return 0;
	}
	
	real Lorentz_Pole::get_b2() const {
		return -get_b0();
	}
	
	/* Deybe pole */
	Deybe_Pole::Deybe_Pole(real rel_perm, real relaxation): epsilon(rel_perm), tau(relaxation) {}
	
	real Deybe_Pole::get_a1() const {
		return (1 - 0.5 * dt / tau) / (1 + 0.5 * dt / tau);
	}
	
	real Deybe_Pole::get_a2() const {
		return 0;
	}
	
	real Deybe_Pole::get_b0() const {
		return epsilon / tau / (1 + 0.5 * dt / tau);
	}
	
	real Deybe_Pole::get_b1() const {
		return -get_b0();
	}
	
	real Deybe_Pole::get_b2() const {
		return 0;
	}
	
	/* Drude pole*/
	Drude_Pole::Drude_Pole(real freq, real inv_relaxation): omega(freq * 2 * pi), gamma(inv_relaxation * 2 * pi) {}
	
	real Drude_Pole::get_a1() const {
		return (1 - gamma * dt / 2) / (1 + gamma * dt / 2);
	}
	
	real Drude_Pole::get_a2() const {
		return 0;
	}
	
	real Drude_Pole::get_b0() const {
		return omega * omega * dt / 2 / (1 + gamma * dt / 2);
	}
	
	real Drude_Pole::get_b1() const {
		return get_b0();
	}
	
	real Drude_Pole::get_b2() const {
		return 0;
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
		res.d1 = 2 * e_inf / dt + sigma_e / e0;
		res.d2 = -2 * e_inf / dt + sigma_e / e0;
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
		res.d1 = 2 * u_inf / dt + sigma_u / u0;
		res.d2 = -2 * u_inf / dt + sigma_u / u0;
		res.d3 = 0;
		
		for(auto& pole : m_poles_ref) {
			res.d1 += pole->b0;
			res.d2 += pole->b1;
			res.d3 += pole->b2;
			res.poles.push_back(pole);
		}
		
		return res;
	}
	
	
	/* medium factory, needs to add garbage collecting for poles*/
	std::vector<std::unique_ptr<Medium>> medium;
	std::vector<std::unique_ptr<Medium_Ref>> e_medium_ref;
	std::vector<std::unique_ptr<Medium_Ref>> m_medium_ref;
	
	void prepare_medium(const real _dt) {
		for(auto& item : medium) {
			item->set_dt(_dt);
			e_medium_ref.push_back(std::make_unique<Medium_Ref>(item->get_e_medium_ref()));
			m_medium_ref.push_back(std::make_unique<Medium_Ref>(item->get_m_medium_ref()));
		}
	}
	
	std::vector<real> get_zero_weights() {
		return std::vector<real>(medium.size(), 0);
	}
	
	Medium_Ref const* get_medium_ref(const bool is_electric_point, const std::vector<real>& weights) {
		
		constexpr real tol = 1e-4;
		int nonzeros = 0;
		real total = 0;
		Medium_Ref* res = nullptr;
		auto& medium_ref = is_electric_point? e_medium_ref : m_medium_ref;
		
		if(weights.size() != medium.size())
			throw std::runtime_error("Mixing numbers have wrong length");
		
		/* rounding down parameters*/
		for(auto& x : weights) {
			if(x > tol) {
				nonzeros ++;
				total += x;
			}
		}
		
		if(nonzeros > 1) {								//for the case of mixing more than 1 material
			medium_ref.push_back(std::unique_ptr<Medium_Ref>{});			//generate a mixing material at the back
			for(int i = 0; i < weights.size(); ++i)
				if (weights[i] > tol){
					*medium_ref.back() += *medium_ref[i] * (weights[i] / total);
				}
			res = medium_ref.back().get();
		}
		else if (nonzeros == 1){						//for the case of only 1 material
			for(int i = 0; i < weights.size(); ++i) {
				if (weights[i] > tol) {
					res = medium_ref[i].get();
					break;
				}
			}
		}
		
		if (res == nullptr)
			throw std::runtime_error("Illegal weights");
		
		return res;
	}
}
