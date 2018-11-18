#include <medium.hpp>

namespace ffip {
	real Medium_Internal::get_omega1() const{
		return 2/ d3;
	}
	
	real Medium_Internal::get_omega2() const{
		if (isinf(d3))
			return 0;
		else
			return 2 * d1 / d3;
	}
	
	real Medium_Internal::get_omega3() const{
		if (isinf(d3))
			return 0;
		else
			return 2 * d2 / d3;
	}
	
	void Medium_Internal::set_d(const real omega1, const real omega2, const real omega3) {
		d1 = omega2 / omega1;
		d2 = omega3 / omega1;
		d3 = 2 / omega1;
	}
	
	real Medium_Internal::update_field(real eh, real eh1, real jmd, Dispersive_Field *f2) const{
		
		if (poles.empty()) {		/* non dispersive material updates*/
			return d1 * eh + d3 * jmd;
		} else {					/* update current contributions*/
			
			real sum_J = 0;
			for(int i = 0; i < poles.size(); ++i) {
				sum_J += f2->jp[i] * poles[i].c1 + f2->jp1[i] * poles[i].c2;
			}
			
			/* buffer E(n + 1) */
			real ehp1 = d1 * eh + d3 * (jmd - sum_J);
			if(d2 > 0) ehp1 += d2 * eh1;					//lorentz material has non-zero d2
			
			/* update currents */
			for(int i = 0; i < poles.size(); ++i) {
				
				if(poles[i].a2 == 0) {			//Deybe or Drude update
					f2->jp[i] = f2->jp[i] * poles[i].a1 + ehp1 * poles[i].b0 + eh * poles[i].b1;
					
				} else {						//Lorentz update
					real tmp = f2->jp[i];		//buffer J(n)
					f2->jp[i] = f2->jp[i] * poles[i].a1 + f2->jp1[i] * poles[i].a2 + ehp1 * poles[i].b0 + eh1 * poles[i].b2;
					f2->jp1[i] = tmp;
				}
			}
			/* update currents*/
			return ehp1;
		}
	}
	
	Medium_Internal& Medium_Internal::operator*=(const real f) {
		d3 /= f;
		for(auto item : poles) {
			item.b0 *= f;
			item.b1 *= f;
			item.b2 *= f;
		}
		
		return *this;
	}
	
	Medium_Internal& Medium_Internal::operator+=(const Medium_Internal &other) {
		poles.insert(poles.end(), other.poles.begin(), other.poles.end());
		set_d(get_omega1() + other.get_omega1(), get_omega2() + other.get_omega2(), get_omega3() + other.get_omega3());
		
		return *this;
	}
	
	Medium_Internal Medium_Internal::operator+(const Medium_Internal& other) const {
		Medium_Internal res{other};
		return res+=other;
	}
	
	Medium_Internal Medium_Internal::operator*(const real f) const{
		Medium_Internal res{*this};
		return res*=f;
	}
	
	/* Pole_base and derived classes*/
	void Pole_Base::set_dt(real _dt) {dt = _dt;}
	void Pole_Base::set_e0(real _e0) {e0 = _e0;}
	int Lorents_Pole::id = 0x100;
	int Deybe_Pole::id = 0x200;
	int Drude_Pole::id = 0x400;
	/* Lorentz Pole */
	
	
	Lorents_Pole::Lorents_Pole(real _epsilon, real _omega, real _delta): epsilon(_epsilon), omega(_omega), delta(_delta) {}
	
	real Lorents_Pole::get_alpha() {
		return (2 - omega * omega * dt * dt) / (1 + delta * dt);
	}
	
	real Lorents_Pole::get_xi() {
		return (delta * dt - 1) / (delta * dt + 1);
	}
	
	real Lorents_Pole::get_gamma() {
		return e0 * epsilon * omega * omega * dt / 2 / (1 + delta * dt);
	}
	
	real Lorents_Pole::get_sigma1() {
		return get_gamma();
	}
	
	real Lorents_Pole::get_sigma2() {
		return 0;
	}
	
	real Lorents_Pole::get_sigma3() {
		return get_gamma();
	}
	
	Pole_Internal Lorents_Pole::get_pole_internal() {
		return Pole_Internal{get_alpha(), get_xi(), get_gamma(), 0, -get_gamma(), 0.5 * (1 + get_alpha()), get_xi()};
	}
	
	/* Deybe pole */
	Deybe_Pole::Deybe_Pole(real _epsilon, real _tau): epsilon(_epsilon), tau(_tau) {}
	
	real Deybe_Pole::get_kappa() {
		return (1 - 0.5 * dt / tau) / (1 + 0.5 * dt / tau);
	}
	
	real Deybe_Pole::get_beta() {
		return e0 * epsilon / tau / (1 + 0.5 * dt / tau);
	}
	
	real Deybe_Pole::get_sigma1() {
		return get_beta();
	}
	
	real Deybe_Pole::get_sigma2() {
		return get_beta();
	}
	
	real Deybe_Pole::get_sigma3() {
		return 0;
	}
	
	Pole_Internal Deybe_Pole::get_pole_internal() {
		return Pole_Internal{get_kappa(), 0, get_beta(), 0, get_beta(), 0.5 * (1 + get_kappa()), 0};
	}
	
	/* Drude pole*/
	Drude_Pole::Drude_Pole(real _omega, real _gamma): omega(_omega), gamma(_gamma) {}
	
	real Drude_Pole::get_kappa() {
		return (1 - gamma * dt / 2) / (1 + gamma * dt / 2);
	}
	
	real Drude_Pole::get_beta() {
		return omega * omega * e0 * dt / 2 / (1 + gamma * dt / 2);
	}
	
	real Drude_Pole::get_sigma1() {
		return get_beta();
	}
	
	real Drude_Pole::get_sigma2() {
		return -get_beta();
	}
	
	real Drude_Pole::get_sigma3() {
		return 0;
	}
	
	Pole_Internal Drude_Pole::get_pole_internal() {
		return Pole_Internal{get_kappa(), 0, get_beta(), get_beta(), 0, 0.5 * (1 + get_kappa()), 0};
	}
	
	/* Medium used for scripting*/
	Medium_Type::Medium_Type(real _epsilon_inf, real _sigma, real _e0): epsilon_inf(_epsilon_inf), sigma(_sigma), e0(_e0) {
		sigma /= abs(e0);
	}
	
	Medium_Type::~Medium_Type() {
		for(auto i : poles) {
			delete i;
		}
	}
	
	void Medium_Type::add_Lorentz_pole(real epsilon, real fp, real delta) {
		poles.push_back(new Lorents_Pole{epsilon, fp * 2 * pi, delta * 2 * pi});
		poles.back()->set_e0(e0);
	}
	
	void Medium_Type::add_Deybe_pole(real epsilon, real tau) {
		poles.push_back(new Deybe_Pole{epsilon, tau / (2 * pi)});
		poles.back()->set_e0(e0);
	}
	
	void Medium_Type::add_Drude_pole(real fp, real gamma) {
		poles.push_back(new Drude_Pole{fp * 2 * pi, gamma * 2 * pi});
		poles.back()->set_e0(e0);
	}
	
	void Medium_Type::set_dt(real _dt) {
		dt = _dt;
		for(auto pole : poles) {
			pole->set_dt(dt);
		}
	}
	
	void Medium_Type::set_id(const int _id) {
		id = _id;
	}
	
	Medium_Internal Medium_Type::get_medium_internal() const {
		Medium_Internal res;
		
		real omega1 = 2 * e0 * epsilon_inf / dt + sigma * e0;
		real omega2 = 2 * e0 * epsilon_inf / dt - sigma * e0;
		real omega3 = 0;
		
		for(auto pole : poles) {
			omega1 += pole->get_sigma1();
			omega2 += pole->get_sigma2();
			omega3 += pole->get_sigma3();
			
			res.poles.push_back(pole->get_pole_internal());
		}
		
		res.set_d(omega1, omega2, omega3);
		
		return res;
	}
	
	
	/* medium factory */
	size_t top = 0;
	std::vector<Medium_Type> electric_medium_type_holder;
	std::vector<Medium_Type> magnetic_medium_type_holder;
	std::vector<Medium_Internal> electric_medium_internal_holder;
	std::vector<Medium_Internal> magnetic_medium_internal_holder;
	
	int make_medium(const real epsilon_inf, const real sigma_epsilon, const real mu_inf, const real sigma_mu) {
		
		electric_medium_type_holder.push_back(Medium_Type(epsilon_inf, sigma_epsilon, e0));
		electric_medium_type_holder.back().set_id(top);
		
		magnetic_medium_type_holder.push_back(Medium_Type(mu_inf, sigma_mu, -u0));
		magnetic_medium_type_holder.back().set_id(top);
		
		return top++;
	}
	
	Medium_Internal* get_medium_internal(const std::vector<real>& weights, const bool is_electric_point) {
		constexpr real tol = 1e-4;
		int nonzeros = 0;
		real total = 0;
		Medium_Internal* res = nullptr;
		
		if(weights.size() != top)
			throw std::runtime_error("Mixing numbers have wrong length");
		
		/* rounding down parameters*/
		for(auto& x : weights) {
			if(x > tol) {
				nonzeros ++;
				total += x;
			}
		}
		
		auto& medium_internal_holder = is_electric_point? electric_medium_internal_holder : magnetic_medium_internal_holder;
		
		if(nonzeros > 1) {								//for the case of mixing more than 1 material
			medium_internal_holder.push_back(Medium_Internal{});			//generate a mixing material at the back
			for(int i = 0; i < weights.size(); ++i)
				if (weights[i] > tol){
					medium_internal_holder.back() += medium_internal_holder[i] * (weights[i] / total);
				}
			res = &medium_internal_holder.back();
		}
		else if (nonzeros == 1){						//for the case of only 1 material
			for(int i = 0; i < weights.size(); ++i) {
				if (weights[i] > tol) {
					res = &medium_internal_holder[i];
					break;
				}
			}
		}
		return res;
	}
	
	std::vector<real> get_zero_weights() {
		return std::vector<real>(top, 0);
	}
	
	void prepare_medium_internal(const real _dt) {
		for(auto& item : electric_medium_type_holder) {
			item.set_dt(_dt);
			electric_medium_internal_holder.push_back(item.get_medium_internal());
		}
		
		for(auto& item : magnetic_medium_type_holder) {
			item.set_dt(_dt);
			magnetic_medium_internal_holder.push_back(item.get_medium_internal());
		}
	}

}
