#pragma once

#include <utility.hpp>

namespace ffip {
	// forward declarations
	class Medium_Type;
	
	/* medium and pole for internal use*/
	struct Pole_Internal {
		real a1, a2, b0, b1, b2, c1, c2;			//J(n+1) = a1 * J(n) + a2 * J(n - 1) + b0 * F(n + 1) + b1 * F(n) + b2 * F(n - 1), sum(J) = c1 * J(n) + c2 * J(n - 1)
	};
	
	class Medium_Internal {
		friend Medium_Type;
	public:
		std::vector<Pole_Internal> poles;
		/* F(n+1) = d1 * F(n) + d2 * F(n - 1) + d3 * (G(n + 0.5) - sum(J))
		 inf means there is no material, sigma = epsilon = 0;
		 */
		real d1{INFINITY}, d2{INFINITY}, d3{INFINITY};
		
	public:
		Medium_Internal() = default;
		Medium_Internal(const Medium_Internal&)= default;				//copy
		Medium_Internal& operator=(const Medium_Internal&) = default;
		Medium_Internal(Medium_Internal&&) = default;					//move
		Medium_Internal& operator=(Medium_Internal&&) = default;
		
		real update_field(real eh, real eh1, real jmd, Dispersive_Field* f2) const;	//F(n) -> F(n+1), retunr F(n + 1)
		Medium_Internal operator+(const Medium_Internal& other) const;
		Medium_Internal operator*(const real f) const;
		
		Medium_Internal& operator+=(const Medium_Internal& other);
		Medium_Internal& operator*=(const real f);
		
		real get_omega1() const;
		real get_omega2() const;
		real get_omega3() const;
		void set_d(const real omega1, const real omega2, const real omega3);	//from omegas to d1,d2,d3
	};
	
	
	
	/* Poles used to derive Pole_Internal*/
	struct Pole_Base {
		virtual real get_sigma1() = 0;			//sigmas used for calculating omegas
		virtual real get_sigma2() = 0;
		virtual real get_sigma3() = 0;
		virtual Pole_Internal get_pole_internal() = 0;
		virtual ~Pole_Base() {}
		
		void set_dt(real _dt);
		real dt;
	};
	
	struct Lorents_Pole : public Pole_Base{
		static int id;
		
		real epsilon, omega, delta;
		
		/* derived members*/
		real get_alpha();
		real get_xi();
		real get_gamma();
		
		Lorents_Pole(real _epsilon, real _omega, real _delta);
		~Lorents_Pole() = default;
		
		/* override functions*/
		real get_sigma1() override;
		real get_sigma2() override;
		real get_sigma3() override;
		Pole_Internal get_pole_internal() override;
	};
	
	struct Deybe_Pole : public Pole_Base{
		static int id;
		
		real epsilon, tau;
		/* derived members */
		real get_kappa();
		real get_beta();
		
		Deybe_Pole(real _epsilon, real _tau);
		~Deybe_Pole() = default;
		
		/* override functions*/
		real get_sigma1() override;
		real get_sigma2() override;
		real get_sigma3() override;
		Pole_Internal get_pole_internal() override;
	};
	
	struct Drude_Pole : public Pole_Base{
		static int id;
		
		real omega, gamma;
		
		/* */
		real get_kappa();
		real get_beta();
		
		/* override functions*/
		Drude_Pole(real _omega, real _gamma);
		~Drude_Pole() = default;
		real get_sigma1() override;
		real get_sigma2() override;
		real get_sigma3() override;
		Pole_Internal get_pole_internal() override;
	};
	
	class Medium_Type {
	private:
		std::vector<Pole_Base*> poles;
		real dt;
		real epsilon_inf{1}, sigma{0};
	
	public:
		int id;					//id of the medium
		
		Medium_Type(real _epsilon_inf, real _simga);
		Medium_Type() = default;
		Medium_Type(const Medium_Type&) = default; //copy
		Medium_Type& operator=(const Medium_Type&) = default;
		Medium_Type(Medium_Type&&) = default; //move
		Medium_Type& operator=(Medium_Type&&) = default;
		~Medium_Type();
		
		//add poles based on given parameters, all parameters are specified in Hz instead of rad/s in accordance with Meep. But internal calculation are all based on rad/s, so there is unit conversion
		void add_Lorentz_pole(real epsilon, real fp, real delta);
		void add_Deybe_pole(real epsilon, real tau);
		void add_Drude_pole(real fp, real gamma);
		
		void set_dt(real _dt);
		void set_id(int _id);
		Medium_Internal get_medium_internal() const;
	};
	
	/* medium factory
	   manage all the medium used in simulation
	   save space for points that share the same material
	 */
	std::vector<Medium_Type> medium_type_holder;
	std::vector<Medium_Internal> medium_internal_holder;
	
	Medium_Type* make_medium(const real epsilon_inf, const real sigma = 0);
	Medium_Internal* get_medium_internal(const std::vector<real>& weights);
	std::vector<real> get_zero_weights();
}

