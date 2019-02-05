#pragma once

#include <utility.hpp>

namespace ffip {
	/* compact Poles structure for computation*/
	struct Pole_Ref {
		real a1, a2, b0, b1, b2;
	};

	/* configuration of medium */
	using Medium_Voxel = typename std::valarray<real>;
	
	/* Poles for scripting*/
	struct Pole_Base {
		virtual real get_a1() const = 0;
		virtual real get_a2() const = 0;
		virtual real get_b0() const = 0;
		virtual real get_b1() const = 0;
		virtual real get_b2() const = 0;
		virtual ~Pole_Base() {}
		
		void set_dt(const real _dt);
		void set_index(const int _index);
		Pole_Ref get_ref() const;
		real dt;
		int index;
	};
	
	struct CP_Pole : public Pole_Base {
		real a0, a1, b0, b1, b2;
		
		CP_Pole() = default;
		CP_Pole(real A, real phi, real Omega, real Gamma);
		
		/* override functions*/
		real get_a1() const override;
		real get_a2() const override;
		real get_b0() const override;
		real get_b1() const override;
		real get_b2() const override;
		real get_cp() const;
	};
	
	struct Lorentz_Pole : public CP_Pole{
		
		Lorentz_Pole(real rel_perm, real freq, real damp);
		~Lorentz_Pole() = default;
	};
	
	struct Deybe_Pole : public CP_Pole{
		
		Deybe_Pole(real rel_perm, real relaxation);
		~Deybe_Pole() = default;
	};
	
	struct Drude_Pole : public CP_Pole{
		
		Drude_Pole(real freq, real inv_relaxation);
		~Drude_Pole() = default;
	};
	
	/* medium reference for use in actual computation
	   it holds a mixed bag of polarization information by pointers
	 it can be added and multiplied to construct new material reference
	 */
	class Medium_Ref {
	public:
		std::vector<Pole_Ref const*> poles;
		std::vector<real> weights;
		
		real d1{0}, d2{0}, d3{0};
	public:
		Medium_Ref() = default;
		//copiable
		Medium_Ref(const Medium_Ref&)= default;
		Medium_Ref& operator=(const Medium_Ref&) = default;
		//movable
		Medium_Ref(Medium_Ref&&) = default;
		Medium_Ref& operator=(Medium_Ref&&) = default;
		
		real update_field(real eh, real eh1, real norm_db, Dispersive_Field* f2) const;	//F(n) -> F(n+1), retunr F(n + 1), norm_jmd = jmd / e0
		Medium_Ref& operator+=(const Medium_Ref& other);
		Medium_Ref& operator*=(const real f);
		
		real get_d1() const;
		real get_d2() const;
		real get_d3() const;
		size_t get_size_poles() const;
	};
	
	Medium_Ref operator+(const Medium_Ref& A, const Medium_Ref& B);
	
	Medium_Ref operator*(const Medium_Ref& A, const real f);
	
	Medium_Ref operator*(const real f, const Medium_Ref& A);
	
	/* Store medium information for */
	class Medium {
	private:
		real e_inf{1}, sigma_e{0};
		real u_inf{1}, sigma_u{0};
		real dt;
		
		std::vector<Pole_Base*> e_poles;
		std::vector<Pole_Base*> m_poles;
		std::vector<Pole_Ref*> e_poles_ref;
		std::vector<Pole_Ref*> m_poles_ref;
	public:
		int index;					//idex of the medium
		~Medium();
		Medium() = default;
		Medium(const Medium&) = delete; //not copiable
		Medium& operator=(const Medium&) = delete;
		Medium(Medium&&) = default; //move
		Medium& operator=(Medium&&) = default;
		
		//constructors
		Medium(const real _e_inf, const real _sigma_e);
		Medium(const real _e_inf, const real _sigma_e, const real _u_inf, const real _sigma_u);
		
		//add e or m poles to the material
		void add_e_poles(Pole_Base* const pole);
		template<typename... Args>
		void add_e_poles(Pole_Base* const pole, Args... args);
		
		void add_m_poles(Pole_Base* const pole);
		template<typename... Args>
		void add_m_poles(Pole_Base* const pole, Args... args);
		
		//for use in simulation class
		Medium_Ref get_e_medium_ref() const;
		Medium_Ref get_m_medium_ref() const;
		void set_index(const int index);
		void set_dt(const real dt);

		real get_e_inf() const;
		real get_u_inf() const;
		real get_c() const;			//light speed
		real get_z() const;			//impedence

		void read_config(const Config& config);
		void init();
	};
	
	template<typename... Args>
	void Medium::add_e_poles(Pole_Base* const pole, Args... args) {
		add_e_poles(pole);
		add_e_poles(args...);
	}
	
	template<typename... Args>
	void Medium::add_m_poles(Pole_Base* const pole, Args... args) {
		add_m_poles(pole);
		add_m_poles(args...);
	}
}

