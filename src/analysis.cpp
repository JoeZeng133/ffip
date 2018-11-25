#include <analysis.hpp>

namespace ffip {
	Probe_Time::Probe_Time(const fVec3& _pos, const int _time_step): pos(_pos), time_step(_time_step) {}
	
	void Probe_Time::update(const Simulation& sim) {
		if (sim.get_step() == time_step) {
			ex = sim.at(pos, Ex);
			ey = sim.at(pos, Ey);
			ez = sim.at(pos, Ez);

			hx = sim.at(pos, Hx) / 2;
			hy = sim.at(pos, Hy) / 2;
			hz = sim.at(pos, Hz) / 2;
		}

		if (sim.get_step() == time_step + 1) {
			hx += sim.at(pos, Hx) / 2;
			hy += sim.at(pos, Hy) / 2;
			hz += sim.at(pos, Hz) / 2;
		}
	}
	
	void Probe_Time::output(std::ostream & o) {
		o << std::scientific;
		o << ex << " " << ey << " " << ez << " " << hx << " " << hy <<  " " << hz << std::endl;
	}
	
	Probe_Frequency::Probe_Frequency(const fVec3& _pos, const real freq): pos(_pos), omega(2 * pi * freq) {}

	void Probe_Frequency::update(const Simulation& sim) {
		dt = sim.get_dt();
		real phase = -sim.get_step() * dt * omega;
		complex_num exp_omega_n = { cos(phase), sin(phase) };

		ex += sim.at(pos, Ex) * exp_omega_n;
		ey += sim.at(pos, Ey) * exp_omega_n;
		ez += sim.at(pos, Ez) * exp_omega_n;

		hx += sim.at(pos, Hx) * exp_omega_n;
		hy += sim.at(pos, Hy) * exp_omega_n;
		hz += sim.at(pos, Hz) * exp_omega_n;
	}

	void Probe_Frequency::output(std::ostream &o) {
		if (!phase_corrected) {
			phase_corrected = true;
			complex_num correction = complex_num{cos(0.5 * dt * omega), sin(0.5 * dt * omega)};		//exp(0.5dtw)
			
			hx *= correction;
			hy *= correction;
			hz *= correction;
		}
		
		o << std::scientific;
		o << ex << " " << ey << " " << ez << " "
		<< hx << " " << hy << " " << hz << std::endl;
	}
}
