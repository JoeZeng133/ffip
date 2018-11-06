#include <analysis.hpp>

namespace ffip {
	void Probe_Time::update(const Simulation& sim) {
		if (sim.get_step() == time_step) {
			ex = sim(pos, Ex);
			ey = sim(pos, Ey);
			ez = sim(pos, Ez);

			hx = sim(pos, Hx) / 2;
			hy = sim(pos, Hy) / 2;
			hz = sim(pos, Hz) / 2;
		}

		if (sim.get_step() == time_step + 1) {
			hx += sim(pos, Hx) / 2;
			hy += sim(pos, Hy) / 2;
			hz += sim(pos, Hz) / 2;
		}
	}

	void Probe_Frequency::update(const Simulation& sim) {
		real phase = -sim.get_step() * sim.get_dt() * omega;
		complex_num exp_omega_n = { cos(phase), sin(phase) };

		ex += sim(pos, Ex) * exp_omega_n;
		ey += sim(pos, Ey) * exp_omega_n;
		ez += sim(pos, Ez) * exp_omega_n;

		hx += sim(pos, Hx) * exp_omega_n;
		hy += sim(pos, Hx) * exp_omega_n;
		hz += sim(pos, Hz) * exp_omega_n;
	}


}
