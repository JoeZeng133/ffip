#pragma once

#include <utility.hpp>
#include <simulation.hpp>
#include <iostream>

namespace ffip {
	class Probe;
	class Simulation;
	
	class Probe {
	public:
		virtual ~Probe() {}
		virtual void update(const Simulation& sim) = 0;
		virtual void output(std::ostream&) = 0;
	};

	class Probe_Time : public Probe {
	private:
		fVec3 pos;
		int time_step;
		
		real ex, ey, ez, hx, hy, hz;
		
	public:
		Probe_Time() = delete;
		Probe_Time(const fVec3& _pos, const int _time_step);
		void update(const Simulation& sim) override;
		void output(std::ostream&) override;
		
	};
	
	class Probe_Frequency : public Probe{
	private:
		fVec3 pos;
		real omega;
		real dt{0};

		complex_num ex{0}, ey{0}, ez{0}, hx{0}, hy{0}, hz{0};
		bool phase_corrected{0};
	public:
		Probe_Frequency() = delete;
		Probe_Frequency(const fVec3& _pos, const real freq);
		void update(const Simulation& sim) override;
		void output(std::ostream&) override;
	};
	
}
