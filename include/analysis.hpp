#pragma once

#include <utility.hpp>
#include <simulation.hpp>

namespace ffip {
	class Probe {
		virtual ~Probe() {}
		virtual void update(const Simulation& sim) = 0;
	};

	class Probe_Time : public Probe {
	private:
		fVec3 pos;
		int time_step;
		
		real ex, ey, ez, hx, hy, hz;
		
	public:
		Probe_Time() = delete;

		void update(const Simulation& sim) override;
	};
	
	class Probe_Frequency : public Probe{
	private:
		fVec3 pos;
		real omega;

		complex_num ex, ey, ez, hx, hy, hz;
	public:
		Probe_Frequency() = delete;

		void update(const Simulation& sim) override;
	};
	
}
