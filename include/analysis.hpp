#pragma once

#include <utility.hpp>

namespace ffip {
	class Probe_Time {
	private:
		std::vector<fVec3> pos_list;
		std::vector<int> time_list;
		Coord_Type ctype;
		
	public:
		Probe_Time(const std::vector<fVec3>& _pos_list, const std::vector<int>& _time_list, Coord_Type _ctype);
	};
	
	class Probe_Frequency {
		
	};
	
}
