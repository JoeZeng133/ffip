//#include <simulation.hpp>
#include <utility.hpp>
//#include <chunk.hpp>
#include <geometry.hpp>
#include <source.hpp>
//#include <medium.hpp>
//#include <analysis.hpp>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;
using namespace ffip;

using D = dir_x_tag;

int main() {
	iVec3 tmp {1, 2, 3};
	
	tmp = rotate_frame(tmp, D{});
	cout << tmp.x << " " << tmp.y << " " << tmp.z << endl;
	
	tmp = rotate_frame(tmp, dir_traits<D>::z{});
	
	
	cout << tmp.x << " " << tmp.y << " " << tmp.z << endl;
}
