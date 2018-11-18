// test geometry
#include <geometry.hpp>
#include <iostream>

using namespace ffip;
using namespace std;


int main(int argc, char const *argv[]) {
	fstream fo{"data.in", ios::out};
	
	auto sph1 = make_sphere(fVec3{0, 0, 0}, 0.5);
	auto sph2 = make_sphere(fVec3{0, -0.6, 0}, 0.5);
	auto box1 = make_box(fVec3{0.1, 0.1, -0.5}, fVec3{0.8, 0.8, 0.5});
	auto box2 = make_box(fVec3{0.26, -0.9, -0.5}, fVec3{0.64, -0.49, 0.5});
	auto sth = sph1 - box1 + sph2 * box2;
	
	
	for(double x = -1; x <= 1; x += 0.01)
		for(double y = -1; y <= 1; y += 0.01)
			fo << sth.is_in_exterior({x, y, 0}) << " ";
	
	fo.close();
    return 0;
}
