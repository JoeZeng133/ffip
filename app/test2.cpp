/* test plane wave projector */

#include <simulation.hpp>

using namespace std;
using namespace ffip;

int main() {
	double courant = 1;
	double dt = 2e-17 / 50;
	double dx = c * dt / courant;
	int n = 30;
	int Np = 40;
	
	auto func = Sinuosuidal_Func(c / (Np * dx));
	
	auto plane_sim = Plane_Wave(dx, dt, n);
	
	plane_sim.set_excitation(func);
	plane_sim.set_PML_neg(PML(X, Low, 0));
	plane_sim.set_PML_pos(PML(X, High, 6, PML::optimal_sigma_max(3, dx)));
	plane_sim.init();
	

	for(int i = 0; i < 200; ++i) {
		plane_sim.advance();
//		for(int j = 0; j <= n; ++j)
//			cout << plane_sim.at({0, 0, dx * j}, Ex) << " ";
//		cout << endl;
	}
//	
	
//	Simulation tmp{1e-9, 1e-9, iVec3{20, 50, 50}};
}
