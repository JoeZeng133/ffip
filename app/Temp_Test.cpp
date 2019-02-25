#include <functional>
#include <iostream>
#include <utility.hpp>
#include <geometry.hpp>
#include <thread>
#include <chrono>
#include <array>

using namespace std;
using namespace ffip;

void func() {
	vector<double> w(8);
	interpn<3> interp{ 2, 2, 2 };
	interp.put(w, 1.0, 0.1, 0.0, 0.0);
	interp.put(w, 1.0, 0.3, 0.9, 1.0);
	for (auto x : w)
		cout << x << " ";
}

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	func();
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
