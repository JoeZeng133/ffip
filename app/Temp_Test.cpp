#include <functional>
#include <iostream>
#include <utility.hpp>
#include <geometry.hpp>
#include <thread>
#include <chrono>

using namespace std;
using namespace ffip;

int func(double x) {
	return 1 / x;
}

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	
	Primitive* ptr = new Disk({0, 0, 0}, 1, 1, Z);
	cout << ptr->is_in_closure({ 0.51, 0.501, 0.51 }) << "\n";
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
