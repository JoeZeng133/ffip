#include <functional>
#include <iostream>
#include <utility.hpp>
#include <thread>
#include <chrono>

using namespace std;
using namespace ffip;

int func(double x) {
	return 1 / x;
}

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	cout << sizeof(Dispersive_Field) << endl;
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
