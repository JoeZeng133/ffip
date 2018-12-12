#include <functional>
#include <iostream>
#include <utility.hpp>
#include <thread>
#include <chrono>

using namespace std;

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	ffip::iVec3 p1{ -11, 3, 20 };
	ffip::iVec3 p2{ 22, 33, 40 };

	int n = 10;
	for (int i = 0; i < n; ++i) {
		auto tmp = ffip::divide_region(p1, p2, i, n);
		cout << tmp.first << " " << tmp.second << endl;
	}

	cout << thread::hardware_concurrency() << endl;
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
