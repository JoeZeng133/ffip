#include <functional>
#include <iostream>
#include <utility.hpp>
#include <thread>
#include <chrono>

using namespace std;

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	ffip::iVec3 p1{ 0, 0, 0 };
	ffip::iVec3 p2{ 4, 4, 4 };

	int n = 3;
	for(int i = 0; i < n; ++i) {
		auto itr = ffip::my_iterator(p1, p2, ffip::Null, i, n);
		for(;!itr.is_end(); itr.advance()) {
			cout << itr.x << " " << itr.y << " " << itr.z << "\n";
		}
	}

	cout << thread::hardware_concurrency() << endl;
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
