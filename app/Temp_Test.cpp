#include <functional>
#include <iostream>
#include <utility.hpp>
#include <thread>

using namespace std;

int  main(int argc, char const *argv[]) {
	ffip::iVec3 p1{ -11, 3, 20 };
	ffip::iVec3 p2{ 22, 33, 40 };

	int n = 10;
	for (int i = 0; i < n; ++i) {
		auto tmp = ffip::divide_region(p1, p2, i, n);
		cout << tmp.first << " " << tmp.second << endl;
	}

	cout << thread::hardware_concurrency() << endl;
}
