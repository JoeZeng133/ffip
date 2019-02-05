#include <functional>
#include <iostream>
#include <utility.hpp>
#include <geometry.hpp>
#include <thread>
#include <chrono>
#include <array>

using namespace std;
using namespace ffip;

template<typename... Args>
int sum(Args&... args) {
	return (... + args.val());
}


int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	cout << sizeof(std::array<char, 3>) << sizeof(std::array<char, 100>) << endl;

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
