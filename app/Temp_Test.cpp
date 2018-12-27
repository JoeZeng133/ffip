#include <functional>
#include <iostream>
#include <utility.hpp>
#include <thread>
#include <chrono>

using namespace std;

int func(double x) {
	return 1 / x;
}

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();

	int num_proc = thread::hardware_concurrency();
	int time_step = 500;

	for (int time = 0; time < 500; ++time) {
		vector<thread> threads;
		for (int i = 1; i < num_proc; ++i) {
			threads.push_back(thread(func, i));
		}
		for (auto& item : threads)
			item.join();

		cout << "Current Threads" << time << "\r";
	}
	

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
