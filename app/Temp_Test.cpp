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

	int num_proc = 32;

	for (int time = 0; time < 1000; ++time) {
		vector<thread> threads;
		for (int i = 1; i < num_proc; ++i) {
			threads.push_back(thread(func, i));
		}
		for (auto& item : threads)
			item.join();

		cout << "Current Threads" << time << "\r";
	}
	

	cout << "Terminated" << endl;
}
