#include <functional>
#include <iostream>
#include <utility.hpp>
#include <geometry.hpp>
#include <thread>
#include <chrono>

using namespace std;
using namespace ffip;

constexpr int N = 10000000;
atomic<size_t> top;
int a[N];

struct A {
	int x;
	A(const int x): x(x) {}
};

void func() {
	size_t cur;
	
	while ((cur = top++) < N) {
		a[cur]++;
	}
}

int  main(int argc, char const *argv[]) {
	auto start = std::chrono::system_clock::now();
	
	A tmp(2);
	cout << tmp.x << endl;
	std::vector<std::future<void>> task_list;
	for(int i = 0; i < 3; ++i) {
		task_list.push_back(std::async(std::launch::async, func));
	}
	for(auto &item : task_list)
		item.get();
	
	bool no = false;
	for (int i = 0; i < N; ++i)
		if (a[i] > 1)
			no = true;
	
	if (no)
		cout << "No" << endl;
	else
		cout << "Ye" << endl;
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
