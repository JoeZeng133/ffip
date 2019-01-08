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
	
	interpn<3> tmp{2, 2, 2};
	vector<double> arr(8);
	vector<double> arr2{2,3,4,5,12,3,1,2};
	double res = 0;
	
	tmp.put(arr, 1.0, 0.1, 0.3, 0.76);
	for(int i = 0; i < 8; ++i)
		res += arr2[i] * arr[i];
	
	cout << res << " " << tmp.get(arr2, 0.1, 0.3, 0.76) << "\n";
	
	
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	cout << "\nRunning Time: " << elapsed_seconds.count() << endl;
}
