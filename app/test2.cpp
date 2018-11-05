//template<int N, typename T, typename... Args>
//class A: public A<N - 1, Args...> {
//protected:
//	using base = A<N - 1, Args...>;
//	int extend;
//
//private:
//	T arr;
//	bool active{0};
//
//
//public:
//	A(const T& _arr, const Args&... args): base(args...), arr(_arr)
//	{
//		if(arr.size() > 1) {
//			active = 1;
//			base::val.resize(base::val.size() * 2);
//			base::jump.resize(base::jump.size() * 2);
//		}
//	}
//
//	int get_coeff(const double x) {
//		if (!active) {
//			extend = base::extend;
//			return extend;
//		}
//
//		int prev = base::get_coeff();
//		if(prev == -1) //capture if point is outside
//			return -1;
//
//		auto itr = lower_bound(arr.begin(), arr.end(), x);
//		if(itr == arr.begin() || itr == arr.end()) //point is outside
//			return -1;
//
//		int idx = distance(itr, arr.begin());
//
//		//update val array
//		for(int i = prev; i < (prev << 1); ++i)
//			base::val[i] = base::val[i - prev] * (arr[idx] - x);
//
//		for(int i = 0; i < prev; ++i)
//			base::val[i] *= (x - arr[idx - 1]);
//
//		extend = prev << 1;
//		return extend;
//	}
////
////	int get_jump() {
////		if(!active)
////			return base::get_jump();
////
////		int res = arr.size() * base::get_jump();
////
////
////
////		return res;
////	}
//};
//
//template<>
//class A<1, vector<double>> {
//protected:
//	vector<double> val;
//	vector<int> jump;
//
//private:
//	vector<double> arr;
//public:
//	A(const vector<double>& _arr): arr(_arr) {}
//};

int main() {
	
}
