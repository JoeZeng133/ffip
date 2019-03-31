#include <utility.hpp>
#include <utility.cpp>


int main() {
    std::valarray<double> a, b;

    auto c = a + b;
    std::cout << c.size() << "\n";
}