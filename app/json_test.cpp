#include <nlohmann/json.hpp>
#include <iostream>
#include <string>

using json = nlohmann::json;

int main() {
    json j;
    j["sth"] = 3.1;

    auto val = j.at("sthl").get<double>();

    std::cout << val << "\n";
}