#include <nlohmann/json.hpp>
#include <iostream>

using json = nlohmann::json;

int main() {
    auto str = R"(
        {
            "medium": {
				"epsilon": 1.0,
				"mu": 1.0,
				"Electric conductivity": 1.0,
				"Magnetic conductivity": 1.0,
				"Electric Susceptibility": [
					{
						"type": "Lorentz",
						"frequency": 1.0,
						"gamma": 1.0,
						"amplitude": 1.0
					},
					{
						"type": "Drude",
						"frequency": 1.0,
						"gamma": 1.0,
						"amplitude": 1.0
					}
				],
				"Magnetic Susceptibility": []
			}
        }
        )";
    
    
    // std::cout << str << "\n";
    json m = json::parse(str);
    std::cout << m["medium"];
}