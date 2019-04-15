#include <simulation.hpp>
#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace ffip;

TEST_CASE("Json Conversions") {
    
    SECTION("Checking json2double_arr") {
        std::vector<double> arr = {0.1, 0.2, 0.3};
        json j = arr;
        auto c_arr = json2double_arr(j);
        REQUIRE(c_arr.size() == arr.size());
        for(int i = 0; i < c_arr.size(); ++i) {
            REQUIRE(c_arr[i] == arr[i]);
        }

        CHECK_NOTHROW(json2double_arr(json{}));
    }

    SECTION("Checking json2ctype") {
        json j = "Hx";
        CHECK(json2ctype(j) == Hx);
        j = "Hy";
        CHECK(json2ctype(j) == Hy);
        j = "Hz";
        CHECK(json2ctype(j) == Hz);
        j = "Dx";
        CHECK(json2ctype(j) == Dx);
        j = "Dy";
        CHECK(json2ctype(j) == Dy);
        j = "Dz";
        CHECK(json2ctype(j) == Dz);

        CHECK_THROWS(json2ctype(json{"adf"}));
    }

    SECTION("Checking medium") {
        json m = R"(
            {
                "medium": {
                    "epsilon": 1.0,
                    "mu": 1.0,
                    "electric_conductivity": 1.0,
                    "magnetic_conductivity": 1.0,
                    "electric_Susceptibility": [
                        {
                            "type": "Lorentz",
                            "frequency": 1.3,
                            "gamma": 1.2,
                            "amplitude": 1.52e9
                        },
                        {
                            "type": "Drude",
                            "frequency": 2.1e2,
                            "gamma": 1.5e-3,
                            "amplitude": 2.0e-5
                        }
                    ],
                    "magnetic_susceptibility": []
			    }
            }
        )"_json;

        auto medium = json2medium(m["medium"]);
        REQUIRE(medium.epsilon == 1.0);
        REQUIRE(medium.mu == 1.0);

        CHECK(medium.e_sus.size() == 3);
        CHECK(medium.e_sus_amp.size() == 3);

        CHECK(medium.m_sus.size() == 1);
        CHECK(medium.m_sus_amp.size() == 1);

        CHECK(medium.e_sus_amp[0] == 1.0);
        CHECK(medium.e_sus_amp[1] == 1.52e9);
        CHECK(medium.e_sus_amp[2] == 2.0e-5);

        CHECK(medium.e_sus[0] == make_conductivity_susceptibility());
        CHECK(medium.e_sus[1] == make_Lorentz_susceptibility(1.3, 1.2));
        CHECK(medium.e_sus[2] == make_Drude_susceptibility(2.1e2, 1.5e-3));
    }

    SECTION("checking json2ivec3") {
        auto j = json::parse(R"(
            [
                1.2, 1.3, 0.5e-4
            ]
        )");

        auto t1 = json2fvec3(j);
        auto t2 = json2ivec3(j);
        CHECK(t1 == fVec3(1.2, 1.3, 0.5e-4));
        CHECK(t2 == iVec3(1, 1, 0));
    }
}

