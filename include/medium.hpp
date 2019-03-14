#pragma once

#include <utility.hpp>

namespace ffip {
	/* compact Poles structure for computation*/
	struct Susceptibility {
		double a0, b0, b1, b2;
	};

    struct Medium  {
        double epsilon, mu, e_cond, h_cond;
        std::valarray<double> e_sus_amp, h_sus_amp;

        Medium& operator+(const Medium& other);
        Medium& operator*(const double mult);
        Medium& operator+=(const Medium& other);
    };


}

