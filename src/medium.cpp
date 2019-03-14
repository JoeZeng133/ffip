
#include <medium.hpp>

namespace ffip {
    Medium& Medium::operator+(const Medium& other) {
        (*this) += other;
        return *this;
    }

    Medium& Medium::operator+=(const Medium& other) {
        epsilon += other.epsilon;
        mu += other.mu;
        e_cond += other.e_cond;
        h_cond += other.h_cond;
        e_sus_amp += other.e_sus_amp;
        h_sus_amp += other.h_sus_amp;
        return *this;
    }

    Medium& Medium::operator*(const double mult) {
        epsilon *= mult;
        mu *= mult;
        e_cond *= mult;
        h_cond *= mult;
        e_sus_amp *= mult;
        h_sus_amp *= mult;
    }    
}