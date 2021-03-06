#pragma once

#include <utility.hpp>

namespace ffip
{
    //Susceptibility in polynomial form
    struct Susceptibility
    {
        double a0, a1, b0, b1, b2;

        //return complex value
		Susceptibility(double a0, double a1, double b0, double b1, double b2);
        std::complex<double> get_val(double frequency) const;
    };

    //Susceptibility comparisons to eliminate duplicates
    bool operator==(const Susceptibility &x, const Susceptibility &y);

    //Abstract Susceptibility for computation
    struct Abstract_Susceptibility
    {
        double c1, c2, c3;
        Abstract_Susceptibility(const Susceptibility &sus, double dt);
    };

    //Medium
    struct Medium
    {
        double epsilon{0}, mu{0};
        std::vector<Susceptibility> e_sus, m_sus;
        std::vector<double> e_sus_amp, m_sus_amp;

        Medium();
        Medium(double epsilon, double mu);
        void add_e_susceptibility(const Susceptibility &sus, double amp);
        void add_m_susceptibility(const Susceptibility &sus, double amp);
		double get_imp() const;
    };

    //Abstract Medium
    //stores only amplitudes of each susceptibility
    struct Abstract_Medium
    {
        double epsilon{0}, mu{0};
        std::valarray<double> e_sus_amp, m_sus_amp;

        // Abstract_Medium& operator=(const Abstract_Medium& other) = default;
        // Abstract_Medium& operator=(Abstract_Medium&& other) = default;

        Abstract_Medium &operator*=(double mult);
        Abstract_Medium &operator+=(const Abstract_Medium &other);
    };

    //return a measure of norm on the medium
    //used for medium averaging to measure error
    double norm(const Abstract_Medium &medium);

    //linear space operations
    Abstract_Medium operator+(const Abstract_Medium &a, const Abstract_Medium &b);
    Abstract_Medium operator*(const Abstract_Medium &a, double b);
    Abstract_Medium operator*(double b, const Abstract_Medium &a);

    //Susceptibility factories
    Susceptibility make_Lorentz_susceptibility(double frequency, double gamma);

    Susceptibility make_Drude_susceptibility(double frequeency, double gamma);

    Susceptibility make_Deybe_susceptibility(double tau);

    Susceptibility make_conductivity_susceptibility();
	

} // namespace ffip
