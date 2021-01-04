
#include <medium.hpp>

namespace ffip
{
    //Susceptibility
	Susceptibility::Susceptibility(double b0, double b1, double a0, double a1, double a2):
	b0(b0), b1(b1), a0(a0), a1(a1), a2(a2) {}
	
    std::complex<double> Susceptibility::get_val(double frequency) const
    {
        double omega = 2 * pi * frequency;
		auto jw = std::complex<double>{0, omega};
		
        return (b0 + b1 * jw) / (a0 + a1 * jw - a2 * omega * omega);
    }

    bool operator==(const Susceptibility &x, const Susceptibility &y)
    {
		return x.b0 == y.b0 && x.a0 == y.a0 && x.a1 == y.a1 && x.a2 == y.a2 && x.b1 == y.b1;
    }

    //Abstract Susceptibility
    Abstract_Susceptibility::Abstract_Susceptibility(const Susceptibility &sus, double dt)
    {
        double a0 = sus.a0 + 2 * sus.a1 / dt + 4 * sus.a2 / dt / dt;
        a1 = 2 * sus.a0 - 8 * sus.a2 / dt / dt;
        a2 = sus.a0 - 2 * sus.a1 / dt + 4 * sus.a2 / dt / dt;
        b0 = sus.b0 + 2 * sus.b1 / dt;
        b1 = 2 * sus.b0;
        b2 = sus.b0 - 2 * sus.b1 / dt;
        //normalize to a0=1
        a1 /= a0;
        a2 /= a0;
        b0 /= a0;
        b1 /= a0;
        b2 /= a0;
    }

    //Medium
    Medium::Medium() : Medium(0, 0) {}
    
    Medium::Medium(double epsilon, double mu) : epsilon(epsilon), mu(mu) {}
	
	double Medium::get_imp() const
	{
		return std::sqrt(mu / epsilon);
	}

    void Medium::add_e_susceptibility(const Susceptibility &sus, double amp)
    {
        e_sus.push_back(sus);
        e_sus_amp.push_back(amp);
    }

    void Medium::add_m_susceptibility(const Susceptibility &sus, double amp)
    {
        m_sus.push_back(sus);
        m_sus_amp.push_back(amp);
    }

    //Abstract_Medium
    Abstract_Medium &Abstract_Medium::operator+=(const Abstract_Medium &other)
    {
        epsilon += other.epsilon;
        mu += other.mu;
        e_sus_amp += other.e_sus_amp;
        m_sus_amp += other.m_sus_amp;
        return *this;
    }

    Abstract_Medium &Abstract_Medium::operator*=(double mult)
    {
        epsilon *= mult;
        mu *= mult;
        e_sus_amp *= mult;
        m_sus_amp *= mult;
        return *this;
    }

    Abstract_Medium operator+(const Abstract_Medium &a, const Abstract_Medium &b)
    {
        Abstract_Medium res{a};
        return res += b;
    }

    Abstract_Medium operator*(const Abstract_Medium &a, double b)
    {
        Abstract_Medium res{a};
        return res *= b;
    }

    Abstract_Medium operator*(double b, const Abstract_Medium &a)
    {
        Abstract_Medium res{a};
        return res *= b;
    }

    double norm(const Abstract_Medium &medium)
    {
        return std::abs(medium.epsilon) + std::abs(medium.mu) +
            std::abs(medium.e_sus_amp).sum() + std::abs(medium.m_sus_amp).sum();
    }

    Susceptibility make_Lorentz_susceptibility(double frequency, double gamma)
    {
        double omega = 2 * pi * frequency;
        return Susceptibility{omega * omega, 0, omega * omega, 2 * pi * gamma, 1};
    }

    Susceptibility make_Drude_susceptibility(double frequency, double gamma)
    {
        double omega = 2 * pi * frequency;
        return Susceptibility{omega * omega, 0, 0, 2 * pi * gamma, 1};
    }

    Susceptibility make_Deybe_susceptibility(double tau)
    {
        return Susceptibility{1, 0, 1, tau, 0};
    }

    Susceptibility make_conductivity_susceptibility()
    {
        return Susceptibility{1, 0, 0, 1, 0};
    }
} // namespace ffip
