#pragma once

#include <utility.hpp>
#include <medium.hpp>
#include <geometry.hpp>

namespace ffip {
    struct Medium_Stepping {
        struct Abstract_Susceptibility {
            double c1, c2, c3;
            Abstract_Susceptibility(const Susceptibility& sus);
        };

        //abstract susceptibility pool, 
        std::vector<Abstract_Susceptibility> sus_pool;
        //susceptibility intensity, num_points x num_sus
        std::vector<double> sus_sigma;
        //polarization currents, num_points x num_sus
        std::vector<double> polarization;
        //polarization currents previous step, num_points x num_ss
        std::vector<double> polarization1;
        
        struct Stepping_Point {
            size_t index;
            bool is_dispersive;
        };

        Medium_Stepping() = default;
        //add a point with Medium
        void add_point(size_t index, const Medium& medium);
        //add susceptibility to pool
        void add_susceptibility(const Susceptibility& sus);
        //re-order to improve performance
        void organize();
        //update
        void step(const std::vector<double>& db, std::vector<double>& eh);
    };

    class Structure {
        //avering method used for material averaging
        enum Average_Method {No_Average, Line_Average, Volume_Average};

        public:
            //add electric susceptibility to the pool
            void add_e_susceptibility(const Susceptibility& sus);

            //add magnetic susceptibility to the pool
            void add_h_susceptibility(const Susceptibility& sus);

            //set material properties from a list of geometry objects
            void set_materials_from_geometry
                (const std::vector<Geometry>& geom_list, const Medium& default_medium, Average_Method method);
            
            //step electric fields
            void step_e(const std::vector<double>& d, std::vector<double>& e, std::vector<double>& e1);

            //step magnetic fields
            void setp_h(const std::vector<double>& b, std::vector<double>& h, std::vector<double>& h1);

            //getters for material properties
            std::complex<double> get_epsilon_x(const fVec3& pos, double frequency);
            std::complex<double> get_epsilon_y(const fVec3& pos, double frequency);
            std::complex<double> get_epsilon_z(const fVec3& pos, double frequency);
            std::complex<double> get_mu_x(const fVec3& pos, double frequency);
            std::complex<double> get_mu_y(const fVec3& pos, double frequency);
            std::complex<double> get_mu_z(const fVec3& pos, double frequency);

        private:
            Grid3 grid;
            std::vector<Susceptibility> e_sus_pool, h_sus_pool;
            std::vector<Medium> medium;
            Medium_Stepping e_medium, h_medium;
    };
}