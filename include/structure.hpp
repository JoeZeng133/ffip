#pragma once

#include <utility.hpp>
#include <medium.hpp>
#include <geometry.hpp>

namespace ffip
{
    struct Medium_Stepping
    {
        struct Stepping_Point
        {
            size_t index;
            double g_inf;
        };

        //time step
        double dt;
        //abstract susceptibility pool,
        std::vector<Abstract_Susceptibility> sus_pool;
        //susceptibility intensity, num_points x num_sus
        std::vector<double> sus_amp;
        //polarization currents, num_points x num_sus
        std::vector<double> polarization;
        //polarization currents previous step, num_points x num_ss
        std::vector<double> polarization1;
        //collection of stepping points
        std::vector<Stepping_Point> points;

        //set dt
        void set_dt(double dt);

        //set material pool
        void set_susceptibility_pool(const std::vector<Abstract_Susceptibility> &sus_pool);
        void set_susceptibility_pool(const std::vector<Susceptibility> &sus_pool);

        //add a point with Medium
        void add_point(size_t index, double g_inf, const std::valarray<double> &sus_amp);

        //re-order to improve performance
        void organize();

        //update
        void step(const std::vector<double> &accdb, std::vector<double> &eh);

        void output_details(std::ostream& os, const Yee3& grid) const;
    };

    class Structure
    {
    private:
        Yee3 grid;
        double dx, dt;

        //geometry list
        std::vector<Geometry> geom_list;

        //Susceptibility pool, both abstract and original
        std::vector<Susceptibility> e_sus_pool, m_sus_pool;
        std::vector<Abstract_Susceptibility> e_ab_sus_pool, m_ab_sus_pool;
        std::vector<Medium> medium;

        //group medium points by their susceptibility non-zeros
        std::unordered_map<size_t, Medium_Stepping> e_stepping, m_stepping;

    public:
        //avering method used for material averaging
        enum Average_Method
        {
            No_Average,
            Line_Average,
            Volume_Average
        };

        // Structure() = default;

        void set_grid(const Yee3 &grid, double dx, double dt);

        //build susceptibility pool for constructing abstract materials
        //need to include all possible materials
        //maximum 64 types of susceptibility
        void add_to_material_pool(const std::vector<Medium> &materials);

        //return abstract medium from medium using a pool of susceptibility
        //used to build geometry objects
        Abstract_Medium get_abstract_medium(const Medium &medium) const;

        //set material properties from a list of geometry objects
        void set_materials_from_geometry(const std::vector<std::reference_wrapper<Geometry>> &geom_list, const Medium &default_medium, Average_Method method);

        //get non-zeros pattern
        size_t get_non_zeros_from_array(const std::valarray<double> &arr, double tol = 1e-4) const;

        //mask susceptibility based on non-zeros
        std::vector<Abstract_Susceptibility> mask_susceptibility_pool(size_t mask, const std::vector<Abstract_Susceptibility> &ab_sus_pool) const;
		
		std::valarray<double> mask_susceptibility_amp(size_t mask, const std::valarray<double> &arr) const;

        //step electric fields
        void step_e(const std::vector<double> &accd, std::vector<double> &e);

        //step magnetic fields
        void step_m(const std::vector<double> &accb, std::vector<double> &h);

        void output_details(std::ostream& os) const;
    };
} // namespace ffip
