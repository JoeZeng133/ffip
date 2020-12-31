#pragma once

#include <utility.hpp>

namespace ffip
{

    //for generating k, b, c
    struct PML
    {
        double sigma_max;
        double k_max;
        //thickness of PML
        double d;
        //polynomial order
        int m;

        //sigma_max = 0 uses recommended sigma_max with Reflection error as 1e-4
        PML(double d, double sigma_max, double k_max = 1, int m = 3);

        //getters
        double get_b(double x, double dt) const;
        double get_c(double x, double dt) const;
        double get_sigma(double x) const;
        double get_k(double x) const;

    };

    //Dipole Sources (Gaussian Types)
    class Gaussian_Dipoles_Stepping
    {
    private:
        struct Stepping_Point
        {
            size_t index;
            double amp;
            double start_time;
            double end_time;
            double width;
        };

        double max_end_time{0};

        std::vector<Stepping_Point> gaussian1_points;
        std::vector<Stepping_Point> gaussian2_points;

    public:
        // Gaussian_Dipoles_Stepping() = default;

        //increase possibility of locality
        void organize();

        //time-0.5dt to time+0.5dt
        void step(double time, std::vector<double> &accdb);

        //add gaussian 1 dipole
        void add_gaussian1(size_t index, double amp, double start_time, double end_time, double width);

        //add gaussian 2 dipole
        void add_gaussian2(size_t index, double amp, double start_time, double end_time, double width);

        //output debug information
        void output_details(std::ostream &os, const Yee3& grid) const;
    };

    //PML Regions curls E, H calculation
    class PML_Stepping
    {
    private:
        struct Stepping_Point
        {
            size_t index;

            double phi1;
            double b1;
            double c1;
            double inv_k1;

            double phi2;
            double b2;
            double c2;
            double inv_k2;
        };

        double inv_dx;

        sVec3 strides;
        //points in x, y, z direction
        std::vector<Stepping_Point> points[3];

    public:
        // PML_Stepping() = default;

        //set dx for use in curl calculation
        void set_dx(double dx);

        //set strides for use in stepping
        void set_strides(sVec3 strides);

        //reorder points to improve spatial locality
        void organize();

        //stepping, updates accdb from eh
        void step(const std::vector<double> &eh, std::vector<double> &accdb);

        //add one PML points
        void add_stepping_point(Direction dir, size_t index, double b1, double c1, double k1, double b2, double c2, double k2);
		
		void output_details(std::ostream& os, const Yee3& grid) const;
    };

    //Non PML Regions
    class Curl_Stepping
    {
    private:
        sVec3 strides;
        std::vector<size_t> points[3];
        double inv_dx;

    public:
        // Curl_Stepping() = default;

        //set dx for calculating curl
        void set_dx(double dx);

        //strides for calculating index offset in 3 direction
        void set_strides(sVec3 strides);

        //reorder to improve locality
        void organize();

        //stepping, updates accdb from eh
        void step(const std::vector<double> &eh, std::vector<double> &accdb);

        //add a stepping point, direction has to be specified
        void add_stepping_point(Direction dir, size_t index);
    };

    enum Boundary_Condition
    {
        PEC,  //same as symmetry with phase -1
        PMC,  //same as symmetry with phase 1
        None, //ghost points forced to be zeros
        Sync,  //ghost points comming from another chunk
        Period, //forcing fields to be the same at two sides
    };

    //Store fields
    //Responsible for updating curl F + C, aware of signs
    class Fields
    {
    private:
        //grid information for the current process
        Yee3 grid;

        double dx, dt;

        PML_Stepping e_pml, m_pml;

        Curl_Stepping e_curl, m_curl;

        Gaussian_Dipoles_Stepping e_gaussian_dipoles, m_gaussian_dipoles;

        Boundary_Condition bc[3][2];

        double turnoff_time{0};

    public:
        //eh = [E, H], accdb = sum(-curl E - M) or sum(curl H - J))
        //d = accdb * dt;
        std::vector<double> eh, accdb;

        // Fields() = default;
		void set_boundary_conditions(Boundary_Condition _bc[3][2]);

        void set_grid(const Yee3 &grid, double dx, double dt);

        //add gaussian1 dipole, flip amplitude
        void add_dipole_source_gaussian1(const fVec3 &pos, Coord_Type ctype, double amp, double start_time, double frequency, double cutoff);

        //add gaussian2 dipole
        void add_dipole_source_gaussian2(const fVec3 &pos, Coord_Type ctype, double amp, double start_time, double frequency, double cutoff);

        //Initializations of curl
        //providing k, b, c from PML and the index extend [p1, p2] covering the whole simulation
        void init(
            const std::array<double_arr, 3> &k, const std::array<double_arr, 3> &b, const std::array<double_arr, 3> &c,
            const iVec3 &p1, const iVec3 &p2);

        //PML layer initialization helper functions
        template <unsigned int F3>
        void PML_init_helper(
            const std::array<double_arr, 3> &k, const std::array<double_arr, 3> &b, const std::array<double_arr, 3> &c,
            const iVec3 &p1, const iVec3 &p2);

        //getter helper functions
        double get_eh_helper(const fVec3 &pos, Coord_Type ctype) const;
        double get_db_helper(const fVec3 &pos, Coord_Type ctype) const;
        //return a e,h,d,b field components
        double get(const fVec3& pos, Coord_Type ctype) const;

        //raw getter functions, no range checking
        double get_eh_raw(const iVec3 &pos) const;
        double get_db_raw(const iVec3 &pos) const;

        //synchronize ghost points of boundaries in a particular direction
        void sync_boundary(MPI_Comm comm, Direction dir);

        //extract fields on a particular face with offset and store it in the buffer and return face size
        size_t extract_fields_face(Direction x3, Side side, int offset, std::vector<double> &buf);

        //symmetry boundary with phase = 1, -1
        size_t extract_symmetry_face(Direction x3, Side side, int phase, std::vector<double> &buf);

        //assign fields to a particular face with offset
        void assign_fields_face(Direction x3, Side side, int offset, const std::vector<double> &buf);
        
        //update accd = curlH - J
        void step_accd(MPI_Comm comm, double time);

        //update accb = -curlE - M
        void step_accb(MPI_Comm comm, double time);

        //return turned off time of source
        double get_source_turnoff_time() const;
		
		//report something
		void report() const;
		void output_details(std::ostream& os) const;
        void output_fields(std::ostream& os) const;
    };

    template <unsigned int F3>
    void Fields::PML_init_helper(const std::array<double_arr, 3> &k, const std::array<double_arr, 3> &b, const std::array<double_arr, 3> &c, const iVec3 &p1, const iVec3 &p2)
    {

        constexpr auto ctype = static_cast<Coord_Type>(F3);
        constexpr auto x3 = get_dir_from_ctype(ctype);
        constexpr auto x1 = static_cast<Direction>((x3 + 1) % 3);
        constexpr auto x2 = static_cast<Direction>((x3 + 2) % 3);

        PML_Stepping &pml = (is_e_point(ctype) ? e_pml : m_pml);
        Curl_Stepping &curl = (is_e_point(ctype) ? e_curl : m_curl);

        for (auto itr = Yee_Iterator(grid.get_grid_p1(), grid.get_grid_p2(), ctype); !itr.is_end(); itr.next())
        {
            iVec3 pos = itr.get_coord();

            size_t index = grid.get_index_from_coord(pos);
            int x1_index = pos.get<x1>() - p1.get<x1>();
            int x2_index = pos.get<x2>() - p1.get<x2>();

            //if it is a PML point
            if (c[x1][x1_index] != 0 || c[x2][x2_index] != 0)
                pml.add_stepping_point(x3, index,
                                    b[x1][x1_index], c[x1][x1_index], k[x1][x1_index],
                                    b[x2][x2_index], c[x2][x2_index], k[x2][x2_index]);
            //if it is not a PML point
            else
                curl.add_stepping_point(x3, index);
        }
    }
} // namespace ffip
