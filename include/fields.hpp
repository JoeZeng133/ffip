#pragma once

#include <utility.hpp>


namespace ffip {
    /* Dipole sources interface */
    class Gaussian_Dipoles_Stepping {
	private:
		struct Stepping_Point {
			size_t index;
            double amp;
			double start_time;
			double end_time;
			double width;
		};

        double max_end_time;
		
		std::vector<Stepping_Point> gaussian1_points;
        std::vector<Stepping_Point> gaussian2_points;
	public:
        Gaussian_Dipoles_Stepping() = default;

        //increase possibility of locality
		void organize();
		//time-0.5dt to time+0.5dt
        void step(double time, std::vector<double>& diffdb);
		//add gaussian 1 dipole
		void add_gaussian1(size_t index, double amp, double start_time, double end_time, double width);
        //add gaussian 2 dipole
        void add_gaussian2(size_t index, double amp, double start_time, double end_time, double width);
		//output debug information
		void output(std::ostream& os);
	};
    
    /* PML sources interface */
    class PML_Stepping {
	private:
		struct Stepping_Point {
            size_t index;

			double phi1;
			double b1;
			double c1;
			double k1;

			double phi2;
			double b2;
			double c2;
			double k2;
		};

		sVec3 strides;
		std::vector<Stepping_Point> points[3];
	public:
		PML_Stepping() = default;

        //set strides
		void set_strides(sVec3 strides);

        //reorder points to improve spatial locality
		void organize();

        //stepping, updates diffdb from eh
        void step(const std::vector<double>& eh, std::vector<double>& diffdb);
        
        //add one PML points
		void add_stepping_Point(Direction dir, size_t index, double b1, double c1, double k1, double b2, double c2, double k2);
	};

    class Fields {
    private:
        std::vector<double> eh, eh1, eh2, diffdb, db;    
        Yee3 cell;
        
        PML_Stepping e_pml, m_pml;
        Gaussian_Dipoles_Stepping e_gaussian_dipoles, m_gaussian_dipoles;

    public:
        //add gaussian1 dipole
        void add_dipole_source_gaussian1(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff);

        //add gaussian2 dipole
        void add_dipole_source_gaussian2(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff);

        //PML layer initializations
        void PML_init(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c);
        
        //PML layer initialization helper functions
        template<unsigned int F3>
        void PML_init_helper(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c);

        //getters inside chunk, return 0 if outside chunk
        double get_ex(const fVec3& pos);
        double get_ey(const fVec3& pos);
        double get_ez(const fVec3& pos);
        double get_hx(const fVec3& pos);
        double get_hy(const fVec3& pos);
        double get_hz(const fVec3& pos);
        double get_dx(const fVec3& pos);
        double get_dy(const fVec3& pos);
        double get_dz(const fVec3& pos);
        double get_bx(const fVec3& pos);
        double get_by(const fVec3& pos);
        double get_bz(const fVec3& pos);

        //getter helper functions
        double get_eh_helper(const fVec3& pos, Coord_Type ctype);
        double get_db_helper(const fVec3& pos, Coord_Type ctype);

        //raw getter functions
        double get_eh_raw(const iVec3& pos);
        double get_db_raw(const iVec3& pos);

        //update jwD = curlE - J
        void update_diffd(double time);
        //update jwB = -curlH - M
        void update_diffb(double time);
        //save fields into hdf5 for reloading 
        void save_hdf5();
    };


    template<unsigned int F3>
    void Fields::PML_init_helper(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c) {
        
        constexpr Coord_Type ctype = F3;
        constexpr int x3 = get_dir_from_ctype(ctype);
        constexpr int x1 = (x3 + 1) % 3;
        constexpr int x2 = (x3 + 2) % 3;

        iVec3 p1 = cell.get_grid_p1() + 1;
        iVec3 p2 = cell.get_grid_p2() - 1;
        PML_Stepping& pml = (is_e_point(ctype)? e_pml : m_pml);

        for (auto itr = Yee_Iterator(p1, p2, ctype); !itr.is_end(); itr.next()) {
            iVec3 pos = itr.get_coord();

            size_t index = cell.get_index_from_coord(pos);
            int x1_index = pos.get<x1>() - p1.get<x1>();
            int x2_index = pos.get<x2>() - p1.get<x2>();
            
            pml.add_stepping_Point(x3, index, 
                b[x1][x1_index], c[x1][x1_index] / dx, 1 / dx / k[x1][x1_index],
                b[x2][x2_index], c[x2][x2_index] / dx, 1 / dx / k[x2][x2_index]);
        }
    }
}