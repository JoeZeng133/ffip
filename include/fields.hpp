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
        void PML_init_helper(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c) {
            constexpr Coord_Type ctype = F3;
            constexpr int x3 = get_dir_from_ctype(ctype);
            constexpr int x1 = (x3 + 1) % 3;
            constexpr int x2 = (x3 + 2) % 3;

            for (auto itr = my_iterator(cell.get_grid_p1(), cell.get_grid_p2(), ctype); !itr.is_end(); itr.advance()) {
                iVec3 pos = itr.get_vec();
                if (Is_Inside_Box(config.phys_p1, config.phys_p2, pos))
                    continue;

                size_t index = get_index_ch(pos);
                int x1_index = pos.get<x1>() - choose<x1>::get(ch_p1);
                int x2_index = choose<x2>::get(pos) - choose<x2>::get(ch_p1);
                
                w_PML->add_update_point<x3>(index, 
                    b[x1][x1_index], c[x1][x1_index] / dx, 1 / dx / k[x1][x1_index],
                    b[x2][x2_index], c[x2][x2_index] / dx, 1 / dx / k[x2][x2_index]);
            }
        }

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
        double get_field(const iVec3& pos); //raw getter
        //getter helper functions
        double get_eh_helper(const fVec3& pos, Coord_Type ctype);
        double get_db_helper(const fVec3& pos, Coord_Type ctype);

        //update jwD = curlE - J
        void update_diffd(double time);
        //update jwB = -curlH - M
        void update_diffb(double time);
        //save fields into hdf5 for reloading 
        void save_hdf5();
    };


    template<int D>
	void Fields::PML_init_helper(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c) {
		using x3 = typename F::dir_base;
		using x1 = typename F::dir_base::x1;
		using x2 = typename F::dir_base::x2;

		CU_PML * w_PML = (is_E_Point(F::ctype)) ? e_PML : m_PML;

		for (auto itr = my_iterator(ch_p1, ch_p2, F::ctype); !itr.is_end(); itr.advance()) {
			iVec3 pos = itr.get_vec();
			if (Is_Inside_Box(config.phys_p1, config.phys_p2, pos))
				continue;

			size_t index = get_index_ch(pos);
			int x1_index = choose<x1>::get(pos) - choose<x1>::get(ch_p1);
			int x2_index = choose<x2>::get(pos) - choose<x2>::get(ch_p1);
			
			w_PML->add_update_point<x3>(index, 
				b[x1::val][x1_index], c[x1::val][x1_index] / dx, 1 / dx / k[x1::val][x1_index],
				b[x2::val][x2_index], c[x2::val][x2_index] / dx, 1 / dx / k[x2::val][x2_index]);
		}
	}
}