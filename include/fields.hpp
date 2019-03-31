#pragma once

#include <utility.hpp>


namespace ffip {
    //Dipole Sources (Gaussian Types)
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
        void step(double time, std::vector<double>& accdb);

		//add gaussian 1 dipole
		void add_gaussian1(size_t index, double amp, double start_time, double end_time, double width);

        //add gaussian 2 dipole
        void add_gaussian2(size_t index, double amp, double start_time, double end_time, double width);

		//output debug information
		void output(std::ostream& os);
	};
    
    //PML Regions curls E, H
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

        //set strides for use in stepping
		void set_strides(sVec3 strides);

        //reorder points to improve spatial locality
		void organize();

        //stepping, updates accdb from eh
        void step(const std::vector<double>& eh, std::vector<double>& accdb);
        
        //add one PML points
		void add_stepping_point(Direction dir, size_t index, double b1, double c1, double k1, double b2, double c2, double k2);
	};

    //Non PML Regions
    class Curl_Stepping {
    private:
        sVec3 strides;
        std::vector<size_t> points[3];
        double dx;

    public:
        Curl_Stepping() = default;

        //set dx for calculating curl
        void set_dx(double dx);

        //strides for calculating index offset in 3 direction
        void set_strides(sVec3 strides);

        //reorder to improve locality
        void organize();

        //stepping, updates accdb from eh
        void step(const std::vector<double>& eh, std::vector<double>& accdb);

        //add a stepping point, direction has to be specified
        void add_stepping_point(Direction dir, size_t index);
        
    };

    enum Boundary_Condition {
        PEC, //same as symmetry with phase -1 
        PMC, //same as symmetry with phase 1
        None, //ghost points forced to be zeros
        Sync //ghost points comming from another chunk
    };


    //Store fields 
    //Responsible for updating curl F + C, aware of signs
    class Fields {
    private:
        //eh = [E, H], accdb = sum(-curl E - M) or sum(curl H - J))
        //d = accdb * dt;
        std::vector<double> eh, accdb;  
        //tentative fields, might not be necessary
        std::vector<double> eh1, eh2, db;
        //cell information for the current process
        Yee3 cell;
        double dx, dt;
        
        PML_Stepping e_pml, m_pml;
        Curl_Stepping e_curl, m_curl;
        Gaussian_Dipoles_Stepping e_gaussian_dipoles, m_gaussian_dipoles;

        Boundary_Condition bc[3][2];

    public:
        Fields() = default;

        //add gaussian1 dipole, flip amplitude
        void add_dipole_source_gaussian1(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff);

        //add gaussian2 dipole
        void add_dipole_source_gaussian2(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff);

        //Initializations of curl
        //providing k, b, c from PML and the index extend [p1, p2] covering the whole simulation
        //provide dx for discretized curl calculation
        void init(
                const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c, 
                const iVec3& p1, const iVec3& p2,
                double dx, double dt, const Yee3& cell);
        
        //PML layer initialization helper functions
        template<unsigned int F3>
        void PML_init_helper(
            const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c, 
            const iVec3& p1, const iVec3& p2);

        //getter helper functions
        double get_eh_helper(const fVec3& pos, Coord_Type ctype);
        double get_db_helper(const fVec3& pos, Coord_Type ctype);

        //raw getter functions, no range checking
        double get_eh_raw(const iVec3& pos);
        double get_db_raw(const iVec3& pos);

        //synchronize ghost points of boundaries in a particular direction
        void sync_boundary(MPI_Comm comm, Direction dir);

        //extract fields on a particular face and store it in the buffer
        template<unsigned int x3, int side>
        void extract_fields_face(std::vector<double>& buf);

        //assign fields to a particular face, reject empty buffer
        template<unsigned int x3, int side>
        void assign_fields_face(const std::vector<double>& buf);

        //symmetry boundary with phase = 1, -1
        template<unsigned int dir, int side>
        void sync_symmetry_boundary(int phase);

        //update accd = curlH - J
        void update_accd(MPI_Comm comm, double time);
        //update accb = -curlE - M
        void update_accb(MPI_Comm comm, double time);
        //save fields into hdf5 for reloading 
        void save_hdf5();
    };

    template<unsigned int x3, int side>
    void Fields::extract_fields_face(std::vector<double>& buf) {

        iVec3 p1 = cell.get_grid_p1();
        iVec3 p2 = cell.get_grid_p2();

        auto face = get_face<x3, side>(p1, p2);
        auto itr = Yee_Iterator(p1, p2, All);
        buf.resize(itr.get_size());

        for(int index = 0; !itr.is_end(); itr.next(), ++index) {
            buf[index] = cell.get_raw_val(eh, itr.get_coord());
        }
    }

    template<unsigned int x3, int side>
    void Fields::assign_fields_face(const std::vector<double>& buf) {
        constexpr iVec3 face_norm = get_norm_vec(x3, side);
        iVec3 p1 = cell.get_grid_p1() + face_norm;
        iVec3 p2 = cell.get_grid_p2() + face_norm;

        auto face = get_face<x3, side>(p1, p2);
        auto itr = Yee_Iterator(p1, p2, All);
        if (buf.size() < itr.get_size()) return;

        for(int index = 0; !itr.is_end(); itr.next(), ++index) {
            eh[cell.get_index_from_coord(itr.get_coord())] = buf[index];
        }
    }

    template<unsigned int x3, int side>
    void Fields::sync_symmetry_boundary(int phase) {
        
        //orientation (x1, x2, x3)
        constexpr int x1 = (x3 + 1) % 3;
        constexpr int x2 = (x3 + 2) % 3;

        //construct face
        constexpr iVec3 face_norm = get_norm_vec(x3, side);
        auto p1 = cell.get_grid_p1() + face_norm;
        auto p2 = cell.get_grid_p2() + face_norm;
        auto face = get_face<x3, side>(p1, p2);

        //index offset in the normal direction
        long long norm_index_offset = cell.get_index_offset(face_norm);
        //actual phase applied
        int mult = (p1.get<x3>() & 1) ? -phase : phase; 

        //loop through points of inner faces
        for(auto itr = Yee_Iterator(face, All); !itr.is_end(); itr.next()) {
        
            //index of the point
            size_t index = cell.get_index_from_coord(itr.get_coord());    
            eh[index] = mult * eh[index - norm_index_offset * 2];
        }

    }

    template<unsigned int F3>
    void Fields::PML_init_helper(const std::array<double_arr, 3>& k, const std::array<double_arr, 3>& b, const std::array<double_arr, 3>& c, const iVec3& p1, const iVec3& p2) {
        
        constexpr Coord_Type ctype = F3;
        constexpr int x3 = get_dir_from_ctype(ctype);
        constexpr int x1 = (x3 + 1) % 3;
        constexpr int x2 = (x3 + 2) % 3;

        PML_Stepping& pml = (is_e_point(ctype)? e_pml : m_pml);
        Curl_Stepping& curl = (is_e_point(ctype)? e_curl : m_curl);

        for (auto itr = Yee_Iterator(cell.get_grid_p1(), cell.get_grid_p2(), ctype); !itr.is_end(); itr.next()) {
            iVec3 pos = itr.get_coord();

            size_t index = cell.get_index_from_coord(pos);
            int x1_index = pos.get<x1>() - p1.get<x1>();
            int x2_index = pos.get<x2>() - p1.get<x2>();
            
            //if it is a PML point
            if (c[x1][x1_index] != 0  || c[x2][x2_index] != 0)
                pml.add_stepping_Point(x3, index, 
                b[x1][x1_index], c[x1][x1_index] / dx, 1 / dx / k[x1][x1_index],
                b[x2][x2_index], c[x2][x2_index] / dx, 1 / dx / k[x2][x2_index]);
            //if it is not a PML point
            else
                curl.add_stepping_point(x3, index);
        }
    }
}