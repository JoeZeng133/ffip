#pragma once

#include <utility.hpp>
#include <fields.hpp>
#include <highfive/H5File.hpp>

namespace ffip {
    using Volume_Metadata = std::tuple<double_arr, double_arr, double_arr>;
    using Yee3_Metadata = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;

    // struct Frequencies_Not_Registered : public std::exception{
    //     std::string str;

    //     const char* what () const throw () {
    //         return str.c_str();
    //     }
    // };

    // struct Volume_Not_Registered : public std::exception{
    //     const char* what () const throw () {
    //         return "Frequencies not registered";
    //     }
    // };

    struct DFT_Unit {
        //start and end point of the unit
        iVec3 p1, p2;
        //one unit only has one coordinate type
        Coord_Type ctype;
        //portions of points in the current process
        iVec3 p1_local, p2_local;
        //sorted and unique list of frequencies
        double_arr freq;
        //stored as MATLAB style [x, y, z, freq]
        double_arr real, imag;

        //check ctype compatibility
        DFT_Unit(const iVec3& p1, const iVec3& p2, Coord_Type ctype);

        //return (offset, count) corresponding to select_real function
        std::pair<std::vector<size_t>, std::vector<size_t>>
        get_selection(const iVec3& s, const iVec3& e, size_t freq_dim) const;

        //direct to get_selection(p1_local, p2_local, freq_dim)
        std::pair<std::vector<size_t>, std::vector<size_t>>
        get_selection(size_t freq_dim) const;

        //return dimensions of the chunk
        std::vector<size_t> get_dimension(size_t freq_dim) const;
        
        //select real parts of the chunk, non local points are zeros, no check
        double_arr select_real(const iVec3& s, const iVec3& e, const double_arr& freqs) const;

        //direct to select_real(p1_loca, p2_local, freqs)
        double_arr select_real(const double_arr& freqs) const;

        //select real parts of the chunk, non local points are zeros, no check
        double_arr select_imag(const iVec3& s, const iVec3& e, const double_arr& freqs) const;

        //direct to select_imag(p1_local, p2_local, freqs)
        double_arr select_imag(const double_arr& freqs) const;

        //select helper
        double_arr select_helper(const iVec3& s, const iVec3& e, const double_arr& freqs, const double_arr& data) const;

        //add frequencies to update
        void add_frequencies(const double_arr& freqs);

        //make frequencies unique
        void unique_frequencies();

        //whether it contains requested frequencies
        bool has_frequencies(const double_arr& freqs) const;

        //whether two unit is equal
        bool operator==(const DFT_Unit& other) const;
    };

    //store discrete fourier transform
    class DFT_Hub {
    private:
        std::vector<DFT_Unit> units;
        Yee3 grid;
        double dx;

    public:
        //update fourier transforms
        void step(const double time, const Fields& fields);

        //register a volume to get fourier transforms of fields in freq_list
        void register_volume_dft(const fVec3& p1, const fVec3& p2, Coord_Type ctype, const double_arr& freqs);

        //paralell output to a hdf5 group
        void output_fields_collective
        (const fVec3& p1, const fVec3& p2, Coord_Type ctype, const double_arr& freqs, MPI_Comm comm, HighFive::Group& group);

        //return fields in a collective function call
        complex_arr get_fields_collective
        (const fVec3& p1, const fVec3& p2, Coord_Type ctype, const double_arr& freqs, MPI_Comm comm);

        //look for existence of a unit
        DFT_Unit* find_unit(const iVec3& p1, const iVec3& p2, Coord_Type ctype);
        DFT_Unit* find_unit(const DFT_Unit& x);

        //return metadata for the volume
        Volume_Metadata get_metadata(const fVec3& p1, const fVec3& p2, Coord_Type ctype);
    public:
        //register a flux monitor with norm specified
        void register_flux_monitor
        (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs);

        //output flux collectively
        void output_flux_collective
        (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs, MPI_Comm comm, HighFive::Group& group);

        //return flux collectively
        double_arr get_flux_collective
        (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs, MPI_Comm comm);

    private:
        double_arr get_flux_local
        (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs);

        //get local portion of fields
        complex_arr get_fields_local
        (const fVec3& p1, const fVec3& p2, Coord_Type ctype, const double_arr& freqs);

        //get metadata in single process
        Volume_Metadata get_metadata_local(const fVec3& p1, const fVec3& p2, Coord_Type ctype);
    };
}