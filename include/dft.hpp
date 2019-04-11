#pragma once

#include <utility.hpp>
#include <fields.hpp>
#include <highfive/H5File.hpp>

namespace ffip
{
    using Volume_Metadata = std::tuple<double_arr, double_arr, double_arr>;

    //For maintaining discrete fourier transform of a particular field in specified volume
    struct DFT_Unit
    {
        //start and end point of the unit
        iVec3 p1, p2;

        //one unit only has one coordinate type
        Coord_Type ctype;

        //portions of points in the current process
        iVec3 p1_local, p2_local;

        //size of the array
        size_t size;

        //dimension of the unit
        sVec3 dim;

        //sorted and unique list of frequencies
        double_arr freqs;

        //c-tyle matrix [freq, z, y, x]
        double_arr real, imag;

        //check ctype compatibility
        DFT_Unit(const iVec3 &p1, const iVec3 &p2, Coord_Type ctype);

        //update fourier transform
        void update(double time, const Fields &fields);

        //return dimensions of the chunk
        std::vector<size_t> get_dimension(size_t freq_dim) const;

        //return local portions of the fields
        double_arr select_local_real(const double_arr &freqs) const;
        double_arr select_local_imag(const double_arr &freqs) const;
        double_arr select_local_helper(const double_arr &freqs, const double_arr &data) const;

        //return selection details to be used on hdf5 partial io
        std::pair<std::vector<size_t>, std::vector<size_t>>
        get_local_selection(size_t freq_dim) const;

        //return size
        size_t get_size() const;

        //return local size
        size_t get_local_size() const;

        //add frequencies to update
        void add_frequencies(const double_arr &freqs);

        //make frequencies unique
        void unique_frequencies();

        //whether it contains requested frequencies
        bool has_frequencies(const double_arr &freqs) const;

        //make it update only the region within the chunk
        void set_local_region(const iVec3 &grid_p1, const iVec3 &grid_p2);

        //whether two unit is equal
        bool operator==(const DFT_Unit &other) const;
    };

    //store discrete fourier transform of volume
    class DFT_Hub
    {
    private:
        std::vector<DFT_Unit> units;
        Yee3 grid;
        double dx, dt;

    public:
        double get_dx() const;
        double get_dt() const;

        //initialization
        void set_grid(const Yee3 &grid, double dx, double dt);

        //ndt for electric fields
        void step_e(const double time, const Fields &fields);

        //(n-0.5)dt for magnetic fields
        void step_m(const double time, const Fields &fields);

        //register a volume to get fourier transforms of fields in freq_list
        void register_volume_dft(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs);

        //paralell output to a hdf5 group
        void output_fields_collective(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs, MPI_Comm comm, HighFive::Group &group);

        //return fields in a volume collectively
        complex_arr get_fields_collective(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs, MPI_Comm comm);

        //look for existence of a unit
        DFT_Unit *find_unit(const iVec3 &p1, const iVec3 &p2, Coord_Type ctype);
        DFT_Unit *find_unit(const DFT_Unit &x);

    public:
        //register a flux monitor with norm specified
        void register_flux_monitor(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs);

        //return flux collectively
        double_arr get_flux_collective(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs, MPI_Comm comm);

    private:
        //return local contribution of flux
        double_arr get_flux_local(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs);

        //return local portion of fields
        double_arr get_fields_local(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs);
    };
} // namespace ffip