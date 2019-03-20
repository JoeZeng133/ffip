#pragma once

#include <utility.hpp>
#include <structure.hpp>

namespace ffip {
    using Volume_Metadata = std::tuple<double_arr, double_arr, double_arr>;
    using Yee3_Metadata = std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>;

    /* store and update Fourier transforms of fields */
    class Discrete_Fourier_Chunk {
    private:
        

    public:
        //update fourier transforms
        void step(const double time, const Structure& structure);

        //register a volume to get fourier transforms of fields in freq_list
        void register_volume_dft(const fVec3& p1, const fVec3& p2, const std::vector<double>& freq_list, const Coord_Type ctype);
        
        //get fields out of the volume, only portion of the current chunk is returned
        std::complex<double> get_volume_dft(const fVec3& p1, const fVec3& p2, const std::vector<double>& freq_list, const Coord_Type ctype) const;

        //parallel version, gather fields to one node, entire fields are returned
        std::complex<double> get_volume_dft_parallel(const fVec3& p1, const fVec3& p2, const double_arr& freq_list, const Coord_Type ctype, MPI_Comm comm, int dest) const;

        //get metadata out of the volume in order to build griddedinterpolation
        Volume_Metadata get_volume_metadata(const fVec3& p1, const fVec3& p2, const Coord_Type ctype) const;

        //get grid coordinates of the metadata in the volume
        Yee3_Metadata get_grid_metadata(const fVec3& p1, const fVec3& p2, const Coord_Type ctype) const;

        //used in get_volume_dft_parallel
        void send_grid_metadata(const fVec3& p1, const fVec3& p2, const Coord_Type ctype, MPI_Comm comm, int dest) const;

        //output to hdf5 file
        void output_hdf5() const;
    };
}