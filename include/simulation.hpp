#pragma once

#include <fields.hpp>
#include <geometry.hpp>
#include <medium.hpp>
#include <structure.hpp>
#include <dft.hpp>
#include <chrono>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

namespace ffip
{
    //json string to coord type
    Coord_Type json2ctype(const json &j);

    //json number array to ivec3
    iVec3 json2ivec3(const json &j);

    //json number array to fvec3
    fVec3 json2fvec3(const json &j);

    //json object to medium
    Medium json2medium(const json &medium_json);

    //json object to susceptibility
    Susceptibility json2susceptibility(const json &sus_json);

    //json number array to vector
    double_arr json2double_arr(const json &j);

    //json string to direction
    Direction json2direction(const json &j);

    //json string to side
    Side json2side(const json& j);

    //flux spectrum around a box
    struct Box_Flux
    {
        //box coordinates
        fVec3 p1, p2;
        //frequencies
        double_arr freqs;
        //output group name string
        std::string dataset_name;

        Box_Flux(const json &config, DFT_Hub &hub);

        void output(MPI_Comm comm, DFT_Hub &hub, HighFive::File &file);
    };

    //gridded fields in a volume
    struct Volume_Fields_DFT
    {
        //box coordinates
        fVec3 p1, p2;
        //frequencies
        double_arr freqs;
        //Coordinate type
        Coord_Type ctype;
        //output group name string
        std::string group_name;

        Volume_Fields_DFT(const json &config, DFT_Hub &hub);

        void output(MPI_Comm comm, DFT_Hub &hub, HighFive::File &file);
    };

    //For converting json configurations into simulation
    class Simulation
    {
    private:
        Fields fields;

        Structure structure;

        DFT_Hub dft_hub;

        double dx, dt, resolution, dx2;

        int np, rank, coords[3];

        //physical size for the simulation
        fVec3 sim_size;

        //adjusted corner points in grid coordinates
        iVec3 sim_p1, sim_p2, sim_span, sim_dim;

        //local grid corner points
        iVec3 grid_p1, grid_p2;

        Yee3 grid;

        Medium bg_medium;

        MPI_Comm cart_comm;

        Boundary_Condition bcs[3][2];

        HighFive::File *input_file{nullptr}, *fields_output_file{nullptr};

        std::vector<Volume_Fields_DFT> volume_fields_dft;

        std::vector<Box_Flux> box_flux;

        int time_step;

        double progress_interval{4};

        //simulation start time, previous report time
        decltype(std::chrono::system_clock::now()) sim_start_time, sim_prev_time;

    public:
        //return field component at a particular point
		double at(const fVec3& pt, const Coord_Type ctype);

        //initialize configurations using json file
        void init(const json &config);

        //set up symmetry(PEC, PMC), sync boundary conditions for the whole simulation
        void set_boundary_conditions(const json &config);

        //round size of cell to be multiple of dx
        void set_size(const json &config);

        //PML regions are calculated based on original size,
        //so it might be either truncated or elongated
        void set_pmls(const json &pmls);

        //decompose grid into chunks, and set up boundary conditions
        void set_grid();

        void set_geometry(const json &geoms);

        void set_source(const json &src);

        void set_volume_source(const json &src);

        void set_scattered_source(const json &src);

        void set_fields_output(const json &j);

        //n to n + 1
        void step_e();

        //(n-0.5) to (n + 0.5)
        void step_m();

        //run until reaching a certain time
        void run_until_time(double time);

        //run until the absolute squared value of field at a particular point decays to decay_by
        void run_until_fields_decayed(double interval, const fVec3 &pt, const Coord_Type ctype, double decay_by);

        //run the simulation with a stop condition
        void run(const json &stop_cond);

        //output everything
        void output();

        //print information every progress interval (seconds)
        void report();

        //json conversions
        std::reference_wrapper<Geometry> build_geometry_from_json(const json &geom_json);

        Abstract_Medium build_abstract_medium_from_json(const json &medium_json);
    };

    //non gridded fields output
    // struct Scattered_Fields_DFT {
    //     double_arr x, y, z, freqs;
    //     std::vector<Coord_Type> ctypes;
    //     std::string group_str;

    //     double_arr real, imag;

    //     //configure from json
    //     Scattered_Fields_DFT(const json& config);

    //     //update each time
    //     void step_e(double time, const Fields& fields);
    //     void step_m(double time, const Fields& fields);

    //     //output collectively to a data group
    //     void output_collective(MPI_Comm comm);
    // };

    // //gridded fields in time domain
    // struct Volume_Fields {
    //     double_arr x, y, z;
    //     double start_time, end_time;
    //     Coord_Type ctype;
    //     std::string group_str;

    //     Volume_Fields(const json& config);

    //     void step_e(double time, const Fields& fields);
    //     void step_m(double time, const Fields& fields);
    // };
} // namespace ffip
