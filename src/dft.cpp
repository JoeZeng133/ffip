#include <dft.hpp>

namespace ffip {
    //DFT_Unit

    //DFT_Hub
    DFT_Unit* DFT_Hub::find_unit(const iVec3& p1, const iVec3& p2, Coord_Type ctype) {
        return find_unit({p1, p2, ctype});
    }

    DFT_Unit* DFT_Hub::find_unit(const DFT_Unit& x) {
        for(auto& unit : units)
            if (unit == x) 
                return &unit;
        
        return nullptr;
    }

    void DFT_Hub::register_volume_dft(const fVec3& p1, const fVec3& p2, Coord_Type ctype, const std::vector<double>& freq_list) {
        auto closure = get_component_closure(p1, p2, ctype);

        auto new_unit = DFT_Unit(closure.first, closure.second, ctype);
        
        //check whether the unit alread exists
        if(auto unit_ptr = find_unit(new_unit); unit_ptr) {
            unit_ptr->add_frequencies(freq_list);
        }
    
        new_unit.add_frequencies(freq_list);
        units.push_back(std::move(new_unit));
    }

    void DFT_Hub::output_fields_collective
    (const fVec3& p1, const fVec3& p2, Coord_Type ctype, const std::vector<double>& freqs, MPI_Comm comm, HighFive::Group& group) {
        
        auto closure = get_component_closure(p1, p2, ctype);
        auto unit_ptr = find_unit(closure.first, closure.second, ctype);

        if (unit_ptr == nullptr)
            throw std::runtime_error("Volume not registered");

        if (!unit_ptr->has_frequencies(freqs))
            throw std::runtime_error("Frequencies not registered");

        auto selection = unit_ptr->get_selection(freqs.size());
        auto dim = unit_ptr->get_dimension(freqs.size());
        
        //dataset creation
        auto re_dataset = group.createDataSet<double>("real", HighFive::DataSpace(dim));
        auto im_dataset = group.createDataSet<double>("imag", HighFive::DataSpace(dim));

        //independent write to different slabs of dataset
        re_dataset.select(selection.first, selection.second).write(unit_ptr->select_real(freqs).data());
        im_dataset.select(selection.first, selection.second).write(unit_ptr->select_imag(freqs).data());

        //metadata creation
        auto x_dataset = group.createDataSet<double>("x", HighFive::DataSpace({dim[3]}));
        auto y_dataset = group.createDataSet<double>("y", HighFive::DataSpace({dim[2]}));
        auto z_dataset = group.createDataSet<double>("z", HighFive::DataSpace({dim[1]}));
        auto f_dataset = group.createDataSet<double>("f", HighFive::DataSpace({dim[0]}));

        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0) {
            x_dataset.write(linspace(closure.first.x, closure.second.x, 2).data());
            y_dataset.write(linspace(closure.first.y, closure.second.y, 2).data());
            z_dataset.write(linspace(closure.first.z, closure.second.z, 2).data());
            f_dataset.write(freqs.data());
        }
    }

    void DFT_Hub::output_flux_collective
    (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs, MPI_Comm comm, HighFive::Group& group) {
        auto flux = get_flux_collective(p1, p2, norm, freqs, comm);
        int rank;
        MPI_Comm_rank(comm, &rank);
        auto dataset = group.createDataSet<double>("flux", HighFive::DataSpace({freqs.size()}));

        if (rank == 0) {
            dataset.write(flux.data());
        }
    }

    double_arr DFT_Hub::get_flux_collective
    (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs, MPI_Comm comm) {
        auto flux_local = get_flux_local(p1, p2, norm, freqs);
        double_arr res(freqs.size());

        MPI_Allreduce(flux_local.data(), res.data(), freqs.size(), MPI_DOUBLE, MPI_SUM, comm);
        return res;
    }

    double_arr DFT_Hub::get_flux_local
    (const fVec3& p1, const fVec3& p2, const iVec3& norm, const double_arr& freqs) {
        
        //check normal vector
        if (std::abs(norm.x) + std::abs(norm.y) + std::abs(norm.z) != 1)
            throw std::runtime_error("Invalid Normal Vector");

        //get x1, x2, x3
        int x3, x2, x1;
        if (norm.get<0>() != 0) x3 = 0;
        if (norm.get<1>() != 0) x3 = 1;
        if (norm.get<2>() != 0) x3 = 2;
        x1 = (x2 + 1) % 3;
        x2 = (x2 + 2) % 3;
        auto E1 = get_e_ctype_from_dir(x1);
        auto H2 = get_m_ctype_from_dir(x2);

        //check whehter the volume is registered
        auto closure_e = get_component_closure(p1, p2, E1);
        auto closure_m = get_component_closure(p1, p2, H2);

        auto unit_ptr_e = find_unit(closure_e.first, closure_e.second, E1);
        auto unit_ptr_m = find_unit(closure_m.first, closure_m.second, H2);
        
        if (unit_ptr_e == nullptr || unit_ptr_m == nullptr)
            throw std::runtime_error("Volume not registered");

        if (!unit_ptr_e->has_frequencies(freqs))
            throw std::runtime_error("Frequencies not registered");

        //retrieve real and imaginary part
        auto grid_e = Grid_3(unit_ptr_e->p1_local, unit_ptr_e->p2_local);
        auto grid_m = Grid_3(unit_ptr_m->p1_local, unit_ptr_m->p2_local);
        
        iVec3 count(((p2 - p1) / 2).round());
        fVec3 dlen = (p2 - p1) / count;
        double w = dlen.prod();
        double_arr res(freqs.size());

        for(int f = 0; f < freqs.size(); ++f) {

            auto real_e = unit_ptr_e->select_real({freqs[f]});
            auto imag_e = unit_ptr_e->select_imag({freqs[f]});
            auto real_m = unit_ptr_m->select_real({freqs[f]});
            auto imag_m = unit_ptr_m->select_imag({freqs[f]});

            for(int k = 0; k <= count.z; ++k)
            for(int j = 0; j <= count.y; ++j)
            for(int i = 0; i <= count.x; ++i) {

                auto p = p1 + dlen * iVec3{i, j, k};
                double trap_w = ((i == 0 || i == count.x)? 0.5 : 1) * 
                                ((j == 0 || j == count.y)? 0.5 : 1) *
                                ((k == 0 || k == count.z)? 0.5 : 1);

                res[f] += trap_w * 0.5 * (
                    grid_e.interp(real_e, p) * grid_m.interp(real_m, p) + 
                    grid_e.interp(imag_e, p) * grid_m.interp(imag_m, p));
            }

            res[f] *= w;
        }
        
        return res;
    }
    
}