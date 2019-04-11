#include <dft.hpp>

namespace ffip
{
    //DFT_Unit
    DFT_Unit::DFT_Unit(const iVec3 &p1, const iVec3 &p2, Coord_Type ctype)
    {
        if (p1.get_type() != p2.get_type() || p1.get_type() != ctype & 0b111)
            throw std::runtime_error("Invalid Coord Type in DFT_Unit constructor");

        if (!leq_vec3(p1, p2))
            throw std::runtime_error("Invali volume in DFT_Unit constructor");

        p1_local = p1;
        p2_local = p2;

        dim = (p2 - p1) / 2 + 1;
        size = dim.prod();
    }

    void DFT_Unit::update(double time, const Fields &fields)
    {
        size_t index = 0;

        for (int f = 0; f < freqs.size(); ++f)
        {
            double cos_wt = std::cos(2 * pi * freqs[f] * time);
            double sin_wt = -std::sin(2 * pi * freqs[f] * time);
            double *real_ptr = real.data() + f * size;
            double *imag_ptr = imag.data() + f * size;

            if (is_eh_point(ctype))
            {
                size_t index = 0;
                for (auto itr = Yee_Iterator(p1_local, p2_local, ctype); !itr.is_end(); itr.next(), ++index)
                {
                    double val = fields.get_eh_raw(itr.get_coord());
                    imag_ptr[index] += sin_wt * val;
                    real_ptr[index] += cos_wt * val;
                }
            }
            else
            {
                size_t index = 0;
                for (auto itr = Yee_Iterator(p1_local, p2_local, ctype); !itr.is_end(); itr.next(), ++index)
                {
                    double val = fields.get_db_raw(itr.get_coord());
                    imag_ptr[index] += sin_wt * val;
                    real_ptr[index] += cos_wt * val;
                }
            }
        }
    }

    std::pair<std::vector<size_t>, std::vector<size_t>>
    DFT_Unit::get_local_selection(size_t freq_dim) const
    {

        auto offset_vec = (p1_local - p1) / 2;
        auto count_vec = (p2_local - p1_local) / 2 + 1;

        return {{0, offset_vec.z, offset_vec.y, offset_vec.x},
                {freq_dim, count_vec.z, count_vec.y, count_vec.x}};
    }

    std::vector<size_t> DFT_Unit::get_dimension(size_t freq_dim) const
    {
        auto count_vec = (p2_local - p1_local) / 2 + 1;
        return {freq_dim, count_vec.z, count_vec.y, count_vec.x};
    }

    double_arr DFT_Unit::select_local_real(const double_arr &freqs) const
    {
        return select_local_helper(freqs, real);
    }

    double_arr DFT_Unit::select_local_imag(const double_arr &freqs) const
    {
        return select_local_helper(freqs, imag);
    }

    double_arr DFT_Unit::select_local_helper(const double_arr &freqs, const double_arr &data) const
    {

        size_t loc_size = get_local_size();
        double_arr res(freqs.size() * loc_size);

        for (int i = 0; i < freqs.size(); ++i)
        {
            auto itr = std::lower_bound(this->freqs.begin(), this->freqs.end(), freqs[i]);

            if (itr == this->freqs.end() || *itr != freqs[i])
                continue;

            size_t index = std::distance(this->freqs.begin(), itr);

            std::copy(data.begin() + index * loc_size, data.begin() + index * (loc_size + 1), res.begin() + i * loc_size);
        }

        return res;
    }

    size_t DFT_Unit::get_size() const
    {
        return ((p2.z - p1.z) / 2 + 1) * ((p2.y - p1.y) / 2 + 1) * ((p2.x - p1.x) / 2 + 1);
    }

    size_t DFT_Unit::get_local_size() const
    {
        return ((p2_local.z - p1_local.z) / 2 + 1) * ((p2_local.y - p1_local.y) / 2 + 1) * ((p2_local.x - p1_local.x) / 2 + 1);
    }

    void DFT_Unit::add_frequencies(const double_arr &freqs)
    {
        for (auto freq : freqs)
            this->freqs.push_back(freq);
    }

    void DFT_Unit::unique_frequencies()
    {
        std::sort(freqs.begin(), freqs.end());
        auto itr = std::unique(freqs.begin(), freqs.end());
        freqs.resize(std::distance(freqs.begin(), itr));
    }

    bool DFT_Unit::has_frequencies(const double_arr &freqs) const
    {
        for (auto freq : freqs)
            if (auto itr = std::lower_bound(this->freqs.begin(), this->freqs.end(), freq);
                itr == this->freqs.end() || *itr != freq)
                return false;

        return true;
    }

    void DFT_Unit::set_local_region(const iVec3 &grid_p1, const iVec3 &grid_p2)
    {
        auto intersect = get_intersection(p1, p2, grid_p1, grid_p2);
        auto tmp = get_component_interior(intersect.first, intersect.second, ctype);

        p1_local = tmp.first;
        p2_local = tmp.second;
        dim = (p2_local - p1_local) / 2 + 1;
        size = dim.prod();
    }

    bool DFT_Unit::operator==(const DFT_Unit &other) const
    {
        return ctype == other.ctype && p1 == other.p1 && p2 == other.p2;
    }

    //DFT_Hub
    double DFT_Hub::get_dx() const
    {
        return dx;
    }

    double DFT_Hub::get_dt() const
    {
        return dt;
    }

    DFT_Unit *DFT_Hub::find_unit(const iVec3 &p1, const iVec3 &p2, Coord_Type ctype)
    {
        return find_unit({p1, p2, ctype});
    }

    DFT_Unit *DFT_Hub::find_unit(const DFT_Unit &x)
    {
        for (auto &unit : units)
            if (unit == x)
                return &unit;

        return nullptr;
    }

    void DFT_Hub::step_e(const double time, const Fields &fields)
    {
        for (auto &unit : units)
            if (is_e_point(unit.ctype))
            {
                unit.update(time, fields);
            }
    }

    void DFT_Hub::step_m(const double time, const Fields &fields)
    {
        for (auto &unit : units)
            if (is_m_point(unit.ctype))
            {
                unit.update(time, fields);
            }
    }

    void DFT_Hub::register_volume_dft(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const std::vector<double> &freqs)
    {

        auto closure = get_component_closure(p1, p2, ctype);
        auto new_unit = DFT_Unit(closure.first, closure.second, ctype);

        //check whether the unit alread exists
        if (auto unit_ptr = find_unit(new_unit); unit_ptr)
        {
            unit_ptr->add_frequencies(freqs);
            return;
        }

        new_unit.set_local_region(grid.get_grid_p1(), grid.get_grid_p2());
        new_unit.add_frequencies(freqs);
        units.push_back(std::move(new_unit));
    }

    void DFT_Hub::output_fields_collective(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const std::vector<double> &freqs, MPI_Comm comm, HighFive::Group &group)
    {

        auto closure = get_component_closure(p1, p2, ctype);
        auto unit_ptr = find_unit(closure.first, closure.second, ctype);

        if (unit_ptr == nullptr)
            throw std::runtime_error("Volume not registered");

        if (!unit_ptr->has_frequencies(freqs))
            throw std::runtime_error("Frequencies not registered");

        auto selection = unit_ptr->get_local_selection(freqs.size());
        auto dim = unit_ptr->get_dimension(freqs.size());

        //dataset creation
        auto re_dataset = group.createDataSet<double>("real", HighFive::DataSpace(dim));
        auto im_dataset = group.createDataSet<double>("imag", HighFive::DataSpace(dim));

        //independent write to different slabs of dataset
        re_dataset.select(selection.first, selection.second).write(unit_ptr->select_local_real(freqs).data());
        im_dataset.select(selection.first, selection.second).write(unit_ptr->select_local_imag(freqs).data());

        //metadata creation
        auto x_dataset = group.createDataSet<double>("x", HighFive::DataSpace({dim[3]}));
        auto y_dataset = group.createDataSet<double>("y", HighFive::DataSpace({dim[2]}));
        auto z_dataset = group.createDataSet<double>("z", HighFive::DataSpace({dim[1]}));
        auto f_dataset = group.createDataSet<double>("f", HighFive::DataSpace({dim[0]}));

        int rank;
        MPI_Comm_rank(comm, &rank);

        if (rank == 0)
        {
            x_dataset.write(linspace(closure.first.x, closure.second.x, 2).data());
            y_dataset.write(linspace(closure.first.y, closure.second.y, 2).data());
            z_dataset.write(linspace(closure.first.z, closure.second.z, 2).data());
            f_dataset.write(freqs.data());
        }
    }

    complex_arr DFT_Hub::get_fields_collective(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs, MPI_Comm comm)
    {
        //yet to be implemented
        return {};
    }

    void DFT_Hub::register_flux_monitor(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs)
    {
        if (!leq_vec3(p1, p2))
            throw std::runtime_error("Invalid face: register_flux_monitor");

        //check normal vector
        if (std::abs(norm.x) + std::abs(norm.y) + std::abs(norm.z) != 1)
            throw std::runtime_error("Invalid Normal Vector: register_flux_monitor");

        //get x1, x2, x3
        int x3, x2, x1;
        if (norm.get<0>() != 0)
            x3 = 0;
        if (norm.get<1>() != 0)
            x3 = 1;
        if (norm.get<2>() != 0)
            x3 = 2;
        x1 = (x2 + 1) % 3;
        x2 = (x2 + 2) % 3;
        auto E1 = get_e_ctype_from_dir(static_cast<Direction>(x1));
        auto H2 = get_m_ctype_from_dir(static_cast<Direction>(x2));

        if (p1[x3] != p2[3])
            throw std::runtime_error("Invalid face: register_flux_monitor");

        register_volume_dft(p1, p2, E1, freqs);
        register_volume_dft(p1, p2, H2, freqs);
    }

    double_arr DFT_Hub::get_flux_collective(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs, MPI_Comm comm)
    {
        auto flux_local = get_flux_local(p1, p2, norm, freqs);
        double_arr res(freqs.size());

        MPI_Allreduce(flux_local.data(), res.data(), freqs.size(), MPI_DOUBLE, MPI_SUM, comm);
        return res;
    }

    double_arr DFT_Hub::get_flux_local(const fVec3 &p1, const fVec3 &p2, const iVec3 &norm, const double_arr &freqs)
    {

        //check normal vector
        if (std::abs(norm.x) + std::abs(norm.y) + std::abs(norm.z) != 1)
            throw std::runtime_error("Invalid Normal Vector");

        //get x1, x2, x3
        int x3, x2, x1;
        if (norm.get<0>() != 0)
            x3 = 0;
        if (norm.get<1>() != 0)
            x3 = 1;
        if (norm.get<2>() != 0)
            x3 = 2;
        x1 = (x2 + 1) % 3;
        x2 = (x2 + 2) % 3;
        auto E1 = get_e_ctype_from_dir(static_cast<Direction>(x1));
        auto H2 = get_m_ctype_from_dir(static_cast<Direction>(x2));

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
        double w = dlen.prod() * (dx * dx * dx);
        double_arr res(freqs.size());

        for (int f = 0; f < freqs.size(); ++f)
        {

            auto real_e = unit_ptr_e->select_local_real({freqs[f]});
            auto imag_e = unit_ptr_e->select_local_imag({freqs[f]});
            auto real_m = unit_ptr_m->select_local_real({freqs[f]});
            auto imag_m = unit_ptr_m->select_local_imag({freqs[f]});

            for (int k = 0; k <= count.z; ++k)
                for (int j = 0; j <= count.y; ++j)
                    for (int i = 0; i <= count.x; ++i)
                    {

                        auto p = p1 + dlen * iVec3{i, j, k};
                        double trap_w = ((i == 0 || i == count.x) ? 0.5 : 1) *
                                        ((j == 0 || j == count.y) ? 0.5 : 1) *
                                        ((k == 0 || k == count.z) ? 0.5 : 1);

                        res[f] += trap_w * 0.5 * (grid_e(real_e, p) * grid_m(real_m, p) + grid_e(imag_e, p) * grid_m(imag_m, p));
                    }

            res[f] *= w;
        }

        return res;
    }

    double_arr DFT_Hub::get_fields_local(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs)
    {
        //yet to be implemented
        return {};
    }

    void DFT_Hub::set_grid(const Yee3 &grid, double dx, double dt)
    {
        this->grid = grid;
        this->dx = dx;
        this->dt = dt;
    }
} // namespace ffip