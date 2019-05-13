#include <dft.hpp>

namespace ffip
{
    //DFT_Unit
	DFT_Unit::DFT_Unit(const iVec3 &p1, const iVec3 &p2, Coord_Type ctype): p1(p1), p2(p2), ctype(ctype)
    {
        if (p1.get_type() != p2.get_type() || p1.get_type() != (ctype & 0b111))
        {
            std::cout << "p1=" << p1 << ",p2=" << p2 << ",ctype=" << ctype << "\n";

            throw std::runtime_error("Invalid Coord Type in DFT_Unit constructor");
            
        }
            

        if (!leq_vec3(p1, p2))
            throw std::runtime_error("Invali volume in DFT_Unit constructor");
		
        p1_local = p1;
        p2_local = p2;

        local_dim = dim = (p2 - p1) / 2 + 1;
        local_size = size = dim.prod();
    }

    DFT_Unit::DFT_Unit(const pair_Vec3<int>& vol, Coord_Type ctype): DFT_Unit(vol.first, vol.second, ctype) {}

    void DFT_Unit::update(double time, const Fields &fields)
    {
        if (local_size == 0) return;

        buf.resize(get_local_size());
        if (is_eh_point(ctype)) {
            size_t index = 0;
            for (auto itr = Yee_Iterator(p1_local, p2_local, ctype); !itr.is_end(); itr.next(), ++index)
            {
                buf[index] = fields.get_eh_raw(itr.get_coord());
            }
        }
        else
        {
            size_t index = 0;
            for (auto itr = Yee_Iterator(p1_local, p2_local, ctype); !itr.is_end(); itr.next(), ++index)
            {
                buf[index] = fields.get_db_raw(itr.get_coord());
            }
        }
        
        for (int f = 0; f < freqs.size(); ++f)
        {
            double cos_wt = std::cos(2 * pi * freqs[f] * time);
            double sin_wt = -std::sin(2 * pi * freqs[f] * time);
            double *real_ptr = real.data() + f * local_size;
            double *imag_ptr = imag.data() + f * local_size;

            for(size_t index = 0; index < local_size; ++index) {    
                imag_ptr[index] += sin_wt * buf[index];
                real_ptr[index] += cos_wt * buf[index];
            }
        }
    }

    std::pair<std::vector<size_t>, std::vector<size_t>>
    DFT_Unit::get_local_selection(size_t freq_dim) const
    {
        Vec3<size_t> offset_vec = (p1_local - p1) / 2;

        return {{0, offset_vec.z, offset_vec.y, offset_vec.x},
                get_local_dimension(freq_dim)};
    }

    std::vector<size_t> DFT_Unit::get_local_dimension(size_t freq_dim) const
    {
        return {freq_dim, (size_t)local_dim.z, (size_t)local_dim.y, (size_t)local_dim.x};
    }

    std::vector<size_t> DFT_Unit::get_dimension(size_t freq_dim) const
    {
        return {freq_dim, (size_t)dim.z, (size_t)dim.y, (size_t)dim.x};
    }

    double_arr DFT_Unit::select_local_real(const double_arr &freqs) const
    {
        return select_local_helper(freqs, real);
    }

    double_arr DFT_Unit::select_local_imag(const double_arr &freqs) const
    {
        return select_local_helper(freqs, imag);
    }

    double DFT_Unit::get_local_norm(const double freq) const
    {
        auto re = select_local_real({freq});
        auto im = select_local_imag({freq});

        double ans = 0;
        for(size_t i = 0; i < re.size(); ++i)
        {
            double val = re[i]*re[i] + im[i] * im[i];
            if (val > ans) ans = val;
        }

        return ans;
    }

    double_arr DFT_Unit::select_local_helper(const double_arr &freqs, const double_arr &data) const
    {
        if (local_size == 0) return {};

        double_arr res(freqs.size() * local_size);

        for (int i = 0; i < freqs.size(); ++i)
        {
            auto itr = std::lower_bound(this->freqs.begin(), this->freqs.end(), freqs[i]);

            if (itr == this->freqs.end() || *itr != freqs[i])
                continue;

            size_t index = std::distance(this->freqs.begin(), itr);

            std::copy(data.begin() + index * local_size, data.begin() + (index + 1) * local_size, res.begin() + i * local_size);
        }

        return res;
    }

    size_t DFT_Unit::get_size() const
    {
        return size;
    }

    size_t DFT_Unit::get_local_size() const
    {
        return local_size;
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
        if (!is_intersect(p1, p2, grid_p1, grid_p2))
        {
            //for degenerate case, make sure the h5 write offset is valid
            local_dim = {0, 0, 0};
            local_size = 0;
            return;
        }

        auto intersect = get_intersection(p1, p2, grid_p1, grid_p2);
        auto tmp = get_component_interior(intersect.first, intersect.second, ctype);

        p1_local = tmp.first;
        p2_local = tmp.second;
        local_dim = (p2_local - p1_local) / 2 + 1;
        local_size = local_dim.prod();
    }

    bool DFT_Unit::operator==(const DFT_Unit &other) const
    {
        return ctype == other.ctype && p1 == other.p1 && p2 == other.p2;
    }
	
	void DFT_Unit::init()
	{
		unique_frequencies();
		real.resize(local_size * freqs.size(), 0);
		imag.resize(local_size * freqs.size(), 0);
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
        if (is_e_point(x.ctype)){
            for (auto &unit : e_units)
                if (unit == x)
                    return &unit;
        }
        else
        {
            for (auto &unit : m_units)
                if (unit == x)
                    return &unit;
        }
        
        return nullptr;
    }

	//DFT_Hub
	void DFT_Hub::init()
	{
		for (auto &unit : e_units)
		{
			unit.init();
		}
		
		for (auto &unit : m_units)
		{
			unit.init();
		}
	}
	
    void DFT_Hub::step_e(const double time, const Fields &fields)
    {
        for (auto &unit : e_units)
            unit.update(time, fields);
    }

    void DFT_Hub::step_m(const double time, const Fields &fields)
    {
        for (auto &unit : m_units)
            unit.update(time, fields);
    }

    void DFT_Hub::register_volume_dft(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const std::vector<double> &freqs)
    {
        auto closure = get_component_closure(p1, p2, ctype);
        auto new_unit = DFT_Unit(closure, ctype);

        //check whether the unit alread exists
        if (auto unit_ptr = find_unit(new_unit); unit_ptr)
        {
            unit_ptr->add_frequencies(freqs);
            return;
        }

        new_unit.set_local_region(grid.get_grid_p1(), grid.get_grid_p2());
        new_unit.add_frequencies(freqs);
        if (is_e_point(ctype))
            e_units.push_back(std::move(new_unit));
        else
            m_units.push_back(std::move(new_unit));
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

        //metadata creation
        auto x_dataset = group.createDataSet<double>("x", HighFive::DataSpace({dim[3]}));
        auto y_dataset = group.createDataSet<double>("y", HighFive::DataSpace({dim[2]}));
        auto z_dataset = group.createDataSet<double>("z", HighFive::DataSpace({dim[1]}));
        auto f_dataset = group.createDataSet<double>("f", HighFive::DataSpace({dim[0]}));

        //independent write to different slabs of dataset, degenerate case considered
        if (unit_ptr->get_local_size())
        {
            auto tmp = unit_ptr->select_local_real(freqs);
            re_dataset.select(selection.first, selection.second).write<double>(tmp.data());
		
		    tmp = unit_ptr->select_local_imag(freqs);
            im_dataset.select(selection.first, selection.second).write<double>(tmp.data());    
        }

        int rank;
        MPI_Comm_rank(comm, &rank);

        if (1)
        {
            auto phys_p1 = closure.first * (dx / 2);
            auto phys_p2 = closure.second * (dx / 2);

            x_dataset.write(linspace(phys_p1.x, phys_p2.x, dim[3]).data());
            y_dataset.write(linspace(phys_p1.y, phys_p2.y, dim[2]).data());
            z_dataset.write(linspace(phys_p1.z, phys_p2.z, dim[1]).data());
            f_dataset.write(freqs.data());
        }
    }

    complex_arr DFT_Hub::get_fields_collective(const fVec3 &p1, const fVec3 &p2, Coord_Type ctype, const double_arr &freqs, MPI_Comm comm)
    {
        //yet to be implemented
        return {};
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
