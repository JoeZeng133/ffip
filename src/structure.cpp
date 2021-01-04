#include <structure.hpp>

namespace ffip
{
    //Medium_Stepping
    void Medium_Stepping::set_dt(double dt)
    {
        this->dt = dt;
    }

    void Medium_Stepping::set_susceptibility_pool(const std::vector<Abstract_Susceptibility> &sus_pool)
    {
        this->sus_pool = sus_pool;
    }

    void Medium_Stepping::set_susceptibility_pool(const std::vector<Susceptibility> &sus_pool)
    {
        for (auto &sus : sus_pool)
        {
            this->sus_pool.push_back({sus, dt});
        }
    }

    void Medium_Stepping::add_point(size_t index, double g_inf, const std::valarray<double> &sus_amp)
    {
        if (sus_amp.size() != sus_pool.size())
            throw std::runtime_error("doesn't match susceptibility pool");

        for (auto x : sus_amp)
            this->sus_amp.push_back(x);

        points.push_back({index, g_inf});
    }

    void Medium_Stepping::organize() {}


    void Medium_Stepping::output_details(std::ostream& os, const Yee3& grid) const
    {
        os << "Reporting susceptibility pool\n";
        os << "b0,b1,b2,a1,a2\n";
        for(auto& sus : sus_pool)
            os << sus.b0
            << "," << sus.b1
            << "," << sus.b2
            << "," << sus.a1
            << "," << sus.a2 << "\n";

        os << "Reporting each point\n";
        for (int i = 0; i < points.size(); ++i)
        {
            os << "coord=" << grid.get_coord_from_index(points[i].index) << ",g_inf=" << points[i].g_inf;
            os << ",sigmas=[";
            for(int j = 0; j <  sus_pool.size(); ++j)
                os << "," << sus_amp[i * sus_pool.size() + j];
            os << "]\n";
        }
            
    }

    void Medium_Stepping::step(
        const std::vector<double> &accdb,
        std::vector<double> &eh,
        std::vector<double> &eh_prev
        )
    {
        const size_t stride = sus_pool.size();
        std::vector<double> q(stride);
        
        if (polarization1.size() == 0)
        {
            polarization1.resize(stride * points.size(), 0);
            polarization.resize(stride * points.size(), 0);
        }

        for (int n = 0; n < points.size(); ++n)
        {
            size_t index = points[n].index;
            double g_inf = points[n].g_inf;
            double sum_b0 = 0;
            double sum_q = 0;
            double eh_1 = eh[index];
            double eh_2 = eh_prev[index];

			double *p_1 = &polarization[n * stride];
			double *p_2 = &polarization1[n * stride];
			double *amp = &sus_amp[n * stride];

            //step 1
            for (int i = 0; i < stride; ++i)
            {
                auto &sus = sus_pool[i];
                q[i] = sus.b1 * eh_1 + sus.b2 * eh_2
                    - sus.a1 * p_1[i] - sus.a2 * p_2[i];
                sum_q += sus_amp[i] * q[i];
                sum_b0 += sus_amp[i] * sus.b0;
            }

            //step 2
            double eh_0 = (accdb[index] * dt - sum_q) / (g_inf + sum_b0);

            //step 3
            for (int i = 0; i < stride; ++i)
            {
                auto &sus = sus_pool[i];
                p_2[i] = p_1[i];
                p_1[i] = sus.b0 * eh_0 + q[i];
            }

            eh[index] = eh_0;
            eh_prev[index] = eh_1;
        }
    }

    //Structure
    void Structure::set_grid(const Yee3 &grid, double dx, double dt)
    {
        this->grid = grid;
        this->dx = dx;
        this->dt = dt;
        this->eh_prev.resize(grid.get_size(), 0);
    }

    void Structure::add_to_material_pool(const std::vector<Medium> &materials)
    {
        auto add_sus_to_pool = [](const Susceptibility &sus, std::vector<Susceptibility> &sus_pool) {
            bool found = 0;
            for (auto &pool_sus : sus_pool)
                if (pool_sus == sus)
                {
                    found = 1;
                    break;
                }

            if (!found)
                sus_pool.push_back(sus);
        };
		
        for (auto &m : materials)
        {
            for (auto &e_sus : m.e_sus)
				add_sus_to_pool(e_sus, e_sus_pool);

            for (auto &m_sus : m.m_sus)
                add_sus_to_pool(m_sus, m_sus_pool);
        }
		
		for(int i = e_ab_sus_pool.size(); i < e_sus_pool.size(); ++i)
			e_ab_sus_pool.push_back({e_sus_pool[i], dt});
		
		for(int i = m_ab_sus_pool.size(); i < m_sus_pool.size(); ++i)
			m_ab_sus_pool.push_back({m_sus_pool[i], dt});
    }

    Abstract_Medium Structure::get_abstract_medium(const Medium &medium) const
    {
        Abstract_Medium res;
        res.epsilon = medium.epsilon;
        res.mu = medium.mu;

        auto build_sus_amp = [](const Susceptibility &sus, double amp, const std::vector<Susceptibility> &sus_pool, std::valarray<double> &sus_amp)
		{
            for (int j = 0; j < sus_pool.size(); ++j)
                if (sus_pool[j] == sus)
                {
                    sus_amp[j] = amp;
                    break;
                }
        };
		
		//resizing valarray reassigns the values
		res.e_sus_amp.resize(e_sus_pool.size());
		res.m_sus_amp.resize(m_sus_pool.size());
		
        for (int i = 0; i < medium.e_sus.size(); ++i)
            build_sus_amp(medium.e_sus[i], medium.e_sus_amp[i], e_sus_pool, res.e_sus_amp);

        for (int i = 0; i < medium.m_sus.size(); ++i)
            build_sus_amp(medium.m_sus[i], medium.m_sus_amp[i], m_sus_pool, res.m_sus_amp);

        return res;
    }

    void Structure::set_materials_from_geometry(const std::vector<std::reference_wrapper<Geometry>> &geom_list, const Medium &default_medium, Average_Method method)
    {
		auto default_ab_medium = get_abstract_medium(default_medium);
		
		auto set_medium = [&](const iVec3 p, const Abstract_Medium& medium) {
			if (is_e_point(p.get_type()))
			{
				//get nonzeros pattern
				size_t nonzeros = get_non_zeros_from_array(medium.e_sus_amp, 1e-4);
				
				//create a new Medium_Stepping corresponding to the pattern
				if (auto itr = e_stepping.find(nonzeros); itr == e_stepping.end())
				{
					e_stepping[nonzeros] = Medium_Stepping();
					e_stepping[nonzeros].set_susceptibility_pool(mask_susceptibility_pool(nonzeros, e_ab_sus_pool));
				}
				
				//add the point to the stepping
				e_stepping[nonzeros].add_point(grid.get_index_from_coord(p), medium.epsilon,
											   mask_susceptibility_amp(nonzeros, medium.e_sus_amp));
			}
			else
			{
				//get nonzeros pattern
				size_t nonzeros = get_non_zeros_from_array(medium.m_sus_amp, 1e-4);
				
				//create a new Medium_Stepping corresponding to the pattern
				if (auto itr = m_stepping.find(nonzeros); itr == m_stepping.end())
				{
					m_stepping[nonzeros] = Medium_Stepping();
					m_stepping[nonzeros].set_susceptibility_pool(mask_susceptibility_pool(nonzeros, m_ab_sus_pool));
				}
				
				//add the point to the stepping
				m_stepping[nonzeros].add_point(grid.get_index_from_coord(p), medium.mu, mask_susceptibility_amp(nonzeros, medium.m_sus_amp));
			}
		};

        switch (method)
        {
        case Line_Average:
        case Volume_Average:
        case No_Average:
        {
            auto p1 = grid.get_grid_p1();
            auto p2 = grid.get_grid_p2();

            geom_map.resize(Yee_Iterator::get_size(p1, p2), -1);
            size_t index = 0;

            for (auto itr = Yee_Iterator(p1, p2); !itr.is_end(); itr.next(), ++index)
            {
                auto grid_p = itr.get_coord();
                fVec3 phys_p{grid_p};
                
                //return if it is not a material point
                Coord_Type ctype = grid_p.get_type();
                //first 4 bits represent coord ctype
                // geom_map[index] = ctype;

                if (!is_eh_point(ctype))
                    continue;

                bool found = 0;
                int num = 0;
                for (auto item : geom_list)
                {
                    if (auto &geom = item.get(); geom.is_inside(phys_p))
                    {
						found = 1;
                        auto medium = geom.get_medium(phys_p);
						set_medium(grid_p, medium);
                        //rest of the bits represent geometry number
                        geom_map[index] = num;
                        break;
                    }
                    num++;
                }
				
				//use default medium
				if (!found)
				{
					set_medium(grid_p, default_ab_medium);
				}
            }
        }
        }
		
		for(auto& item : e_stepping)
			item.second.set_dt(dt);
		
		for(auto& item : m_stepping)
			item.second.set_dt(dt);
		
    }

    const std::vector<int>& Structure::get_geom_map() const
    {
        return geom_map;
    }

    size_t Structure::get_non_zeros_from_array(const std::valarray<double> &arr, double tol) const
    {
        size_t res = 0;
        for (int i = 0; i < arr.size(); ++i)
            if (std::abs(arr[i]) > tol)
                res |= (1 << i);

        return res;
    }

    std::vector<Abstract_Susceptibility>
    Structure::mask_susceptibility_pool(size_t mask, const std::vector<Abstract_Susceptibility> &ab_sus_pool) const
    {
        std::vector<Abstract_Susceptibility> res;
        for (int i = 0; mask > 0; ++i)
        {
            if (mask & 1)
            {
                res.push_back(ab_sus_pool[i]);
            }

            mask >>= 1;
        }

        return res;
    }
	
	std::valarray<double>
	Structure::mask_susceptibility_amp(size_t mask, const std::valarray<double> &arr) const
	{
		std::vector<double> res;
		
		for (int i = 0; mask > 0; ++i)
		{
			if (mask & 1)
			{
				res.push_back(arr[i]);
			}
			
			mask >>= 1;
		}
		
		return {res.data(), res.size()};
	}

    void Structure::step_e(const std::vector<double> &accd, std::vector<double> &e)
    {
        for (auto &item : e_stepping)
        {
            item.second.step(accd, e, eh_prev);
        }
    }

    void Structure::step_m(const std::vector<double> &accb, std::vector<double> &h)
    {
        for (auto &item : m_stepping)
        {
            item.second.step(accb, h, eh_prev);
        }
    }

    void Structure::output_details(std::ostream& os) const
    {
        
        size_t index=0;
        for(auto& item : e_stepping)
        {
            os << "Reporting the " << index++ << "th e medium\n";
            item.second.output_details(os, grid);
        }

        index = 0;
        for(auto& item : m_stepping)
        {
            os << "Reporting the " << index++ << "th m medium\n";
            item.second.output_details(os, grid);
        }
    }

} // namespace ffip
