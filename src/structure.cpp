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

    void Medium_Stepping::step(const std::vector<double> &accdb, std::vector<double> &eh)
    {
        const size_t stride = sus_pool.size();

        for (int i = 0; i < points.size(); ++i)
        {
            size_t index = points[i].index;
            double g_inf = points[i].g_inf;
			double *p = &polarization[i * stride];
			double *p1 = &polarization1[i * stride];

            for (int n = 0; n < stride; ++n)
            {
                auto &sus = sus_pool[n];

                p[n] = sus.c1 * p[n] + sus.c2 * p1[n] + sus.c3 * eh[index];
            }

            eh[index] = (accdb[index] * dt - std::accumulate(p, p + stride, 0)) / g_inf;
        }
    }

    //Structure
    void Structure::set_grid(const Yee3 &grid, double dx, double dt)
    {
        this->grid = grid;
        this->dx = dx;
        this->dt = dt;
    }

    void Structure::build_material_pool(const std::vector<Medium> &materials)
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
    }

    Abstract_Medium Structure::get_abstract_medium(const Medium &medium) const
    {
        Abstract_Medium res;
        res.epsilon = medium.epsilon;
        res.mu = medium.mu;

        auto build_sus_amp = [](const Susceptibility &sus, double amp, const std::vector<Susceptibility> &sus_pool, std::valarray<double> &sus_amp) {
            sus_amp.resize(sus_pool.size());
            for (int j = 0; j < sus_pool.size(); ++j)
                if (sus_pool[j] == sus)
                {
                    sus_amp[j] = amp;
                    break;
                }
        };

        for (int i = 0; i < medium.e_sus.size(); ++i)
            build_sus_amp(medium.e_sus[i], medium.e_sus_amp[i], e_sus_pool, res.e_sus_amp);

        for (int i = 0; i < medium.m_sus.size(); ++i)
            build_sus_amp(medium.m_sus[i], medium.m_sus_amp[i], m_sus_pool, res.m_sus_amp);

        return res;
    }

    void Structure::set_materials_from_geometry(const std::vector<std::reference_wrapper<Geometry>> &geom_list, const Medium &default_medium, Average_Method method)
    {
		auto default_ab_medium = get_abstract_medium(default_medium);

        switch (method)
        {
        case Line_Average:
        case Volume_Average:
        case No_Average:
        {
            auto p1 = grid.get_grid_p1();
            auto p2 = grid.get_grid_p2();

            for (auto itr = Yee_Iterator(p1, p2); !itr.is_end(); itr.next())
            {
                auto grid_p = itr.get_coord();
                fVec3 phys_p{grid_p};
                
                //return if it is not a material point
                if (!is_eh_point(grid_p.get_type()))
                    continue;

                bool found = 0;
                for (auto item : geom_list)
                {
                    if (auto &geom = item.get(); geom.is_inside(phys_p))
                    {
						found = 1;
                        auto medium = geom.get_medium(phys_p);
                        //electrical point
                        if (is_e_point(grid_p.get_type()))
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
                            e_stepping[nonzeros].add_point(grid.get_index_from_coord(grid_p), medium.epsilon, medium.e_sus_amp);
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
                            m_stepping[nonzeros].add_point(grid.get_index_from_coord(grid_p), medium.mu, medium.m_sus_amp);
                        }

                        break;
                    }
                }
				
				//use default medium
				if (!found)
				{
					auto& medium = default_ab_medium;
					//electrical point
					if (is_e_point(grid_p.get_type()))
					{
						//create a new Medium_Stepping corresponding to the pattern
						if (auto itr = e_stepping.find(0); itr == e_stepping.end())
						{
							e_stepping[0] = Medium_Stepping();
						}
						
						//add the point to the stepping
						e_stepping[0].add_point(grid.get_index_from_coord(grid_p), medium.epsilon, medium.e_sus_amp);
					}
					else
					{
						//create a new Medium_Stepping corresponding to the pattern
						if (auto itr = m_stepping.find(0); itr == m_stepping.end())
						{
							m_stepping[0] = Medium_Stepping();
						}
						
						//add the point to the stepping
						m_stepping[0].add_point(grid.get_index_from_coord(grid_p), medium.mu, medium.m_sus_amp);
					}
				}
            }
        }
        }
		
		for(auto& item : e_stepping)
			item.second.set_dt(dt);
		
		for(auto& item : m_stepping)
			item.second.set_dt(dt);
		
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
        for (int i = 0; mask >= 0; ++i)
        {
            if (mask & 1)
            {
                res.push_back(ab_sus_pool[i]);
            }

            mask >>= 1;
        }

        return res;
    }

    void Structure::step_e(const std::vector<double> &accd, std::vector<double> &e, std::vector<double> &e1)
    {
        for (auto &item : e_stepping)
        {
            item.second.step(accd, e);
        }
    }

    void Structure::step_m(const std::vector<double> &accb, std::vector<double> &h, std::vector<double> &h1)
    {
        for (auto &item : m_stepping)
        {
            item.second.step(accb, h);
        }
    }

} // namespace ffip
