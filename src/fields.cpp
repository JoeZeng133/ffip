#include <fields.hpp>

namespace ffip
{
    //Gaussian Dipoles
    void Gaussian_Dipoles_Stepping::organize()
    {
        sort(gaussian1_points.begin(), gaussian1_points.end(), [](const Stepping_Point &a, const Stepping_Point &b) {
            return a.index < b.index;
        });

        sort(gaussian2_points.begin(), gaussian2_points.end(), [](const Stepping_Point &a, const Stepping_Point &b) {
            return a.index < b.index;
        });
    }

    void Gaussian_Dipoles_Stepping::step(double time, std::vector<double> &accdb)
    {
        if (time > max_end_time)
            return;
		
        for (auto &point : gaussian1_points)
        {
            if (time > point.start_time && time < point.end_time)
                accdb[point.index] += point.amp * Gaussian1(time - (point.start_time + point.end_time) / 2, point.width);
        }

        for (auto &point : gaussian2_points)
        {
            if (time > point.start_time && time < point.end_time)
                accdb[point.index] += point.amp * Gaussian2(time - (point.start_time + point.end_time) / 2, point.width);
        }
    }

    void Gaussian_Dipoles_Stepping::add_gaussian1(size_t index, double amp, double start_time, double end_time, double width)
    {
        gaussian1_points.push_back({index, amp, start_time, end_time, width});
        max_end_time = std::max(end_time, max_end_time);
    }

    void Gaussian_Dipoles_Stepping::add_gaussian2(size_t index, double amp, double start_time, double end_time, double width)
    {
        gaussian2_points.push_back({index, amp, start_time, end_time, width});
        max_end_time = std::max(end_time, max_end_time);
    }

    void Gaussian_Dipoles_Stepping::output_details(std::ostream &os, const Yee3& grid) const
	{
		os << "Reporting gaussian1 points\n";
		for(auto &point : gaussian1_points)
            os << "coord=" <<  grid.get_coord_from_index(point.index) << ",amp=" << point.amp << ",st=" << point.start_time << ",et=" << point.end_time << ",w=" << point.width << "\n";
		
		os << "Reporting gaussian2 points\n";
		for(auto &point : gaussian2_points)
			os << "coord=" <<  grid.get_coord_from_index(point.index) << ",amp=" << point.amp << ",st=" << point.start_time << ",et=" << point.end_time << ",w=" << point.width << "\n";
	}

    //PML Stepping
    void PML_Stepping::set_strides(sVec3 strides)
    {
        this->strides = strides;
    }

    void PML_Stepping::set_dx(double dx) 
    {
        this->inv_dx = 1 / dx;
    }

    void PML_Stepping::organize()
    {
        for (int i = 0; i < 3; ++i)
            sort(points[i].begin(), points[i].end(), [](const Stepping_Point &a, const Stepping_Point &b) {
                return a.index < b.index;
            });
    }
	
	void PML_Stepping::output_details(std::ostream& os, const Yee3& grid) const
	{
		for(int i = 0; i < 3; ++i)
		{
			os << "Reporting PML" << i << "direction points\n";
			
			for(auto &point : points[i])
				os << "coord=" << grid.get_coord_from_index(point.index) <<
				",b1=" << point.b1 << ",b2=" << point.b2 << ",invk1=" << point.inv_k1 << ",invk2=" << point.inv_k2 << "\n";
		}
	}

    void PML_Stepping::step(const std::vector<double> &eh, std::vector<double> &accdb)
    {
        for (int i = 0; i < 3; ++i)
        {
            //stride in x2
            long long stride1 = strides[(i + 1) % 3];
            //stride in x3
            long long stride2 = strides[(i + 2) % 3];

            for (auto &point : points[i])
            {
                double curl1 = (eh[point.index + stride1] - eh[point.index - stride1]) * inv_dx;
                double curl2 = (eh[point.index + stride2] - eh[point.index - stride2]) * inv_dx;

                point.phi1 = point.b1 * point.phi1 + point.c1 * curl1;
                point.phi2 = point.b2 * point.phi2 + point.c2 * curl2;

                accdb[point.index] += point.phi1 - point.phi2 + curl1 * point.inv_k1 - curl2 * point.inv_k2;
            }
        }
    }

    void PML_Stepping::add_stepping_point(Direction dir, size_t index, double b1, double c1, double k1, double b2, double c2, double k2)
    {
        points[dir].push_back({index, 0, b1, c1, 1 / k1, 0, b2, c2, 1 / k2});
    }

    //Curl_Steppping
    void Curl_Stepping::set_dx(double dx)
    {
        this->inv_dx = 1 / dx;
    }

    void Curl_Stepping::set_strides(sVec3 strides)
    {
        this->strides = strides;
    }

    void Curl_Stepping::organize()
    {
        for (int i = 0; i < 3; ++i)
            sort(points[i].begin(), points[i].end());
    }

    void Curl_Stepping::step(const std::vector<double> &eh, std::vector<double> &accdb)
    {
        for (int i = 0; i < 3; ++i)
        {
            long long stride1 = strides[(i + 1) % 3];
            long long stride2 = strides[(i + 2) % 3];

            for (auto point : points[i])
            {
                double curl1 = (eh[point + stride1] - eh[point - stride1]) * inv_dx;
                double curl2 = (eh[point + stride2] - eh[point - stride2]) * inv_dx;

                accdb[point] += (curl1 - curl2);
            }
        }
    }

    void Curl_Stepping::add_stepping_point(Direction dir, size_t index)
    {
        points[dir].push_back(index);
    }

    //Fields
	void Fields::set_boundary_conditions(Boundary_Condition _bc[3][2])
	{
		for(int i = 0; i < 3; ++i)
			for(int j = 0; j < 2; ++j)
				bc[i][j] = _bc[i][j];
	}
	
    void Fields::set_grid(const Yee3 &grid, double dx, double dt)
    {
        this->grid = grid;
        this->dx = dx;
        this->dt = dt;
		
		e_pml.set_strides(grid.get_stride());
		e_pml.set_dx(dx);
		e_curl.set_strides(grid.get_stride());
		e_curl.set_dx(dx);
		
		m_pml.set_strides(-grid.get_stride());
		m_pml.set_dx(dx);
		m_curl.set_strides(-grid.get_stride());
		m_curl.set_dx(dx);
    }

    void Fields::add_dipole_source_gaussian1(const fVec3 &pos, Coord_Type ctype, double amp, double start_time, double frequency, double cutoff)
    {
        amp = -amp;
        auto weights = grid.get_interp_weights(pos, ctype);
        auto base = grid.get_base_point(pos, ctype);

        auto &gaussian_dipoles = (is_e_point(ctype) ? e_gaussian_dipoles : m_gaussian_dipoles);
        double width = 1 / (2 * pi * frequency);
        double end_time = start_time + 2 * width * cutoff;

        if (end_time > turnoff_time)
            turnoff_time = end_time;
		
		size_t index = 0;
        for (int k = 0; k < 4; k += 2)
            for (int j = 0; j < 4; j += 2)
                for (int i = 0; i < 4; i += 2)
                {
                    auto p = base + iVec3{i, j, k};
                    //if not inside, don't care
                    if (grid.is_inside(p))
                        gaussian_dipoles.add_gaussian1(grid.get_index_from_coord(p), amp * weights[index], start_time, end_time, width);
                    ++index;
                }
    }

    void Fields::add_dipole_source_gaussian2(const fVec3 &pos, Coord_Type ctype, double amp, double start_time, double frequency, double cutoff)
    {
        amp = -amp;
        auto weights = grid.get_interp_weights(pos, ctype);
        auto base = grid.get_base_point(pos, ctype);

        auto &gaussian_dipoles = (is_e_point(ctype) ? e_gaussian_dipoles : m_gaussian_dipoles);
        double width = std::sqrt(2.0) / (2 * pi * frequency);
        double end_time = start_time + 2 * width * cutoff;

        if (end_time > turnoff_time)
            turnoff_time = end_time;

        for (int k = 0, index = 0; k < 4; k += 2)
            for (int j = 0; j < 4; j += 2)
                for (int i = 0; i < 4; i += 2)
                {
                    auto p = base + iVec3{i, j, k};
                    //if not inside, don't care
                    if (grid.is_inside(p))
                        gaussian_dipoles.add_gaussian2(grid.get_index_from_coord(p), amp * weights[index++], start_time, end_time, width);
                    ++index;
                }
    }

    void Fields::init(
        const std::array<double_arr, 3> &k, const std::array<double_arr, 3> &b, const std::array<double_arr, 3> &c,
        const iVec3 &p1, const iVec3 &p2)
    {
        eh.resize(grid.get_size());
        accdb.resize(grid.get_size());

        PML_init_helper<Ex>(k, b, c, p1, p2);
        PML_init_helper<Ey>(k, b, c, p1, p2);
        PML_init_helper<Ez>(k, b, c, p1, p2);
        PML_init_helper<Hx>(k, b, c, p1, p2);
        PML_init_helper<Hy>(k, b, c, p1, p2);
        PML_init_helper<Hz>(k, b, c, p1, p2);
    }

    double Fields::get_eh_helper(const fVec3 &pos, Coord_Type ctype) const
    {
        return grid.interp(eh, pos, ctype);
    }

    double Fields::get_db_helper(const fVec3 &pos, Coord_Type ctype) const
    {
        return grid.interp(accdb, pos, ctype) * dt;
    }

    double Fields::get_eh_raw(const iVec3 &pos) const
    {
        return grid.get_raw_val(eh, pos);
    }

    double Fields::get_db_raw(const iVec3 &pos) const
    {
        return grid.get_raw_val(accdb, pos) * dt;
    }

    double Fields::get(const fVec3& pos, Coord_Type ctype) const {
        if (is_eh_point(ctype))
            return get_eh_helper(pos, ctype);
        return get_db_helper(pos, ctype) * dt;
    }
	
	void Fields::report() const
	{
		double tmp=0;
		for(auto x : accdb)
			tmp += std::abs(x);
		std::cout << "accdb=" << tmp;
		
		tmp = 0;
		for(auto x : eh)
			tmp += std::abs(x);
		std::cout << ",eh=" << tmp << "\n";
	}

    void Fields::output_fields(std::ostream& os) const
    {
        os << "Reporting eh fields:\n";
        auto p1 = iVec3(-5, -5, -5);
        auto p2 = iVec3(5, 5, 5);

		for(auto itr = Yee_Iterator(p1, p2); !itr.is_end(); itr.next())
        {
			size_t index = grid.get_index_from_coord(itr.get_coord());
			if (eh[index] != 0 || accdb[index] != 0) {
				os << "index=" << index << ",coord=" <<  itr.get_coord() << ",eh=" << eh[index] << ",accdb=" << accdb[index] << "\n";
			}
        }
//        os << "\n";
    }

    void Fields::output_details(std::ostream& os) const
    {
        os << "e_pml\n";
        e_pml.output_details(os, grid);
        os << "m_pml\n";
        m_pml.output_details(os, grid);
        os << "e gaussian dipoles\n";
        e_gaussian_dipoles.output_details(os, grid);
        os << "m gaussian dipoles\n";
        m_gaussian_dipoles.output_details(os, grid);
    }

    void Fields::step_accd(MPI_Comm comm, double time)
    {
        e_pml.step(eh, accdb);
        e_curl.step(eh, accdb);
        e_gaussian_dipoles.step(time, accdb);
    }

    void Fields::step_accb(MPI_Comm comm, double time)
    {
        m_pml.step(eh, accdb);
        m_curl.step(eh, accdb);
        m_gaussian_dipoles.step(time, accdb);
    }

    void Fields::sync_boundary(MPI_Comm comm, Direction dir)
    {
        //send, receive for the positive face
        static std::vector<double> pos_send_buf, pos_receive_buf;
        //send, receive for the negative face
        static std::vector<double> neg_send_buf, neg_receive_buf;

        int rank, pos_rank, neg_rank;
        size_t face_size;
        size_t pos_send_size, pos_receive_size, neg_send_size, neg_receive_size;
        Boundary_Condition pos_bc = bc[dir][1], neg_bc = bc[dir][0];

        MPI_Comm_rank(comm, &rank);
        MPI_Cart_shift(comm, (int)dir, 1, &neg_rank, &pos_rank);

        //get face size
        iVec3 p1 = grid.get_grid_p1();
        iVec3 p2 = grid.get_grid_p2();
        auto face = get_face(p1, p2, dir, Positive);
        face_size = Yee_Iterator(face).get_size();
		
		//resizing buffer when necessary
		if (face_size > pos_send_buf.size()) {
			pos_send_buf.resize(face_size);
			neg_send_buf.resize(face_size);
			pos_receive_buf.resize(face_size);
			neg_receive_buf.resize(face_size);
		}

        //extract fields depending on particular boundary conditiosn
        auto extraction = [&]
        (Boundary_Condition bc, Direction dir, Side side, std::vector<double>& buf) -> size_t
        {
            switch(bc)
            {
            case None:
                return 0;
                break;

            case PEC:
                extract_symmetry_face(dir, side, -1, buf);
                return 0;
                break;

            case PMC:
                extract_symmetry_face(dir, side, 1, buf);
                return 0;
                break;

            case Period:
                return extract_fields_face(dir, side, -1, buf);
                break;

            case Sync:
                return extract_fields_face(dir, side, 0, buf);
                break;
            }
        };

        //assigning fields
        auto assignment = [&]
        (Boundary_Condition bc, Direction dir, Side side, const std::vector<double>& send_buf, const std::vector<double>& receive_buf)
        {
            switch(bc)
            {
                case None:
					break;
					
                case PMC:
                case PEC:
					assign_fields_face(dir, side, 1, send_buf);
					break;
					
                case Sync:
                case Period:
					assign_fields_face(dir, side, 1, receive_buf);
					break;
            }
        };

        neg_receive_size = neg_send_size = extraction(neg_bc, dir, Negative, neg_send_buf);
        pos_receive_size = pos_send_size = extraction(pos_bc, dir, Positive, pos_send_buf);

        MPI_Status status;

        //send in positive direction
        MPI_Sendrecv(
            pos_send_buf.data(), pos_send_size, MPI_DOUBLE, pos_rank, dir,
            neg_receive_buf.data(), neg_receive_size, MPI_DOUBLE, neg_rank, dir,
            comm, &status);

        //send in negative direction
        MPI_Sendrecv(
            neg_send_buf.data(), neg_send_size, MPI_DOUBLE, neg_rank, dir,
            pos_receive_buf.data(), pos_receive_size, MPI_DOUBLE, pos_rank, dir,
            comm, &status);

        assignment(neg_bc, dir, Negative, neg_send_buf, neg_receive_buf);
        assignment(pos_bc, dir, Positive, pos_send_buf, pos_receive_buf);
    }

    double Fields::get_source_turnoff_time() const
    {
        return turnoff_time;
    }

    size_t Fields::extract_fields_face(Direction x3, Side side, int offset, std::vector<double> &buf)
    {
        iVec3 p1 = grid.get_grid_p1();
        iVec3 p2 = grid.get_grid_p2();

        auto norm = get_norm_vec(x3, side);
        auto face = get_face(p1 + norm * offset, p2 + norm * offset, x3, side);
        auto itr = Yee_Iterator(face);

        for (int index = 0; !itr.is_end(); itr.next(), ++index)
        {
            buf[index] = grid.get_raw_val(eh, itr.get_coord());
        }

        return itr.get_size();
    }

    void Fields::assign_fields_face(Direction x3, Side side, int offset, const std::vector<double> &buf)
    {
        iVec3 p1 = grid.get_grid_p1();
        iVec3 p2 = grid.get_grid_p2();

        auto norm = get_norm_vec(x3, side);
        auto face = get_face(p1 + norm * offset, p2 + norm * offset, x3, side);
        auto itr = Yee_Iterator(face);

        for (int index = 0; !itr.is_end(); itr.next(), ++index)
        {
            eh[grid.get_index_from_coord(itr.get_coord())] = buf[index];
        }
    }

    size_t Fields::extract_symmetry_face(Direction x3, Side side, int phase, std::vector<double> &buf)
    {
        iVec3 p1 = grid.get_grid_p1();
		iVec3 p2 = grid.get_grid_p2();
		int pos;
		if (side == Negative)
			pos = p1[x3];
		else
			pos = p2[x3];
		
        //orientation (x1, x2, x3)
        int mult = (pos & 1) ? phase : -phase;
        size_t len = extract_fields_face(x3, side, -1, buf);
        for(size_t i = 0; i < len; ++i)
            buf[i] *= mult;
        
        return len;
    }

    //PML class
    PML::PML(double d, double sigma_max, double k_max, int m) : d(d), sigma_max(sigma_max), k_max(k_max), m(m)
    {}

    double PML::get_b(double x, double dt) const
    {
        return std::exp(-get_sigma(x) / get_k(x) * dt);
    }

    double PML::get_c(double x, double dt) const
    {
        return 1 / get_k(x) * (get_b(x, dt) - 1);
    }

    double PML::get_sigma(double x) const
    {
        return std::pow(x / d, m) * sigma_max;
    }

    double PML::get_k(double x) const
    {
        return 1 + (k_max - 1) * std::pow(x / d, m);
    }
} // namespace ffip
