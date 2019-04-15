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

    void Gaussian_Dipoles_Stepping::output(std::ostream &os) {}

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
            long stride1 = strides[(i + 1) % 3];
            long stride2 = strides[(i + 2) % 3];

            for (auto point : points[i])
            {
                double curl1 = (eh[point + stride1] - eh[point - stride1]);
                double curl2 = (eh[point + stride2] - eh[point - stride2]);

                accdb[point] += (curl1 - curl2) * inv_dx;
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
		
		e_pml.set_dx(dx);
		e_curl.set_dx(dx);
		
		m_curl.set_dx(dx);
		m_pml.set_dx(dx);
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
        return get_db_helper(pos, ctype);
    }

    void Fields::step_accd(MPI_Comm comm, double time)
    {
        e_pml.step(eh, accdb);
        e_curl.step(eh, accdb);
        e_gaussian_dipoles.step(time, accdb);

        sync_boundary(comm, X);
        sync_boundary(comm, Y);
        sync_boundary(comm, Z);
    }

    void Fields::step_accb(MPI_Comm comm, double time)
    {
        m_pml.step(eh, accdb);
        m_curl.step(eh, accdb);
        m_gaussian_dipoles.step(time, accdb);

        sync_boundary(comm, X);
        sync_boundary(comm, Y);
        sync_boundary(comm, Z);
    }

    void Fields::sync_boundary(MPI_Comm comm, Direction dir)
    {

        //send, receive for the positive face
        static std::vector<double> pos_send_buf, pos_receive_buf;
        //send, receive for the negative face
        static std::vector<double> neg_send_buf, neg_receive_buf;

        int rank, pos_rank, neg_rank;
        int face_size;
        int pos_send_size, pos_receive_size, neg_send_size, neg_receive_size;
        int phase;

        MPI_Comm_rank(comm, &rank);
        MPI_Cart_shift(comm, (int)dir, 1, &neg_rank, &pos_rank);

        //get face size
        iVec3 p1 = grid.get_grid_p1();
        iVec3 p2 = grid.get_grid_p2();
        auto face = get_face(p1, p2, dir, Positive);
        face_size = Yee_Iterator(face).get_size();

        //positive face
        phase = 1;
        switch (bc[dir][1])
        {
        //PEC
        case PEC:
            phase = -1;

        //PEC, PMC
        case PMC:
            switch (dir)
            {
            case X:
                sync_symmetry_boundary<X, Positive>(phase);
                break;

            case Y:
                sync_symmetry_boundary<Y, Positive>(phase);
                break;

            case Z:
                sync_symmetry_boundary<Z, Positive>(phase);
                break;
            }

        //PEC, PMC, None
        case None:
            pos_send_size = pos_receive_size = 0;
            pos_send_buf.resize(0);
            pos_receive_buf.resize(0);
            break;

        case Sync:
            pos_send_size = pos_receive_size = face_size;
            pos_receive_buf.resize(face_size);
            pos_send_buf.resize(face_size);
            switch (dir)
            {
            case X:
                extract_fields_face<X, Positive>(pos_send_buf);
                break;

            case Y:
                extract_fields_face<Y, Positive>(pos_send_buf);
                break;

            case Z:
                extract_fields_face<Z, Positive>(pos_send_buf);
                break;
            }
            break;
        }

        //negative face
        phase = 1;
        switch (bc[dir][0])
        {
        //PEC
        case PEC:
            phase = -1;

        //PEC, PMC
        case PMC:
            switch (dir)
            {
            case X:
                sync_symmetry_boundary<X, Negative>(phase);
                break;

            case Y:
                sync_symmetry_boundary<Y, Negative>(phase);
                break;

            case Z:
                sync_symmetry_boundary<Z, Negative>(phase);
                break;
            }

        //PEC, PMC, None
        case None:
            neg_send_size = neg_receive_size = 0;
            neg_send_buf.resize(0);
            neg_receive_buf.resize(0);
            break;

        case Sync:
            neg_send_size = neg_receive_size = face_size;
            neg_receive_buf.resize(face_size);
            neg_send_buf.resize(face_size);
            switch (dir)
            {
            case X:
                extract_fields_face<X, Negative>(neg_send_buf);
                break;

            case Y:
                extract_fields_face<Y, Negative>(neg_send_buf);
                break;

            case Z:
                extract_fields_face<Z, Negative>(neg_send_buf);
                break;
            }
            break;
        }

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

        //assigning fields
        switch (dir)
        {
        case X:
            assign_fields_face<X, Positive>(pos_receive_buf);
            assign_fields_face<X, Negative>(neg_receive_buf);
            break;

        case Y:
            assign_fields_face<Y, Positive>(pos_receive_buf);
            assign_fields_face<Y, Negative>(neg_receive_buf);
            break;

        case Z:
            assign_fields_face<Z, Positive>(pos_receive_buf);
            assign_fields_face<Z, Negative>(neg_receive_buf);
            break;
        }
    }

    double Fields::get_source_turnoff_time() const
    {
        return turnoff_time;
    }

    //PML class
    PML::PML(double d, double sigma_max, double k_max, int m) : d(d), sigma_max(sigma_max), k_max(k_max), m(m)
    {
		if (std::abs(sigma_max) < 1e-9)
            this->sigma_max = get_sigma_max(d, 1e-4, m, 1);
    }

    double PML::get_sigma_max(double d, double R, double m, double imp)
    {
        return -(m + 1) * std::log(R) / (imp * 2 * d);
    }

    double PML::get_b(double x, double dt)
    {
        return std::exp(-get_sigma(x) / get_k(x) * dt);
    }

    double PML::get_c(double x, double dt)
    {
        return 1 / get_k(x) * (get_b(x, dt) - 1);
    }

    double PML::get_sigma(double x)
    {
        return std::pow(x / d, m) * sigma_max;
    }

    double PML::get_k(double x)
    {
        return 1 + (k_max - 1) * std::pow(x / d, m);
    }
} // namespace ffip
