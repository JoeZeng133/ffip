#include <simulation.hpp>

namespace ffip
{
    using namespace HighFive;

    std::map<Boundary_Condition, std::string> bc_map = {{PEC, "PEC"}, {PMC, "PMC"}, {None, "None"}, {Sync, "Sync"}, {Period, "Period"}};
    std::map<Direction, std::string> dir_map = {{X, "x"}, {Y, "y"}, {Z, "z"}};
    std::map<Side, std::string> side_map = {{Positive, "+"}, {Negative, "-"}};

    void wait_until_attach()
    {
        int i = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while(i==0) sleep(5);
        std::cout << "Process Attached\n";
    }

    Side json2side(const json &j)
    {
        static std::unordered_map<std::string, Side> map = {
            {"positive", Positive}, {"negative", Negative}};

        return map.at(j.get<std::string>());
    }

    Coord_Type json2ctype(const json &j)
    {
        static std::unordered_map<std::string, Coord_Type> map = {
            {"Ex", Ex}, {"Ey", Ey}, {"Ez", Ez}, {"Hx", Hx}, {"Hy", Hy}, {"Hz", Hz}, {"Dx", Dx}, {"Dy", Dy}, {"Dz", Dz}, {"Bx", Bx}, {"By", By}, {"Bz", Bz}};

        std::string str = j.get<std::string>();
        return map.at(str);
    }

    Direction json2direction(const json &j)
    {
        static std::unordered_map<std::string, Direction> map = {
            {"x", Direction::X}, {"y", Direction::Y}, {"z", Direction::Z}};

        std::string str = j.get<std::string>();
        return map.at(str);
    }

    fVec3 json2fvec3(const json &j)
    {
        return {j.at(0).get<double>(), j.at(1).get<double>(), j.at(2).get<double>()};
    }

    iVec3 json2ivec3(const json &j)
    {
        return {j.at(0).get<int>(), j.at(1).get<int>(), j.at(2).get<int>()};
    }

    Medium json2medium(const json &medium_json)
    {
        Medium res(medium_json.at("epsilon").get<double>(), medium_json.at("mu").get<double>());

        if (auto e_cond = medium_json.at("electric conductivity").get<double>(); e_cond != 0)
            res.add_e_susceptibility(make_conductivity_susceptibility(), e_cond);

        if (auto m_cond = medium_json.at("magnetic conductivity").get<double>(); m_cond != 0)
            res.add_m_susceptibility(make_conductivity_susceptibility(), m_cond);

        if (auto e_sus = medium_json.find("electric susceptibility"); e_sus != medium_json.end())
            for (auto itr = e_sus->begin(); itr != e_sus->end(); itr++)
            {
                res.add_e_susceptibility(json2susceptibility(*itr), itr->at("amplitude").get<double>());
            }

        if (auto m_sus = medium_json.find("magnetic susceptibility"); m_sus != medium_json.end())
            for (auto itr = m_sus->begin(); itr != m_sus->end(); itr++)
            {
                res.add_m_susceptibility(json2susceptibility(*itr), itr->at("amplitude").get<double>());
            }

        return res;
    }

    double_arr json2double_arr(const json &j)
    {
        return double_arr(j.begin(), j.end());
    }

    Susceptibility json2susceptibility(const json &sus_json)
    {
        auto type_str = sus_json.at("type").get<std::string>();

        if (type_str == "Lorentz")
            return make_Lorentz_susceptibility(sus_json.at("frequency").get<double>(), sus_json.at("gamma").get<double>());
        else if (type_str == "Drude")
            return make_Drude_susceptibility(sus_json.at("frequency").get<double>(), sus_json.at("gamma").get<double>());
        else if (type_str == "Deybe")
            return make_Deybe_susceptibility(sus_json.at("tau").get<double>());

        throw std::runtime_error("Unknonw susceptibility");
    }

    //Volume_Fields_DFT
    Volume_Fields_DFT::Volume_Fields_DFT(const json &config, DFT_Hub &hub)
    {
        double dx = hub.get_dx();
        auto center = json2fvec3(config.at("center")) / (dx / 2);
        auto size = json2fvec3(config.at("size")) / (dx / 2);
        freqs = json2double_arr(config.at("frequency"));
        ctype = json2ctype(config.at("field component"));
        config.at("output group").get_to(group_name);

        p1 = center - size / 2;
        p2 = center + size / 2;

        // wait_until_attach();

        hub.register_volume_dft(p1, p2, ctype, freqs);
    }

    void Volume_Fields_DFT::output(MPI_Comm comm, DFT_Hub &hub, HighFive::File &file)
    {
        auto group = file.createGroup(group_name);
        hub.output_fields_collective(p1, p2, ctype, freqs, comm, group);
    }

    //Simulation
    void Simulation::init(const json &config)
    {
        double sc;

        
        //dx dt
        if (auto itr = config.find("courant number"); itr != config.end())
            itr->get_to(sc);
        else
            sc = 0.5;

        config.at("resolution").get_to(resolution);
        dx = 1 / resolution;
        dt = dx * sc;
        dx2 = dx / 2;

        set_size(config);
        read_boundary_conditions(config);
        set_grid();

        input_file = new File(config.at("input file").get<std::string>(),
                            File::ReadOnly | File::Create,
                            MPIOFileDriver(cart_comm, MPI_INFO_NULL));

        fields_output_file = new File(config.at("fields output file").get<std::string>(),
                                    File::ReadWrite | File::Create | File::Truncate,
                                    MPIOFileDriver(cart_comm, MPI_INFO_NULL));

        if (auto itr = config.find("default material"); itr != config.end())
            bg_medium = json2medium(*itr);
        else
            bg_medium = Medium(1, 1);

        //add background medium as well
        structure.add_to_material_pool({bg_medium});

        if (auto itr = config.find("PML"); itr != config.end())
            set_pmls(*itr);
		else
			set_pmls(json{});

        if (auto itr = config.find("geometry"); itr != config.end())
            set_geometry(*itr);
        else
            set_geometry(json{});

        if (auto itr = config.find("geometry info request"); itr != config.end())
            geometry_info_request(*itr);

        if (auto itr = config.find("sources"); itr != config.end())
            set_source(*itr);

        if (auto itr = config.find("fields output"); itr != config.end())
            set_fields_output(*itr);

        if (auto itr = config.find("progress interval"); itr != config.end())
            itr->get_to(progress_interval);
    }

    void Simulation::geometry_info_request(const json &requests)
    {
        for(auto itr = requests.begin(); itr != requests.end(); ++itr)
        {
            std::string group_name = itr->at("output group").get<std::string>();
            Coord_Type ctype = json2ctype(itr->at("field component"));
            int geom_index = itr->at("geometry index").get<int>();
            auto&geom_map = structure.get_geom_map();
            auto strides = sVec3((grid_p2 - grid_p1 + 1).to_strides());

            std::vector<int> x, y, z;

            for(auto itr = Yee_Iterator(grid_p1, grid_p2, ctype); !itr.is_end(); itr.next())
            {
                auto coord = itr.get_coord();
                size_t index = strides.dot(coord - grid_p1);

                if (geom_map[index] == geom_index)
                {
                    x.push_back(coord.x*dx2);
                    y.push_back(coord.y*dx2);
                    z.push_back(coord.z*dx2);
                }
            }

            int total_size, size, end;
            size = x.size();

            MPI_Allreduce(&size, &total_size, 1, MPI_INT, MPI_SUM, cart_comm);
            MPI_Scan(&size, &end, 1, MPI_INT, MPI_SUM, cart_comm);

            auto group = fields_output_file->createGroup(group_name);
            auto x_dataset = group.createDataSet<double>("x", HighFive::DataSpace({size_t(total_size)}));
            auto y_dataset = group.createDataSet<double>("y", HighFive::DataSpace({size_t(total_size)}));
            auto z_dataset = group.createDataSet<double>("z", HighFive::DataSpace({size_t(total_size)}));

            if (size > 0)
            {
                x_dataset.select({size_t(end - size)}, {size_t(size)}).write(x);
                y_dataset.select({size_t(end - size)}, {size_t(size)}).write(y);
                z_dataset.select({size_t(end - size)}, {size_t(size)}).write(z);
            }
        }
    }

    void Simulation::set_volume_source(const json &src)
    {
        //read information
        auto center = json2fvec3(src.at("center")) / dx2;
        auto size = json2fvec3(src.at("size")) / dx2;
        double amp = src.at("amplitude").get<double>();
        auto ctype = json2ctype(src.at("field component"));
        auto func_json = src.at("function");
        auto func_type_str = func_json.at("type").get<std::string>();
        auto freq = func_json.at("frequency").get<double>();
        auto cutoff = func_json.at("cutoff").get<double>();
        auto stime = func_json.at("start time").get<double>();

        auto p1 = center - size / 2;
        iVec3 count = (size / 2).round();
        fVec3 dlen;

        //rescale amplitude, reduction of dimension is considered
        for (int i = 0; i < 3; ++i)
        {
            if (count[i] != 0)
            {
                dlen[i] = size[i] / count[i];
                amp *= dlen[i] / 2;
            }
            else
                dlen[i] = 0;
        }

        if (func_type_str == "Gaussian1")
        {
            for (auto itr = Yee_Iterator({0, 0, 0}, count); !itr.is_end(); itr.next())
            {
                auto p = p1 + itr.get_coord() * dlen;
                fields.add_dipole_source_gaussian1(p, ctype, amp, stime, freq, cutoff);
            }
        }
        else if (func_type_str == "Gaussian2")
        {
            for (auto itr = Yee_Iterator({0, 0, 0}, count); !itr.is_end(); itr.next())
            {
                auto p = p1 + itr.get_coord() * dlen;
                fields.add_dipole_source_gaussian2(p, ctype, amp, stime, freq, cutoff);
            }
        }
        else
            std::cout << "Source function type not supported\n";
    }

    void Simulation::set_inhom_source(const json &src)
    {
        //read information
        auto center = json2fvec3(src.at("center")) / dx2;
        auto size = json2fvec3(src.at("size")) / dx2;
        auto dim = json2ivec3(src.at("dimension"));
        auto amp_dataset_name = src.at("amplitude dataset").get<std::string>();
        auto phase_dataset_name = src.at("phase dataset").get<std::string>();
        auto ctype = json2ctype(src.at("field component"));
        auto fcen = src.at("frequency").get<double>();

        auto func_json = src.at("function");
        auto func_type_str = func_json.at("type").get<std::string>();
        auto freq = func_json.at("frequency").get<double>();
        auto cutoff = func_json.at("cutoff").get<double>();
        auto stime = func_json.at("start time").get<double>();

        auto amp_dataset = input_file->getDataSet(amp_dataset_name);
        auto phase_dataset = input_file->getDataSet(phase_dataset_name);

        std::vector<double> amp, phase;
        amp_dataset.read(amp);
        phase_dataset.read(phase);

        auto p1 = center - size / 2;
        fVec3 dlen{};
        if (dim.x > 1) dlen.x = size.x / (dim.x - 1);
        if (dim.y > 1) dlen.y = size.y / (dim.y - 1);
        if (dim.z > 1) dlen.z = size.z / (dim.z - 1);

        if (func_type_str == "Gaussian1")
        {
            size_t index = 0;
            for (auto itr = Yee_Iterator({0, 0, 0}, dim - 1); !itr.is_end(); itr.next())
            {
                auto p = p1 + itr.get_coord() * dlen;
                fields.add_dipole_source_gaussian1(p, ctype, amp[index], stime - phase[index] / (2 * pi * fcen), freq, cutoff);
                ++index;
            }
        }
        else if (func_type_str == "Gaussian2")
        {
            size_t index = 0;
            for (auto itr = Yee_Iterator({0, 0, 0}, dim - 1); !itr.is_end(); itr.next())
            {
                auto p = p1 + itr.get_coord() * dlen;
                fields.add_dipole_source_gaussian2(p, ctype, amp[index], stime - phase[index] / (2 * pi * fcen), freq, cutoff);
                ++index;
            }
        }
        else
            std::cout << "Source function type not supported\n";
    }

    void Simulation::set_source(const json &src)
    {
        for (auto itr = src.begin(); itr != src.end(); ++itr)
        {
            auto type_str = itr->at("type").get<std::string>();
            
            if (rank == 0)
                std::cout << "setting up source " << type_str << "\n";

            if (type_str == "volume source")
                set_volume_source(*itr);

            else if (type_str == "inhomogeneous source")
                set_inhom_source(*itr);
            else
                std::cout << "Unknown souce type\n";
        }

        // std::cout << "Process " << rank << " sources initialized\n";
    }

    void Simulation::set_geometry(const json &geoms)
    {
        std::vector<Medium> materials;
        std::vector<std::reference_wrapper<Geometry>> geom_list;

        for (auto itr = geoms.begin(); itr != geoms.end(); itr++)
        {
            auto type_str = itr->at("type").get<std::string>();
            if (rank == 0)
            {
                std::cout << "Setting up geometry " << type_str << "\n";
            }

            if (type_str != "two medium box")
            {
                auto medium_json = itr->at("medium");
                materials.push_back(json2medium(medium_json));
            }
            else
            {
                auto medium1_json = itr->at("medium1");
                auto medium2_json = itr->at("medium2");

                materials.push_back(json2medium(medium1_json));
                materials.push_back(json2medium(medium2_json));
            }
        }

        structure.add_to_material_pool(materials);

        for (auto itr = geoms.begin(); itr != geoms.end(); itr++)
        {
            geom_list.push_back(build_geometry_from_json(*itr));
        }

        //for the moment, no average is used
        Structure::Average_Method method = Structure::No_Average;

        structure.set_materials_from_geometry(geom_list, bg_medium, method);

        // std::cout << "Rank " << rank << " geometry initialized\n";
    }

    std::reference_wrapper<Geometry> Simulation::build_geometry_from_json(const json &geom_json)
    {
        Geometry *res = nullptr;
        std::string type_str = geom_json.at("type").get<std::string>();

        if (type_str == "box")
        {
            auto medium = build_abstract_medium_from_json(geom_json.at("medium"));
            fVec3 center = json2fvec3(geom_json.at("center"));
            fVec3 size = json2fvec3(geom_json.at("size"));

            res = new Box(center / dx2, size / dx2, medium);
        }
        else if (type_str == "sphere")
        {
            auto medium = build_abstract_medium_from_json(geom_json.at("medium"));
            double radius = geom_json.at("radius").get<double>();
            fVec3 center = json2fvec3(geom_json.at("center"));

            res = new Sphere(center / dx2, radius / dx2, medium);
        }
        else if (type_str == "two medium box")
        {
            auto medium1 = build_abstract_medium_from_json(geom_json.at("medium1"));
            auto medium2 = build_abstract_medium_from_json(geom_json.at("medium2"));
            auto center = json2fvec3(geom_json.at("center"));
            auto size = json2fvec3(geom_json.at("size"));
            auto dim = json2ivec3(geom_json.at("dimension"));

            std::string dataset_name = geom_json.at("density dataset").get<std::string>();

            //read rho
            std::vector<double> rho;
            auto dataset = input_file->getDataSet(dataset_name);
            dataset.read(rho);

            res = new Two_Medium_Box(center / dx2, size / dx2, dim, medium1, medium2, rho);
        }

        return std::ref(*res);
    }

    Abstract_Medium Simulation::build_abstract_medium_from_json(const json &medium_json)
    {
        return structure.get_abstract_medium(json2medium(medium_json));
    }

    void Simulation::set_grid()
    {
        MPI_Comm_size(MPI_COMM_WORLD, &np);
        sim_dim = sim_p2 - sim_p1 + 1;

        auto num = decompose_domain(sim_dim, np);
        int periods[3] = {bcs[0][0] == Period, bcs[1][0] == Period, bcs[2][0] == Period};

        MPI_Cart_create(MPI_COMM_WORLD, 3, num.to_array().data(), periods, 1, &cart_comm);
        MPI_Comm_rank(cart_comm, &rank);
        MPI_Cart_coords(cart_comm, rank, 3, coords);

        auto ch = get_chunk_from_coords(sim_dim, num, {coords[0], coords[1], coords[2]});

        grid_p1 = sim_p1 + ch.first;
        grid_p2 = sim_p1 + ch.second - 1;

        grid = Yee3(grid_p1 - 1, grid_p2 + 1);

        structure.set_grid(grid, dx, dt);
        fields.set_grid(grid, dx, dt);
        dft_hub.set_grid(grid, dx, dt);

        //set up boundary conditions for the current chunk
        for (unsigned int dir = 0; dir < 3; ++dir)
        {
            if (coords[dir] != 0)
                bcs[dir][0] = Sync;
            if (coords[dir] != num[dir] - 1)
                bcs[dir][1] = Sync;
        }

        if (rank == 0)
            std::cout << "decomposition=" << num << "\n";

        // std::cout << "Process " << rank << " :the grid spans from " << grid_p1 << " to " << grid_p2 << "\n";
        // std::cout << "Process " << rank << " Direction x boundaries:" << bc_map[bcs[0][0]] << "," << bc_map[bcs[0][1]] << "\n";
        // std::cout << "Process " << rank << " Direction y boundaries:" << bc_map[bcs[1][0]] << "," << bc_map[bcs[1][1]] << "\n";
        // std::cout << "Process " << rank << " Direction z boundaries:" << bc_map[bcs[2][0]] << "," << bc_map[bcs[2][1]] << "\n";
        fields.set_boundary_conditions(bcs);
    }

    void Simulation::set_size(const json &config)
    {
        json size_config = config["size"];

        size_config.at(0).get_to(sim_size.x);
        size_config.at(1).get_to(sim_size.y);
        size_config.at(2).get_to(sim_size.z);

        sim_span = (sim_size / dx).round() * 2;
        sim_p1 = -sim_span / 2;
        sim_p2 = sim_span / 2;
    }

    void Simulation::read_boundary_conditions(const json &config)
    {
        //Initialize to default boundary conditions
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 2; ++j)
                bcs[i][j] = None;

        //set up symmetry
        if (auto sym = config.find("symmetry"); sym != config.end())
        {

            for (auto itr = sym->begin(); itr != sym->end(); ++itr)
            {
                std::string dir;
                int phase;

                try
                {
                    itr->at("direction").get_to(dir);
                    itr->at("phase factor").get_to(phase);

                    if (dir == "x")
                    {
                        sim_p1.x = 0;
                        bcs[0][0] = (phase == 1) ? PMC : PEC;
                    }
                    else if (dir == "y")
                    {
                        sim_p1.y = 0;
                        bcs[1][0] = (phase == 1) ? PMC : PEC;
                    }
                    else if (dir == "z")
                    {
                        sim_p1.z = 0;
                        bcs[2][0] = (phase == 1) ? PMC : PEC;
                    }

                    std::cout << dir + " symmetry set up\n";
                }
                catch (...)
                {
                    std::cout << "Failed to read symmetry conditions\n";
                }
            }
        }

        if (auto pd = config.find("periodic"); pd != config.end())
        {
            try
            {
                for (int i = 0; i < 3; ++i)
                {
                    int c;
                    pd->at(i).get_to(c);
                    if (c)
                    {
                        //if one direction is both symmetry and periodic
                        if (bcs[i][0] == PMC || bcs[i][0] == PEC)
                            bcs[i][1] = bcs[i][0];
                        //if one direction is only peridoic
                        else
                            bcs[i][0] = bcs[i][1] = Period;
                    }
                }
            }
            catch (...)
            {
                std::cout << "Failed to read periodic boundary conditions\n";
            }
        }
    }

    void Simulation::set_pmls(const json &pmls)
    {
        std::array<double_arr, 3> k, b, c;
        k[0].resize(sim_dim.x, 1);
        k[1].resize(sim_dim.y, 1);
        k[2].resize(sim_dim.z, 1);

        b[0].resize(sim_dim.x);
        b[1].resize(sim_dim.y);
        b[2].resize(sim_dim.z);

        c[0].resize(sim_dim.x);
        c[1].resize(sim_dim.y);
        c[2].resize(sim_dim.z);

        auto one_side_setup = [&](const Direction dir, const Side side, double thickness, const PML &pml) {
            for (int i = sim_p1[dir]; i <= sim_p2[dir]; ++i)
                if (side == Positive && (i * dx2 >= sim_size[dir] / 2 - thickness))
                {
                    double x = i * dx2 - (sim_size[dir] / 2 - thickness);
                    int index = i - sim_p1[dir];

                    k[dir][index] = pml.get_k(x);
                    b[dir][index] = pml.get_b(x, dt);
                    c[dir][index] = pml.get_c(x, dt);
                }
                else if (side == Negative && (i * dx2 <= -sim_size[dir] / 2 + thickness))
                {
                    double x = -sim_size[dir] / 2 + thickness - i * dx2;
                    int index = i - sim_p1[dir];

                    k[dir][index] = pml.get_k(x);
                    b[dir][index] = pml.get_b(x, dt);
                    c[dir][index] = pml.get_c(x, dt);
                }
        };

        for (auto itr = pmls.begin(); itr != pmls.end(); ++itr)
        {
            double thickness = itr->at("thickness").get<double>();
            double k_max = itr->at("k max").get<double>();
            double sigma_max = itr->at("sigma max").get<double>();
            int order = itr->at("order").get<int>();

            PML pml(thickness, sigma_max, k_max, order);

            int dir0 = 0, dir1 = 2;
            int side0 = -1, side1 = 1;

            if (auto find = itr->find("direction"); find != itr->end())
            {
                dir0 = dir1 = json2direction(itr->at("direction"));
            }

            if (auto find = itr->find("side"); find != itr->end())
                side0 = side1 = json2side(itr->at("side"));

            for (int dir = dir0; dir <= dir1; ++dir)
                for (int side = side0; side <= side1; side += 2)
                {
                    if (rank == 0)
                    {
                        std::cout << "Setting up pml, direction=" << dir_map[static_cast<Direction>(dir)] <<
                        ",side=" << side_map[static_cast<Side>(side)] << ",thickness=" << thickness << "\n";
                    }
                    one_side_setup(static_cast<Direction>(dir), static_cast<Side>(side), thickness, pml);
                }
        }

        fields.init(k, b, c, sim_p1, sim_p2);
        // std::cout << "Process " << rank << " PML initiliazed\n";
    }

    void Simulation::set_fields_output(const json &j)
    {
        // wait_until_attach();
        for (auto itr = j.begin(); itr != j.end(); ++itr)
        {
            std::string type_str = itr->at("type").get<std::string>();

            if (type_str == "volume fields dft")
            {
                volume_fields_dft.push_back(Volume_Fields_DFT(*itr, dft_hub));
            }
            else
                std::cout << "Unable to read fields_output configuration\n";
        }

        dft_hub.init();

        // std::cout << "Process " << rank << " Fields output initialized\n";
    }

    void Simulation::step_e()
    {
        fields.step_accd(cart_comm, (time_step + 0.5) * dt);
        structure.step_e(fields.accdb, fields.eh);
        fields.sync_boundary(cart_comm, X);
        fields.sync_boundary(cart_comm, Y);
        fields.sync_boundary(cart_comm, Z);
        dft_hub.step_e((time_step + 1) * dt, fields);
    }

    void Simulation::step_m()
    {
        fields.step_accb(cart_comm, (time_step)*dt);
        structure.step_m(fields.accdb, fields.eh);
        fields.sync_boundary(cart_comm, X);
        fields.sync_boundary(cart_comm, Y);
        fields.sync_boundary(cart_comm, Z);
        dft_hub.step_m((time_step + 0.5) * dt, fields);
    }

    void Simulation::run_until_time(double time, std::ostream &os)
    {
        sim_start_time = std::chrono::system_clock::now();
        sim_prev_time = sim_start_time;

        while (time_step * dt < time)
        {
            step_m();
            step_e();
            ++time_step;
            report();
        }
    }

    void Simulation::run_until_dft(const fVec3& p1, const fVec3& p2, const Coord_Type ctype, double frequency, double interval, double var, std::ostream &os)
    {
        sim_start_time = std::chrono::system_clock::now();
        sim_prev_time = sim_start_time;
        double src_turnoff_time = fields.get_source_turnoff_time();
        auto closure = get_component_closure(p1, p2, ctype);

        auto ptr = dft_hub.find_unit(closure.first, closure.second, ctype);

        if (ptr == nullptr)
            throw std::runtime_error("dft is not registered");
        else
        {
            if (!ptr->has_frequencies({frequency}))
                throw std::runtime_error("frequency is not registered");
        }
        

        while (time_step * dt < src_turnoff_time)
        {
            step_m();
            step_e();
            ++time_step;
            report();
        }

        int extension = 1;
        double max_norm;
        double min_norm;
        double cur_var;

        do
        {
            max_norm = std::numeric_limits<double>::min();
            min_norm = std::numeric_limits<double>::max();

            while (time_step * dt < src_turnoff_time + extension * interval)
            {
                step_m();
                step_e();
                ++time_step;
                report();

                double local_norm = ptr->get_local_norm(frequency);
                double glob_norm = 0;

                MPI_Allreduce(&local_norm, &glob_norm, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
                
                if (glob_norm > max_norm) max_norm = glob_norm;
                if (glob_norm < min_norm) min_norm = glob_norm;
            }
            extension++;

            cur_var = (max_norm - min_norm) / max_norm;

            if (rank == 0)
                os << "max norm=" << max_norm << ",min norm=" << min_norm << ",var=" << cur_var << std::endl;

        } while (cur_var > var);

    }

    void Simulation::run_until_fields_decayed(double interval, const fVec3 &pt, const Coord_Type ctype, double decay_by, std::ostream &os)
    {
        sim_start_time = std::chrono::system_clock::now();
        sim_prev_time = sim_start_time;
        double src_turnoff_time = fields.get_source_turnoff_time();
        double max_field_sqr = 0;

        while (time_step * dt < src_turnoff_time)
        {
            step_m();
            step_e();
            ++time_step;
            report();

            double local_field = fields.get_eh_helper(pt, ctype);
            double global_field = 0;

            MPI_Allreduce(&local_field, &global_field, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
            global_field *= global_field;


            if (global_field > max_field_sqr)
                max_field_sqr = global_field;
        }

        int extension = 1;
        double cur_max_field_sqr;
        double decay;
        do
        {
            cur_max_field_sqr = 0;
            while (time_step * dt < src_turnoff_time + extension * interval)
            {
                step_m();
                step_e();
                ++time_step;
                report();

                double local_field = fields.get_eh_helper(pt, ctype);
                double global_field = 0;

                MPI_Allreduce(&local_field, &global_field, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
                global_field *= global_field;

                if (global_field > cur_max_field_sqr)
                    cur_max_field_sqr = global_field;
            }
            extension++;

            if (cur_max_field_sqr > max_field_sqr)
                max_field_sqr = cur_max_field_sqr;

            decay = cur_max_field_sqr / max_field_sqr;
            if (rank == 0)
                os << "Field Decay = " << cur_max_field_sqr  << "/" <<  max_field_sqr << "=" << decay << std::endl;

        } while (decay > decay_by);
    }

    void Simulation::run(const json &stop_cond, std::ostream& os)
    {
        time_step = 0;
        std::string cond_str;
        stop_cond.at("type").get_to(cond_str);

        //run until a certain reaches
        if (cond_str == "time")
        {
            double stop_time = stop_cond.at("time").get<double>();
            
            if (rank == 0)
                os << "Simulation stops at time=" << stop_time << "\n";

            run_until_time(stop_time, os);
        }
        else if (cond_str == "decay")
        {
            //run until the field decays at a particular point
            auto pt = json2fvec3(stop_cond.at("position"));
            Coord_Type ctype = json2ctype(stop_cond.at("field component"));
            auto ctype_str = stop_cond.at("field component").get<std::string>();
            double time_interval = stop_cond.at("time interval examined").get<double>();
            double decayed_by = stop_cond.at("decayed by").get<double>();

            if (rank == 0)
                os << "Simulation stops when pt=" << pt << "(" << ctype_str << ")"
                    << " decalys " << decayed_by << "\n";

            run_until_fields_decayed(time_interval, pt / dx2, ctype, decayed_by, os);
        }
        else
        {
            fVec3 center = json2fvec3(stop_cond.at("center"));
            fVec3 size = json2fvec3(stop_cond.at("size"));
            Coord_Type ctype = json2ctype(stop_cond.at("field component"));
            auto ctype_str = stop_cond.at("field component").get<std::string>();
            double interval = stop_cond.at("time interval examined").get<double>();
            double var = stop_cond.at("variation").get<double>();
            double freq = stop_cond.at("frequency").get<double>();

            if (rank == 0)
                os << "Simulation stops when region=(center=" << center << ",size=" << size << ",component=" << ctype_str 
                << ") has steady dft with variation=" << var << "\n"; 

            run_until_dft((center-size/2) / dx2, (center + size/2) / dx2, ctype, freq, interval, var, os);
        }

        if (rank == 0)
        {
            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - sim_start_time;
            os << "Simulation complete, Total time spent = " << elapsed_seconds.count() << " s\n";
            os << "Last time step = " << time_step << "\n";
        }
    }

    void Simulation::output_details(std::ostream& os) const
    {
        os<<"Report simulation details\n";
        fields.output_details(os);
        structure.output_details(os);
    }

    void Simulation::output()
    {
        for (auto &item : volume_fields_dft)
        {
            item.output(cart_comm, dft_hub, *fields_output_file);
        }

        MPI_Barrier(cart_comm);
        delete fields_output_file;
        delete input_file;
    }

    void Simulation::report()
    {
        if (rank != 0)
            return;

        auto end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - sim_prev_time;

        if (elapsed_seconds.count() > progress_interval)
        {
            sim_prev_time = end;
            std::cout << "Time Step = " << time_step << ", time = " << time_step * dt << std::endl;
        }
    }

    double Simulation::at(const fVec3 &pt, const Coord_Type ctype)
    {
        return fields.get(pt / dx2, ctype);
    }
} // namespace ffip
