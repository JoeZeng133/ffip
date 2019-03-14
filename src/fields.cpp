#include <fields.hpp>

namespace ffip {
    /*
        Serialized Gaussian Dipoles 
    */
    void Gaussian_Dipoles_Stepping::organize() {
        sort(gaussian1_points.begin(), gaussian1_points.end(), [](const Stepping_Point& a, const Stepping_Point& b){
            return a.index < b.index;
        });

        sort(gaussian2_points.begin(), gaussian2_points.end(), [](const Stepping_Point& a, const Stepping_Point& b){
            return a.index < b.index;
        });
    }

    void Gaussian_Dipoles_Stepping::step(double time, std::vector<double>& diffdb) {
        if (time > max_end_time) return;

        for (auto& point : gaussian1_points) {
            if (time > point.start_time && time < point.end_time)
                diffdb[point.index] += point.amp * Gaussian1(time - point.start_time, point.width);
        }

        for (auto& point : gaussian2_points) {
            if (time > point.start_time && time < point.end_time)
                diffdb[point.index] += point.amp * Gaussian2(time - point.start_time, point.width);
        }
    }

    void Gaussian_Dipoles_Stepping::add_gaussian1(size_t index, double amp, double start_time, double end_time, double width) {
        gaussian1_points.push_back({index, amp, start_time, end_time, width});
    }

    void Gaussian_Dipoles_Stepping::add_gaussian2(size_t index, double amp, double start_time, double end_time, double width) {
        gaussian2_points.push_back({index, amp, start_time, end_time, width});
    }

    void Gaussian_Dipoles_Stepping::output(std::ostream& os) {
    }

    /*
        Serialize PML steppping
    */
    void PML_Stepping::set_strides(sVec3 strides) {
        this->strides = strides;
    }

    void PML_Stepping::organize() {
        for(int i = 0; i < 3; ++i)
            sort(points[i].begin(), points[i].end(), [](const Stepping_Point& a, const Stepping_Point& b){
                return a.index < b.index;
            });
    }

    void PML_Stepping::step(const std::vector<double>& eh, std::vector<double>& diffdb) {
        for(int i = 0; i < 3; ++i) {
            
			size_t stride1 = strides[(i + 1) % 3];
			size_t stride2 = strides[(i + 2) % 3];
			auto& rpoints = points[i];

			for (auto& point : points[i]) {
				double curl1 = (eh[point.index + stride1] - eh[point.index - stride1]);
				double curl2 = (eh[point.index + stride2] - eh[point.index - stride2]);

				point.phi1 = point.b1 * point.phi1 + point.c1 * curl1;
				point.phi2 = point.b2 * point.phi2 + point.c2 * curl2;

				diffdb[point.index] = point.phi1 - point.phi2 + curl1 * point.k1 - curl2 * point.k2;
			}
        }
    }

    void PML_Stepping::add_stepping_Point(Direction dir, size_t index, double b1, double c1, double k1, double b2, double c2, double k2) {
        points[dir].push_back({index, 0, b1, c1, k1, 0, b2, c2, k2});
    }

    /*
        fields
    */
    void Fields::add_dipole_source_gaussian1(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff) {

        auto interp_info = cells.get_interp_info(0, pos, ctype);
        auto& gaussian_dipoles = (is_E_Point(ctype)? e_gaussian_dipoles : m_gaussian_dipoles);
        double width = std::sqrt(2.0) / (2 * pi * frequency);
        double actual_end_time = std::min(end_time, start_time + 2 * width * cutoff);

        for(auto& point : interp_info) {
            gaussian_dipoles.add_gaussian1(point.first, amp, start_time, actual_end_time, width);
        }
    }

    void Fields::add_dipole_source_gaussian2(fVec3 pos, Coord_Type ctype, double amp, double start_time, double end_time, double frequency, double cutoff) {
        auto interp_info = cells.get_interp_info(0, pos, ctype);
        auto& gaussian_dipoles = (is_E_Point(ctype)? e_gaussian_dipoles : m_gaussian_dipoles);
        double width = std::sqrt(2.0) / (2 * pi * frequency);
        double actual_end_time = std::min(end_time, start_time + 2 * width * cutoff);

        for(auto& point : interp_info) {
            gaussian_dipoles.add_gaussian2(point.first, amp, start_time, actual_end_time, width);
        }
    }

    

}