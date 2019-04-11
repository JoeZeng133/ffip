#include <geometry.hpp>

namespace ffip {
    Point_3 vec3_to_point3(const fVec3& p) {
        return {p.x, p.y, p.z};
    }

    fVec3 point3_to_vec3(const Point_3& p) {
        return {p.x(), p.y(), p.z()};
    }

    //Sphere
    bool Sphere::is_inside(const fVec3& pos) const {
        return geom.has_on_bounded_side({pos.x, pos.y, pos.z});
    }

    bool Sphere::is_homogeneous() const {
        return 1;
    }

    Abstract_Medium Sphere::get_medium(const fVec3& pos) const {
        return medium;
    }

    Sphere::Sphere(const fVec3& center, const double radius, const Abstract_Medium& medium):
    geom({center.x, center.y, center.z}, radius * radius) ,medium(medium), center(center), radius(radius) {}

    //Block
    bool Block::is_inside(const fVec3& pos) const {
        return geom.has_on_bounded_side({pos.x, pos.y, pos.z});
    }

    bool Block::is_homogeneous() const {
        return 1;
    }

    Abstract_Medium Block::get_medium(const fVec3& pos) const {
        return medium;
    }

    Block::Block(const fVec3& center, const fVec3& size, const Abstract_Medium& medium):
    geom(vec3_to_point3(center - size / 2), vec3_to_point3(center + size / 2)), medium(medium)
    {
        p1 = center - size / 2;
        p2 = center + size / 2;
    }

    //Mixe2
    Mixed2::Mixed2
    (const fVec3& center, const fVec3& size, const iVec3& dim, const Abstract_Medium& m1, const Abstract_Medium& m2, const std::vector<double>& rho):
    geom(vec3_to_point3(center - size / 2), vec3_to_point3(center + size / 2)), m1(m1), m2(m2), dim(dim)
    {
        interp = interpn<3>{dim.z, dim.y, dim.x};
        p1 = center - size / 2;
        p2 = center + size / 2;
    }

    bool Mixed2::is_inside(const fVec3& pos) const {
        return geom.has_on_bounded_side({pos.x, pos.y, pos.z});
    }

    bool Mixed2::is_homogeneous() const {
        return 0;
    }

    Abstract_Medium Mixed2::get_medium(const fVec3& p) const {
        double val = interp(rho, 
        (p.z - p1.z) / (p2.z - p1.z) * (dim.z - 1), 
        (p.y - p1.y) / (p2.y - p1.y) * (dim.y - 1),
        (p.x - p1.x) / (p2.x - p1.x) * (dim.x - 1));

        return m1 * val + m2 * (1 - val);
    }
}