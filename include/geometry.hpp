#pragma once

#include <medium.hpp>
#include <utility.hpp>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>

namespace ffip {
    using Spherical_K = CGAL::Exact_spherical_kernel_3;
    using Point_3 = CGAL::Point_3<Spherical_K>;
    using Sphere_3 = CGAL::Sphere_3<Spherical_K>;
    using Segment_3 = CGAL::Segment_3<Spherical_K>;
    using Line_3 = CGAL::Line_3<Spherical_K>;
    using ArcPoint_3 = CGAL::Circular_arc_point_3<Spherical_K>;
    using LineArc_3 = CGAL::Line_arc_3<Spherical_K>;
    using Cuboid_3 = CGAL::Iso_cuboid_3<Spherical_K>;

    Point_3 vec3_to_point3(const fVec3& p);

    fVec3 point3_to_vec3(const Point_3& p);

    //geometry class wrapper
    //alNegative material function
    class Geometry {
    public:
        //return whether a given physical coordinate is inside the geometry
        virtual bool is_inside(const fVec3& pos) const = 0;

        //return whether the geometry is homogeneous
        //most geometry should be homogeneous
        virtual bool is_homogeneous() const = 0;

        //return medium at particular position, if homogeneous, pos ignored
        virtual Abstract_Medium get_medium(const fVec3& pos) const = 0;
    };

    //Sphere geometry
    class Sphere : public Geometry {
    public:
        bool is_inside(const fVec3& pos) const override;
        
        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3& pos) const override;

        Sphere(const fVec3& center, const double radius, const Abstract_Medium& medium);
    private:

        Sphere_3 geom;
        Abstract_Medium medium;
        fVec3 center;
        double radius;
    };

    //Block geometry
    class Block : public Geometry {
    public:
        bool is_inside(const fVec3& pos) const override;
        
        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3& pos) const override;

        Block(const fVec3& center, const fVec3& size, const Abstract_Medium& medium);
    private:

        Cuboid_3 geom;
        Abstract_Medium medium;
        fVec3 p1, p2;
    };

    //Two material mixing in a block
    class Mixed2 : public Geometry {
    public:
        bool is_inside(const fVec3& pos) const override;

        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3& pos) const override;

        Mixed2(const fVec3& center, const fVec3& size, const iVec3& dim, const Abstract_Medium& m1, const Abstract_Medium& m2, const std::vector<double>& rho);
    private:

        Cuboid_3 geom;
        Abstract_Medium m1, m2;

        fVec3 p1, p2;
        iVec3 dim;
        std::vector<double> rho;
        interpn<3> interp;
    };


}