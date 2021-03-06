#pragma once

#include <medium.hpp>
#include <utility.hpp>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

namespace ffip
{
    using Spherical_K = CGAL::Exact_spherical_kernel_3;
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

    using Point_3 = CGAL::Point_3<Kernel>;
    using Sphere_3 = CGAL::Sphere_3<Kernel>;
    using Line_3 = CGAL::Line_3<Kernel>;
    using Cuboid_3 = CGAL::Iso_cuboid_3<Kernel>;

    Point_3 vec3_to_point3(const fVec3 &p);

    fVec3 point3_to_vec3(const Point_3 &p);

    //geometry class wrapper
    class Geometry
    {
    public:
        //return whether a given physical coordinate is inside the geometry
        virtual bool is_inside(const fVec3 &pos) const = 0;

        //return whether the geometry is homogeneous
        //most geometry should be homogeneous
        virtual bool is_homogeneous() const = 0;

        //return medium at particular position, if homogeneous, pos ignored
        virtual Abstract_Medium get_medium(const fVec3 &pos) const = 0;
		
		virtual ~Geometry() {};
    };

    //Sphere geometry
    class Sphere : public Geometry
    {
    public:
        bool is_inside(const fVec3 &pos) const override;

        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3 &pos) const override;

        Sphere(const fVec3 &center, const double radius, const Abstract_Medium &medium);

    private:
        Sphere_3 geom;
        Abstract_Medium medium;
        fVec3 center;
        double radius;
    };

    //Box geometry
    class Box : public Geometry
    {
    public:
        bool is_inside(const fVec3 &pos) const override;

        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3 &pos) const override;

        Box(const fVec3 &center, const fVec3 &size, const Abstract_Medium &medium);

    private:
        Cuboid_3 geom;
        Abstract_Medium medium;
        fVec3 p1, p2;
    };

    //Two material mixing in a block
    class Two_Medium_Box : public Geometry
    {
    public:
        bool is_inside(const fVec3 &pos) const override;

        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3 &pos) const override;

        Two_Medium_Box(const fVec3 &center, const fVec3 &size, const iVec3 &dim, const Abstract_Medium &m1, const Abstract_Medium &m2, const std::vector<double> &rho);

    private:
        Cuboid_3 geom;
        Abstract_Medium m1, m2;

        fVec3 p1, p2;
        iVec3 dim;
        std::vector<double> rho;
        interpn<3> interp;
    };

    //independent pole intesnity and epsilon
    class General_Medium_Box : public Geometry
    {
    public:
        bool is_inside(const fVec3 &pos) const override;

        bool is_homogeneous() const override;

        Abstract_Medium get_medium(const fVec3 &pos) const override;

        General_Medium_Box(const fVec3 &center, const fVec3 &size, const iVec3 &dim, const std::vector<double> &rho, const std::vector<double> &rho_fun, const std::vector<Abstract_Medium> &medium_fun);

    private:
        Cuboid_3 geom;

        fVec3 p1, p2;
        iVec3 dim;
        std::vector<double> rho;
        interpn<3> interp;

        //nonlinear dependence material properties
        interpn<1> medium_interp;
        std::vector<double> rho_fun;
        std::vector<Abstract_Medium> medium_fun;
    };
} // namespace ffip
