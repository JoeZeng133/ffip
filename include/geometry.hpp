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

    //geometry class wrapper
    //allow material function with a discriptor
    class Geometry {
        public:
            //return whether a given physical coordinate is inside the geometry
            virtual bool is_inside(const fVec3& pos) const = 0;

            //return whether the geometry is homogeneous, i.e., does it rely on a discriptor
            //most geometry should be homogeneous
            virtual bool is_homogeneous() const = 0;

            //return discriptor at a given phsyical coordinate, 
            //if homogeneous, should return empty discriptor
            //if outside, 
            virtual std::valarray<double> get_discriptor_from_coord(const fVec3& phys_coord) const = 0;

            //return Medium from a given discrptor
            //if homogeneous, should return the same medium no matter what
            //if discriptor does not match, throw error
            virtual Medium get_medium_from_discriptor(const std::valarray<double>& disc) const = 0;
    };


    class Sphere : Geometry {
        public:
            bool is_inside(const fVec3& pos) const override;
            
            bool is_homogeneous() const override;

            std::valarray<double> get_discriptor_from_coord(const fVec3& phys_coord) const = 0;

            Medium get_medium_from_discriptor(const std::valarray<double>& disc) const = 0;

        private:
            Sphere_3 geom;

    };

    class Block : Geometry {
    public:
        bool is_inside(const fVec3& pos) const override;
        
        bool is_homogeneous() const override;

        std::valarray<double> get_discriptor_from_coord(const fVec3& phys_coord) const = 0;

        Medium get_medium_from_discriptor(const std::valarray<double>& disc) const = 0;

    private:
        

    };
}