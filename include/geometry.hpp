#pragma once

#include <medium.hpp>
#include <utility.hpp>

#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>

using Spherical_K = CGAL::Exact_spherical_kernel_3;
using Point_3 = CGAL::Point_3<Spherical_K>;
using Sphere_3 = CGAL::Sphere_3<Spherical_K>;
using Segment_3 = CGAL::Segment_3<Spherical_K>;
using Line_3 = CGAL::Line_3<Spherical_K>;
using ArcPoint_3 = CGAL::Circular_arc_point_3<Spherical_K>;
using LineArc_3 = CGAL::Line_arc_3<Spherical_K>;

namespace ffip {
    class Geometry {
        public:
        private:
    };
}