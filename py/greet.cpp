#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <string>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <CGAL/Spherical_kernel_intersections.h>
#include <vector>
#include <iostream>


namespace bp = boost::python;
namespace np = boost::python::numpy;

// a newly creaed pyObject has reference count = 1
// a handle will decrease the reference count by 1 when destructed
// a borrow will increase the reference count by 1 when constructed
// a borrowed reference won't increase the count so use handle directly would be wrong
// a new reference increases count for you however, so you don't need to use borrow

using object = boost::python::object;

std::string greet() {
    if constexpr(1) {
        return "Oh Yeah";
    } else {
        return "Fuck Yeah";
    }
    
}

// int greet3(double r, double x1, double y1, double z1, double x2, double y2, double z2) {
//     auto p1 = Point_3(x1, y1, z1);
//     auto p2 = Point_3(x2, y2, z2);

//     auto arc = LineArc_3(Line_3(p1, p2), ArcPoint_3(p1), ArcPoint_3(p2));

//     auto ce = Point_3(0.0, 0.0, 0.0);
//     auto sh = Sphere_3(ce, r*r);

//     std::vector< CGAL::Object > intersecs;
// 	CGAL::intersection(arc, sh,std::back_inserter(intersecs));

//     // return CGAL::do_intersect(sh, arc);
//     return intersecs.size();
// }

std::string greet2() {
    return "Oh Yeah";
}

object add(object a, object b) {
    return a + b;
}

struct Medium {
    double sigma;
    Medium(double sigma): sigma(sigma) {}
};

struct Vec3 {
    double x, y, z;
    Vec3(double x, double y, double z): x(x), y(y), z(z) {}
};

// a simple work around to pass function from python to c++ is to 
// use only c++ exposed objects in the input arguments and return arguments
// speed is OK, 1 Million calls takes 1 second
// extract<A>(func(args)) or call<Medium>(func.ptr(), args)

void sth(bp::object func, Vec3 vec, int num) {
    for(int i = 0; i < num; ++i) {
        bp::call<Medium>(func.ptr(), vec);
    }

    Medium res = bp::call<Medium>(func.ptr(), vec);
    // Medium res = bp::extract<Medium>(func(vec));
    std::cout << "The Medium is " << res.sigma << "\n";
}

void spit(bp::list& vec_list) {
    int len = bp::len(vec_list);
    for(int i = 0; i < len; ++i) {
        Vec3 item = bp::extract<Vec3>(vec_list[i]);
        std::cout << item.x << " " << item.y << " " << item.z << "\n";
    }
}

void set(bp::list obj) {
	obj[0] = 100;
}

//able to pass reference properly
void see(Vec3& p) {
    p.x = p.y = p.z = 0;
}

BOOST_PYTHON_MODULE(greet)
{
    using namespace boost::python;
    def("greet", greet);
    def("greet2", greet2);
    // def("greet3", greet3);
    def("add", add);
    def("set", set);
    def("see", see);

    class_<Vec3>("Vec3", init<double, double, double>())
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z);

    class_<Medium>("Medium", init<double>())
        .def_readwrite("sigma", &Medium::sigma);

    def("spit", spit);
    def("sth", sth);
}
