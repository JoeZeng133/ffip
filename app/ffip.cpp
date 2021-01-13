#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <iostream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_2.h>
#include <iostream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_2;
using Polygon_2 = CGAL::Polygon_2<K>;
namespace py = pybind11;

template <int N>
class interpn : public interpn<N - 1>
{
public:
    using base_class = interpn<N - 1>;
    size_t dimn{1}, stride{1};

    interpn() = default;

    template <typename... Args>
    interpn(const size_t dimn, Args... args) : base_class(args...), dimn(dimn)
    {
        if (dimn < 1)
            throw std::runtime_error("Invalid Dimension");

        stride = base_class::get_size();
    }

    //get the array size
    size_t get_size() const
    {
        return dimn * base_class::get_size();
    }

    //interpolate xn, xn-1, ..., x1
    template <typename T, typename... Args>
    T operator()(T const *data, double xq, Args &&... args) const
    {
        //ignore this dimension if it is 1
        if (dimn == 1)
            return base_class::operator()(data, args...);

        //0 padding extrapolation
        if (xq <= 0)
        {
            if (xq <= -1)
                return T{};
            return (xq + 1) * base_class::operator()(data, args...);
        }

        if (xq >= dimn - 1)
        {
            if (xq >= dimn)
                return T{};
            return (dimn - xq) * base_class::operator()(data + stride * (dimn - 1), args...);
        }

        size_t index = xq;
        double tx = xq - index;

        return tx * base_class::operator()(data + stride * (index + 1), args...) +
               (1 - tx) * base_class::operator()(data + stride * index, args...);
    }

    //interpolate xn, xn-1, ..., x1
    template <typename T, typename... Args>
    T operator()(const std::vector<T> &vec, Args &&... args) const
    {
        return operator()(vec.data(), args...);
    }

    //transpose_helper
    template <typename T, typename... Args>
    void transpose(T *data, const T &val, double xq, Args &&... args) const
    {

        //ignore 1 dimension
        if (dimn == 1)
        {
            base_class::transpose(data, val, args...);
            return;
        }

        //0 padding extrapolation
        if (xq <= 0)
        {
            if (xq <= -1)
                return;
            base_class::transpose(data, val * (1 + xq), args...);
            return;
        }

        if (xq >= dimn - 1)
        {
            if (xq >= dimn)
                return;
            base_class::transpose(data + stride * (dimn - 1), val * (dimn - xq), args...);
            return;
        }

        size_t index = xq;
        double tx = xq - index;

        base_class::transpose(data + stride * (index + 1), val * tx, args...);
        base_class::transpose(data + stride * index, val * (1 - tx), args...);

        return;
    }

    //transpose val additively back to original grid
    template <typename T, typename... Args>
    void transpose(std::vector<T> &vec, const T &val, Args &&... args) const
    {
        transpose(vec.data(), val, args...);
    }
};

template <>
class interpn<1>
{
public:
    size_t dimn{1}, stride{1};

    interpn() = default;

    interpn(const size_t _dimn) : dimn(_dimn)
    {

        if (dimn < 1)
            throw std::runtime_error("Invalid Dimension");

        stride = 1;
    }

    size_t get_size() const
    {
        return dimn;
    }

    template <typename T>
    T operator()(const T *data, double xq) const
    {

        //ignore this dimension if it is 1
        if (dimn == 1)
            return data[0];

        //0 padding extrapolation
        if (xq <= 0)
        {
            if (xq <= -1)
                return T{};
            return (xq + 1) * data[0];
        }

        if (xq >= dimn - 1)
        {
            if (xq >= dimn)
                return T{};
            return (dimn - xq) * data[dimn - 1];
        }

        size_t index = xq;
        double tx = xq - index;

        return tx * data[index + 1] + (1 - tx) * data[index];
    }

    template <typename T>
    T operator()(const std::vector<T> &vec, double xq) const
    {
        return operator()(vec.data(), xq);
    }

    template <typename T>
    void transpose(std::vector<T> &vec, const T &val, double xq) const
    {
        return transpose(vec.data(), val, xq);
    }

    template <typename T>
    void transpose(T *data, const T &val, double xq) const
    {

        //ignore this dimension if it is 1
        if (dimn == 1)
        {
            data[0] += val;
            return;
        }

        //0 padding extrapolation
        if (xq <= 0)
        {
            if (xq <= -1)
                return;
            data[0] += val * (1 + xq);
            return;
        }

        if (xq >= dimn - 1)
        {
            if (xq >= dimn)
                return;
            data[dimn - 1] += val * (dimn - xq);
            return;
        }

        size_t index = xq;
        double tx = xq - index;

        data[index + 1] += val * tx;
        data[index] += val * (1 - tx);
    }
};

py::array TO_convolve(py::array_t<double> in1, py::array_t<double> in2)
{
    if (in1.ndim() != 2 || in2.ndim() != 2)
        throw py::key_error("dimensions are not 2");
    
    if (in2.shape(0) % 2 == 0 || in2.shape(1) % 2 == 0)
        throw py::key_error("filter shape is not odd");
    
    auto p1 = in1.unchecked<2>();
    auto p2 = in2.unchecked<2>();

    auto res = py::array_t<double>();
    res.resize({p1.shape(0), p1.shape(1)});

    auto pr = res.mutable_unchecked<2>();

    int r0 = in2.shape(0) / 2;
    int r1 = in2.shape(1) / 2;

    for(ssize_t i = 0; i < in1.shape(0); ++i)
        for(ssize_t j = 0; j < in1.shape(1); ++j)
        {
            double w = 0;
            pr(i, j) = 0;
            ssize_t l0 = std::max(i - r0, 0L);
            ssize_t l1 = std::max(j - r1, 0L);
            ssize_t u0 = std::min(i + r0 + 1, in1.shape(0));
            ssize_t u1 = std::min(j + r1 + 1, in1.shape(1)); 

            for(ssize_t m = l0; m < u0; ++m)
                for(ssize_t n =l1; n < u1; ++n)
                {
                    w += p2(m - i + r0, n - j + r0);
                    pr(i, j) += p1(m, n) * p2(m - i + r0, n - j + r0);
                }
            
            pr(i, j) /= w;
        }

    return res;
}

py::array TO_convolve_transpose(py::array_t<double> in1, py::array_t<double> in2)
{
    if (in1.ndim() != 2 || in2.ndim() != 2)
        throw py::key_error("dimensions are not 2");
    
    if (in2.shape(0) % 2 == 0 || in2.shape(1) % 2 == 0)
        throw py::key_error("filter shape is not odd");
    
    auto p1 = in1.unchecked<2>();
    auto p2 = in2.unchecked<2>();

    auto res = py::array_t<double>();
    res.resize({p1.shape(0), p1.shape(1)});

    auto pr = res.mutable_unchecked<2>();

    int r0 = in2.shape(0) / 2;
    int r1 = in2.shape(1) / 2;

    for(ssize_t i = 0; i < in1.shape(0); ++i)
        for(ssize_t j = 0; j < in1.shape(1); ++j)
        {
            double w = 0;
            ssize_t l0 = std::max(i - r0, 0L);
            ssize_t l1 = std::max(j - r1, 0L);
            ssize_t u0 = std::min(i + r0 + 1, in1.shape(0));
            ssize_t u1 = std::min(j + r1 + 1, in1.shape(1)); 

            for(ssize_t m = l0; m < u0; ++m)
                for(ssize_t n =l1; n < u1; ++n)
                {
                    w += p2(m - i + r0, n - j + r0);
                }
            
            for(ssize_t m = l0; m < u0; ++m)
                for(ssize_t n =l1; n < u1; ++n)
                {
                    pr(m, n) += p1(i, j) * p2(m - i + r0, n - j + r0) / w;
                }
        }

    return res;
}

py::array transpose(py::array_t<double> gx, py::array_t<double> gy, py::array_t<double> gz,
                    py::array_t<double> pt)
{
    if (gx.ndim() != 1 || gy.ndim() != 1 || gz.ndim() != 1 || pt.ndim() != 2)
        throw py::key_error("some dimensions are not right");
    
    if (gx.size() < 1 || gy.size() < 1 || gz.size() < 1)
        throw py::key_error("some grid size is 0");

    if (pt.shape(1) != 4)
        throw std::range_error("shape(1) of pts is not 4");

    auto gx_ptr = gx.unchecked<1>();
    auto gy_ptr = gy.unchecked<1>();
    auto gz_ptr = gz.unchecked<1>();
    auto pt_ptr = pt.unchecked<2>();

    double dx = 1, dy = 1, dz = 1;
    if (gx_ptr.shape(0) > 1)
        dx = gx_ptr(1) - gx_ptr(0);
    if (gy_ptr.shape(0) > 1)
        dy = gy_ptr(1) - gy_ptr(0);
    if (gz_ptr.shape(0) > 1)
        dz = gz_ptr(1) - gz_ptr(0);

    interpn<3> interp(gz_ptr.shape(0), gy_ptr.shape(0), gx_ptr.shape(0));
    auto res = py::array_t<double>();
    res.resize({gz_ptr.shape(0), gy_ptr.shape(0), gx_ptr.shape(0)});
    double* data = res.mutable_data();

    // std::cout << "number of points to transpose=" << pt.size() << "\n";

    // double sum = 0;
    
    for (ssize_t i = 0; i < pt_ptr.shape(0); ++i)
    {
        double sx = (pt_ptr(i, 0) - gx_ptr(0)) / dx;
        double sy = (pt_ptr(i, 1) - gy_ptr(0)) / dy;
        double sz = (pt_ptr(i, 2) - gz_ptr(0)) / dz;

        interp.transpose(data, pt_ptr(i, 3), sz, sy, sx);
    }

    return res;
}

py::array check_inside(py::array_t<double> req_pts, py::array_t<double> poly_pts)
{
    if (req_pts.ndim() < 2 || poly_pts.ndim() < 2)
        throw py::value_error("ndim is not bigger than 1");

    if (req_pts.shape(req_pts.ndim() - 1) != 2)
        throw py::value_error("req_pts.shape[-1] is not 2");
    
    if (poly_pts.shape(poly_pts.ndim() - 1) != 2)
        throw py::value_error("poly_pts.shape[-1] is not 2");

    double const *poly_pts_ptr = poly_pts.data();
    double const *req_ptr_ptr = req_pts.data();

    auto res = py::array_t<bool>();

    res.resize({req_pts.shape(), req_pts.shape() + req_pts.ndim() - 1});
    bool* res_ptr = res.mutable_data();

    auto poly_pts_vec = std::vector<Point>{};
    
    for(ssize_t i = 0; i < poly_pts.size(); i += 2)
    {
        poly_pts_vec.push_back({poly_pts_ptr[i], poly_pts_ptr[i + 1]});
    }

    // auto pgn = Polygon_2(poly_pts_vec.begin(), poly_pts_vec.end());

    for(ssize_t i = 0; i < req_pts.size(); i+=2)
    {
        res_ptr[i / 2] = bounded_side_2(
            poly_pts_vec.begin(), 
            poly_pts_vec.end(), 
            {req_ptr_ptr[i], req_ptr_ptr[i + 1]}) == CGAL::ON_BOUNDED_SIDE;
        // res_ptr[i / 2] = (pgn.bounded_side({req_ptr_ptr[i], req_ptr_ptr[i + 1]}) == CGAL::ON_BOUNDED_SIDE);
    }

    return res;
}

PYBIND11_MODULE(_ffip, m)
{
    m.doc() = "ffip c++ libraries";

    m.def("check_inside", &check_inside,
        py::arg("req_pts"),
        py::arg("poly_pts"),
        "check inside 2 D polygon");

    m.def("transpose", &transpose,
        py::arg("gx"),
        py::arg("gy"),
        py::arg("gz"),
        py::arg("pt"),
        "shape(pt) = ... x 4"
        );
    
    m.def("TO_convolve", &TO_convolve,
        py::arg("in1"),
        py::arg("in2"),
        "convolution in topology optimization");
    
    m.def("TO_convolve_transpose", &TO_convolve_transpose,
        py::arg("in1"),
        py::arg("in2"),
        "transpose of convoluation in topology optimization");
}