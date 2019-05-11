#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

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

void add_arrays(py::array_t<double> &input1, py::array_t<double> input2)
{
    py::buffer_info buf1 = input1.request(), buf2 = input2.request();

    if (buf1.ndim != 1 || buf2.ndim != 1)
        throw std::runtime_error("Number of dimensions must be one");

    if (buf1.size != buf2.size)
        throw std::runtime_error("Input shapes must match");

    /* No pointer is passed, so NumPy will allocate the buffer */
    auto result = py::array_t<double>(buf1.size);

    py::buffer_info buf3 = result.request();

    double *ptr1 = (double *)buf1.ptr,
           *ptr2 = (double *)buf2.ptr,
           *ptr3 = (double *)buf3.ptr;

    for (size_t idx = 0; idx < buf1.shape[0]; idx++)
        ptr1[idx] = ptr1[idx] + ptr2[idx];

    // return result;
}

void add1_array(py::array_t<double> input)
{
    auto r = input.mutable_unchecked<1>();
    for (ssize_t i = 0; i < r.shape(0); ++i)
        r(i) += 1;
}

py::array zeros(py::array a)
{
    auto res = py::array(
        a.dtype(), {}, {});

    res.resize({a.shape(), a.shape() + a.ndim()});
    return res;
}



template<typename T>
py::array transpose(py::array_t<double> gx, py::array_t<double> gy, py::array_t<double> gz,
                    py::array_t<double> x, py::array_t<double> y, py::array_t<double> z, py::array_t<T> v)
{
    if (gx.ndim() != 1 || gy.ndim() != 1 || gz.ndim() != 1 ||
        x.ndim() != 1 || y.ndim() != 1 || z.ndim() != 1 || v.ndim() != 1)
        throw py::key_error("some dimension is not 1");
    
    if (gx.size() < 1 || gy.size() < 1 || gz.size() < 1)
        throw py::key_error("some grid size is 0");

    if (x.size() != y.size() || y.size() != z.size() || z.size() != v.size())
        throw std::range_error("size of points are not consistent");

    auto gx_ptr = gx.unchecked<1>();
    auto gy_ptr = gy.unchecked<1>();
    auto gz_ptr = gz.unchecked<1>();

    auto x_ptr = x.unchecked<1>();
    auto y_ptr = y.unchecked<1>();
    auto z_ptr = z.unchecked<1>();
    const T* v_ptr = v.data();

    double dx = 1, dy = 1, dz = 1;
    if (gx_ptr.shape(0) > 1)
        dx = gx_ptr(1) - gx_ptr(0);
    if (gy_ptr.shape(0) > 1)
        dy = gy_ptr(1) - gy_ptr(0);
    if (gz_ptr.shape(0) > 1)
        dz = gz_ptr(1) - gz_ptr(0);

    interpn<3> interp(gz_ptr.shape(0), gy_ptr.shape(0), gx_ptr.shape(0));
    auto res = py::array_t<T>();
    res.resize({gz_ptr.shape(0), gy_ptr.shape(0), gx_ptr.shape(0)});
    T* data = res.mutable_data();

    for (ssize_t i = 0; i < x_ptr.size(); ++i)
    {
        double sx = (x_ptr(i) - gx_ptr(0)) / dx;
        double sy = (y_ptr(i) - gy_ptr(0)) / dy;
        double sz = (z_ptr(i) - gz_ptr(0)) / dz;

        interp.transpose(data, v_ptr[i], sz, sy, sx);
    }

    return res;
}

py::array linearfilter2D(py::array_t<double> input, py::array_t<int> nbhd, py::array_t<double> weight)
{
    if (input.ndim() != 2)
        throw std::invalid_argument("input data is not 2D");

    if (nbhd.ndim() != 2 || nbhd.shape(1) != input.ndim())
        throw std::invalid_argument("neigborhood shape is invalid");
    
    if (weight.ndim() != 1 || weight.shape(0) != nbhd.shape(0))
        throw std::invalid_argument("weight shape is not compatible with neighborhood");

    py::array_t<double> res;
    res.resize({input.shape(), input.shape() + input.ndim()});

    auto in = input.unchecked<2>();
    auto ou = res.mutable_unchecked<2>();
    auto nb = nbhd.unchecked<2>();
    auto wt = weight.unchecked<1>();

    for(int i = 0; i < ou.shape(0); ++i)
    for(int j = 0; j < ou.shape(1); ++j)
    for(int k = 0; k < nb.shape(0); ++k)
    {
        int ioffset = i + nb(k, 0);
        int joffset = j + nb(k, 1);

        if (ioffset >= 0 && ioffset < ou.shape(0) && joffset >= 0 && joffset < ou.shape(1))
            ou(i, j) += wt(k) * in(ioffset, joffset);
    }

    return res;
}

PYBIND11_MODULE(pybind_practice, m)
{
    m.doc() = "examples";
    m.def("add_arrays", &add_arrays, "Add two NumPy arrays");
    m.def("add1_array", &add1_array, "Add 1 to numpy array");
    m.def("zeros", &zeros, "return zeros of the same size");

    m.def("linearfilter2D", &linearfilter2D, 
        py::arg("input"),
        py::arg("nbhd"),
        py::arg("weight"),
        "linear 2d filter providing neighborhood and weight");

    m.def("transpose", &transpose<std::complex<double>>, 
        py::arg("gx"),
        py::arg("gy"),
        py::arg("gz"),
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("v").noconvert(),
        "linear transpose");

    m.def("transpose", &transpose<double>, 
        py::arg("gx"),
        py::arg("gy"),
        py::arg("gz"),
        py::arg("x"),
        py::arg("y"),
        py::arg("z"),
        py::arg("v").noconvert(),
        "linear transpose");
}