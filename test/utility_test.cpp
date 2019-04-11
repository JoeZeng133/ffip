#include <utility.hpp>
#include <utility.cpp>

using namespace ffip;

void Yee3_test() {

    Yee3 grid{{-2, -2, -3}, {3, 5, 2}};
    double_arr arr(grid.get_size(), 1);

    if (grid.get_size() != 288 || !(grid.get_grid_p1() == iVec3(-1, -1, -2)) || !(grid.get_grid_p2() == iVec3(2, 4, 1)) ||
        !(grid.get_stride() == sVec3(1, 6, 48)))
        std::cout << "Yee3 failed at getter functions\n";
    
    if (!(grid.get_base_point(fVec3{1.2, 3, 4}, Ex) == iVec3(1, 2, 4)) || !(grid.get_base_point(fVec3{1.2,4,1.1}, Hx) == iVec3(0, 3, 1)))
        std::cout << "get_base_point failed\n";

    if (!(grid.get_index_from_coord(1, 1, 1) == 213))
        std::cout << "get_index_from_coord failed\n";

    std::cout << grid.interp(arr, {-1, -1.2, 0.9}, Ex) << " interp at {-1, -1.2, 0.9}, Ex\n";
    {
        auto tmp = grid.get_interp_weights({-1, -1.2, 0.9}, Ex);
        for(auto x : tmp)
            std::cout << x << " ";
        std::cout << "interpolation weights\n";
    }

    std::cout << grid.interp(arr, {-5, -4, -6}, Ex) << " interp at {-5, -4, -6}, Ex\n";
    {
        auto tmp = grid.get_interp_weights({-2, -1.1, -2.2}, Ex);
        for(auto x : tmp)
            std::cout << x << " ";
        std::cout << "interpolation weights\n";
    }

    std::cout << grid.interp(arr, {0, 0, 0}, Hy) << " interp at {0, 0, 0}, Hy\n";
    {
        auto tmp = grid.get_interp_weights({0, 0, 0}, Hy);
        for(auto x : tmp)
            std::cout << x << " ";
        std::cout << "interpolation weights\n";
    }
}

void basic_functions_test() {

    if (count_hex_bits(3) != 2 || count_hex_bits(14) != 3 || count_hex_bits(7) != 3)
        std::cout << "count_hex_bits failed\n";

    if (get_dir_from_ctype(Ex) != X || get_dir_from_ctype(Hy) != Y || get_dir_from_ctype(Dz) != Z)
        std::cout << "get_dir_from_ctype failed\n";

    if (!is_e_point(Ex) || is_e_point(Hy) || is_e_point(Bx) || !is_e_point(Dz))
        std::cout << "is_e_point failed\n";

    if (is_m_point(Ex) || !is_m_point(Hy) || !is_m_point(Bx) || is_m_point(Dz))
        std::cout << "is_m_point failed\n";

    if (!is_eh_point(Ex) || !is_eh_point(Hz) || is_eh_point(Bx) || is_eh_point(Dz))
        std::cout << "is_eh_point failed\n";

    if (is_db_point(Ex) || is_db_point(Hz) || !is_db_point(Bx) || !is_db_point(Dz))
        std::cout << "is_db_point failed\n";

    if (get_e_ctype_from_dir(X) != Ex || get_e_ctype_from_dir(Y) != Ey || get_e_ctype_from_dir(Z) != Ez)
        std::cout << "get_e_ctype_from_dir failed\n";
    
    if (get_m_ctype_from_dir(X) != Hx || get_m_ctype_from_dir(Y) != Hy || get_m_ctype_from_dir(Z) != Hz)
        std::cout << "get_m_ctype_from_dir failed\n";

    if (!(iVec3(1, 1, 1) * iVec3(2, 3, 4) == iVec3(2, 3, 4)) ||
        !(iVec3(1, 1, 1) + iVec3(2, 3, 4) == iVec3(3, 4, 5)) ||
        !(iVec3(3, 6, 7) / iVec3(2, 2, 2) == iVec3(1, 3, 3)) ||
        !(iVec3(3, 2, 1) - iVec3(4, 2, 1) == iVec3(-1, 0, 0)))
        std::cout << "vec3 arithmic wrong\n";
}

void Grid_3_test() {
    Grid_3 grid3({-10, -10, -11}, {22, 22, 15});
    double_arr arr(grid3.get_size(), 1);

    std::cout << grid3(arr, {-13, 0, 0}) << " at (-13, 0, 0)\n";
    std::cout << grid3(arr, {-12, 0, 0}) << " at (-12, 0, 0)\n";
    std::cout << grid3(arr, {-11, -11, -12}) << " at (-11, -11, -12)\n";
    std::cout << grid3(arr, {0, 0, 0}) << " at (0, 0, 0)\n";
}

void decomposition_test() {
    iVec3 dim = {49, 25, 9};
    auto num = decompose_domain(dim, 105);
    std::cout << "num = " << num << "\n";

    auto max_chunk = get_max_size_chunk(dim, num);
    std::cout << "max_chunk = " << max_chunk << "\n";
}

int main() {
    init();
    basic_functions_test();
    Yee3_test();
    Grid_3_test();
    decomposition_test();
    std::cout << "Program ended\n";

}