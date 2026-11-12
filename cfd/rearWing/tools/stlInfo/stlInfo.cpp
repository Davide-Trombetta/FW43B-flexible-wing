// ============================================================================
//  stlInfo.cpp  --  geometric summary of an STL surface
//
//  Reads an ASCII (multi-solid) or binary STL file and prints, for every
//  solid (region) and for the whole surface:
//    - number of triangles and wetted area
//    - bounding box and dimensions (length, width, height)
//    - frontal area projected on the y-z plane (seen by a flow along x)
//    - planform area projected on the x-y plane (seen from above)
//
//  Why a dedicated tool: OpenFOAM's surfaceCheck prints only the global
//  bounding box, in a format meant for people. stlInfo also gives the
//  bounding box of each region and the exact projected frontal area
//  (needed for the blockage ratio), and with -shell it prints KEY=value
//  lines that a bash script can load with eval.
//  Plain C++17, no OpenFOAM library: the same source builds with any
//  OpenFOAM distribution (the OpenFOAM-API version of the Ahmed-body case
//  used argList::optionFound and triFace::area of OpenFOAM.org, which
//  ESI v2306 does not have).
//
//  Build:  g++ -O2 -std=c++17 -o stlInfo stlInfo.cpp
//  Usage:  stlInfo file.stl            human-readable table
//          stlInfo -shell file.stl     KEY=value lines (for bash)
// ============================================================================

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using Vec = std::array<double, 3>;

struct Region
{
    explicit Region(std::string n) : name(std::move(n)) {}

    std::string name;
    Vec min{1e300, 1e300, 1e300};
    Vec max{-1e300, -1e300, -1e300};
    std::size_t nTri = 0;
    double area = 0.0;          // wetted area
    double frontalArea = 0.0;   // area projected on the y-z plane
    double planformArea = 0.0;  // area projected on the x-y plane

    void add(const Vec& a, const Vec& b, const Vec& c)
    {
        for (const Vec* p : {&a, &b, &c})
            for (int i = 0; i < 3; ++i)
            {
                min[i] = std::min(min[i], (*p)[i]);
                max[i] = std::max(max[i], (*p)[i]);
            }
        // triangle area vector = 0.5 * (b - a) x (c - a)
        const Vec u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const Vec v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const Vec n{u[1] * v[2] - u[2] * v[1],
                    u[2] * v[0] - u[0] * v[2],
                    u[0] * v[1] - u[1] * v[0]};
        area += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        // For a surface enclosing a volume, every line along x crosses it
        // twice: the projected (frontal) area is half the sum of |A_x|.
        // Faces lying on a symmetry plane (normal along y) do not contribute,
        // so the formula also holds for the half-model.
        frontalArea += 0.25 * std::fabs(n[0]);
        planformArea += 0.25 * std::fabs(n[2]);   // same rule along z
        ++nTri;
    }

    void merge(const Region& r)
    {
        for (int i = 0; i < 3; ++i)
        {
            min[i] = std::min(min[i], r.min[i]);
            max[i] = std::max(max[i], r.max[i]);
        }
        nTri += r.nTri;
        area += r.area;
        frontalArea += r.frontalArea;
        planformArea += r.planformArea;
    }
};

std::vector<Region> readSTL(const std::string& file)
{
    std::ifstream in(file, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open " + file);

    in.seekg(0, std::ios::end);
    const auto size = static_cast<std::uint64_t>(in.tellg());
    in.seekg(0);

    // Binary STL: 80-byte header, uint32 triangle count, 50 bytes/triangle
    if (size >= 84)
    {
        char header[80];
        std::uint32_t n = 0;
        in.read(header, 80);
        in.read(reinterpret_cast<char*>(&n), 4);
        if (84 + 50ull * n == size)
        {
            std::vector<Region> regions{Region("solid")};
            for (std::uint32_t t = 0; t < n; ++t)
            {
                float f[12];
                std::uint16_t attribute;
                in.read(reinterpret_cast<char*>(f), 48);
                in.read(reinterpret_cast<char*>(&attribute), 2);
                regions[0].add({f[3], f[4], f[5]}, {f[6], f[7], f[8]},
                               {f[9], f[10], f[11]});
            }
            return regions;
        }
        in.seekg(0);
    }

    // ASCII STL, possibly with several named solids
    std::vector<Region> regions;
    std::vector<Vec> pts;
    std::string line, word;
    while (std::getline(in, line))
    {
        std::istringstream is(line);
        if (!(is >> word)) continue;
        if (word == "solid")
        {
            std::string name = "solid" + std::to_string(regions.size());
            is >> name;
            regions.emplace_back(name);
        }
        else if (word == "vertex")
        {
            Vec p;
            is >> p[0] >> p[1] >> p[2];
            pts.push_back(p);
            if (pts.size() == 3)
            {
                if (regions.empty()) regions.emplace_back("solid0");
                regions.back().add(pts[0], pts[1], pts[2]);
                pts.clear();
            }
        }
    }
    if (regions.empty()) throw std::runtime_error("no triangles in " + file);
    return regions;
}

void printTable(const std::vector<Region>& regions, const Region& all)
{
    std::printf("%-10s %9s %10s %10s   %-33s %-33s\n", "region", "triangles",
                "area[m2]", "Afront[m2]", "min (x y z) [m]", "max (x y z) [m]");
    auto row = [](const Region& r) {
        std::printf("%-10s %9zu %10.5f %10.5f   (%9.5f %9.5f %9.5f)   (%9.5f %9.5f %9.5f)\n",
                    r.name.c_str(), r.nTri, r.area, r.frontalArea,
                    r.min[0], r.min[1], r.min[2], r.max[0], r.max[1], r.max[2]);
    };
    for (const auto& r : regions) row(r);
    row(all);
    std::printf("\nlength L (x) = %.5f m   width (y) = %.5f m   height (z) = %.5f m\n",
                all.max[0] - all.min[0], all.max[1] - all.min[1],
                all.max[2] - all.min[2]);
}

void printShell(const std::vector<Region>& regions, const Region& all)
{
    auto out = [](const std::string& prefix, const Region& r) {
        const char* c = "XYZ";
        for (int i = 0; i < 3; ++i)
        {
            std::printf("%s%cMIN=%.8g\n", prefix.c_str(), c[i], r.min[i]);
            std::printf("%s%cMAX=%.8g\n", prefix.c_str(), c[i], r.max[i]);
        }
        std::printf("%sAREA=%.8g\n%sAFRONT=%.8g\n%sAPLAN=%.8g\n", prefix.c_str(),
                    r.area, prefix.c_str(), r.frontalArea, prefix.c_str(),
                    r.planformArea);
    };
    std::printf("STL_REGIONS=\"");
    for (std::size_t i = 0; i < regions.size(); ++i)
        std::printf("%s%s", i ? " " : "", regions[i].name.c_str());
    std::printf("\"\n");
    out("STL_", all);
    std::printf("STL_L=%.8g\nSTL_W=%.8g\nSTL_H=%.8g\n", all.max[0] - all.min[0],
                all.max[1] - all.min[1], all.max[2] - all.min[2]);
    for (const auto& r : regions) out("STL_" + r.name + "_", r);
}

int main(int argc, char* argv[])
{
    bool shell = false;
    std::string file;
    for (int i = 1; i < argc; ++i)
    {
        const std::string a = argv[i];
        if (a == "-shell") shell = true;
        else if (a == "-h" || a == "-help")
        {
            std::cout << "Usage: stlInfo [-shell] file.stl\n";
            return 0;
        }
        else file = a;
    }
    if (file.empty())
    {
        std::cerr << "Usage: stlInfo [-shell] file.stl\n";
        return 1;
    }

    try
    {
        const auto regions = readSTL(file);
        Region all("TOTAL");
        for (const auto& r : regions) all.merge(r);
        if (shell) printShell(regions, all);
        else printTable(regions, all);
    }
    catch (const std::exception& e)
    {
        std::cerr << "stlInfo: " << e.what() << "\n";
        return 1;
    }
    return 0;
}


