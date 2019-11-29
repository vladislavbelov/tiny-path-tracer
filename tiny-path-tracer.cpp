#include <cstdlib>
#include <fstream>
#include <vector>

struct Vector3 {
    float x, y, z; // Coordinates in R3.
};

template<typename T>
T clamp(const T& val, const T& min, const T& max) {
    return std::min(max, std::max(min, val));
}

void OutputLinearColor(std::ostream& out, const Vector3& v) {
    out << clamp<int>(v.x * 255, 0, 255) << " "
        << clamp<int>(v.y * 255, 0, 255) << " "
        << clamp<int>(v.z * 255, 0, 255) << "\n";
}

int main(int argc, char* argv[]) {
    int width = 640, height = 480;

    std::ofstream out("result.ppm");
    if (!out)
        return EXIT_FAILURE;
    out << "P3\n" << width << " " << height << "\n" << 255 << "\n";
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            OutputLinearColor(out, {1.0, 1.0, 1.0});
    return EXIT_SUCCESS;
}

