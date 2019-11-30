#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>

template<typename T>
T clamp(const T& val, const T& min, const T& max) {
    return std::min(max, std::max(min, val));
}

struct Vector3 {
    float x, y, z; // Coordinates in R3.
};

struct Ray3 {
    Vector3 o, d; // Origin and direction respectively.
};

Vector3 operator+(const Vector3& v1, const Vector3& v2) {
    return {v1.x + v2.x, v1.y + v2.y, v1.z + v2.z};
}
Vector3 operator-(const Vector3& v1, const Vector3& v2) {
    return {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
}
Vector3 operator-(const Vector3& v) {
    return {-v.x, -v.y, -v.z};
}
Vector3 operator*(const Vector3& v, const float k) {
    return {v.x * k, v.y * k, v.z * k};
}
Vector3 operator/(const Vector3& v, const float k) {
    return {v.x / k, v.y / k, v.z / k};
}
Vector3 operator*(const Vector3& v1, const Vector3& v2) {
    return {v1.x * v2.x, v1.y * v2.y, v1.z * v2.z};
}
float dot(const Vector3& v1, const Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
float length(const Vector3& v) {
    return std::sqrt(dot(v, v));
}
Vector3 normalize(const Vector3& v) {
    const float len = length(v);
    return len > 0.0 ? v / len : v;
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

