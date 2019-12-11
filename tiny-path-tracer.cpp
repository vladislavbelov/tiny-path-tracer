#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
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

struct Sphere {
    Vector3 p; // Position.
    float r; // Radius.
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
Vector3& operator+=(Vector3& v1, const Vector3& v2) {
    v1.x += v2.x; v1.y += v2.y; v1.z += v2.z;
    return v1;
}
Vector3& operator*=(Vector3& v, const float k) {
    v.x *= k; v.y *= k; v.z *= k;
    return v;
}
float dot(const Vector3& v1, const Vector3& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
float length(const Vector3& v) {
    return std::sqrt(dot(v, v));
}
Vector3 normalize(const Vector3& v) {
    const float len = length(v);
    return len > 0.0f ? v / len : v;
}

std::istream& operator>>(std::istream& in, Vector3& v) {
    return in >> v.x >> v.y >> v.z;
}
std::ostream& operator<<(std::ostream& out, const Vector3& v) {
    return out << "(" << v.x << "," << v.y << "," << v.z << ")";
}
std::istream& operator>>(std::istream& in, Sphere& sphere) {
    return in >> sphere.p >> sphere.r;
}
template<typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& vector) {
    size_t size;
    in >> size;
    if (!in || size == std::numeric_limits<size_t>::max())
        return in;
    vector.resize(size);
    for (T& item : vector)
        in >> item;
    return in;
}

bool Intersect(const Ray3& ray, const Sphere& sphere, Vector3& x, Vector3& normal) {
    const Vector3 k = ray.o - sphere.p;
    const float b = dot(k, ray.d);
    const float c = dot(k, k) - sphere.r * sphere.r;
    if (b > 0.0f && c > 0.0f)
        return false; // The origin is outside of the sphere and pointing away.
    const float d = b * b - c;
    if (d < 0.0f)
        return false; // The ray doesn't intersect the sphere.
    const float ds = std::sqrt(d);
    const float t0 = -b - ds;
    const float t1 = -b + ds;
    x = ray.o + ray.d * (t0 < 0.0f ? t1 : t0);
    normal = normalize(x - sphere.p);
    return true;
}

std::vector<Sphere> spheres;

bool FindIntersection(const Ray3& ray, Vector3& x, Vector3& normal) {
    Vector3 tempX, tempNormal;
    bool found = false;
    for (const Sphere& sphere : spheres) {
        if (!Intersect(ray, sphere, tempX, tempNormal))
            continue;
        if (dot(ray.d, tempX) >= dot(ray.d, x))
            continue;
        x = tempX;
        normal = tempNormal;
        found = true;
    }
    return found;
}

Vector3 Trace(const Ray3& ray) {
    Vector3 x = ray.o + ray.d * 1e7f, normal;
    if (!FindIntersection(ray, x, normal))
        return Vector3{0, 0, 0};
    return normal;
}

void Render(int x0, int y0, int x1, int y1, int width, int height, int spp,
            std::vector<Vector3>& out) {
    const float aspect_ratio = 1.0f * width / height;
    for (int y = y0; y < y1; ++y)
        for (int x = x0, index = x0 + y * width; x < x1; ++x, ++index) {
            const float u = 1.0f * x / width - 0.5f;
            const float v = (1.0f * y / height - 0.5f) / aspect_ratio;
            Ray3 ray{Vector3{0, 0, -1}, normalize(Vector3{u, v, 1})};
            for (int sample = 0; sample < spp; ++sample)
                out[index] += (Trace(ray) + Vector3{1, 1, 1}) / 2;
            out[index] *= 1.0f / spp;
        }
}

void OutputColor(std::ostream& out, const Vector3& v) {
    out << clamp<int>(v.x * 255, 0, 255) << " "
        << clamp<int>(v.y * 255, 0, 255) << " "
        << clamp<int>(v.z * 255, 0, 255) << "\n";
}

int main(int argc, char* argv[]) {
    std::string scene_path = argc > 1 ? argv[1] : "scenes/default.txt";
    std::ifstream scene_raw_in(scene_path); // Read scene from the file.
    if (!scene_raw_in)
        return EXIT_FAILURE;
    std::stringstream scene;
    for (std::string line; getline(scene_raw_in, line);) {
        // Trim leading spaces and skip comments.
        size_t start_pos = line.find_first_not_of(" \t\r");
        if (start_pos != std::string::npos)
            line = line.substr(start_pos);
        if (!line.empty() && line.front() != '#' && line.front() != ';')
            scene << line << '\n';
    }
    scene_raw_in.close();

    int width = 640, height = 480, spp = 1;
    scene >> width >> height >> spp;
    scene >> spheres;

    std::vector<Vector3> pixels(width * height);
    Render(0, 0, width, height, width, height, spp, pixels);

    std::ofstream out("result.ppm");
    if (!out)
        return EXIT_FAILURE;
    out << "P3\n" << width << " " << height << "\n" << 255 << "\n";
    for (const Vector3& pixel_color : pixels)
        OutputColor(out, pixel_color);
    return EXIT_SUCCESS;
}

