#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <limits>
#include <random>
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

struct Material {
    Vector3 c; // Color.
    Vector3 e; // Emission.
};

struct Sphere {
    Vector3 p; // Position.
    float r; // Radius.
    size_t m; // Material id.
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
Vector3 reflect(const Vector3& v, const Vector3& n) {
    return v - n * (dot(n, v) * 2);
}

std::istream& operator>>(std::istream& in, Vector3& v) {
    return in >> v.x >> v.y >> v.z;
}
std::ostream& operator<<(std::ostream& out, const Vector3& v) {
    return out << "(" << v.x << "," << v.y << "," << v.z << ")";
}
std::istream& operator>>(std::istream& in, Material& material) {
    return in >> material.c >> material.e;
}
std::istream& operator>>(std::istream& in, Sphere& sphere) {
    return in >> sphere.p >> sphere.r >> sphere.m;
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

std::vector<Material> materials;
std::vector<Sphere> spheres;

bool FindIntersection(const Ray3& ray, Vector3& x, Vector3& normal, size_t& material_id) {
    Vector3 tempX, tempNormal;
    bool found = false;
    for (const Sphere& sphere : spheres) {
        if (!Intersect(ray, sphere, tempX, tempNormal))
            continue;
        if (dot(ray.d, tempX) >= dot(ray.d, x))
            continue;
        x = tempX;
        normal = tempNormal;
        material_id = sphere.m;
        found = true;
    }
    return found;
}

const float PI = acosf(-1.0f);

const Vector3 GetRandomDirectionOnSphere() {
    static std::mt19937 gen;
    static std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    const float theta = 2.0f * PI * dist(gen);
    const float phi = acosf(1.0f - 2.0f * dist(gen));
    const float sin_phi = sinf(phi);
    return Vector3{sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta)};
}

Vector3 Trace(const Ray3& ray, int depth = 0) {
    Vector3 x = ray.o + ray.d * 1e7f, normal;
    size_t material_id;
    if (depth > 4 || !FindIntersection(ray, x, normal, material_id))
        return Vector3{0, 0, 0};
    Vector3 dir = GetRandomDirectionOnSphere();
    float cos_theta = dot(dir, normal);
    if (cos_theta < 0.0f) {
        cos_theta -= cos_theta;
        dir = reflect(dir, normal);
    }
    const Vector3 incoming_light = Trace(Ray3{x + dir * 1e-4f, dir}, depth + 1);
    const float pdf = 2.0f * PI;
    return materials[material_id].e +
               materials[material_id].c * incoming_light * cos_theta / pdf;
}

void Render(int x0, int y0, int x1, int y1, int width, int height, int spp,
            std::vector<Vector3>& out) {
    const float aspect_ratio = 1.0f * width / height;
    for (int y = y0; y < y1; ++y)
        for (int x = x0, index = x0 + y * width; x < x1; ++x, ++index) {
            const float u = 1.0f * x / width - 0.5f;
            const float v = (0.5f - 1.0f * y / height) / aspect_ratio;
            const Ray3 ray{Vector3{0, 0, -1}, normalize(Vector3{u, v, 1})};
            for (int sample = 0; sample < spp; ++sample)
                out[index] += Trace(ray);
            out[index] *= 1.0f / spp;
        }
}

void OutputColor(std::ostream& out, const Vector3& color) {
    out << clamp<int>(color.x * 255, 0, 255) << " "
        << clamp<int>(color.y * 255, 0, 255) << " "
        << clamp<int>(color.z * 255, 0, 255) << "\n";
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
    if (width < 1 || height < 1 || spp < 1)
        return EXIT_FAILURE;
    scene >> materials >> spheres;

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

