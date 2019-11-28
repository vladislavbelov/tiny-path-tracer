#include <cstdlib>
#include <fstream>

int main(int argc, char* argv[]) {
    int width = 640, height = 480;

    std::ofstream out("result.ppm");
    if (!out)
        return EXIT_FAILURE;
    out << "P3\n" << width << " " << height << "\n" << 255 << "\n";
    for (int y = 0; y < height; ++y)
        for (int x = 0; x < width; ++x)
            out << 255 << " " << 255 << " " << 255 << "\n";
    return EXIT_SUCCESS;
}

