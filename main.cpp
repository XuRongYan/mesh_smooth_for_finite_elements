#include <iostream>
#include <SurfaceMesh/SurfaceMesh.h>
#include "algorithm/MeshSmooth.h"
#include "common/string_utils.h"

int main(int argc, char *argv[]) {
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    std::string file_name = argv[1];
    SurfaceMesh mesh;
    mesh.read(DATA_SHARED_PATH"/models/" + file_name);
    alg::MeshSmooth smoother(mesh);
    auto res = smoother.smooth();
    std::string prefix = common::findPrefix(file_name);
    std::string suffix = common::findSuffix(file_name);
    res.write(DATA_SHARED_PATH"/output/" + prefix + "_res" + suffix);
    return 0;
}
