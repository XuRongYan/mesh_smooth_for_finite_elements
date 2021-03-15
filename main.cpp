#include <iostream>
#include <fstream>
#include <filesystem>
#include <SurfaceMesh/SurfaceMesh.h>
#include <spdlog/spdlog.h>
#include "algorithm/MeshSmooth.h"
#include "common/string_utils.h"
#include "common/timer.h"

std::vector<std::string> folder_ordered_files(const std::string &path) {
    std::vector<std::string> files;
    if(!std::filesystem::exists(path)) {
        std::cerr << "path does not exists." << std::endl;
        return files;
    }
    std::filesystem::directory_entry entry(path);
    for(const auto &fs: std::filesystem::directory_iterator(entry)) {
        const auto file = fs.path().string();
        const auto ext = common::findSuffix(file);
        if(ext != ".obj"&&ext != ".off"&&ext != ".vtk"&&ext != ".mesh") continue;
        files.push_back(file);
    }
    std::sort(files.begin(), files.end());
    return files;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        exit(EXIT_FAILURE);
    }
    auto files = folder_ordered_files(DATA_SHARED_PATH"/models");
    spdlog::set_level(spdlog::level::debug);
    Timer<> timer;
    std::ofstream ofs(DATA_SHARED_PATH"/output/time.txt");
    for (const auto &file : files) {
        timer.beginStage("START SMOOTHING");
        SurfaceMesh mesh;
        mesh.read(file);
        //std::cout << "input:" << file << std::endl;
        alg::MeshSmooth meshSmooth(mesh);
        auto res = meshSmooth.smooth();
        std::string prefix = common::findPrefix(file);
        std::string suffix = common::findSuffix(file);
        std::string out_path = DATA_SHARED_PATH"/output/res_fem_" + prefix;
        res.write(out_path);
        timer.endStage("END SMOOTHING");
        ofs << prefix << " " << timer.value() << std::endl;
        //std::cout << "output:" << out_path << std::endl;
    }
    ofs.close();

//    std::string file_name = argv[1];
//    SurfaceMesh mesh;
//    mesh.read(DATA_SHARED_PATH"/models/" + file_name);
//
//    alg::MeshSmooth smoother(mesh);
//    auto res = smoother.smooth();
//    std::string prefix = common::findPrefix(file_name);
//    std::string suffix = common::findSuffix(file_name);
//    res.write(DATA_SHARED_PATH"/output/" + prefix + "_res" + suffix);
    return 0;
}
