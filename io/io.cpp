
#include "io.h"
#include <fstream>
#include <algorithm>
//#include <spdlog/spdlog.h>
#include "../common/string_utils.h"
#include "../common/surface_mesh_utils.h"
#include "../common/util.h"
#include "vtk_utils.h"
#include "vtk.h"

namespace io {
//    void curve2vtk(const std::string &filename, const discrete_curve &curve) {
//        if (!checkSuffix(filename, ".vtk"))
//            return;
//        std::ofstream ofs(filename);
//        if (!ofs.is_open()) {
//            //spdlog::error("open output file failed. path = {}", filename);
//            return;
//        }
//        Eigen::Matrix3Xd V(3, curve.dim());
//        for (size_t i = 0; i < curve.dim(); i++) {
//            auto pos = curve[i].pos_;
//            V.col(i) << pos[0], pos[1], pos[2];
//        }
//        writeLines(ofs, V, curve.loop_);
//        //spdlog::debug("write file to {}", filename);
//    }

    void points2vtk(const std::string &filename,
                   const Eigen::Matrix3Xd &pts) {
        if (!checkSuffix(filename, ".vtk"))
            return;
        std::ofstream ofs(filename);
        if (!ofs.is_open()) {
            //spdlog::error("open output file failed. path = {}", filename);
            return;
        }
        std::vector<int> ids(pts.cols());
        for (int i = 0; i < pts.cols(); i++) {
            ids[i] = i;
        }
        point2vtk(ofs, pts.data(), pts.cols(), ids.data(), ids.size());
    }

    void spts2vtk(const std::string &filename,
                  const std::vector<surface_point> &spts) {
        Eigen::Matrix3Xd pts(3, spts.size());
        for (size_t i = 0; i < spts.size(); i++) {
            const auto &pos = spts[i].pos_;
            pts.col(i) << pos.x(), pos.y(), pos.z();
        }
        points2vtk(filename, pts);
    }

    void readSpt(const std::string &filename,
                 const SurfaceMesh &mesh,
                 std::vector<surface_point> &spts) {
        spts.clear();
        if (!checkSuffix(filename, ".spt"))
            return;
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            return;
        }
        size_t n = 0;
        ifs >> n;
        spts.reserve(n);
        for (size_t i = 0; i < n; i++) {
            int fid;
            float x, y;
            // face_id barycentric
            ifs >> fid >> x >> y;
            Eigen::Vector3f bary(x, y, 1 - x - y);
            if (!isBaryValid(bary)) {
                abort();
            }
            spts.emplace_back(mesh, x, y, fid);
        }
    }

    void readGpt(const std::string &filename,
                     const SurfaceMesh &mesh,
                     std::vector<surface_point> &gpts) {
        gpts.clear();
        if (!checkSuffix(filename, ".gpt"))
            return;
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            return;
        }
        size_t n = 0;
        ifs >> n;
        gpts.reserve(n);
        for (size_t i = 0; i < n; i++) {
            size_t fid;
            float x, y, z;
            ifs >> fid >> x >> y >> z;
            Point p(x, y, z);
            gpts.emplace_back(mesh, p, fid);
        }
    }

    int readSurfaceMesh(SurfaceMesh &mesh,
                        const Eigen::Matrix3Xf &pts,
                        const Eigen::Matrix3Xi &tris) {
        return buildSurfaceMesh(mesh, pts, tris);
    }

    int writeSurfaceMesh(const std::string &filename,
                         const Eigen::Matrix3Xf &pts,
                         const Eigen::Matrix3Xi &tris) {
        SurfaceMesh mesh;
        bool success = buildSurfaceMesh(mesh, pts, tris);
        if (success) {
            return mesh.write(filename);
        }
        return 0;
    }

    void readIds(const std::string &filename, std::vector<size_t> &ids) {
        if (!checkSuffix(filename, ".pid"))
            return;
        std::ifstream ifs(filename);
        if (!ifs.is_open()) {
            //spdlog::error("open input file failed. path = {}", filename);
            return;
        }
        while (!ifs.eof()) {
            long long id;
            ifs >> id;
            if (id < 0) {
                //spdlog::error("invalid id: {}", id);
                return;
            }
            ids.emplace_back((size_t) id);
        }
    }

    int checkSuffix(const std::string &filename, const std::string &suffix) {
        std::string real_suffix = findSuffix(filename);
        std::string suffix1(real_suffix), suffix2(real_suffix);
        std::transform(suffix.begin(), suffix.end(), suffix1.begin(), ::toupper);
        std::transform(suffix.begin(), suffix.end(), suffix2.begin(), ::tolower);
        if (real_suffix != suffix1 && real_suffix != suffix2) {
            //spdlog::error("the file is not a {} file. path = {}", suffix1, filename);
            return 0;
        }
        return 1;
    }
}
