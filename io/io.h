#pragma once

#include <string>
#include <vector>
#include "../common/types.h"

using namespace common;
namespace io {
//    void curve2vtk(const std::string &filename, const discrete_curve &curve);

    void points2vtk(const std::string &filename, const Eigen::Matrix3Xd &pts);

    void spts2vtk(const std::string &filename, const std::vector<surface_point> &spts);

    /**
     * read surface points from .spt file
     * @param filename
     * @param mesh
     * @param spt face_id barycentric_x barycentric_y
     */
    void readSpt(const std::string &filename,
                 const SurfaceMesh &mesh,
                 std::vector<surface_point> &spt);

    /**
     * read surface points from .gpt file
     * @param filename
     * @param gpts face_id x y z
     */
    void readGpt(const std::string &filename,
                 const SurfaceMesh &mesh,
                 std::vector<surface_point> &gpts);

    /**
     * read surface mesh from matrix
     * @param pts
     * @param tris
     * @param mesh
     * @return 1 success 0 failed
     */
    int readSurfaceMesh(SurfaceMesh &mesh,
                        const Eigen::Matrix3Xf &pts,
                        const Eigen::Matrix3Xi &tris);

    /**
     * write surface mesh from matrix
     * @param filename
     * @param pts
     * @param tris
     * @return 1 success 0 failed
     */
    int writeSurfaceMesh(const std::string &filename,
                         const Eigen::Matrix3Xf &pts,
                         const Eigen::Matrix3Xi &tris);

    void readIds(const std::string &filename, std::vector<size_t> &ids);

    int checkSuffix(const std::string &filename, const std::string &suffix);
}
