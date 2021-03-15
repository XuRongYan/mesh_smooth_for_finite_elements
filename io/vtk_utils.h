//
// Created by 徐溶延 on 2020/11/18.
//

#ifndef MESH_CUTTING_VTK_UTILS_H
#define MESH_CUTTING_VTK_UTILS_H
#include <Eigen/Dense>
#include <vector>
#include "vtk.h"

namespace io {
template<typename OS>
void writeLines(OS& os, const Eigen::Matrix3Xd &V, bool loop) {
    std::vector<int> line_idx;
    long seg = V.cols() - 1;
    if (loop) seg++;
    line_idx.reserve(seg);
    for (long i = 0; i < seg; i++) {
        line_idx.emplace_back(i);
        line_idx.emplace_back((i + 1) % V.cols());
    }
    line2vtk(os, V.data(), V.cols(), line_idx.data(), seg);
}

template<typename OS>
void writePoints(OS& os, const Eigen::Matrix3Xd &V) {
    std::vector<int> points(V.cols());
    for (int i = 0; i < points.size(); i++) {
        points[i] = i;
    }
    point2vtk(os, V.data(), V.cols(), points.data(), points.size());
}
}


#endif //MESH_CUTTING_VTK_UTILS_H
