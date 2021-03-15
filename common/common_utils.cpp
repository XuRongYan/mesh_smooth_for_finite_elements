//
// Created by 徐溶延 on 2020/11/25.
//

#include "common_utils.h"

namespace common {
Eigen::MatrixX2i getLineIdx(const Eigen::MatrixXd &V, bool loop) {
    assert(V.rows() > 1);
    size_t segments = loop ? V.rows() : V.rows() - 1;
    Eigen::MatrixX2i ids(segments, 2);
    for (size_t i = 0; i < segments; i++) {
        ids.row(i) << i, (i + 1) % V.rows();
    }
    return ids;
}

bool insidePolygon(const std::vector<Eigen::Vector3f> &points,
                   const Eigen::Vector3f &test_point) {
    bool inside = false;
    long n = points.size();
    for (long i = 0, j = n - 1; i < n; j = i++) {
        if (((points[i].y() > test_point.y()) != (points[j].y() > test_point.y()))
                && (test_point.x() < (points[j].x() - points[i].x()) * (test_point.y() - points[i].y()) /
                    (points[j].y() - points[i].y()) + points[i].x())) {
            inside = !inside;
        }
    }
    return inside;
}
} // namespace common


