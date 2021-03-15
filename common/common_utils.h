//
// Created by 徐溶延 on 2020/11/25.
//

#ifndef MESH_CUTTING_COMMON_UTILS_H
#define MESH_CUTTING_COMMON_UTILS_H

#include <vector>
#include <map>
#include <Eigen/Dense>

namespace common {
    Eigen::MatrixX2i getLineIdx(const Eigen::MatrixXd &V, bool loop);

    bool insidePolygon(const std::vector<Eigen::Vector3f> &points,
                       const Eigen::Vector3f &test_point);

    template<typename K, typename V>
    std::map<V, K> inverseMap(const std::map<K, V> &mp) {
        std::map<V, K> res;
        for (auto it = mp.begin(); it != mp.end(); it++) {
            res[it->second] = it->first;
        }
        return res;
    }
} // namespace common


#endif //MESH_CUTTING_COMMON_UTILS_H
