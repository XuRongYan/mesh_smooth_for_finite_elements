//
// Created by 徐溶延local on 2021/3/13.
//

#include "eigen_utils.h"
#include <Eigen/Dense>

namespace common {
    VecXf mat2vec(RowMatf &mat) {
        return Eigen::Map<VecXf>(mat.data(), 1, mat.size());
    }
}