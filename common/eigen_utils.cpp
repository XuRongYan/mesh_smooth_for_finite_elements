//
// Created by 徐溶延local on 2021/3/13.
//

#include "eigen_utils.h"
#include <Eigen/Dense>

namespace common {
    RowVecXf mat2vec(RowMatf &mat) {
        return Eigen::Map<RowVecXf>(mat.data(), 1, mat.size());
    }
}