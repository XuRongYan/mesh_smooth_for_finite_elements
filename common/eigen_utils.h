//
// Created by 徐溶延local on 2021/3/13.
//

#ifndef MESH_SMOOTH_FOR_FINITE_ELEMENTS_EIGEN_UTILS_H
#define MESH_SMOOTH_FOR_FINITE_ELEMENTS_EIGEN_UTILS_H

#include "types.h"

namespace common {
    RowVecXf mat2vec(RowMatf &mat);

    inline RowVecXf col2row(const VecXf &vec) {
        return vec.transpose();
    }

    inline VecXf row2col(const RowVecXf &vec) {
        return vec.transpose();
    }
}


#endif //MESH_SMOOTH_FOR_FINITE_ELEMENTS_EIGEN_UTILS_H
