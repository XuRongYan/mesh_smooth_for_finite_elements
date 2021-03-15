//
// Created by 徐溶延local on 2021/3/14.
//

#ifndef MESH_SMOOTH_FOR_FINITE_ELEMENTS_TEST_UTILS_H
#define MESH_SMOOTH_FOR_FINITE_ELEMENTS_TEST_UTILS_H

#endif //MESH_SMOOTH_FOR_FINITE_ELEMENTS_TEST_UTILS_H
#include <gtest/gtest.h>
#include "common/types.h"

inline void EXPECT_MAT_EQ(const RowMatf &mat1, const RowMatf &mat2) {
    EXPECT_EQ(mat1.rows(), mat2.rows());
    EXPECT_EQ(mat1.cols(), mat2.cols());
    for (size_t i = 0; i < mat1.rows(); i++)
        for (size_t j = 0; j < mat1.cols(); j++)
            EXPECT_NEAR(mat1(i, j), mat2(i, j), 1e-8);
}