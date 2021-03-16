//
// Created by 徐溶延local on 2021/3/14.
//

#include <gtest/gtest.h>
#include <Eigen/Sparse>
#include "test_utils.h"
#include "common/math_utils.h"
#include "algorithm/common.h"

using namespace common::math;
class CommonTest : public ::testing::Test {
protected:
    void SetUp() override {
        tri1 << 0, 0,
                1, 0,
                0.5, sqrt(3) / 2;
        tri2 << 0, 0,
                1, 0,
                1, 1;
        tri3 << 0, 0,
                6, 0,
                3, 0.1;
        tri4 << 5, 5,
                5, 4,
                6, 5;
        tri5 << 0, 0,
                3, 0,
                0, 3;
    }

    RowMat32f tri1, tri2, tri3, tri4, tri5;
};

TEST_F(CommonTest, localCoord2dTest) {
    RowMat32f ref4;
    ref4 << 0, 0,
            1, 0,
            0, 1;
    auto Cl1 = localCoord2d(tri1);
    auto Cl2 = localCoord2d(tri2);
    auto Cl3 = localCoord2d(tri3);
    auto Cl4 = localCoord2d(tri4);
    EXPECT_MAT_EQ(Cl1, tri1);
    EXPECT_MAT_EQ(Cl2, tri2);
    EXPECT_MAT_EQ(Cl3, tri3);
    EXPECT_MAT_EQ(Cl4, ref4);
}

TEST_F(CommonTest, avgCoordTest) {
    auto avg_coord = computeAvgCoord(tri5);
    RowVec2f ref_coord;
    ref_coord << 1, 1;
    EXPECT_MAT_EQ(avg_coord, ref_coord);
}

TEST_F(CommonTest, getStrain2StressProjectMatrixTest) {
    auto D = getStrain2StressProjectMatrix(1, 0);
    RowMat3f D_ref;
    D_ref << 1, 0, 0,
             0, 1, 0,
             0, 0, 0.5;
    EXPECT_MAT_EQ(D_ref, D);
}

TEST_F(CommonTest, SparseMatrixTest) {
    std::vector<Eigen::Triplet<float >> S;
    Eigen::SparseMatrix<double > G(5, 5);
    for (size_t i = 0; i < 5; i++)
        S.emplace_back(0, 0, 1);
    G.setFromTriplets(S.begin(), S.end());
    std::cout << G << std::endl;
}