//
// Created by 徐溶延local on 2021/3/13.
//
#include "common.h"
#include "common/math_utils.h"

using namespace common::math;

VecXf computeAvgCoord(const RowMatf &C) {
    const size_t d = C.cols();
    const size_t n = C.rows();
    VecXf avg_coord(d);
    avg_coord.setZero();
    for (size_t i = 0; i < n; i++) {
        avg_coord += C.row(i);
    }
    avg_coord /= n;
    return avg_coord;
}

RowMat32f computeD(const RowMat32f &C, const Vec2f &C_avg) {
    RowMat32f D = C;
    for (size_t i = 0; i < D.rows(); i++) {
        D.row(i) -= C_avg;
    }
    return D;
}

RowMat2f computeR(const RowMat2f &H) {
    RowMat2f R;
    R.setZero();
    Eigen::JacobiSVD<RowMat2f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    R = svd.matrixV() * svd.matrixU().transpose();
    if (fabs(svd.singularValues()(1)) < 1e-8) {
        R.setIdentity();
    }
    if (R.determinant() < 0) {
        auto U = svd.matrixU();
        U.col(1) *= -1;
        R = svd.matrixV() * U.transpose();
    }
    return R;
}

RowMat32f computeTm(const Vec2f &T) {
    RowMat32f Tm;
    for (size_t i = 0; i < 3; i++) {
        Tm.row(i) = T;
    }
    return Tm;
}

RowMatf computeKe(const RowMat32f &Ce, double t, double E, double v) {
    double d = detJ(Ce);
    double A = area(Ce);
    auto B = computeB(Ce);
    auto D = getStrain2StressProjectMatrix(E, v);
    return t * A * B.transpose() * D * B;
}

RowMat3f getStrain2StressProjectMatrix(double E, double v) {
    RowMat3f D;
    D << 1, v, 0,
         v, 1, 0,
         0, 0, (1 - v) * 0.5;
    return (E / (1 - v * v)) * D;
}

RowMat2f JacobianOf2dTriangle(const RowMat32f &Ce) {
    RowMat2f J;
    double x13 = Ce(0, 0) - Ce(2, 0);
    double x23 = Ce(1, 0) - Ce(2, 0);
    double y13 = Ce(0, 1) - Ce(2, 1);
    double y23 = Ce(1, 1) - Ce(2, 1);
    J << x13, y13,
         x23, y23;

    return J;
}

double detJ(const RowMat32f &Ce) {
    double x13 = Ce(0, 0) - Ce(2, 0);
    double x23 = Ce(1, 0) - Ce(2, 0);
    double y13 = Ce(0, 1) - Ce(2, 1);
    double y23 = Ce(1, 1) - Ce(2, 1);
    return x13 * y23 - y13 * x23;
}

RowMatf computeB(const RowMat32f &Ce) {
    double x13 = Ce(0, 0) - Ce(2, 0);
    double x21 = Ce(1, 0) - Ce(0, 0);
    double x23 = Ce(1, 0) - Ce(2, 0);
    double x32 = Ce(2, 0) - Ce(1, 0);
    double y12 = Ce(0, 1) - Ce(1, 1);
    double y13 = Ce(0, 1) - Ce(2, 1);
    double y23 = Ce(1, 1) - Ce(2, 1);
    double y31 = Ce(2, 1) - Ce(0, 1);
    double det_j = detJ(Ce);

    RowMatf B(3, 6);
    B << y23, 0, y31, 0, y12, 0,
         0, x32, 0, x13, 0, x21,
         x32, y23, x13, y31, x21, y12;
    return (1.0 / det_j) * B;
}


