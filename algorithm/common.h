//
// Created by 徐溶延local on 2021/3/13.
//

#ifndef MESH_SMOOTH_FOR_FINITE_ELEMENTS_COMMON_H
#define MESH_SMOOTH_FOR_FINITE_ELEMENTS_COMMON_H

#endif //MESH_SMOOTH_FOR_FINITE_ELEMENTS_COMMON_H
#include <cmath>
#include "common/types.h"
#include "common/eigen_utils.h"

/**
 * compute scale coefficient
 * reference to Eq.22 in the paper
 * @param A triangle area
 * @return
 */
inline double compute2dS(double A) {
    return sqrt(A);
}

/**
 * reference to Eq.21, Eq.36
 * @param Cl local coordinate of an element
 * @param s  scale coefficient compute by Eq.22
 * @param alpha
 * @return
 */
inline RowMat32f computeRefLocalCoord(const RowMat32f &Cl, double s, double alpha = 1.0) {
    return alpha * s * Cl;
}

/**
 * reference to Eq.29
 * @param Cr local coordinate of an reference element
 * @param R  rotate matrix
 * @param Tm translate matrix
 * @return
 */
inline RowMat32f rotateAndTranslate(const RowMat32f &Cr, const RowMatf &R, const RowMat32f &Tm) {
    assert(R.rows() == 2 && R.cols() == 2);
    return Cr * R.transpose() + Tm;
}

/**
 * @param C     compute by Eq.29
 * @param C0    original triangle coordinate
 * @return Ue
 */
inline RowMat32f computeLocalDisplacement(const RowMat32f &C, const RowMat32f &C0) {
    return C - C0;
}

/**
 * reference to Eq.26 & Eq.27
 * @param C
 * @return
 */
RowVecXf computeAvgCoord(const RowMatf &C);

/**
 * reference to Eq.24 & Eq.25
 * @param C
 * @param C_avg
 * @return
 */
RowMat32f computeD(const RowMat32f &C, const RowVec2f &C_avg);

/**
 * reference to Eq.23
 * @param D
 * @param D_prime
 * @return
 */
inline RowMat2f computeH(const RowMat32f &D, const RowMat32f &D_prime) {
    return D.transpose() * D_prime;
}

/**
 * reference to Eq.14
 * @param H
 * @return
 */
RowMat2f computeR(const RowMat2f &H);

/**
 * reference to Eq.28
 * @param R
 * @param Ce_avg
 * @param Cr_avg
 * @return
 */
inline RowVec2f computeT(const RowMat2f &R, const RowVec2f &Ce_avg, const RowVec2f &Cr_avg) {
    return (Ce_avg.transpose() - R * Cr_avg.transpose()).transpose();
}

/**
 * reference to Eq.30
 * @param T
 * @return
 */
RowMat32f computeTm(const RowVec2f &T);

/**
 * reference to Eq.29
 * @param Cr
 * @param R
 * @param T
 * @return
 */
inline RowMat32f computeC(const RowMat32f &Cr, const RowMat2f &R, const RowVec2f &T) {
    auto Tm = computeTm(T);
    return Cr * R.transpose() + Tm;
}

/**
 * compute required displacements
 * @return
 */
inline RowMat32f computeUe(const RowMat32f &C, const RowMat32f &C0) {
    return C - C0;
}

/**
 * stiffness matrix
 * @param Ce
 * @param t
 * @param E
 * @param v
 * @return
 */
RowMatf computeKe(const RowMat32f &Ce, double t, double E, double v);

/**
 * reference to Eq.31
 * @param Ke
 * @param Ue
 * @param extend_approach
 * @return
 */
inline RowVecXf computeFe(const RowMatf &Ke, const RowMat32f &Ue) {
    assert(Ke.rows() == 6 && Ke.cols() == 6);
    RowMatf Utmp = Ue;
    RowVecXf U = common::mat2vec(Utmp);
    return (Ke * U.transpose()).transpose();
}

/**
 * reference to Eq.36
 * @param Ke
 * @param Ue
 * @param q
 * @param q_min
 * @return
 */
inline RowVecXf computeFe(const RowMat3f &Ke, const RowMat32f &Ue, double q, double q_min) {
    return (1 - q) / (q - q_min) * computeFe(Ke, Ue);
}

/**
 * transform strain to stress
 * reference to (4.14) in 2D Triangular Elements
 * @param E Young's modules
 * @param v Poisson's ratio
 * @return
 */
RowMat3f getStrain2StressProjectMatrix(double E, double v);

/**
 * Triangle Jacobian
 * @param Ce
 * @return
 */
RowMat2f JacobianOf2dTriangle(const RowMat32f &Ce);

/**
 * determinant of Jacobian
 * @param Ce
 * @return
 */
double detJ(const RowMat32f &Ce);

/**
 * reference to (4.61) in 2D Triangular Elements
 * @param Ce
 * @return
 */
RowMatf computeB(const RowMat32f &Ce);

inline bool updatable(double p_after, double p_before, double gamma = 0.8) {
    return p_after > gamma * p_before;
}








