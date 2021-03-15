//
// Created by 徐溶延 on 2020/11/16.
//

#ifndef MESH_CUTTING_MATH_UTILS_H
#define MESH_CUTTING_MATH_UTILS_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <Eigen/Dense>
#include "types.h"

namespace common {
    namespace math {
        double angle(const Vec3f &A, const Vec3f &B);

        double angle2D(const Vec2f &A, const Vec2f &B);

        double angle(const Vec3f &A, const Vec3f &B, const Vec3f &C);

        double area(const RowMat32f &e);

        double perimeter(const RowMat32f &e);

        RowMat32f localCoord2d(const RowMat32f &e);

        double cosOf2Vec(const Vec3f &A, const Vec3f &B);

        Vec3f cross2D(const Vec2f &a, const Vec2f &b);

        bool intersect(const Eigen::Vector3f &A,
                       const Eigen::Vector3f &B,
                       const Eigen::Vector3f &C,
                       const Eigen::Vector3f &D);

        double distanceToEdge(const Eigen::Vector3f &p, const Eigen::Vector3f &A, const Eigen::Vector3f &B);

        long mod(long i, long size);

        inline size_t prev(size_t i, size_t n) {
            return i == 0 ? n - 1 : i - 1;
        }

        inline size_t next(size_t i, size_t n) {
            return (i + 1) % n;
        }

        inline size_t prev(size_t i, size_t step, size_t n) {
            return step == 1 ? prev(i, n) : prev(prev(i, n), step - 1, n);
        }

        inline size_t next(size_t i, size_t step, size_t n) {
            return step == 1 ? next(i, n) : next(next(i, n), step - 1, n);
        }
    } // namespace math
} // namespace common




#endif //MESH_CUTTING_MATH_UTILS_H
