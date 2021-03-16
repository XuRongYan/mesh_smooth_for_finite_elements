//
// Created by 徐溶延 on 2020/11/16.
//

#include <cmath>
#include "math_utils.h"

namespace common {
    namespace math {
        double angle(const RowVec3f &A, const RowVec3f &B) {
            return atan2(A.cross(B).norm(), A.dot(B));
        }

        double angle2D(const RowVec2f &A, const RowVec2f &B) {
            return atan2(cross2D(A, B).norm(), A.dot(B));
        }

        double angle(const RowVec3f &A, const RowVec3f &B, const RowVec3f &C) {
            return angle(B - A, C - A);
        }

        double area(const RowMat32f &e) {
            auto v1 = e.row(1) - e.row(0);
            auto v2 = e.row(2) - e.row(0);
            return 0.5 * (cross2D(v1, v2)).norm();
        }

        double perimeter(const RowMat32f &e) {
            double p = 0;
            for (size_t i = 0; i < 3; i++) {
                p += (e.row((i + 1) % 3) - e.row(i)).norm();
            }
            return p;
        }

        RowMat32f localCoord2d(const RowMat32f &e) {
            RowMat32f Cl;
            Cl.row(0) << 0, 0;
            const double l1 = (e.row(1) - e.row(0)).norm();
            Cl.row(1) << l1, 0;
            double theta = angle2D(e.row(1) - e.row(0), e.row(2) - e.row(0));
            const double l2 = (e.row(2) - e.row(0)).norm();
            Cl.row(2) << l2 * cos(theta), l2 * sin(theta);
            return Cl;
        }

        double cosOf2Vec(const RowVec3f &A, const RowVec3f &B) {
            return A.dot(B) / (A.norm() * B.norm());
        }

        RowVec3f cross2D(const RowVec2f &a, const RowVec2f &b) {
            RowVec3f v;
            v.setZero();
            v[2] = a.x() * b.y() - a.y() * b.x();
            return v;
        }

        int compare(const Eigen::Vector3f &A, const Eigen::Vector3f &B) {
            if (A.x() != B.x()) {
                return A.x() < B.x();
            }
            return A.y() < B.y();
        }


        bool
        intersect(const Eigen::Vector3f &A,
                  const Eigen::Vector3f &B,
                  const Eigen::Vector3f &C,
                  const Eigen::Vector3f &D) {
            const auto ab = B - A;
            const auto ad = D - A;
            const auto ac = C - A;
            const auto cd = D - C;
            const auto ca = A - C;
            const auto cb = B - C;
            const auto acab = ac.cross(ab);
            const auto adab = ad.cross(ab);
            const auto cacd = ca.cross(cd);
            const auto cbcd = cb.cross(cd);
            const auto abcd = ab.cross(cd);
            float side_cd = acab.dot(adab);
            float side_ab = cacd.dot(cbcd);

            // normal condition
            if (side_ab < 0 && side_cd < 0) {
                return true;
            }
            // only one point intersected
            if ((side_ab < 0 && side_cd == 0) || (side_cd < 0 && side_ab == 0)) {
                return true;
            }
            // overlap
            if (abcd.norm() == 0) {
                const auto abac = ab.cross(ac);
                if (abac.norm() == 0) {
                    Eigen::Vector3f s1, s2, e1, e2;
                    s1 = A, e1 = B, s2 = C, e2 = D;
                    if (!compare(A, B))
                        std::swap(s1, e1);
                    if (!compare(C, D))
                        std::swap(s2, e2);
                    if (compare(s1, s2) && compare(e2, e1)) {
                        return true;
                    }
                }
            }
            return false;
        }

        double distanceToEdge(const Eigen::Vector3f &p,
                              const Eigen::Vector3f &A,
                              const Eigen::Vector3f &B) {
            Eigen::Vector3f ap = p - A;
            Eigen::Vector3f ab = B - A;
            double d = ap.cross(ab).norm() / ab.norm();
            return d;
        }

        long mod(long i, long size) {
            if (i >= 0) return i % size;
            return size + i;
        }
    } // namespace math
} // namespace common

