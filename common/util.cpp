#include <iostream>
#include "util.h"


namespace common {
    int barycentric(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                    const Eigen::Vector3f &c, const Eigen::Vector3f &p,
                    Eigen::Vector3f &bary) {
        const double BARY_EPS = 1e-6;
        Eigen::Vector3f u = b - a;
        Eigen::Vector3f v = c - a;
        Eigen::Vector3f w = p - a;

        Eigen::Vector3f vw = v.cross(w);
        Eigen::Vector3f vu = v.cross(u);
        double vu_norm = vu.norm();

        if (vu_norm < BARY_EPS) {
            std::cerr << "degenerated triangle " << vu_norm << std::endl;
            return 0;
        }
        vu /= vu_norm;

        // test sign of r
        if (vw.dot(vu) < -1e-4) {
            std::cerr << "barycentric error vw.dot(vu): " << vw.dot(vu) << std::endl;
            return 0;
        }

        Eigen::Vector3f uw = u.cross(w);
        Eigen::Vector3f uv = -vu;

        // test sign of t
        if (uw.dot(uv) < -1e-4) {
            std::cerr << "barycentric error uw.dot(uv): " << uw.dot(uv) << std::endl;
            return 0;
        }

        bary[1] = vw.norm() / vu_norm;
        bary[2] = uw.norm() / vu_norm;
        bary[0] = 1 - bary[1] - bary[2];

        return 1;
        //if (bary[0] >= -10 * 1e-6) {
        //    return 1;
        //}
        //else {
        //    std::cerr << "barycentric error: " << bary[0] << ' ' << bary[1] << ' ' << bary[2] << std::endl;
        //    return 0;
        //}
    }

    bool isBaryValid(Eigen::Vector3f &bary) {
        if (bary[0] >= -10 * eps && bary[1] >= -10 * eps && bary[2] >= -10 * eps
            && bary[0] <= 1 + 10 * eps && bary[2] <= 1 + 10 * eps && bary[2] <= 1 + 10 * eps) {
            if (fabs(bary.sum() - 1) > 10 * eps)
                return false;
            // refine boundary values
            for (size_t i = 0; i < bary.size(); i++) {
                if (bary[i] <= 10 * eps) {
                    bary[i] = 0;
                } else if (fabs(bary[i] - 1) <= 10 * eps) {
                    bary[i] = 1;
                }
            }
            return true;
        }
        return false;
    }
} // namespace common
