#ifndef COMMON_UTIL_H
#define COMMON_UTIL_H

#include <Eigen/Dense>

namespace common {
static const float eps = 1.1e-6;

int barycentric(const Eigen::Vector3f &a, const Eigen::Vector3f &b,
                const Eigen::Vector3f &c, const Eigen::Vector3f &p,
                Eigen::Vector3f &bary);

bool isBaryValid(Eigen::Vector3f &bary);

} // namespace common

#endif //COMMON_UTIL_H
