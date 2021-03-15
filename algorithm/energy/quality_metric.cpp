//
// Created by 徐溶延local on 2021/3/13.
//

#include "quality_metric.h"
#include "common/math_utils.h"

using namespace common::math;
const double beta = 2;
namespace energy {
    /**
     * reference to Eq.1
     * @param e
     * @return
     */
    double tri2dMetric(const RowMat32f &e) {
        double A = area(e);
        double p = perimeter(e);
        double pr = 6 * sqrt(A / sqrt(3));
        double q2d = pow(pr / p, beta);
        return q2d;
    }
}