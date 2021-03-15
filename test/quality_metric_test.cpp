//
// Created by 徐溶延local on 2021/3/14.
//
#include <gtest/gtest.h>
#include "algorithm/energy/quality_metric.h"

using namespace energy;
class QualityMetricTest : public ::testing::Test {
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

    }

    RowMat32f tri1, tri2, tri3;
};

TEST_F(QualityMetricTest, qualityMetric) {
    double q1 = tri2dMetric(tri1);
    double q2 = tri2dMetric(tri2);
    double q3 = tri2dMetric(tri3);

    EXPECT_FLOAT_EQ(q1, 1.0);
    EXPECT_FLOAT_EQ(q2, 0.89151883);
    EXPECT_FLOAT_EQ(q3, 0.04327723);
}
