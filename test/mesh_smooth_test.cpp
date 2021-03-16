//
// Created by 徐溶延local on 2021/3/14.
//
#include <gtest/gtest.h>
#include <string>
#include "algorithm/energy/quality_metric.h"
#include "common/surface_mesh_utils.h"
#include "common/math_utils.h"

#define private public
#define protected public

#include "algorithm/MeshSmooth.h"

#undef private
#undef protected
using namespace alg;
using namespace common;
using namespace common::math;
using namespace energy;

class MeshSmoothTest : public ::testing::Test {
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

TEST_F(MeshSmoothTest, localCheckTest) {
    SurfaceMesh mesh;
    mesh.read(DATA_SHARED_PATH"/models/etri.obj");
    MeshSmooth smoother(mesh);
    for (const auto &f : mesh.faces()) {
        auto vid = vertexAroundFace(mesh, f);
        auto C0 = get2DFaceCoord(mesh, f);
        auto mesh_c0 = triangle2mesh(C0);
        double A = area(C0);
        auto Cl = smoother.Cl_;
        double s = compute2dS(A);
        auto Cr = computeRefLocalCoord(Cl, s, smoother.alpha_);
        RowMat2f R;
        RowVec2f T;
        smoother.computeRandTMatrices(C0, Cr, R, T);
        auto C = computeC(Cr, R, T);                    // Eq.29
        auto mesh_c = triangle2mesh(C);
        mesh_c0.write(DATA_SHARED_PATH"/output/C0_" + std::to_string(f.idx()) + ".obj");
        mesh_c.write(DATA_SHARED_PATH"/output/C_" + std::to_string(f.idx()) + ".obj");
        double q = tri2dMetric(C);
        double Ac = area(C);
        EXPECT_NEAR(Ac, A, 1e-3);
        EXPECT_FLOAT_EQ(q, 1);
    }
}



