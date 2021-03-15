//
// Created by 徐溶延local on 2021/3/13.
//

#ifndef MESH_SMOOTH_FOR_FINITE_ELEMENTS_MESHSMOOTH_H
#define MESH_SMOOTH_FOR_FINITE_ELEMENTS_MESHSMOOTH_H

#include <SurfaceMesh/SurfaceMesh.h>
#include <Eigen/Sparse>
#include "common.h"

using namespace Surface_Mesh;
using LLT = Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>;
using LDLT = Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>>;
using TripletSet = std::vector<Eigen::Triplet<float>>;

namespace alg {
    class MeshSmooth {
    public:
        MeshSmooth(const SurfaceMesh &mesh,
                   int maxIter = 10,
                   double epsAvg = 1e-3, double epsMin = 1e-3,
                   double alpha = 1.0);

        MeshSmooth(const MeshSmooth &smooth);

        MeshSmooth operator=(const MeshSmooth& rhs);

        SurfaceMesh smooth();

    private:
        std::vector<double> computeQs();

        double computeQAvg();

        double computeQMin();

        std::vector<double> computePatchMinVec(const std::vector<double> &qs);

        double computePatchMin(const std::vector<double> &qs, size_t vid);

        void smooth_pipeline();

        /**
         * Compute best fit rotation and translation matrices
         * @param Ce
         * @param Cr
         * @param R
         * @param T
         * @return
         */
        int computeRandTMatrices(const RowMat32f &Ce, const RowMat32f &Cr,
                                 RowMat2f &R, Vec2f &T);

        void assembleK(const RowMatf &Ke, const VecXi &vid, size_t fid);

        void assembleF(const VecXf &Fe, const VecXi &vid, size_t fid);

        void assembleU(const RowMat32f &Ue, const VecXi &vid);

        void assembleA();

        void assembleB();

        int buildSolveSystem(const TripletSet &Kset, const TripletSet &Aset);

        VecXf buildRhs(const VecXf &F, const VecXf &b);

        RowMatf solve(const VecXf &rhs);

        int updateMesh(const RowMatf &U);

        bool isConverge();

        void reset();

        void printOptInfo(int i);

    private:
        const double t_ = 1;
        const double E_ = 1;
        const double v_ = 0;
        const double l_ = 2 / (sqrt(sqrt(3)));
        SurfaceMesh mesh_;
        int max_iter_{1};
        double eps_avg_{1e-3};
        double eps_min_{1e-3};
        double alpha_{1.0};
        std::vector<double> qs_;
        std::vector<double> q_patch_min_;
        size_t a_rows_{0};
        double q_min_;
        double q_avg_;
        RowMat32f Cl_;
        Eigen::SparseMatrix<float> lhs_;
        TripletSet K_set_;
        TripletSet A_set_;
        LLT llt_;
        LDLT ldlt_;
        VecXf F_;
        VecXf U_;
        VecXf b_;
    };
}


#endif //MESH_SMOOTH_FOR_FINITE_ELEMENTS_MESHSMOOTH_H
