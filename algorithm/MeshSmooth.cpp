//
// Created by 徐溶延local on 2021/3/13.
//

#include <spdlog/spdlog.h>
#include "MeshSmooth.h"
#include "energy/quality_metric.h"
#include "common/surface_mesh_utils.h"
#include "common/math_utils.h"


using namespace common;
using namespace common::math;
using namespace energy;

namespace alg {
    MeshSmooth::MeshSmooth(const SurfaceMesh &mesh,
                           int maxIter,
                           double epsAvg, double epsMin,
                           double alpha) : mesh_(mesh),
                                           max_iter_(maxIter),
                                           eps_avg_(epsAvg), eps_min_(epsMin),
                                           alpha_(alpha) {
        // initial quality calculations
        F_.resize(2 * mesh_.n_vertices());
        b_.resize(2 * mesh_.n_vertices());
        qs_ = computeQs();
        q_patch_min_ = computePatchMinVec(qs_);
        q_min_ = computeQMin();
        q_avg_ = computeQAvg();
        Cl_ << 0, 0,
               l_, 0,
               0.5 * l_, 0.5 * sqrt(3) * l_;

        printOptInfo(-1);
    }

    MeshSmooth::MeshSmooth(const MeshSmooth &smooth) {
        this->mesh_ = smooth.mesh_;
        this->max_iter_ = smooth.max_iter_;
        this->eps_avg_ = smooth.eps_avg_;
        this->eps_min_ = smooth.eps_min_;
        this->alpha_ = smooth.alpha_;
        this->qs_ = smooth.qs_;
        this->q_patch_min_ = smooth.q_patch_min_;
        this->q_min_ = smooth.q_min_;
        this->q_avg_ = smooth.q_avg_;
    }

    MeshSmooth MeshSmooth::operator=(const MeshSmooth &rhs) {
        MeshSmooth smoother(rhs.mesh_, rhs.max_iter_, rhs.eps_avg_, rhs.eps_min_, rhs.alpha_);
        return smoother;
    }


    SurfaceMesh MeshSmooth::smooth() {
        smooth_pipeline();
        return mesh_;
    }

    std::vector<double> MeshSmooth::computeQs() {
        std::vector<double> qs;
        qs.reserve(mesh_.n_faces());
        for (const auto &f : mesh_.faces()) {
            auto Ce = get2DFaceCoord(mesh_, f);
            double q = tri2dMetric(Ce);
            qs.emplace_back(q);
        }
        return qs;
    }

    double MeshSmooth::computeQAvg() {
        assert(!qs_.empty());
        double q_avg = 0;
        for (const auto &q : qs_) {
            q_avg += q;
        }
        return q_avg / qs_.size();
    }

    double MeshSmooth::computeQMin() {
        return *std::min_element(qs_.begin(), qs_.end());
    }

    std::vector<double> MeshSmooth::computePatchMinVec(const std::vector<double> &qs) {
        std::vector<double> patch_min;
        patch_min.reserve(qs.size());
        for (size_t i = 0; i < mesh_.n_vertices(); i++) {
            patch_min.emplace_back(computePatchMin(i));
        }
        return patch_min;
    }

    double MeshSmooth::computePatchMin(size_t vid) {
        const SurfaceMesh::Vertex v(vid);
        double q_min = 1.0;
        for (const auto &f : mesh_.faces(v)) {
            double q = qs_[f.idx()];
            q_min = std::min(q_min, q);
        }
        return q_min;
    }

    void MeshSmooth::smooth_pipeline() {
        // Smoothing algorithm iterative process
        for (size_t iter = 0; iter != max_iter_; iter++) {
            for (const auto &f : mesh_.faces()) {
                auto vid = vertexAroundFace(mesh_, f);
                auto C0 = get2DFaceCoord(mesh_, f);
                double A = area(C0);
                assert(A > 0);
                double s = compute2dS(A);
                auto Cr = computeRefLocalCoord(Cl_, s, alpha_);
                RowMat2f R;
                Vec2f T;
                computeRandTMatrices(C0, Cr, R, T);
                auto C = computeC(Cr, R, T);                    // Eq.29
                auto Ue = computeUe(C, C0);
                // compute Ke
                auto Ke = computeKe(C0, t_, E_, v_);
                // compute Fe
                auto Fe = computeFe(Ke, Ue);
                // Assemble Ke into global K matrix
                assembleK(Ke, vid, f.idx());
                // Assemble Fe into global F vector
                assembleF(Fe, vid, f.idx());
            }
            // mount global system using Lagrange multipliers
            assembleA();
            assembleB();
            int decompose_state = buildSolveSystem(K_set_, A_set_);
            if (!decompose_state) {
                spdlog::error("decompose failed");
            }
            auto rhs = buildRhs(F_, b_);
            // solve
            auto U = solve(rhs);
            // update mesh coords
            updateMesh(U);
            qs_ = computeQs();
            q_patch_min_ = computePatchMinVec(qs_);
            reset();
            bool is_converge = isConverge();
            printOptInfo(iter);
            if (is_converge)
                break;
        }
    }

    int MeshSmooth::computeRandTMatrices(const RowMat32f &Ce, const RowMat32f &Cr,
                                         RowMat2f &R, Vec2f &T) {
        auto Ce_avg = computeAvgCoord(Ce);      // Eq.26
        auto Cr_avg = computeAvgCoord(Cr);      // Eq.27
        assert(Ce_avg.size() == 2);
        assert(Cr_avg.size() == 2);
        auto D_prime = computeD(Ce, Ce_avg);    // Eq.24
        auto D = computeD(Cr, Cr_avg);          // Eq.25
        auto H = computeH(D, D_prime);          // Eq.23
        R = computeR(H);                        // Eq.14
        T = computeT(R, Ce_avg, Cr_avg);        // Eq.28
        return 1;
    }

    void MeshSmooth::assembleK(const RowMatf &Ke, const VecXi &vid, size_t fid) {
        for (size_t i = 0; i < Ke.rows(); i++) {
            for (size_t j = 0; j < Ke.cols(); j++) {
                int row = 2 * vid[i / 2] + (i % 2);
                int col = 2 * vid[j / 2] + (j % 2);
                K_set_.emplace_back(row, col, Ke(i, j));
            }
        }
    }

    void MeshSmooth::assembleF(const VecXf &Fe, const VecXi &vid, size_t fid) {
        F_.setZero();
        for (size_t i = 0; i < Fe.size(); i++) {
            size_t idx = 2 * vid[i / 2] + (i % 2);
            F_[idx] += Fe[i];
        }
    }

    void MeshSmooth::assembleA() {
        size_t x = 2 * mesh_.n_vertices();
        size_t y = 2 * mesh_.n_vertices();
        a_rows_ = 0;
        for (const auto &v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) {
                // A
                A_set_.emplace_back(x + 2 * a_rows_, 2 * v.idx(), 1);
                A_set_.emplace_back(x + 2 * a_rows_ + 1, 2 * v.idx() + 1, 1);
                // At
                A_set_.emplace_back(2 * v.idx(), y + 2 * a_rows_, 1);
                A_set_.emplace_back(2 * v.idx() + 1, y + 2 * a_rows_ + 1, 1);
                a_rows_++;
            }
        }
    }

    void MeshSmooth::assembleB(const RowMat32f &Ue) {
        b_.resize(2 * a_rows_);
        b_.setZero();
        int i = 0;
        for (const auto &v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) {


            }
        }
    }

    int MeshSmooth::buildSolveSystem(const TripletSet &Kset, const TripletSet &Aset) {
        const size_t k_rows = 2 * mesh_.n_vertices();
        const size_t k_cols = 2 * mesh_.n_vertices();
        const size_t a_rows = 2 * a_rows_;
        const size_t a_cols = 2 * mesh_.n_vertices();
        TripletSet lhs_set;
        lhs_set.insert(lhs_set.end(), Kset.begin(), Kset.end());
        lhs_set.insert(lhs_set.end(), Aset.begin(), Aset.end());

        lhs_.resize(k_rows + a_rows, k_cols + a_rows);
        lhs_.setFromTriplets(lhs_set.begin(), lhs_set.end());
        ldlt_.compute(lhs_);
        return ldlt_.info() == Eigen::Success;
    }

    VecXf MeshSmooth::buildRhs(const VecXf &F, const VecXf &b) {
        VecXf rhs(F.size() + b.size());
        for (size_t i = 0; i < F.size(); i++) {
            rhs[i] = F[i];
        }
        for (size_t i = F.size(); i < rhs.size(); i++) {
            rhs[i] = b[i - F.size()];
        }
//        std::cout << rhs.transpose() << std::endl;
        return rhs;
    }

    RowMatf MeshSmooth::solve(const VecXf &rhs) {
        const size_t n = 2 * mesh_.n_vertices();
        Eigen::VectorXf b = rhs.transpose();
        Eigen::VectorXf U_lambda = ldlt_.solve(b);

        if (ldlt_.info() != Eigen::Success) {
            spdlog::error("solve failed");
        }
        std::cout << (lhs_ * U_lambda - b).norm() << std::endl;
//        std::cout << U_lambda << std::endl;
        Eigen::VectorXf Utmp = U_lambda.block(0, 0, n, 1);
//        std::cout << Utmp << std::endl;
        return Eigen::Map<RowMatf>(Utmp.data(), Utmp.size() / 2, 2);
    }

    int MeshSmooth::updateMesh(const RowMatf &U) {
        assert(U.rows() == mesh_.n_vertices());
//        std::cout << U << std::endl;
        for (size_t i = 0; i < U.rows(); i++) {
            auto &p = mesh_.position(SurfaceMesh::Vertex(i));
            Point displacement;
            displacement << U.row(i).x(), U.row(i).y(), 0;
            p += displacement;
        }
        auto new_qs = computeQs();
        auto new_patch_min = computePatchMinVec(new_qs);
//        for (size_t i = 0; i < U.rows(); i++) {
//            if (!updatable(new_patch_min[i], q_patch_min_[i])) {
//                auto &p = mesh_.position(SurfaceMesh::Vertex(i));
//                Point displacement;
//                displacement << U.row(i).x(), U.row(i).y(), 0;
//                p -= displacement;
//            }
//        }
        return 1;
    }

    bool MeshSmooth::isConverge() {
        double q_min_curr = computeQMin();
        double q_avg_curr = computeQAvg();
        double delta_min = fabs(q_min_curr - q_min_);
        double delta_avg = fabs(q_avg_curr - q_avg_);
        bool is_converge = delta_min < eps_min_ && delta_avg < eps_avg_;
        q_min_ = q_min_curr;
        q_avg_ = q_avg_curr;
        return is_converge;
    }

    void MeshSmooth::reset() {
        K_set_.clear();
        A_set_.clear();
        lhs_.setZero();

    }

    void MeshSmooth::printOptInfo(int i) {
        if (i < 0) {
            spdlog::info("init quality: q_avg = {}, q_min = {}", q_avg_, q_min_);
        } else {
            spdlog::info("iter {}: q_avg = {}, q_min = {}", i, q_avg_, q_min_);
        }
    }


} // namespace alg
