//
// Created by Xu Rongyan on 2020/11/17.
//

#include <queue>
#include "surface_mesh_utils.h"

#include "math_utils.h"
#include "util.h"


namespace common {
    Eigen::MatrixX3d getPointsMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        Eigen::MatrixX3d V(mesh.n_vertices(), 3);
        for (const auto &v : mesh.vertices()) {
            const auto &p = mesh.position(v);
            V.row(v.idx()) << p[0], p[1], p[2];
        }
        return V;
    }

    Eigen::MatrixX3i getFaceMatrix(const Surface_Mesh::SurfaceMesh &mesh) {
        Eigen::MatrixX3i F(mesh.n_faces(), 3);
        for (const auto &f : mesh.faces()) {
            size_t i = 0;
            for (const auto &v : mesh.vertices(f)) {
                F(f.idx(), i++) = v.idx();
            }
        }
        return F;
    }

    RowMat32f get2DFaceCoord(const SurfaceMesh &mesh, const SurfaceMesh::Face &f) {
        RowMat32f mat;
        size_t i = 0;
        for (const auto &v : mesh.vertices(f)) {
            const auto& p = mesh.position(v);
            mat.row(i++) << p[0], p[1];
        }
        return mat;
    }

    RowMat3f get3DFaceCoord(const SurfaceMesh &mesh, const SurfaceMesh::Face &f) {
        RowMat3f mat;
        size_t i = 0;
        for (const auto &v : mesh.vertices(f)) {
            const auto& p = mesh.position(v);
            mat.row(i++) << p[0], p[1], p[2];
        }
        return mat;
    }

    Eigen::Vector3f getVector(const SurfaceMesh &mesh,
                              const SurfaceMesh::Halfedge &he) {
        Eigen::Vector3f res;
        if (!mesh.is_valid(he)) {
            return res;
        }
        auto start = mesh.from_vertex(he);
        auto end = mesh.to_vertex(he);
        res = mesh.position(end) - mesh.position(start);
        return res;
    }

    int
    buildSurfaceMesh(SurfaceMesh &mesh,
                     const Eigen::Matrix3Xf &pts,
                     const Eigen::Matrix3Xi &tris) {
        mesh.clear();
        for (size_t i = 0; i < pts.cols(); i++) {
            Point p = pts.col(i);
            mesh.add_vertex(p);
        }
        for (size_t i = 0; i < tris.cols(); i++) {
            const auto &face_ids = tris.col(i);
            std::vector<SurfaceMesh::Vertex> vs(3);
            for (size_t j = 0; j < 3; j++) {
                vs[j] = SurfaceMesh::Vertex(face_ids[j]);
                if (!mesh.is_valid(vs[j])) {
                    return 0;
                }
            }
            mesh.add_face(vs);
        }
        return 1;
    }

    int buildSurfaceMesh(SurfaceMesh &mesh, const Eigen::MatrixXd &pts, const Eigen::MatrixXi &tris) {
        mesh.clear();
        for (size_t i = 0; i < pts.cols(); i++) {
            Point p(pts(0, i), pts(1, i), pts(2, i));
            mesh.add_vertex(p);
        }
        for (size_t i = 0; i < tris.cols(); i++) {
            const auto &face_ids = tris.col(i);
            std::vector<SurfaceMesh::Vertex> vs(3);
            for (size_t j = 0; j < 3; j++) {
                vs[j] = SurfaceMesh::Vertex(face_ids[j]);
                if (!mesh.is_valid(vs[j])) {
                    return __LINE__;
                }
            }
            mesh.add_face(vs);
        }
        return 0;
    }

    int writeSurfaceMesh(const SurfaceMesh &mesh,
                         Eigen::MatrixXd &pts, Eigen::MatrixXi &tris) {
        pts.resize(3, mesh.n_vertices());
        for (const auto &v : mesh.vertices()) {
            const auto &p = mesh.position(v);
            pts.col(v.idx()) << p[0], p[1], p[2];
        }

        tris.resize(3, mesh.n_faces());
        for (const auto &f : mesh.faces()) {
            size_t i = 0;
            for (const auto &v : mesh.vertices(f)) {
                tris(i++, f.idx()) = v.idx();
            }
        }
        return 0;
    }

    double area(const SurfaceMesh &mesh,
                const std::vector<SurfaceMesh::Face> &faces) {
        double res = 0;
        for (const auto &f : faces) {
            assert(mesh.is_valid(f));
            res += area(mesh, f);
        }
        return res;
    }

    double area(const SurfaceMesh &mesh, const SurfaceMesh::Face &f) {
        assert(mesh.is_valid(f));
        std::vector<Point> ps;
        for (const auto &v : mesh.vertices(f)) {
            ps.emplace_back(mesh.position(v));
        }
        auto v1 = ps[1] - ps[0];
        auto v2 = ps[2] - ps[0];
        return (v1.cross(v2)).norm();
    }

    double angle(const Surface_Mesh::SurfaceMesh &mesh,
                 const Surface_Mesh::SurfaceMesh::Halfedge &he1,
                 const Surface_Mesh::SurfaceMesh::Halfedge &he2) {
        using namespace math;
        Eigen::Vector3f e1 = mesh.position(mesh.to_vertex(he1)) - mesh.position(mesh.from_vertex(he1));
        Eigen::Vector3f e2 = mesh.position(mesh.to_vertex(he2)) - mesh.position(mesh.from_vertex(he2));
        return angle(e1, e2);
    }

    Eigen::Vector3i vertexAroundFace(const Surface_Mesh::SurfaceMesh &mesh,
                                     const Surface_Mesh::SurfaceMesh::Face &f) {
        assert(mesh.is_valid(f));
        Eigen::Vector3i vertices;
        size_t i = 0;
        for (const auto &v : mesh.vertices(f)) {
            vertices[i++] = v.idx();
        }
        return vertices;
    }

    void decideFaceVertexOrder(const SurfaceMesh &basic_mesh,
                               const SurfaceMesh &sub_mesh,
                               const SurfaceMesh::Face &base_f,
                               std::vector<size_t> &f) {
        assert(basic_mesh.is_valid(base_f));
        auto basic_n = basic_mesh.compute_face_normal(base_f);
        const Point &p1 = sub_mesh.position(SurfaceMesh::Vertex(f[0]));
        const Point &p2 = sub_mesh.position(SurfaceMesh::Vertex(f[1]));
        const Point &p3 = sub_mesh.position(SurfaceMesh::Vertex(f[2]));
        Eigen::Vector3f n = (p2 - p1).cross(p3 - p1);
        n.normalize();
        double direct = basic_n.dot(n);
        if (direct < 0) {
            std::swap(f[0], f[1]);
        }
    }

    void decideFaceVertexOrder(const SurfaceMesh &basic_mesh,
                               const SurfaceMesh &sub_mesh,
                               const SurfaceMesh::Face &base_f,
                               std::vector<SurfaceMesh::Vertex> &fvs) {
        assert(basic_mesh.is_valid(base_f));
        auto basic_n = basic_mesh.compute_face_normal(base_f);
        const Point &p1 = sub_mesh.position(fvs[0]);
        const Point &p2 = sub_mesh.position(fvs[1]);
        const Point &p3 = sub_mesh.position(fvs[2]);
        Eigen::Vector3f n = (p2 - p1).cross(p3 - p1);
        n.normalize();
        double direct = basic_n.dot(n);
        if (direct < 0) {
            std::swap(fvs[0], fvs[1]);
        }
    }

    double ccwSearch(const SurfaceMesh &mesh, const SurfaceMesh::Halfedge &he,
                     const std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> &backward,
                     std::vector<SurfaceMesh::Halfedge> &ccw_seq) {
        ccw_seq.clear();
        SurfaceMesh::Halfedge hei = he;
        SurfaceMesh::Halfedge he_t;
        double theta = 0;
        bool cut = false;
        if (backward.size() == 1) {
            he_t = backward[0].first;
            cut = true;
        } else {
            he_t = backward[0].second == CW ? backward[0].first : backward[1].first;
        }
        while (hei != he_t) {
            ccw_seq.emplace_back(hei);
            auto he_next = mesh.ccw_rotated_halfedge(hei);
            theta += angle(mesh, hei, he_next);
            hei = he_next;
        }
        ccw_seq.emplace_back(he_t);
        if (cut)
            cutTail(mesh, ccw_seq, theta);
        return theta;
    }

    double cwSearch(const SurfaceMesh &mesh, const SurfaceMesh::Halfedge &he,
                    const std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> &backward,
                    std::vector<SurfaceMesh::Halfedge> &cw_seq) {
        cw_seq.clear();
        SurfaceMesh::Halfedge hei = he;
        SurfaceMesh::Halfedge he_t;
        double theta = 0;
        bool cut = false;
        if (backward.size() == 1) {
            he_t = backward[0].first;
            cut = true;
        } else {
            he_t = backward[0].second == CCW ? backward[0].first : backward[1].first;
        }
        while (hei != he_t) {
            cw_seq.emplace_back(hei);
            auto he_next = mesh.cw_rotated_halfedge(hei);
            theta += angle(mesh, hei, he_next);
            hei = he_next;
        }
        cw_seq.emplace_back(he_t);
        if (cut)
            cutTail(mesh, cw_seq, theta);
        return theta;
    }

    void cutHead(const SurfaceMesh &mesh, std::vector<SurfaceMesh::Halfedge> &seq, double &theta) {
        if (seq.size() < 2) return;
        //		theta -= angle(mesh, seq[0], seq[1]);
        seq.erase(seq.begin());
    }

    void cutTail(const SurfaceMesh &mesh, std::vector<SurfaceMesh::Halfedge> &seq, double &theta) {
        if (seq.size() < 2) return;
        //		theta -= angle(mesh, seq.back(), seq[seq.size() - 2]);
        seq.pop_back();
    }

    int
    barycentric(const SurfaceMesh &mesh,
                const SurfaceMesh::Face &f,
                const Eigen::Vector3f &p,
                std::vector<Eigen::Vector3f> &points,
                Eigen::Vector3f &bary) {
        //std::vector<Eigen::Vector3f> points;
        for (const auto &v : mesh.vertices(f)) {
            auto pos = mesh.position(v);
            points.emplace_back(pos);
        }
        return barycentric(points[0], points[1], points[2], p, bary);
    }

    int
    barycentric(const SurfaceMesh &mesh,
                const SurfaceMesh::Face &f,
                const Eigen::Vector3f &p,
                Eigen::Vector3f &bary) {
        std::vector<Eigen::Vector3f> points;
        return barycentric(mesh, f, p, points, bary);
    }

    SurfaceMesh::Vertex findABoundaryVertex(const SurfaceMesh &mesh) {
        for (const auto &v : mesh.vertices()) {
            if (mesh.is_boundary(v)) {
                return v;
            }
        }
        return SurfaceMesh::Vertex();
    }

    std::vector<SurfaceMesh::Face> neighborFaces(const SurfaceMesh &mesh,
                                                 const SurfaceMesh::Face &f) {
        std::vector<SurfaceMesh::Face> fs;
        for (const auto &he : mesh.halfedges(f)) {
            auto nei_f = mesh.face(mesh.opposite_halfedge(he));
            // some closed mesh might get invalid neighbor face
            if (mesh.is_valid(nei_f)) {
                fs.emplace_back(nei_f);
            }
        }
        return fs;
    }

    std::vector<std::vector<SurfaceMesh::Face>> bfsFace(const SurfaceMesh &mesh) {
        std::vector<std::vector<SurfaceMesh::Face>> res;
        std::vector<bool> visit(mesh.n_faces(), false);
        int visited = 0;

        auto visitFace = [&](const SurfaceMesh::Face &f,
                             std::queue<SurfaceMesh::Face> &q) {
            visit[f.idx()] = true;
            visited++;
            res.back().emplace_back(f);
            q.push(f);
        };

        while (visited != mesh.n_faces()) {
            size_t fid = 0;
            while (fid < visit.size() && visit[fid])
                fid++;
            res.emplace_back(std::vector<SurfaceMesh::Face>());
            SurfaceMesh::Face start_f(fid);
            std::queue<SurfaceMesh::Face> q;
            visitFace(start_f, q);
            while (!q.empty()) {
                auto f = q.front();
                q.pop();
                auto nei_fs = neighborFaces(mesh, f);
                for (const auto nei_f : nei_fs)
                    if (!visit[nei_f.idx()])
                        visitFace(nei_f, q);
            }
        }
        return res;
    }

    std::vector<std::vector<SurfaceMesh::Vertex>> bfsVertex(const SurfaceMesh &mesh) {
        std::vector<std::vector<SurfaceMesh::Vertex>> res;
        std::vector<bool> visit(mesh.n_vertices(), false);
        int visited = 0;
        auto visitFace = [&](const SurfaceMesh::Vertex &v,
                             std::queue<SurfaceMesh::Vertex> &q) {
            visit[v.idx()] = true;
            visited++;
            res.back().emplace_back(v);
            q.push(v);
        };

        while (visited != mesh.n_faces()) {
            size_t vid = 0;
            while (vid < visit.size() && visit[vid])
                vid++;
            res.emplace_back(std::vector<SurfaceMesh::Vertex>());
            SurfaceMesh::Vertex start_v(vid);
            std::queue<SurfaceMesh::Vertex> q;
            visitFace(start_v, q);
            while (!q.empty()) {
                auto v = q.front();
                q.pop();
                for (const auto nei_v : mesh.vertices(v))
                    if (!visit[nei_v.idx()])
                        visitFace(nei_v, q);
            }
        }
        return res;
    }

    std::vector<SurfaceMesh::Vertex> getVertexSet(const SurfaceMesh &mesh,
                                                  const std::vector<SurfaceMesh::Face> &faces) {
        std::vector<bool> visit(mesh.n_vertices(), false);
        std::vector<SurfaceMesh::Vertex> vertices;
        vertices.reserve(3 * faces.size());
        for (const auto &f : faces) {
            for (const auto &v : mesh.vertices(f)) {
                if (!visit[v.idx()]) {
                    visit[v.idx()] = true;
                    vertices.emplace_back(v);
                }
            }
        }
        return vertices;
    }

    void deleteFaces(SurfaceMesh &mesh,
                     std::vector<SurfaceMesh::Face> &faces,
                     bool collection) {
        for (const auto &f : faces) {
            assert(mesh.is_valid(f));
            mesh.delete_face(f);
        }
        if (collection) {
            mesh.garbage_collection();
        }
    }
}
