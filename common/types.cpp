#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include "types.h"
#include <vector>
#include<math.h>
#include <utility>
//#include <spdlog/spdlog.h>
#include "math_utils.h"
#include "surface_mesh_utils.h"



#include "surface_mesh_utils.h"
#include "util.h"

#define EPS2 1e-6

namespace common {
    surface_point::surface_point()
            : flag_(0), type_(UNDEFINED_POINT), idx_(-1) {

    }

    surface_point::surface_point(const Eigen::Vector3f &pt, point_type type, int idx, int flag)
            : pos_(pt), type_(type), idx_(idx), flag_(flag) {
    }

    surface_point::surface_point(const SurfaceMesh &mesh,
                                 const Eigen::Vector3f &pt,
                                 size_t face_idx) : flag_(0) {
        Eigen::Vector3f bary;
        const SurfaceMesh::Face f(face_idx);
        assert(mesh.is_valid(f));
        if (!barycentric(mesh, f, pt, bary)) {
            std::cout << "input point isn't in the face" << std::endl;
            exit(-1);
        }
        if (!isBaryValid(bary)) {
            std::cout << "invalid barycentric" << std::endl;
            exit(-1);
        }
        set_surface_point(mesh, bary, face_idx);
        calcPos(mesh, bary, face_idx);
    }

    surface_point::surface_point(const SurfaceMesh &mesh,
                                 const float bary_x, const float bary_y, int face_idx)
            : flag_(0) {
        Eigen::Vector3f bary(bary_x, bary_y, 1 - bary_x - bary_y);
        set_surface_point(mesh, bary, face_idx);
        calcPos(mesh, bary, face_idx);
    }

    int surface_point::set_surface_point(const SurfaceMesh &mesh,
                                         const Eigen::Vector3f &bary,
                                         int face_idx) {

        if (bary[0] > EPS2 && bary[1] > EPS2
            && bary[2] > EPS2) {//inside the face
            idx_ = face_idx;
            type_ = FACE_POINT;
            return 2;
        }

        const SurfaceMesh::Face face(face_idx);
        size_t i(0);
        std::vector<SurfaceMesh::Vertex> vs;
        for (const auto v : mesh.vertices(face)) {
            if (bary[i++] > EPS2) {
                vs.emplace_back(v);
            }
            if (vs.size() == 2) {
                idx_ = mesh.edge(mesh.find_halfedge(vs[0], vs[1])).idx();
                type_ = EDGE_POINT;
                return 1;
            }
        }
        //		i = 0;
        //		for (const auto &he: mesh.halfedges(face)) {
        //			if (bary[i] > EPS2 && bary[(i + 1) % 3] > EPS2) {//on the edge
        //				const SurfaceMesh::Vertex v0 = mesh.from_vertex(he);
        //				const SurfaceMesh::Vertex v1 = mesh.to_vertex(he);//to check
        //				idx_ = mesh.edge(mesh.find_halfedge(v0, v1)).idx();
        //				type_ = EDGE_POINT;
        //				return 1;
        //			}
        //			++i;
        //		}

        i = 0;
        for (const auto &vit : mesh.vertices(face)) {
            if (bary[i++] > EPS2) {//on the vertex
                idx_ = vit.idx();
                type_ = VERTEX_POINT;
                return 0;
            }
        }
        return -1;
    }

    void surface_point::calcPos(const SurfaceMesh &mesh,
                                const Eigen::Vector3f &bary,
                                int face_idx) {
        std::vector<Eigen::Vector3f> points;
        points.reserve(3);
        SurfaceMesh::Face f(face_idx);
        for (const auto &v : mesh.vertices(f)) {
            points.emplace_back(mesh.position(v));
        }
        pos_ = bary[0] * points[0] + bary[1] * points[1] + bary[2] * points[2];
    }

    bool surface_point::operator==(const surface_point &rhs) const {
        return type_ == rhs.type_ &&
               flag_ == rhs.flag_ &&
               idx_ == rhs.idx_;
    }

    bool surface_point::operator!=(const surface_point &rhs) const {
        return !(rhs == *this);
    }

    Snaxel::Snaxel(Eigen::Vector3f &p, double t, size_t id1, size_t id2, bool isVertex) :
            t(t), id1(id1), id2(id2) {
        regularize_t();
        pos_ = p;
        type_ = EDGE_POINT;
    }

    Snaxel::Snaxel(const SurfaceMesh &mesh, const SurfaceMesh::Vertex &v) : Snaxel(mesh, v.idx()) {}

    Snaxel::Snaxel(const SurfaceMesh &mesh, size_t id) {
        SurfaceMesh::Vertex v(id);
        assert(mesh.is_valid(v));
        typedef SurfaceMesh::Halfedge Halfedge;
        typedef SurfaceMesh::Vertex Vertex;
        t = 0;
        const Eigen::Vector3f q = mesh.position(v);
        const Halfedge he = mesh.halfedge(v);
        if (!he.is_valid()) {
            pos_ = q;
        } else {
            id1 = mesh.from_vertex(he).idx();
            id2 = mesh.to_vertex(he).idx();
            assert(id1 == v.idx());
            update(mesh);
        }
    }

    Snaxel::Snaxel(const SurfaceMesh &mesh, const SurfaceMesh::Halfedge &he, double t) {
        setPos(mesh, he, t);
    }

    Snaxel::Snaxel(const surface_point &p) {
        pos_ = p.pos_;
        type_ = p.type_;
        idx_ = p.idx_;
        fixed = true;

    }

    int Snaxel::setPos(const SurfaceMesh &mesh, double t) {
        if (fixed) return 0;
        SurfaceMesh::Vertex va(id1), vb(id2);
        this->t = t;
        regularize_t();
        update(mesh);
        return 1;
    }

    int
    Snaxel::setPos(const SurfaceMesh &mesh,
                   const SurfaceMesh::Vertex &va,
                   const SurfaceMesh::Vertex &vb,
                   double t) {
        assert(mesh.is_valid(va) && mesh.is_valid(vb));
        if (!mesh.is_valid(va) || !mesh.is_valid(vb)) {
            //spdlog::debug("invalid vertex");
            return 0;
        }
        id1 = va.idx();
        id2 = vb.idx();
        setPos(mesh, t);
        return 1;
    }

    int Snaxel::setPos(const SurfaceMesh &mesh, const SurfaceMesh::Halfedge &he, double t) {
        assert(mesh.is_valid(he));
        if (!mesh.is_valid(he)) {
            //spdlog::debug("invalid halfedge");
            return 0;
        }
        id1 = mesh.from_vertex(he).idx();
        id2 = mesh.to_vertex(he).idx();
        setPos(mesh, t);
        return 1;
    }

    Eigen::Vector3f Snaxel::calcPos(const Point &p1, const Point &p2, double t) {
        regularize_t();
        pos_ = p1 + t * (p2 - p1);
        return pos_;
    }

    Eigen::Vector3f Snaxel::calcPos(const SurfaceMesh &mesh) {
        const SurfaceMesh::Vertex va(id1), vb(id2);
        const Point &p1 = mesh.position(va);
        const Point &p2 = mesh.position(vb);
        regularize_t();
        return calcPos(p1, p2, t);
    }

    void Snaxel::update(const SurfaceMesh &mesh) {
        calcPos(mesh);
        updatePointType(mesh);
    }

    void Snaxel::regularize_t() {
        if (t > 1) t = 1;
        if (t < 0) t = 0;
    }

    void Snaxel::cling(const SurfaceMesh &mesh) {
        if (fixed) return;
        if (t < EPS3) {
            t = 0;
        }
        if (t > 1 - EPS3) {
            t = 1;
        }
        update(mesh);
    }

    size_t Snaxel::getVertexId() const {
        if (type_ != VERTEX_POINT) return -1;
        if (fixed) return idx_;
        return t > 0.5 ? id2 : id1;
    }

    int Snaxel::findEndPoints(const SurfaceMesh &mesh,
                              SurfaceMesh::Vertex &v1,
                              SurfaceMesh::Vertex &v2) {
        if (type_ != EDGE_POINT) return 0;
        if (fixed)
            return findEdgePoints(mesh, *this, v1, v2);
        v1 = SurfaceMesh::Vertex(id1);
        v2 = SurfaceMesh::Vertex(id2);
        return 1;
    }

    int Snaxel::getNeighborFaces(const SurfaceMesh &mesh, std::set<size_t> &faces) const {
        if (type_ == VERTEX_POINT) {
            size_t id = getVertexId();
            const SurfaceMesh::Vertex v(id);
            assert(mesh.is_valid(v));
            for (const auto &f : mesh.faces(v)) {
                faces.insert(f.idx());
            }
        } else if (!fixed) {
            const SurfaceMesh::Vertex v1(id1), v2(id2);
            const auto &e = mesh.find_edge(v1, v2);
            assert(mesh.is_valid(e));
            const auto &f1 = mesh.face(e, 0);
            const auto &f2 = mesh.face(e, 1);
            if (mesh.is_valid(f1)) faces.insert(f1.idx());
            if (mesh.is_valid(f2)) faces.insert(f2.idx());
        } else {
            if (type_ == EDGE_POINT) {
                SurfaceMesh::Vertex v1, v2;
                SurfaceMesh::Edge e(idx_);
                assert(mesh.is_valid(e));
                const auto &f1 = mesh.face(e, 0);
                const auto &f2 = mesh.face(e, 1);
                if (mesh.is_valid(f1)) faces.insert(f1.idx());
                if (mesh.is_valid(f2)) faces.insert(f2.idx());
            } else if (type_ == FACE_POINT) {
                SurfaceMesh::Face f(idx_);
                assert(mesh.is_valid(f));
                faces.insert(idx_);
            }
        }
        return 1;
    }

    Eigen::Vector3f Snaxel::direction(const SurfaceMesh &mesh) const {
        SurfaceMesh::Vertex v1(id1), v2(id2);
        return mesh.position(v2) - mesh.position(v1);
    }

    int Snaxel::updatePointType(const SurfaceMesh &mesh) {
        int tmp = type_;
        const auto &p1 = mesh.position(SurfaceMesh::Vertex(id1));
        const auto &p2 = mesh.position(SurfaceMesh::Vertex(id2));
        double dist1 = (p1 - pos_).norm();
        double dist2 = (p2 - pos_).norm();
        if (dist1 <= EPS2) {
            type_ = VERTEX_POINT;
            t = 0;
            calcPos(mesh);
        } else if (dist2 <= EPS2) {
            type_ = VERTEX_POINT;
            t = 0;
            std::swap(id1, id2);
            assert(mesh.is_valid(mesh.find_halfedge(SurfaceMesh::Vertex(id1), SurfaceMesh::Vertex(id2))));
            calcPos(mesh);
        } else {
            type_ = EDGE_POINT;
        }
        return tmp != type_;
    }

    bool Snaxel::isVertex() const {
        return type_ == VERTEX_POINT;
    }

    int get_vertex_id(const SurfaceMesh &mesh, const face_point &p) {
        const SurfaceMesh::Face face(p.idx_);
        size_t i(0);
        for (const auto &vit : mesh.vertices(face)) {
            if (fabs(1 - p.bary_[i++]) < 1e-6) {
                return vit.idx();
            }
        }
        return -1;
    }

    int get_vertex_id(const SurfaceMesh &mesh, const surface_point &p) {
        if (p.type_ == VERTEX_POINT) {
            return p.idx_;
        } else if (p.type_ == EDGE_POINT) {
            SurfaceMesh::Edge e(p.idx_);
            SurfaceMesh::Halfedge he = mesh.halfedge(e, 0);
            const Eigen::Vector3f &p1(mesh.position(mesh.from_vertex(he)));
            const Eigen::Vector3f &p2(mesh.position(mesh.to_vertex(he)));
            float r = (p1 - p.pos_).squaredNorm() / (p1 - p2).squaredNorm();
            const float eps = 1e-5;
            if (fabs(r - 1.0) < eps) {
                return mesh.to_vertex(he).idx();
            } else if (fabs(r) < eps)
                return mesh.from_vertex(he).idx();
        } else {
            SurfaceMesh::Face f(p.idx_);
            Eigen::Vector3f bary;
            if (barycentric(mesh, f, p.pos_, bary)) {
                size_t i(0);
                for (const auto &vit : mesh.vertices(f)) {
                    if (fabs(1 - bary[i++]) < 1e-3) {
                        return vit.idx();
                    }
                }
            }
        }
        return -1;
    }

    std::vector<surface_point> ids2surface_points(const SurfaceMesh &mesh,
                                                  const std::vector<size_t> &ids) {
        std::vector<surface_point> points;
        points.reserve(ids.size());
        for (auto id : ids) {
            SurfaceMesh::Vertex v(id);
            assert(mesh.is_valid(v));
            const Point p = mesh.position(v);
            points.emplace_back(p, VERTEX_POINT, id);
        }
        return points;
    }

    int findEdgePoints(const SurfaceMesh &mesh,
                       const surface_point &p,
                       SurfaceMesh::Vertex &v1,
                       SurfaceMesh::Vertex &v2) {
        if (p.type_ != EDGE_POINT) return 0;
        SurfaceMesh::Edge e(p.idx_);
        v1 = mesh.vertex(e, 0);
        v2 = mesh.vertex(e, 1);
        return 1;
    }

    int findEndPoints(const SurfaceMesh &mesh, const Snaxel &p, SurfaceMesh::Vertex &v1, SurfaceMesh::Vertex &v2) {
        if (p.type_ == VERTEX_POINT) {
            v1 = SurfaceMesh::Vertex(p.getVertexId());
            v2 = SurfaceMesh::Vertex();
            return 1;
        }
        if (p.type_ != EDGE_POINT) return 0;
        if (p.fixed) {
            return findEdgePoints(mesh, p, v1, v2);
        }
        v1 = SurfaceMesh::Vertex(p.id1);
        v2 = SurfaceMesh::Vertex(p.id2);
        return 1;
    }


    std::vector<size_t> faceOfSnaxelEdge(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2) {
        std::vector<size_t> fs;
        std::set<size_t> st1, st2;
        if (!s1.isVertex() && !s2.isVertex() && !s1.fixed && !s2.fixed) {
            s1.getNeighborFaces(mesh, st1);
            s2.getNeighborFaces(mesh, st2);
            std::set_intersection(st1.begin(), st1.end(), st2.begin(), st2.end(), std::back_inserter(fs));
        } else if (s1.isVertex() && s2.isVertex() && !s1.fixed && !s2.fixed) {
            size_t id1 = s1.getVertexId();
            size_t id2 = s2.getVertexId();
            const SurfaceMesh::Vertex v1(id1), v2(id2);
            const auto &e = mesh.find_edge(v1, v2);
            if (!mesh.is_valid(e)) {
                return fs;
            }
            const auto &f1 = mesh.face(e, 0);
            const auto &f2 = mesh.face(e, 1);
            if (mesh.is_valid(f1)) fs.push_back(f1.idx());
            if (mesh.is_valid(f2)) fs.push_back(f2.idx());
        } else {
            s2.getNeighborFaces(mesh, st1);
            s1.getNeighborFaces(mesh, st2);
            std::set_intersection(st1.begin(), st1.end(), st2.begin(), st2.end(), std::back_inserter(fs));
        }
        return fs;
    }

    SurfaceMesh::Face findFaceBySnaxel(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2) {
        std::vector<size_t> faces = faceOfSnaxelEdge(mesh, s1, s2);
        if (faces.size() < 1) {
            //spdlog::error("s2 is not a neighbor of s1.");
        } else if (faces.size() == 1) {
            return SurfaceMesh::Face(faces.front());
        } else {
            //spdlog::error("too many line in one face");
        }
        return SurfaceMesh::Face();
    }

    bool snaxelOptAble(const Snaxel &s1, const Snaxel &s2, const Snaxel &s3) {
        const auto p1 = s1.pos_;
        const auto p2 = s2.pos_;
        const auto p3 = s3.pos_;
        double theta = math::angle(p2, p1, p3);
        return fabs(theta - M_PI) > 1e-6;
    }

    bool isSnaxelEdgeBoundary(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2) {
        if (!s1.isVertex() || !s2.isVertex()) {
            return false;
        }
        size_t id1 = s1.getVertexId();
        size_t id2 = s2.getVertexId();
        SurfaceMesh::Vertex v1(id1), v2(id2);
        SurfaceMesh::Edge e = mesh.find_edge(v1, v2);
        if (mesh.is_valid(e)) {
            return mesh.is_boundary(e);
        }
        return false;
    }

    std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>>
    findWedgeEdgeDirection(const SurfaceMesh &mesh, const Snaxel &q, const Snaxel &s) {
        std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> res;
        if (s.isVertex()) {
            const auto &he = mesh.find_halfedge(SurfaceMesh::Vertex(q.id1), SurfaceMesh::Vertex(s.id1));
            if (!mesh.is_valid(he)) {
                return res;
            }
            res.emplace_back(he, BILATERAL);
        } else if (s.id1 == q.id1 || s.id2 == q.id1) {
            size_t id2 = s.id1 == q.id1 ? s.id2 : s.id1;
            const auto &he = mesh.find_halfedge(SurfaceMesh::Vertex(q.id1), SurfaceMesh::Vertex(id2));
            if (!mesh.is_valid(he)) {
                return res;
            }
            res.emplace_back(he, BILATERAL);
        } else {
            const auto &he1 = mesh.find_halfedge(SurfaceMesh::Vertex(q.id1), SurfaceMesh::Vertex(s.id1));
            const auto &he2 = mesh.find_halfedge(SurfaceMesh::Vertex(q.id1), SurfaceMesh::Vertex(s.id2));
            DIRECTION dir1 = CW;
            DIRECTION dir2 = CCW;
            if (!mesh.is_valid(he1) || !mesh.is_valid(he2)) {
                return res;
            }
            if (mesh.ccw_rotated_halfedge(he2) == he1) {
                std::swap(dir1, dir2);
            }
            res.emplace_back(he1, dir1);
            res.emplace_back(he2, dir2);
        }
        return res;
    }

    std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>>
    findInterWedgeDirection(const SurfaceMesh &mesh,
                            const Snaxel &s,
                            const surface_point &q) {
        std::vector<SurfaceMesh::Halfedge> hes;
        std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> res;
        assert(s.isVertex());
        SurfaceMesh::Vertex v_center(s.getVertexId());
        if (q.type_ == FACE_POINT) {
            SurfaceMesh::Face f(q.idx_);
            for (auto v : mesh.vertices(f)) {
                if (v != v_center) {
                    auto he = mesh.find_halfedge(v_center, v);
                    if (!mesh.is_valid(he)) {
                        return res;
                    }
                    hes.emplace_back(he);
                }
            }
            DIRECTION dir1 = CW;
            DIRECTION dir2 = CCW;
            if (mesh.ccw_rotated_halfedge(hes[1]) == hes[0]) {
                std::swap(dir1, dir2);
            }
            res.emplace_back(hes[0], dir1);
            res.emplace_back(hes[1], dir2);
        } else if (q.type_ == EDGE_POINT) {
            auto e = SurfaceMesh::Edge(q.idx_);
            assert(mesh.is_valid(e));
            auto v1 = mesh.vertex(e, 0);
            auto v2 = mesh.vertex(e, 1);
            if (v1 == v_center) {
                auto he = mesh.find_halfedge(v_center, v2);
                if (!mesh.is_valid(he)) {
                    return res;
                }
                res.emplace_back(he, BILATERAL);
            } else if (v2 == v_center) {
                auto he = mesh.find_halfedge(v_center, v1);
                if (!mesh.is_valid(he)) {
                    return res;
                }
                res.emplace_back(he, BILATERAL);
            } else {
                const auto &he1 = mesh.find_halfedge(v_center, v1);
                const auto &he2 = mesh.find_halfedge(v_center, v2);
                DIRECTION dir1 = CW;
                DIRECTION dir2 = CCW;
                if (!mesh.is_valid(he1) || !mesh.is_valid(he2)) {
                    return res;
                }
                if (mesh.ccw_rotated_halfedge(he2) == he1) {
                    std::swap(dir1, dir2);
                }
                res.emplace_back(he1, dir1);
                res.emplace_back(he2, dir2);
            }
        } else if (q.type_ == VERTEX_POINT) { // VERTEX_POINT
            SurfaceMesh::Vertex v(q.idx_);
            auto he = mesh.find_halfedge(v_center, v);
            if (!mesh.is_valid(he)) {
                return res;
            }
            res.emplace_back(he, BILATERAL);
        }
        return res;
    }

//    std::vector<Snaxel> getSnaxelsByWedge(const SurfaceMesh &mesh, const wedge &w, double t) {
//        return getSnaxelsByHalfEdge(mesh, w.halfedges_, t);
//    }

    std::vector<Snaxel>
    getSnaxelsByHalfEdge(const SurfaceMesh &mesh, const std::vector<SurfaceMesh::Halfedge> &halfEdges, double t) {
        std::vector<Snaxel> seq;
        seq.reserve(halfEdges.size());
        for (const auto &he : halfEdges) {
            Snaxel s(mesh, he, t);
            seq.emplace_back(s);
        }
        return seq;
    }

    std::vector<Snaxel> computeSnaxelsByIds(const SurfaceMesh &mesh,
                                            const std::vector<size_t> &ids) {
        std::vector<Snaxel> seq;
        for (auto id : ids) {
            seq.emplace_back(mesh, id);
        }
        return seq;
    }

    std::vector<Snaxel> computeSnaxelsByIds(const SurfaceMesh &mesh,
                                            const std::vector<size_t> &ids,
                                            const std::vector<size_t> &fixed_ids,
                                            const std::vector<surface_point> &points) {
        size_t i = 0;
        std::vector<Snaxel> seq;
        for (auto id : ids) {
            if (i < fixed_ids.size() && seq.size() == fixed_ids[i]) {
                seq.emplace_back(points[i++]);
            } else {
                seq.emplace_back(mesh, id);
            }
        }
        return seq;
    }

    int replaceSnaxelSeq(std::vector<Snaxel> &seq,
                         size_t src_s, size_t src_e,
                         const std::vector<Snaxel> &new_seq) {
        return replaceSnaxelSeq(seq, src_s, src_e, new_seq, 0, new_seq.size());
    }

    int replaceSnaxelSeq(std::vector<Snaxel> &seq,
                         size_t src_s, size_t src_e,
                         const std::vector<Snaxel> &new_seq,
                         size_t dst_s, size_t dst_e) {
        assert(src_s >= 0 && src_e <= seq.size());
        size_t j = src_s;
        for (size_t i = dst_s; i < dst_e; i++) {
            if (j < src_e) { // if the replaced sequence is not used completely then replace
                seq[j] = new_seq[i];
                j++;
            } else { // or replace all
                seq.insert(seq.begin() + j, new_seq.begin() + i, new_seq.begin() + dst_e);
                break;
            }
        }
        if (j < src_e) { // delete redundant sequence when new sequence smaller than order sequence
            seq.erase(seq.begin() + j, seq.begin() + src_e);
        }
        return 1;
    }

    bool overlay(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2) {
        if (!s1.isVertex() || !s2.isVertex()) return false;
        size_t id_s1 = s1.getVertexId();
        size_t id_s2 = s2.getVertexId();
        return id_s1 == id_s2;
    }

    bool isMultiLineInSameFace(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2, const Snaxel &s3) {
        std::vector<size_t> st;
        const auto fs1 = faceOfSnaxelEdge(mesh, s1, s2);
        const auto fs2 = faceOfSnaxelEdge(mesh, s2, s3);
        std::set_intersection(fs1.begin(), fs1.end(), fs2.begin(), fs2.end(), std::back_inserter(st));
        return !st.empty();
    }

    surface_edge::surface_edge() {}

    surface_edge::surface_edge(size_t idx1, size_t idx2) : idx1(idx1), idx2(idx2) {}

    bool surface_edge::operator==(const surface_edge &rhs) const {
        return idx1 == rhs.idx1 &&
               idx2 == rhs.idx2;
    }

    bool surface_edge::operator!=(const surface_edge &rhs) const {
        return !(rhs == *this);
    }

    bool surface_edge::operator<(const surface_edge &rhs) const {
        if (idx1 < rhs.idx1)
            return true;
        if (rhs.idx1 < idx1)
            return false;
        return idx2 < rhs.idx2;
    }

    bool surface_edge::operator>(const surface_edge &rhs) const {
        return rhs < *this;
    }

    bool surface_edge::operator<=(const surface_edge &rhs) const {
        return !(rhs < *this);
    }

    bool surface_edge::operator>=(const surface_edge &rhs) const {
        return !(*this < rhs);
    }

//    wedge chooseLowAngle(const wedge &cw_wedge, const wedge &ccw_wedge) {
//        double theta_cw = cw_wedge.angle_;
//        double theta_ccw = ccw_wedge.angle_;
//        if (fabs(theta_cw - theta_ccw) < 1e-6) {
//            return wedge();
//        }
//        return theta_cw < theta_ccw ? cw_wedge : ccw_wedge;
//    }
//
//    wedge chooseLowEnergy(const wedge &cw_wedge,
//                          const wedge &ccw_wedge,
//                          const std::vector<Snaxel> &seq,
//                          std::function<double(const std::vector<Snaxel> &,
//                                               const SurfaceMesh::Halfedge &,
//                                               size_t)> energy,
//                          std::function<double(const SurfaceMesh::Halfedge &he)> step) {
//        if (!cw_wedge.pMesh || !ccw_wedge.pMesh) {
//            return chooseLowAngle(cw_wedge, ccw_wedge);
//        }
//        double val_cw = computeWedgeEnergy(cw_wedge, seq, energy, step);
//        double val_ccw = computeWedgeEnergy(ccw_wedge, seq, energy, step);
//        return val_cw < val_ccw ? cw_wedge : ccw_wedge;
//    }
//
//    double computeWedgeEnergy(const wedge &w,
//                              std::vector<Snaxel> seq,
//                              std::function<double(const std::vector<Snaxel>&, const SurfaceMesh::Halfedge&, size_t)> energy,
//                              std::function<double(const SurfaceMesh::Halfedge &he)> step) {
//        auto mid = w.getMidEdge();
//        double t = step(mid);
//        seq[w.center_id].setPos(*(w.pMesh), mid, t);
//        return energy(seq, mid, w.center_id);
//    }

    double computeSeqEnergy(const std::vector<Snaxel> seq,
                            std::function<double(const std::vector<Snaxel> &seq)> energy) {
        return energy(seq);
    }
} // namespace common
