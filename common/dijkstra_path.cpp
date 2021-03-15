#include "dijkstra_path.h"

#include <unordered_map>

using namespace common;
//using namespace Surface_Mesh;
using Surface_Mesh::SurfaceMesh;
using Surface_Mesh::Point;

#define USING_UNORDER_MAP 1;

namespace alg {
dijkstra_path::dijkstra_path(const SurfaceMesh &mesh, const std::vector<float> &edge_lengths)
    : mesh_(mesh), edge_lengths_(edge_lengths) {
}

int dijkstra_path::get_neighboring_verts(const surface_point &p, std::vector<size_t> &vs) const {
    if (p.type_ == 2) {
        vs.reserve(3);
        const SurfaceMesh::Face face(p.idx_);
        for (const auto &vit : mesh_.vertices(face)) {
            vs.push_back(vit.idx());
        }
        return 2;
    }

    if (p.type_ == 1) {
        vs.reserve(4);
        const SurfaceMesh::Edge edge(p.idx_);
        for (int i = 0; i < 2; ++i) {
            SurfaceMesh::Halfedge he = mesh_.halfedge(edge, i);
            vs.push_back(mesh_.to_vertex(he).idx());
            he = mesh_.next_halfedge(he);
            vs.push_back(mesh_.to_vertex(he).idx());
        }
        return 1;
    }

    if (p.type_ == 0) {
        const SurfaceMesh::Vertex vert(p.idx_);
        for (const auto &he : mesh_.halfedges(vert)) {
            vs.push_back(mesh_.to_vertex(he).idx());
        }
        return 0;
    }
    return -1;
}

int dijkstra_path::compute_path(const std::vector<surface_point> &points,
                                std::vector<size_t> &path,
                                std::vector<size_t> &fixed_ids,
                                bool loop) {
    path.clear();
    size_t end = loop ? points.size() : points.size() - 1;
    size_t n = points.size();
    fixed_ids.emplace_back(0);
    for (size_t i = 0; i < end; i++) {
        std::vector<size_t> sub_path;
        int connect = compute_path(points[i], points[(i + 1) % n], sub_path);
        if (connect == -1) {
            return 0;
        }
        path.insert(path.end(), sub_path.begin(), sub_path.end() - 1);
        if (i != end - 1) {
            fixed_ids.emplace_back(path.size());
        }
        if (!loop && i == end - 1) {
            fixed_ids.emplace_back(path.size());
            path.emplace_back(sub_path.back());
        }
    }
    return 1;
}

int dijkstra_path::compute_path(const surface_point &s,
                                const surface_point &t,
                                std::vector<size_t> &path) const {
    std::vector<surface_point> ts;
    ts.push_back(t);

#ifdef   USING_UNORDER_MAP
    std::unordered_map<size_t, float> min_dist;
    std::unordered_map<size_t, size_t> prev;
    std::unordered_map<size_t, size_t>::iterator prev_iter;
#elif
    Eigen::VectorXf min_dist;
    Eigen::VectorXi prev;
#endif
    const int u = compute_paths(s, ts, min_dist, prev);
    if (u != -1) {
        path.clear();
        int vid = u;
#ifdef   USING_UNORDER_MAP
        while (1) {
            path.push_back(vid);
            prev_iter = prev.find(vid);
            if (prev_iter != prev.end())
                vid = prev_iter->second;
            else
                break;
        }
#elif
        for (; vid != -1; vid = prev[vid]){
            path.push_back(vid);
        }
#endif
        std::reverse(path.begin(), path.end());
        return path.size();
    }
    return -1;
}

//TODO:unfinished
int dijkstra_path::get_nearest_vertex_half_path(const surface_point &s,
                                                const surface_point &t,
                                                const std::vector<size_t> &path) const {
    const size_t num = path.size();
    if (num > 3) {
        Eigen::VectorXf lengths(num - 1);
        if (path[0] >= mesh_.n_vertices()) {
            const auto &p(mesh_.position(SurfaceMesh::Vertex(path[1])));
            lengths[0] = (s.pos_ - p).norm();
        } else {
            const SurfaceMesh::Edge edge = mesh_.find_edge(
                        SurfaceMesh::Vertex(path[0]), SurfaceMesh::Vertex(path[1]));
            lengths[0] = edge_lengths_[edge.idx()];
        }
        for (size_t i = 1; i < num - 1; ++i) {

        }
        if (path[num - 1] >= mesh_.n_vertices()) {
            const auto &p(mesh_.position(SurfaceMesh::Vertex(path[1])));
            lengths[num - 2] = (s.pos_ - p).norm();
        } else {

        }
    } else if (num == 3) {
        return path[1];
    } else if (num == 0)
        return -1;
    return path[0];
}

int dijkstra_path::update_distance(size_t u, size_t v, size_t si, float dist_through_u,
                                   std::unordered_map<size_t, float> &min_dist,
                                   std::unordered_map<size_t, size_t> &prev,
                                   std::set<std::pair<float, size_t>> &vertex_queue) const {
    float &md = min_dist[v];
    if (abs(md) < 1e-6 && v != si) {
        md = std::numeric_limits<float>::max();
    }
    if (dist_through_u < md) {
        vertex_queue.erase(std::make_pair(md, v));
        md = dist_through_u;
        prev[v] = u;
        vertex_queue.insert(std::make_pair(md, v));
        return 1;
    }
    return 0;
}

int dijkstra_path::compute_paths(const surface_point &s,
                                 const std::vector<surface_point> &ts,
                                 Eigen::VectorXf &min_dist,
                                 Eigen::VectorXi &prev) const {
    const size_t numV = mesh_.n_vertices();
    const size_t numVis = numV + 1 + ts.size();
    Eigen::Matrix3Xf pts(3, 1 + ts.size());
    int si = get_vertex_id(mesh_, s);
    if (si < 0) {
        si = numV;
        pts.col(0) = s.pos_;
    }
    min_dist.setConstant(numVis, std::numeric_limits<float>::infinity());
    min_dist[si] = 0;
    prev.setConstant(numVis, -1);
    std::unordered_map<int, int> vv;
    std::set<std::pair<float, size_t> > vertex_queue;
    vertex_queue.insert(std::make_pair(min_dist[si], si));

    std::set<size_t> tset;
    std::vector<size_t> vs;
    vs.reserve(ts.size());
    size_t vi(numV);
    for (const auto &t : ts) {
        ++vi;
        int ti = get_vertex_id(mesh_, t);
        if (ti < 0) {
            ti = vi;
            pts.col(ti - numV) = t.pos_;
            if (t.type_ == FACE_POINT) {
                const SurfaceMesh::Face face(t.idx_);
                for (const auto &vit : mesh_.vertices(face)) {
                    vv.insert(std::make_pair(vit.idx(), ti));
                }
            } else {
                const SurfaceMesh::Edge e(t.idx_);
                for (size_t i = 0; i < 2; ++i) {
                    SurfaceMesh::Halfedge he = mesh_.halfedge(e, i);
                    vv.insert(std::make_pair(mesh_.to_vertex(he).idx(), ti));
                    if (mesh_.face(he).is_valid()) {
                        vv.insert(std::make_pair(mesh_.to_vertex(mesh_.next_halfedge(he)).idx(), ti));
                    }
                }
            }
        }
        tset.insert(ti);
        vs.push_back(ti);
    }
    while (!vertex_queue.empty() && !tset.empty()) {
        const float dist = vertex_queue.begin()->first;
        const size_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());
        const auto it = tset.find(u);
        if (it != tset.end()) {
            tset.erase(it);
        }
        // Visit each edge exiting u
        if (u < numV) {
            const SurfaceMesh::Vertex vi(u);
            for (const auto &he : mesh_.halfedges(vi)) {
                const size_t v = mesh_.to_vertex(he).idx();
                const float dist_through_u = dist + edge_lengths_[mesh_.edge(he).idx()];
                float &md = min_dist[v];
                if (dist_through_u < md) {
                    vertex_queue.erase(std::make_pair(md, v));
                    md = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(md, v));
                }
            }

            const std::unordered_map<int, int>::const_iterator vit = vv.find(u);
            if (vit != vv.end()) {
                int v = vit->second;
                const Surface_Mesh::Point &pt = mesh_.position(vi);

                const float dist_through_u = dist + (pt - pts.col(v - numV)).norm();
                float &md = min_dist[v];
                if (dist_through_u < md) {//TODO: need further optimization
                    vertex_queue.erase(std::make_pair(md, v));
                    md = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(md, v));
                }
            }
        } else {
            const surface_point &fp = u == numV ? s : ts[u - numV - 1];
            std::vector<size_t> vs;
            if (fp.type_ == FACE_POINT) {
                const SurfaceMesh::Face f(fp.idx_);
                for (const auto &vit : mesh_.vertices(f)) {
                    vs.push_back(vit.idx());
                }
            } else {
                const SurfaceMesh::Edge e(fp.idx_);
                for (size_t i = 0; i < 2; ++i) {
                    SurfaceMesh::Halfedge he = mesh_.halfedge(e, i);
                    vs.push_back(mesh_.to_vertex(he).idx());
                    if (mesh_.face(he).is_valid()) {
                        vs.push_back(mesh_.to_vertex(mesh_.next_halfedge(he)).idx());
                    }
                }
            }
            for (auto v:vs) {
                const SurfaceMesh::Vertex vi(v);
                const Surface_Mesh::Point &pos = mesh_.position(vi);
                const float dist_through_u = dist + (pos - fp.pos_).norm();
                if (dist_through_u < min_dist[v]) {
                    vertex_queue.erase(std::make_pair(min_dist[v], v));
                    min_dist[v] = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(min_dist[v], v));
                }
            }
        }
    }
    int u = -1;
    float md = std::numeric_limits<float>::max();
    for (const auto &v : vs) {
        if (md > min_dist[v]) {
            md = min_dist[v];
            u = v;
        }
    }
    return u;
}

int dijkstra_path::compute_paths(const surface_point &s,
                                 const std::vector<surface_point> &ts,
                                 std::unordered_map<size_t, float> &min_dist,
                                 std::unordered_map<size_t, size_t> &prev) const {
    const size_t numV = mesh_.n_vertices();
    Eigen::Matrix3Xf pts(3, 1 + ts.size());
    int si = get_vertex_id(mesh_, s);
    if (si < 0) {
        si = numV;
        pts.col(0) = s.pos_;
    }
    std::unordered_map<int, int> vv;
    std::set<std::pair<float, size_t> > vertex_queue;
    vertex_queue.insert(std::make_pair(min_dist[si] = 1e-6, si));

    std::set<size_t> tset;
    std::vector<size_t> vs;
    vs.reserve(ts.size());
    size_t vi(numV);
    for (const auto &t : ts) {
        ++vi;
        int ti = get_vertex_id(mesh_, t);
        if (ti < 0) {
            ti = vi;
            pts.col(ti - numV) = t.pos_;
            if (t.type_ == FACE_POINT) {
                const SurfaceMesh::Face face(t.idx_);
                for (const auto &vit : mesh_.vertices(face)) {
                    vv.insert(std::make_pair(vit.idx(), ti));
                }
            } else {
                const SurfaceMesh::Edge e(t.idx_);
                for (size_t i = 0; i < 2; ++i) {
                    SurfaceMesh::Halfedge he = mesh_.halfedge(e, i);
                    vv.insert(std::make_pair(mesh_.to_vertex(he).idx(), ti));
                    const SurfaceMesh::Face f = mesh_.face(he);
                    if (f.is_valid()) {
                        vv.insert(std::make_pair(mesh_.to_vertex(mesh_.next_halfedge(he)).idx(), ti));
                    }
                }
            }
        }
        tset.insert(ti);
        vs.push_back(ti);
    }
    while (!vertex_queue.empty() && !tset.empty()) {
        const float dist = vertex_queue.begin()->first;
        const size_t u = vertex_queue.begin()->second;
        vertex_queue.erase(vertex_queue.begin());
        const auto it = tset.find(u);
        if (it != tset.end()) {
            tset.erase(it);
        }
        // Visit each edge exiting u
        if (u < numV) {
            const SurfaceMesh::Vertex vi(u);
            for (const auto &he : mesh_.halfedges(vi)) {
                const size_t v = mesh_.to_vertex(he).idx();
                const size_t e = mesh_.edge(he).idx();
                const float dist_through_u = dist + edge_lengths_[e];
                float &md = min_dist[v];
                if (md == 0 || dist_through_u < md) {
                    vertex_queue.erase(std::make_pair(md, v));
                    md = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(md, v));
                }
            }

            const std::unordered_map<int, int>::const_iterator vit = vv.find(u);
            if (vit != vv.end()) {
                int v = vit->second;
                const Surface_Mesh::Point &pt = mesh_.position(vi);
                const float dist_through_u = dist + (pt - pts.col(v - numV)).norm();
                float &md = min_dist[v];
                if (md == 0 || dist_through_u < md) {//TODO: need further optimization
                    vertex_queue.erase(std::make_pair(md, v));
                    md = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(md, v));
                }
            }
        } else {
            const surface_point &fp = u == numV ? s : ts[u - numV - 1];
            std::vector<size_t> vs;
            if (fp.type_ == FACE_POINT) {
                const SurfaceMesh::Face f(fp.idx_);
                for (const auto &vit : mesh_.vertices(f)) {
                    vs.push_back(vit.idx());
                }
            } else {
                const SurfaceMesh::Edge e(fp.idx_);
                for (size_t i = 0; i < 2; ++i) {
                    SurfaceMesh::Halfedge he = mesh_.halfedge(e, i);
                    vs.push_back(mesh_.to_vertex(he).idx());
                    const SurfaceMesh::Face f = mesh_.face(he);
                    if (f.is_valid()) {
                        vs.push_back(mesh_.to_vertex(mesh_.next_halfedge(he)).idx());
                    }
                }
            }
            for (auto v:vs) {
                const SurfaceMesh::Vertex vi(v);
                const float dist_through_u = dist + (mesh_.position(vi) - fp.pos_).norm();
                float &md = min_dist[v];
                if (md == 0 || dist_through_u < md) {
                    vertex_queue.erase(std::make_pair(md, v));
                    md = dist_through_u;
                    prev[v] = u;
                    vertex_queue.insert(std::make_pair(md, v));
                }
            }
        }
    }
    int u = -1;
    float min_dist_tmp = std::numeric_limits<float>::max();
    for (const auto &v : vs) {
        float &md = min_dist[v];
        if (md > 0 && min_dist_tmp > md) {
            min_dist_tmp = md;
            u = v;
        }
    }
    return u;
}


}
