#ifndef DISKTRA_PATH_H
#define DISKTRA_PATH_H

#include <vector>
#include <memory>
#include <unordered_map>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <SurfaceMesh/SurfaceMesh.h>

#include "types.h"
#include "surface_mesh_utils.h"

using namespace common;
namespace alg {
class dijkstra_path {
public:
    dijkstra_path(const Surface_Mesh::SurfaceMesh &mesh,
                  const std::vector<float> &edge_lengths);

    int compute_path(const std::vector<surface_point> &points,
                     std::vector<size_t> &path,
                     std::vector<size_t> &fixed_ids,
                     bool loop = true);

    int compute_path(const common::surface_point &s,
                     const common::surface_point &t,
                     std::vector<size_t> &path) const;

    int get_nearest_vertex_half_path(const common::surface_point &s,
                                     const common::surface_point &t,
                                     const std::vector<size_t> &path) const;

protected:
    int get_neighboring_verts(const common::surface_point &p,
                              std::vector<size_t> &vs) const;

    int update_distance(size_t u, size_t v, size_t si, float dist_through_u,
                        std::unordered_map<size_t, float> &min_dist,
                        std::unordered_map<size_t, size_t> &prev,
                        std::set<std::pair<float, size_t>> &vertex_queue) const;

    int compute_paths(const common::surface_point &s,
                      const std::vector<common::surface_point> &ts,
                      std::unordered_map<size_t, float> &min_dist,
                      std::unordered_map<size_t, size_t> &prev) const;

    int compute_paths(const common::surface_point &s,
                      const std::vector<common::surface_point> &ts,
                      Eigen::VectorXf &min_dist,
                      Eigen::VectorXi &prev) const;

    int compute_paths(const size_t s,
                      const std::set<size_t> &ts,
                      Eigen::VectorXd &min_dist,
                      Eigen::VectorXi &prev) const;

private:
    const Surface_Mesh::SurfaceMesh &mesh_;
    const std::vector<float> &edge_lengths_;
};
}

#endif
