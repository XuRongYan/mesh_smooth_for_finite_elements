//
// Created by Xu Rongyan on 2020/11/17.
//

#ifndef MESH_CUTTING_SURFACE_MESH_UTILS_H
#define MESH_CUTTING_SURFACE_MESH_UTILS_H

#include <SurfaceMesh/SurfaceMesh.h>
#include "types.h"

//using namespace Surface_Mesh;
using Surface_Mesh::SurfaceMesh;
using Surface_Mesh::Point;
namespace common {
    Eigen::MatrixX3d getPointsMatrix(const SurfaceMesh &mesh);

    Eigen::MatrixX3i getFaceMatrix(const SurfaceMesh &mesh);

    RowMat32f get2DFaceCoord(const SurfaceMesh &mesh,
                             const SurfaceMesh::Face &f);

    RowMat3f get3DFaceCoord(const SurfaceMesh &mesh,
                            const SurfaceMesh::Face &f);

    Eigen::Vector3f getVector(const SurfaceMesh &mesh,
                              const SurfaceMesh::Halfedge &he);

    int buildSurfaceMesh(SurfaceMesh &mesh,
                         const Eigen::Matrix3Xf &pts,
                         const Eigen::Matrix3Xi &tris);

    int buildSurfaceMesh(SurfaceMesh &mesh, const Eigen::MatrixXd &pts, const Eigen::MatrixXi &tris);

    int writeSurfaceMesh(const SurfaceMesh &mesh,
                         Eigen::MatrixXd &pts, Eigen::MatrixXi &tris);

    template<typename Scalar>
    std::vector<Scalar> edgeLengths(const SurfaceMesh &mesh) {
        std::vector<Scalar> vec_lens;
        vec_lens.reserve(mesh.n_edges());
        for (const auto &e : mesh.edges()) {
            vec_lens.emplace_back(mesh.edge_length(e));
        }
        return vec_lens;
    }

    double area(const SurfaceMesh &mesh,
                const std::vector<SurfaceMesh::Face> &faces);

    /**
     * compute unsigned area
     * @param mesh
     * @param f
     * @return
     */
    double area(const SurfaceMesh &mesh, const SurfaceMesh::Face &f);

    double angle(const SurfaceMesh &mesh,
                 const SurfaceMesh::Halfedge &he1,
                 const SurfaceMesh::Halfedge &he2);

    Eigen::Vector3i vertexAroundFace(const SurfaceMesh &mesh,
                                     const SurfaceMesh::Face &f);

    void decideFaceVertexOrder(const SurfaceMesh &basic_mesh,
                               const SurfaceMesh &sub_mesh,
                               const SurfaceMesh::Face &base_f,
                               std::vector<size_t> &f);

    void decideFaceVertexOrder(const SurfaceMesh &basic_mesh,
                               const SurfaceMesh &sub_mesh,
                               const SurfaceMesh::Face &base_f,
                               std::vector<SurfaceMesh::Vertex> &fvs);

    double ccwSearch(const SurfaceMesh &mesh,
                     const SurfaceMesh::Halfedge &he,
                     const std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> &backward,
                     std::vector<SurfaceMesh::Halfedge> &ccw_seq);

    double cwSearch(const SurfaceMesh &mesh,
                    const SurfaceMesh::Halfedge &he,
                    const std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> &backward,
                    std::vector<SurfaceMesh::Halfedge> &cw_seq);

    void cutHead(const SurfaceMesh &mesh, std::vector<SurfaceMesh::Halfedge> &seq, double &theta);

    void cutTail(const SurfaceMesh &mesh, std::vector<SurfaceMesh::Halfedge> &seq, double &theta);

    int barycentric(const SurfaceMesh &mesh,
                    const SurfaceMesh::Face &f,
                    const Eigen::Vector3f &p,
                    Eigen::Vector3f &bary);

    int barycentric(const SurfaceMesh &mesh,
                    const SurfaceMesh::Face &f,
                    const Eigen::Vector3f &p,
                    std::vector<Eigen::Vector3f> &points,
                    Eigen::Vector3f &bary);

    SurfaceMesh::Vertex findABoundaryVertex(const SurfaceMesh &mesh);

    std::vector<SurfaceMesh::Face> neighborFaces(const SurfaceMesh &mesh,
                                                 const SurfaceMesh::Face &f);

    std::vector<std::vector<SurfaceMesh::Face>> bfsFace(const SurfaceMesh &mesh);

    std::vector<std::vector<SurfaceMesh::Vertex>> bfsVertex(const SurfaceMesh &mesh);

    std::vector<SurfaceMesh::Vertex> getVertexSet(const SurfaceMesh &mesh,
                                                  const std::vector<SurfaceMesh::Face> &faces);

    void deleteFaces(SurfaceMesh &mesh,
                     std::vector<SurfaceMesh::Face> &faces,
                     bool collection);

    SurfaceMesh triangle2mesh(const RowMat32f &C);

}


#endif //MESH_CUTTING_SURFACE_MESH_UTILS_H
