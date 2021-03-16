#pragma once

#include <vector>
#include <set>
#include <ostream>

#include <Eigen/Dense>
#include <SurfaceMesh/SurfaceMesh.h>

//using namespace Surface_Mesh;
using Surface_Mesh::SurfaceMesh;
using Surface_Mesh::Point;

using RowMatf = Eigen::Matrix<float, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using RowMat3f = Eigen::Matrix<float, 3, 3, Eigen::RowMajor>;
using RowMat2f = Eigen::Matrix<float, 2, 2, Eigen::RowMajor>;
using RowMat32f = Eigen::Matrix<float, 3, 2, Eigen::RowMajor>;
using RowMat34i = Eigen::Matrix<size_t, 3, 4, Eigen::RowMajor>;
using RowMat34f = Eigen::Matrix<float, 3, 4, Eigen::RowMajor>;
using RowMat43f = Eigen::Matrix<float, 4, 3, Eigen::RowMajor>;
using ColMatd = Eigen::Matrix<float, -1, -1, Eigen::ColMajor>;
using ColMati = Eigen::Matrix<int, -1, -1, Eigen::ColMajor>;
using ColMat3f = Eigen::Matrix<float, 3, 3, Eigen::ColMajor>;
using ColMat34i = Eigen::Matrix<size_t, 3, 4, Eigen::ColMajor>;
using ColMat34f = Eigen::Matrix<float, 3, 4, Eigen::ColMajor>;
using ColMat43f = Eigen::Matrix<float, 4, 3, Eigen::ColMajor>;
using RowVec2f = Eigen::RowVector2f;
using RowVec3f = Eigen::RowVector3f;
using RowVec3i = Eigen::RowVector3i;
using RowVec4i = Eigen::Matrix<int, 4, 1, Eigen::RowMajor>;
using RowVecXf = Eigen::RowVectorXf;
using RowVecXi = Eigen::RowVectorXi;
using Vec2f = Eigen::Vector2f;
using Vec3f = Eigen::Vector3f;
using Vec3i = Eigen::Vector3i;
using Vec4i = Eigen::Vector4i;
using VecXf = Eigen::VectorXf;
using VecXi = Eigen::VectorXi;


namespace common {
    const double EPS3 = 1e-2;
    enum DIRECTION {
        BILATERAL, // bi-direction
        CW,  // clockwise
        CCW //anti-clockwise
    };


    enum LABEL_TYPE {
        TYPE_NONE = -1,             // inter-point
        TYPE_OUT = 0,               // boundary 1
        TYPE_IN = 1,                    // boundary 0
        TYPE_BOTH = 999,            // boundary belongs to both side of curve
        TYPE_BOUND = 998,           // face or vertex on curve segment
        TYPE_MULTI_LINE = 997,      // two segments in one face
        TYPE_UNCONNECT = 996,       // not belong
    };


    struct argument {
        std::string in_path;     // input mesh's file path
        std::string out_path;    // output polyline's file path
        std::string out_init_path;  // output polyline's initial path
        std::string data_path;    // control point's file path
        std::string topo_path;      // script of controlling topological operations
        bool loop;         // is closed curve
    };

    enum point_type {
        UNDEFINED_POINT = -1,
        VERTEX_POINT,
        EDGE_POINT,
        FACE_POINT
    };

    class surface_point {
    public:
        surface_point();

        surface_point(const Eigen::Vector3f &pt, point_type type, int idx, int flag = 0);

        /**
         * be careful about the type of face_idx
         * @param mesh
         * @param pt
         * @param face_idx
         */
        surface_point(const SurfaceMesh &mesh, const Eigen::Vector3f &pt, size_t face_idx);

        /**
         * br careful about the type of face_idx
         * @param mesh
         * @param bary
         * @param face_idx
         */
        surface_point(const SurfaceMesh &mesh, const float bary_x, const float bary_y, int face_idx);

        bool operator==(const surface_point &rhs) const;

        bool operator!=(const surface_point &rhs) const;

        ~surface_point() {}

    public:
        int set_surface_point(const SurfaceMesh &mesh,
                              const Eigen::Vector3f &bary,
                              int face_idx);

        void calcPos(const SurfaceMesh &mesh,
                     const Eigen::Vector3f &bary,
                     int face_idx);

    public:
        point_type type_;
        int flag_;
        int idx_;
        Eigen::Vector3f pos_;

    };

    class surface_edge {
    public:
        surface_edge();

        surface_edge(size_t idx1, size_t idx2);

        bool operator<(const surface_edge &rhs) const;

        bool operator>(const surface_edge &rhs) const;

        bool operator<=(const surface_edge &rhs) const;

        bool operator>=(const surface_edge &rhs) const;

        bool operator==(const surface_edge &rhs) const;

        bool operator!=(const surface_edge &rhs) const;

    public:
        size_t idx1, idx2;
    };


    class face_point : public surface_point {
    public:
        Eigen::Vector3f bary_;
    };

    class Snaxel : public surface_point {
    public:
        Snaxel() = default;

        Snaxel(Eigen::Vector3f &p, double t, size_t id1, size_t id2, bool isVertex);

        Snaxel(const SurfaceMesh &mesh, const SurfaceMesh::Vertex &v);

        Snaxel(const SurfaceMesh &mesh, size_t id);

        Snaxel(const SurfaceMesh &mesh, const SurfaceMesh::Halfedge &he, double t);

        Snaxel(const surface_point &p);

        int setPos(const SurfaceMesh &mesh, double t);

        int setPos(const SurfaceMesh &mesh,
                   const SurfaceMesh::Vertex &va,
                   const SurfaceMesh::Vertex &vb,
                   double t);

        int setPos(const SurfaceMesh &mesh,
                   const SurfaceMesh::Halfedge &he,
                   double t);

        Eigen::Vector3f calcPos(const Point &p1, const Point &p2, double t);

        Eigen::Vector3f calcPos(const SurfaceMesh &mesh);

        void update(const SurfaceMesh &mesh);

        void regularize_t();

        void cling(const SurfaceMesh &mesh);

        size_t getVertexId() const;

        int findEndPoints(const SurfaceMesh &mesh,
                          SurfaceMesh::Vertex &v1,
                          SurfaceMesh::Vertex &v2);

        int getNeighborFaces(const SurfaceMesh &mesh, std::set<size_t> &faces) const;

        Eigen::Vector3f direction(const SurfaceMesh &mesh) const;

        /**
        *
        * @param mesh
        * @return is point type changed
        */
        int updatePointType(const SurfaceMesh &mesh);

        bool isVertex() const;

    public:
        double t{};    // paramter t, for vertex position interpolation
        size_t id1{std::numeric_limits<size_t>::max()}, id2{
                std::numeric_limits<size_t>::max()};  //vertex coordinate
        bool fixed{false};    // is constraint point
    };

    inline void surface_points2eigen_points(const std::vector<surface_point> &pts, Eigen::Matrix3Xf &Q) {
        Q.resize(3, pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            Q.col(i) = pts[i].pos_;
        }
    }

// surface_point utils
    int get_vertex_id(const SurfaceMesh &mesh, const face_point &p);

    int get_vertex_id(const SurfaceMesh &mesh, const surface_point &p);

    std::vector<surface_point> ids2surface_points(const SurfaceMesh &mesh,
                                                  const std::vector<size_t> &ids);

    int findEdgePoints(const SurfaceMesh &mesh,
                       const surface_point &p,
                       SurfaceMesh::Vertex &v1,
                       SurfaceMesh::Vertex &v2);

    int findEndPoints(const SurfaceMesh &mesh,
                      const Snaxel &p,
                      SurfaceMesh::Vertex &v1,
                      SurfaceMesh::Vertex &v2);



// snaxel utils

    /**
    * snaxel the face of edge
    * @return
    */
    std::vector<size_t> faceOfSnaxelEdge(const SurfaceMesh &mesh,
                                         const Snaxel &s1,
                                         const Snaxel &s2);


    SurfaceMesh::Face findFaceBySnaxel(const SurfaceMesh &mesh,
                                       const Snaxel &s1,
                                       const Snaxel &s2);

    bool snaxelOptAble(const Snaxel &s1,
                       const Snaxel &s2,
                       const Snaxel &s3);

    bool isSnaxelEdgeBoundary(const SurfaceMesh &mesh,
                              const Snaxel &s1,
                              const Snaxel &s2);

    /**
    * looking for two edges that center at q, angle of qs; when s is the center, return qs
    * @param mesh
    * @param s		snaxel s
    * @param q		snaxel q
    * @return
    */
    std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> findWedgeEdgeDirection(const SurfaceMesh &mesh,
                                                                                    const Snaxel &q,
                                                                                    const Snaxel &s);

    std::vector<std::pair<SurfaceMesh::Halfedge, DIRECTION>> findInterWedgeDirection(const SurfaceMesh &mesh,
                                                                                     const Snaxel &s,
                                                                                     const surface_point &q);

//    std::vector<Snaxel> getSnaxelsByWedge(const SurfaceMesh &mesh,
//                                          const wedge &w,
//                                          double t = 0.5);

    std::vector<Snaxel> getSnaxelsByHalfEdge(const SurfaceMesh &mesh,
                                             const std::vector<SurfaceMesh::Halfedge> &halfEdges,
                                             double t = 0.5);

    std::vector<Snaxel> computeSnaxelsByIds(const SurfaceMesh &mesh,
                                            const std::vector<size_t> &ids);

    std::vector<Snaxel> computeSnaxelsByIds(const SurfaceMesh &mesh,
                                            const std::vector<size_t> &ids,
                                            const std::vector<size_t> &fixed_ids,
                                            const std::vector<surface_point> &points);

    int replaceSnaxelSeq(std::vector<Snaxel> &seq,
                         size_t src_s, size_t src_e,
                         const std::vector<Snaxel> &new_seq);

    int replaceSnaxelSeq(std::vector<Snaxel> &seq,
                         size_t src_s, size_t src_e,
                         const std::vector<Snaxel> &new_seq,
                         size_t dst_s, size_t dst_e);

    bool overlay(const SurfaceMesh &mesh, const Snaxel &s1, const Snaxel &s2);

    bool isMultiLineInSameFace(const SurfaceMesh &mesh,
                               const Snaxel &s1, const Snaxel &s2, const Snaxel &s3);

//    wedge chooseLowAngle(const wedge &cw_wedge, const wedge &ccw_wedge);
//
//    wedge chooseLowEnergy(const wedge &cw_wedge,
//                          const wedge &ccw_wedge,
//                          const std::vector<Snaxel> &seq,
//                          std::function<double(const std::vector<Snaxel> &, const SurfaceMesh::Halfedge &,
//                                               size_t)> energy,
//                          std::function<double(const SurfaceMesh::Halfedge &he)> step);
//
//    double computeWedgeEnergy(const wedge &w,
//                              std::vector<Snaxel> seq,
//                              std::function<double(const std::vector<Snaxel> &, const SurfaceMesh::Halfedge &,
//                                                   size_t)> energy,
//                              std::function<double(const SurfaceMesh::Halfedge &he)> step);

    double computeSeqEnergy(const std::vector<Snaxel> seq,
                            std::function<double(const std::vector<Snaxel> &seq)> energy);

    inline void snaxels2eigen_points(const std::vector<Snaxel> &pts, Eigen::Matrix3Xf &Q) {
        Q.resize(3, pts.size());
        for (size_t i = 0; i < pts.size(); ++i) {
            Q.col(i) = pts[i].pos_;
        }
    }


} // namespace common
