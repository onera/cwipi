#ifndef CONSERVATIVEMESH_HXX_
#define CONSERVATIVEMESH_HXX_

#include <string>
#include <vector>
#include <cmath>

#include "mesh.hxx"

#include "cwipi.h"

namespace cwipi
{
  class ConservativeMesh
  {
  public : 

    ///
    /// \brief Constructor
    /// 
    ///  @param [in] localComm   
    ///  @param [in] sourceMesh   source mesh
    ///  @param [in] targetMesh   target mesh
    ///  @param [in] tolerance    tolerance used for the distance around the vertices
    ///

    ConservativeMesh(const MPI_Comm &localComm,                    
                     Mesh& sourceMesh,
                     Mesh& targetMesh, 
                     const double tolerance
                     );


    ///
    /// \brief Compute the intersection between two edges
    /// 
    ///  @param [in] distMinP1SM       Minimal distance of the first vertex of the source edge
    ///  @param [in] distMinP2SM       Minimal distance of the second vertex of the source edge
    ///  @param [in] distMinP1TM       Minimal distance of the first vertex of the target edge
    ///  @param [in] distMinP2TM       Minimal distance of the second vertex of the target edge
    ///  @param [in] p1SourEdge        Coordinates of the first vertex of the source edge
    ///  @param [in] p2SourEdge        Coordinates of the second vertex of the source edge
    ///  @param [in] p1TarEdge         Coordinates of the first vertex of the target edge
    ///  @param [in] p2TarEdge         Coordinates of the second vertex of the target edge
    ///  @param [in] vois              Tolerance for the local coordinates s and t
    ///  @param [out] locCoordInters   Local coordinates s and t of the intersection points    
    ///  @return                       number of intersections (0,1 or 2)
    
    int intersectionEdge(
                          const double distMinP1SM,
                          const double distMinP2SM,
                          const double distMinP1TM,
                          const double distMinP2TM,
                          const double* p1SourEdge,
                          const double* p2SourEdge,
                          const double* p1TarEdge,
                          const double* p2TarEdge,
                          const double vois,
                          double* locCoordInters);
    


    ///
    /// \brief Compute the Element/Edge connectivity, the Edge/Vertex connectivity tab
    /// \brief the Vertex/Edge connectivity and index tabs
    /// 
    ///  @param [in]   mesh                     Mesh (source or target)
    ///  @param [out]  eltEdgeConnectivity      Element/Edge connectivity tab
    ///  @param [out]  edgeVertexConnectivity   Edge/Vertex connectivity tab
    ///  @param [out]  vertexEdgeConnectivity   Vertex/Edge connectivity tab
    ///  @param [out]  vertexEdgeIndex          Vertex/Edge index tab
    ///  @param [out]  nEdge                    Number of edges of the mesh
    /// 
    

    void compute_eltEdge_edgeVert_connectivity(const Mesh & mesh,
                                               std::vector<int> &eltEdgeConnectivity,
                                               std::vector<int> &edgeVertexConnectivity,
                                               std::vector<int> &vertexEdgeConnectivity,
                                               std::vector<int> &vertexEdgeIndex,
                                               int& nEdge);


    ///
    /// \brief  Compute the minimal distance of each vertex of the target and source mesh
    ///
    /// @param [in]  edgeVertexConnectivitySM  Edge/Vertex connectivity tab of the source mesh
    /// @param [in]  vertexEdgeConnectivitySM  Vertex/Edge connectivity tab of the source mesh
    /// @param [in]  vertexEdgeIndexSM         Vertex/Edge index tab of the source mesh
    /// @param [in]  edgeVertexConnectivityTM  Edge/Vertex connectivity tab of the target mesh
    /// @param [in]  vertexEdgeConnectivityTM  Vertex/Edge connectivity tab of the target mesh
    /// @param [in]  vertexEdgeIndexTM         Vertex/Edge index tab of the target mesh
    /// @param [out] distMinVertex             minimal distance tab
    ///      


    void computeMinDist(
                        const std::vector<int>& edgeVertexConnectivitySM,
                        const std::vector<int>& vertexEdgeConnectivitySM,
                        const std::vector<int>& vertexEdgeIndexSM,
                        const std::vector<int>& edgeVertexConnectivityTM,
                        const std::vector<int>& vertexEdgeConnectivityTM,
                        const std::vector<int>& vertexEdgeIndexTM,
                        std::vector<double>& distMinVertex );
    

    ///
    /// \brief  Compute all intersections between source edges and target edges and
    /// \brief  link the old edges and new edges tab used for the building
    ///
    /// @param [in]  edgeVertexConnectivitySM       Edge/Vertex connectivity tab of the source mesh
    /// @param [in]  vertexEdgeConnectivitySM       Vertex/Edge connectivity tab of the source mesh
    /// @param [in]  vertexEdgeIndexSM              Vertex/Edge index tab of the source mesh
    /// @param [in]  nEdgeSM                        Number of edges of the source mesh
    /// @param [in]  edgeVertexConnectivityTM       Edge/Vertex connectivity tab of the target mesh
    /// @param [in]  vertexEdgeConnectivityTM       Vertex/Edge connectivity tab of the target mesh
    /// @param [in]  vertexEdgeIndexTM              Vertex/Edge index tab of the target mesh
    /// @param [in]  nEdgeTM                        Number of edges of the target mesh
    /// @param [out]  edgeVertexConnectivityIM      Edge/Vertex connectivity tab of the intersection mesh
    /// @param [out]  vertexEdgeConnectivityIM      Vertex/Edge connectivity tab of the intersection mesh
    /// @param [out]  vertexEdgeIndexIM             Vertex/Edge index tab of the intersection mesh
    /// @param [out]  oldEdgeNewEdgeConnectivityIM  Old Edge/New edge connectivity tab
    /// @param [out]  oldEdgeNewEdgeIndexIM         Old Edge/New edge index tab
    /// @param [out]  edgeMeshTag                   Father mesh : 1 (source mesh) 2 (target mesh) 3 (the both)
    /// @param [out]  nEdge                         Number of edges of the intersection mesh
    /// @param [out]  nVertex                       Number of vertices of the intersection mesh
    /// @return                                     True if there is no cutting problem
    /// 

    bool cutMesh(
                 const std::vector<int>& edgeVertexConnectivitySM,
                 const std::vector<int>& vertexEdgeConnectivitySM,
                 const std::vector<int>& vertexEdgeIndexSM,
                 const int nEdgeSM,
                 const std::vector<int>& edgeVertexConnectivityTM,
                 const std::vector<int>& vertexEdgeConnectivityTM,
                 const std::vector<int>& vertexEdgeIndexTM,
                 const int nEdgeTM,
                 std::vector<int>& edgeVertexConnectivityIM,
                 std::vector<int>& vertexEdgeConnectivityIM,
                 std::vector<int>& vertexEdgeIndexIM,
                 std::vector<int>& oldEdgeNewEdgeConnectivityIM,
                 std::vector<int>& oldEdgeNewEdgeIndexIM,
                 std::vector<int>& edgeMeshTag,
                 int& nEdge,
                 int& nVertex);
    


    /// 
    /// \brief  Build the intersection mesh
    /// 
    /// @param [in]  eltEdgeConnectivitySM         Element/Edge connectivity tab of the source mesh
    /// @param [in]  eltEdgeConnectivityTM         Element/Edge connectivity tab of the target mesh
    /// @param [in]  nEdgeTM                       Number of edges of the target mesh                   
    /// @param [in]  edgeVertexConnectivityIM      Edge/Vertex connectivity tab of the intersection mesh
    /// @param [in]  vertexEdgeConnectivityIM      Vertex/Edge connectivity tab of the intersection mesh
    /// @param [in]  vertexEdgeIndexIM             Vertex/Edge index tab of the intersection mesh
    /// @param [in]  oldEdgeNewEdgeConnectivityIM  Old Edge/New edge connectivity tab
    /// @param [in]  oldEdgeNewEdgeIndexIM         Old Edge/New edge index tab
    /// @param [in]  edgeMeshTag                   Father mesh : 1 (source mesh) 2 (target mesh) 3 (the both)
    /// @param [in]  nEdge                         Number of edges of the intersection mesh
    /// @param [out] eltVertConnectivityIM         Element/Vertex connectivity tab of the intersection mesh
    /// @param [out] eltVertIndexIM                Element/Vertex index tab of the intersection mesh
    /// @param [out] nElts                         Number of elements of the intersection mesh
    /// @return                                    True if there is no building problem
    ///

    bool buildIntersectionMesh(
                              const std::vector<int>& eltEdgeConnectivitySM,
                              const std::vector<int>& eltEdgeConnectivityTM,
                              const int nEdgeTM,
                              const std::vector<int>& edgeVertexConnectivityIM,
                              const std::vector<int>& vertexEdgeConnectivityIM,
                              const std::vector<int>& vertexEdgeIndexIM,
                              const std::vector<int>& oldEdgeNewEdgeConnectivityIM,
                              const std::vector<int>& oldEdgeNewEdgeIndexIM,
                              const std::vector<int>& edgeMeshTag,                                 
                              const int nEdge,
                              std::vector<int>& eltVertConnectivityIM,
                              std::vector<int>& eltVertIndexIM,
                              int& nElts);    

    ///
    /// \brief Compute the dot product between two vectors
    /// 
    /// @param [in]  vect1  Coordinates of the first vector
    /// @param [in]  vect2  Coordinates of the second vector
    /// @return             value of the dot product
    ///

    inline double computeDotProduct(const double* vect1 , const double*  vect2);

    ///
    /// \brief Compute the dot product between two vectors
    /// 
    /// @param [in]  p1V1  Coordinates of the first point of the first vector
    /// @param [in]  p2V1  Coordinates of the second point of the first vector
    /// @param [in]  p1V2  Coordinates of the first point of the second vector
    /// @param [in]  p2V2  Coordinates of the second point of the second vector
    /// @return            value of the dot product
    ///

    inline double computeDotProduct(const double* p1V1 , 
                                    const double* p2V1,
                                    const double* p1V2 , 
                                    const double* p2V2);


    
    ///
    /// \brief Compute the cosinus between two vectors
    /// 
    /// @param [in]  vect1  Coordinates of the first vector
    /// @param [in]  vect2  Coordinates of the second vector
    /// @return             value of the cosinus
    ///

    inline double computeCos(const double* vect1 , const double*  vect2);


    ///
    /// \brief Compute the cosinus between two vectors
    /// 
    /// @param [in]  p1V1  Coordinates of the first point of the first vector
    /// @param [in]  p2V1  Coordinates of the second point of the first vector
    /// @param [in]  p1V2  Coordinates of the first point of the second vector
    /// @param [in]  p2V2  Coordinates of the second point of the second vector
    /// @return            value of the cosinus
    ///

    inline double computeCos(const double* p1V1 , 
                             const double* p2V1,
                             const double* p1V2 , 
                             const double* p2V2);


    ///
    /// \brief Compute the norm of a vector
    /// 
    /// @param [in]  vect1  Coordinates of the first vector
    /// @return             value of the norm
    ///


    inline double norme(const double* vect1);


    ///
    /// \brief Compute the norm of a vector
    /// 
    /// @param [in]  p1V1  Coordinates of the first point of the first vector
    /// @param [in]  p2V1  Coordinates of the second point of the first vector
    /// @return             value of the norm
    ///

    inline double norme(const double* p1V1 , 
                        const double* p2V1);

    
    ///
    /// \brief Compute the normal vector of two vectors
    /// 
    /// @param [in]  vect1          Coordinates of the first vector
    /// @param [in]  vect2          Coordinates of the second vector
    /// @param [out] crossProduct   Coordinates of the normal vector
    ///

    inline void computeCrossProduct(const double* vect1 , 
                                    const double* vect2,
                                    double* crossProduct);


    ///
    /// \brief Compute if two vectors in the same direction
    /// 
    /// @param [in]  vect1   Coordinates of the first vector
    /// @param [in]  vect2   Coordinates of the second vector
    /// @return              Return 1 if the two vectors are in the same direction, else -1
    ///

    inline int computeDirection(const double* vect1 ,
                                const double* vect2);

    ///
    /// \brief Return the intersection mesh
    ///

    inline Mesh* getIntersectionMesh();


  private :

    ConservativeMesh();

    ConservativeMesh(const ConservativeMesh&);

    ConservativeMesh& operator=(const ConservativeMesh&);



  private :

    Mesh *_intersectionMesh;                    ///< Intersection mesh
    std::vector<double> _coordsIM;              ///< Coordinates of the intersection mesh
    std::vector<int> _eltVertConnectivityIM;    ///< Element/Vertex connectivity tab of the intersection mesh
    std::vector<int> _eltVertIndexIM;           ///< Element/Vertex index tab of the intersection mesh

    Mesh & _sourceMesh;                         ///< Source mesh
    Mesh & _targetMesh;                         ///< Target mesh

    double _tolerance;                          ///< Tolerance used as geometric epsilon for the distance of the vertices

    std::vector<int> _newEltOldElt;             ///< New Element/Old Element tab used for the computing of value from one mesh to another. For a new Element i, _newEltOldElt[2*i] is the source father element and _newEltOldElt[2*i + 1] is the target father element
    

  }; /// class Conservative Mesh
    
  ///
  /// \brief Return the intersection mesh
  ///
  
  Mesh* ConservativeMesh::getIntersectionMesh(){
    return _intersectionMesh;
  }

} /// Namespace cwipi

#endif /*  CONSERVATIVE_HXX_ */
