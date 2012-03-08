/*
  This file is part of the CWIPI library. 
  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include "conservativeMesh.hxx"

#include <mpi.h>
#include <algorithm>
#include <cmath>
#include <cassert>

#include <bftc_printf.h>

#include <fvmc_writer.h>

#include "locationToLocalMesh.hxx"
#include "locationToDistantMesh.hxx"
#include "applicationProperties.hxx"
#include "vectorUtilities.hxx"

namespace cwipi
{

  ConservativeMesh::ConservativeMesh(const MPI_Comm &localComm,
                                     Mesh & sourceMesh,
                                     Mesh & targetMesh, 
                                     const double tolerance) :

    _sourceMesh(sourceMesh), _targetMesh(targetMesh),_intersectionMesh(NULL),
    _tolerance(tolerance)
  {
    
    std::vector<int> eltEdgeConnectivitySM;          //tableau de connectivité Element/Arete du maillage source
    std::vector<int> edgeVertexConnectivitySM;       //tableau de connectivité Arete/Sommet du maillage source
    std::vector<int> vertexEdgeConnectivitySM;       //tableau de connectivité Sommet/Arete du maillage source
    std::vector<int> vertexEdgeIndexSM;              //tableau d'index Sommet/Arete du maillage source

    std::vector<int> eltEdgeConnectivityTM;          //tableau de connectivité Element/Arete du maillage cible
    std::vector<int> edgeVertexConnectivityTM;       //tableau de connectivité Arete/Sommet du maillage cible
    std::vector<int> vertexEdgeConnectivityTM;       //tableau de connectivité Sommet/Arete du maillage cible
    std::vector<int> vertexEdgeIndexTM;              //tableau d'index Sommet/Arete du maillage cible

    int nEdgeSM;                                     //nombre d'aretes du maillage source
    int nEdgeTM;                                     //nombre d'aretes du maillage cible
    int nEdge;                                       //nombre d'aretes du maillage intersecté
    int nVertex;                                     //nombre de sommets du maillage intersecté
    int nElts;                                       //nombre d'elements du maillage intersecté

    std::vector<int> edgeVertexConnectivityIM;       //tableau de connectivité Arete/Sommet du maillage intersecté
    std::vector<int> vertexEdgeConnectivityIM;       //tableau de connectivité Sommet/Arete du maillage intersecté
    std::vector<int> vertexEdgeIndexIM;              //tableau d'index Sommet/Arete du maillage intersecté
    std::vector<int> oldEdgeNewEdgeConnectivityIM;   //tableau de connectivité Ancienne arete/Nouvelle arete 
                                                     //(l'ancienne arete est decomposee, on a d'abord les aretes du maillages source et ensuite celles du maillage cible)

    std::vector<int> oldEdgeNewEdgeIndexIM;          //tableau d'index Ancienne arete/Nouvelle arete
    std::vector<int> edgeMeshTag;
    bool done;                                       //verifie que le decoupage et la construction s'est bien passee



    /**********   Construction des tableaux de connectivite   **********/

    compute_eltEdge_edgeVert_connectivity(_sourceMesh,
                                          eltEdgeConnectivitySM,
                                          edgeVertexConnectivitySM,
                                          vertexEdgeConnectivitySM,
                                          vertexEdgeIndexSM,
                                          nEdgeSM);

    compute_eltEdge_edgeVert_connectivity(_targetMesh,
                                          eltEdgeConnectivityTM,
                                          edgeVertexConnectivityTM,
                                          vertexEdgeConnectivityTM,
                                          vertexEdgeIndexTM,
                                          nEdgeTM);

    printf("Connectivity done \n");


    /**********   Decoupage   **********/

    printf("nEdgeTM %d \n",nEdgeTM);
    printf("nEdgeSM %d \n",nEdgeSM);

    
    done = cutMesh(
                   edgeVertexConnectivitySM,
                   vertexEdgeConnectivitySM,
                   vertexEdgeIndexSM,
                   nEdgeSM,
                   edgeVertexConnectivityTM,
                   vertexEdgeConnectivityTM,
                   vertexEdgeIndexTM,
                   nEdgeTM,
                   edgeVertexConnectivityIM,
                   vertexEdgeConnectivityIM,
                   vertexEdgeIndexIM,
                   oldEdgeNewEdgeConnectivityIM,
                   oldEdgeNewEdgeIndexIM,
                   edgeMeshTag,
                   nEdge,
                   nVertex);
    
    
    if(done)
      printf("Mesh cut done \n");


    /**********   Construction du maillage intersecte a partir des aretes (verification)   **********/

    std::vector<int> edgeIndex(nEdge + 1,2);

    Mesh* edgeMesh = new Mesh(MPI_COMM_WORLD,
                              1 ,
                              nVertex,
                              nEdge,
                              &_coordsIM[0],
                              &edgeIndex[0],
                              &edgeVertexConnectivityIM[0]);

    fvmc_writer_t *fvmWriter = fvmc_writer_init("EdgeMesh",
                                                          "/home/bfrisull/workspace/cwipi/src/compil/meshes",
                                                          "Ensight Gold",
                                                          "text",
                                                          FVMC_WRITER_FIXED_MESH);
    
    fvmc_writer_export_nodal(fvmWriter, 
                                 &(edgeMesh->getFvmNodal()));
    
    fvmc_writer_finalize(fvmWriter);
    


    /**********   Construction du maillage intersecte   **********/

    done = buildIntersectionMesh(eltEdgeConnectivitySM,
                                 eltEdgeConnectivityTM,
                                 nEdgeTM,
                                 edgeVertexConnectivityIM,
                                 vertexEdgeConnectivityIM,
                                 vertexEdgeIndexIM,
                                 oldEdgeNewEdgeConnectivityIM,
                                 oldEdgeNewEdgeIndexIM,
                                 edgeMeshTag,                                 
                                 nEdge,
                                 _eltVertConnectivityIM,
                                 _eltVertIndexIM,
                                 nElts);
      
    if(done)
      printf("Mesh Built done \n");

    _intersectionMesh = new  Mesh(localComm,
                                  2,
                                  nVertex,
                                  nElts,
                                  &_coordsIM[0],
                                  &_eltVertIndexIM[0],
                                  &_eltVertConnectivityIM[0]);
    

    /**********   Liberation de la memoire   **********/

    eltEdgeConnectivitySM.clear();
    edgeVertexConnectivitySM.clear();
    vertexEdgeConnectivitySM.clear();
    vertexEdgeIndexSM.clear();

    eltEdgeConnectivityTM.clear();
    edgeVertexConnectivityTM.clear();
    vertexEdgeConnectivityTM.clear();
    vertexEdgeIndexTM.clear();


    edgeVertexConnectivityIM.clear();
    edgeVertexConnectivityIM.clear();
    vertexEdgeIndexIM.clear();
    oldEdgeNewEdgeConnectivityIM.clear();
    oldEdgeNewEdgeIndexIM.clear();
    edgeMeshTag.clear();
  }

  int ConservativeMesh::intersectionEdge(
                                          const double distMinP1SM,
                                          const double distMinP2SM,
                                          const double distMinP1TM,
                                          const double distMinP2TM,
                                          const double* p1SourEdge,
                                          const double* p2SourEdge,
                                          const double* p1TarEdge,
                                          const double* p2TarEdge,
                                          const double vois,
                                          double* locCoordInters)
  {


   /*                  ->
   *                   v0(SM)     
   *         P1SM x-------------->P2SM
   *               \
   *                \               P2TM
   *                 \             ^
   *                  \           /
   *             ->    \         /
   *             v2     \       / ->
   *                     \     /  v1(TM)
   *                      \   /
   *                       \ /
   *                        x
   *                       P1TM
   */


    std::vector<double> v0 (3);                          //Arete source
    std::vector<double> v1 (3);                          //Arete cible
    std::vector<double> v2 (3);                          //Vecteur P1SM P1TM

    v0[0] = p2SourEdge[0] - p1SourEdge[0];
    v0[1] = p2SourEdge[1] - p1SourEdge[1];
    v0[2] = p2SourEdge[2] - p1SourEdge[2];

    v1[0] = p2TarEdge[0] - p1TarEdge[0];
    v1[1] = p2TarEdge[1] - p1TarEdge[1];
    v1[2] = p2TarEdge[2] - p1TarEdge[2];
    
    v2[0] = p1TarEdge[0] - p1SourEdge[0];
    v2[1] = p1TarEdge[1] - p1SourEdge[1];
    v2[2] = p1TarEdge[2] - p1SourEdge[2];
    

    double normV0 = dotProduct(&v0[0],&v0[0]);    //norme du vecteur V0
    double normV1 = dotProduct(&v1[0],&v1[0]);    //norme du vecteur V1
    double normV2 = dotProduct(&v2[0],&v2[0]);    //norme du vecteur V2

    double psV0V1 = dotProduct(&v0[0],&v1[0]);    //produit scalaire du (V0,V1)
    double psV0V2 = dotProduct(&v0[0],&v2[0]);    //produit scalaire du (V0,V2)
    double psV1V2 = dotProduct(&v1[0],&v2[0]);    //produit scalaire du (V1,V2)

    v0.clear();    
    v1.clear();
    v2.clear();
    
    bool intP1TM = false;                                //true si intersection avec le premier sommet de l'arete cible
    bool intP2TM = false;                                //true si intersection avec le second sommet de l'arete cible
    bool intersectionFound;                              //true si intersection sommet/sommet

    double distPointArete;                               //distance sommet arete
    double distMinSM;                                    //distance minimale pour un point de l'arete source
    double distMinTM;                                    //distance minimale pour un point de l'arete cible
    
    int nInters = 0;                                     //nombre d'intersections
    
    double t;                                            //valeur de la coordonne locale de l'intersection sur l'arete cible

    /********   Algorithme    ***********/
    /*          ----------
      P1SM = premier sommet de l'arete source
      P2SM = second sommet de l'arete source
      P1TM = premier sommet de l'arete cible
      P2TM = second sommet de l'arete cible

      -On projette P1SM sur l'arete cible et on divise par la norme pour trouver une 
      coordonnee locale.
      -On verifie que la coordonnee locale appartient a [0,1].
      -A partir de la formule  || v2 + (v1,v2)v1 ||² ) = dist²(P1SM,v1) , on en deduit si la distance
      de P1SM a l'arete est inferieur a la distance minimale de P1SM pour laquelle on considere
      qu'il y a intersection.
      -On verifie egalement que la distance de P1SM a l'arete est inferieur a la distance minimale
      de la projection de P1SM sur l'arete.
      -On sait a partir de ce moment qu'il y a intersection entre P1SM et l'arete cible ;
      il reste a determiner sa nature:
      
        -Si la coordonnee locale t est dans un voisinage de 0 et que la distance de P1SM a P1TM est
      inferieur aux distances minimales de P1SM et P1TM alors on en deduit qu'il y a intersection
      entre P1SM et P1TM.
        -Si la coordonnee locale t est dans un voisinage de 1 et que la distance de P1SM a P2TM est
      inferieur aux distances minimales de P1SM et P2TM alors on en deduit qu'il y a intersection
      entre P1SM et P2TM.
        -Sinon si 0 < t < 1 on en deduit une intersection entre P1SM et l'arete cible.
     */

    t= - psV1V2 / normV1;

    if( (t > -vois) && (t < 1 + vois) ){
     
      distPointArete = t*(normV1 * t + 2.0 * psV1V2 ) + normV2; 

      if(distPointArete < distMinP1SM * distMinP1SM){
        
        distMinTM = t*distMinP2TM  + (1.0 - t) * distMinP1TM;

        if(distPointArete < distMinTM * distMinTM ){
          intersectionFound = false;
          if ( t <= vois + 0.01){
            distPointArete = normV2;
            
            if(distPointArete < distMinP1SM * distMinP1SM &&
               distPointArete < distMinP1TM * distMinP1TM){

              locCoordInters[2*nInters] = -0.01;
              locCoordInters[2*nInters + 1] = -0.01;

              intP1TM = true;
              intersectionFound = true;

              nInters++;

            }
          }
          
          else if ( t >= .99 - vois){
            
            distPointArete = normV1 + 2.*psV1V2 + normV2;
            
            if(distPointArete < distMinP1SM * distMinP1SM &&
               distPointArete < distMinP2TM * distMinP2TM){
                          
              locCoordInters[2*nInters ] = -0.01;
              locCoordInters[2*nInters + 1] = 1.01;

              intP2TM = true;
              intersectionFound = true;
          
              nInters++;
            }            
          }
          
          if( !intersectionFound && (t < 1.0) && (t > 0.0) ){
            
            locCoordInters[2*nInters] = -0.01;
            locCoordInters[2*nInters + 1] = t;

            nInters++;
            
          }
        }
      }
      
    }
    
    /********   Algorithme    ***********/
    /*          ----------
      Intersection avec P2SM, idem que pour P1SM avec cette fois si la formule de distance 
      || - v0 + v2 + t*v1 ||² ) = dist²(P2e1,e2) et on travaille sur un voisinage
      de 1 poour la variable locale t.
    */


    t = (-psV1V2 + psV0V1) / normV1;

    if( (t > -vois) && (t < 1 + vois) ){

      distPointArete = normV0 + t*t*normV1 + normV2 - 2*t*psV0V1 - 2*psV0V2 + 2*t*psV1V2;

      if(distPointArete < distMinP2SM * distMinP2SM){
        
        distMinTM = t*distMinP2TM   + (1.0 - t) * distMinP1TM;

        if(distPointArete < distMinTM*distMinTM ){

          intersectionFound = false;

          if ( t <= vois + 0.01){
            distPointArete = normV0 - 2*psV0V2 + normV2;
            
            if(distPointArete < distMinP2SM * distMinP2SM &&
               distPointArete < distMinP1TM * distMinP1TM){
              
              locCoordInters[2*nInters] = 1.01;
              locCoordInters[2*nInters + 1] = -0.01;
              
              intP1TM = true;
              intersectionFound = true;

              nInters++;
            }
          }
          
          else if ( t >= .99 - vois){
            distPointArete = normV0 + normV1 + normV2 - 2*psV0V1 - 2*psV0V2 + 2*psV1V2;
            
            if(distPointArete < distMinP2SM * distMinP2SM &&
               distPointArete < distMinP2TM * distMinP2TM){

              locCoordInters[2*nInters] = 1.01;
              locCoordInters[2*nInters + 1] = 1.01;

              intP2TM = true;
              intersectionFound = true;

              nInters++;              
            }            
          }
          
          if( !intersectionFound && (t < 1.0) && (t > 0.0) ){

            locCoordInters[2*nInters] = 1.01;
            locCoordInters[2*nInters + 1] = t;

            nInters++;                        
          }
        }
      }      
    }


    if(nInters == 2)
      return nInters;

    
    double s;

    /********   Algorithme    ***********/
    /*          ----------
      Nous travaillons cette fois ci sur les sommets de l'arete cible. Les intersections
      Sommet/Sommet ont du etre detecte auparavant donc on verifie seulement s'il y a une 
      intersection Arete/Sommet.
     */

    if (!intP1TM){
      
      s = psV0V2 / normV0 ;
      
      if(s > 0. && s < 1.0){

        distPointArete = s*s*normV0 - 2*s*psV0V2 + normV2;

        if(distPointArete < distMinP1TM*distMinP1TM){

          distMinSM = s*distMinP2SM + ( 1.0 - s)*distMinP1SM;
          
          if(distPointArete < distMinSM*distMinSM){

            locCoordInters[2*nInters] = s;
            locCoordInters[2*nInters + 1] = -0.1;

            nInters ++;

            if(nInters == 2)
              return nInters;

          }       
        }
      }      
    }
    
    if(!intP2TM){


      s = (psV0V2 + psV0V1)/normV0;
      
      if(s > 0. && s < 1.0){

        distPointArete = s*s*normV0 + normV1 + normV2 +
          2*(psV1V2 - s*(psV0V2+psV0V1));
        
        if(distPointArete < distMinP1TM*distMinP1TM){

          distMinSM = s*distMinP2SM + ( 1.0 - s)*distMinP1SM ;
          
          if(distPointArete < distMinSM*distMinSM){

            locCoordInters[2*nInters] = s;
            locCoordInters[2*nInters + 1] = 1.01;

            nInters ++;

            if(nInters == 2)
              return nInters;
          }       
        }
      }            
    }

    if(nInters > 0)
      return nInters;


    /********   Algorithme *********/
    /*          ----------

    Les intersections avec un ou plusieurs sommets ont ete detecte, il reste a verifier s'il y a
    une intersection Arete/Arete.
    On utilise l'algorithme suivant :
    
    Si on considere s et t les coordonnees locales de l'arete source et cible respectivement, on a :

    dist = || v2 + t*v1 - s*v0||
          = s*s*norm(v0) + t*t*norm(v1)+ norm(v2) + 2*t*<v1,v2> - 2*s*<v0,v2> - 2*s*t*<v0,v1>

    On en deduit que :

    min dist <=> grad (dist) = 0
             <=> norm(v0)*s - <v0,v1>*t - <v0,v2> = 0
                 norm(v1)*t - <v0,v1>*s + <v1,v2> = 0
    
    Ce qui nous permet de conclure :

    solution <=> det = norm(v0)*norm(v1) - <v0,v1> * <v0,v1> != 0 
             <=> s = (<v0,v2>*norm(v1) - <v0,v1>*<v1,v2>)/det 
                 t = (<v0,v1>*<v0,v2> - <v1,v2> * norm(v2)) / det
    */

    double det = normV0*normV1 - psV0V1*psV0V1;

    if(det < 0.001*0.001*normV1*normV0) //pas de solutions
      return nInters;

    s = (psV0V2*normV1 - psV0V1*psV1V2)/det;
    t = (psV0V1*psV0V2 - psV1V2*normV0)/det;

    if((s > 0.0) && (t > 0.0) && (s < 1.0) && (t < 1.0)){
      
      distPointArete = s*s*normV0 + t*t*normV1 + normV2 + 2*t*psV1V2 - 2*s*psV0V2 - 2*s*t*psV0V1;

      distMinSM = s*distMinP1SM + (1. - s) * distMinP2SM;
      distMinTM = t*distMinP1TM + (1. - t) * distMinP2TM;
      
      if(distPointArete < distMinSM * distMinSM && distPointArete < distMinTM * distMinTM){
       
        locCoordInters[2*nInters] = s;
        locCoordInters[2*nInters + 1] = t;

        nInters++;
      }
        
    }
    return nInters;
     
  }
  
  void ConservativeMesh::compute_eltEdge_edgeVert_connectivity(
                                                               const Mesh& mesh,
                                                               std::vector<int> &eltEdgeConnectivity,
                                                               std::vector<int> &edgeVertexConnectivity,
                                                               std::vector<int> &vertexEdgeConnectivity,
                                                               std::vector<int> &vertexEdgeIndex,
                                                               int& nEdge)
  {

    const int *eltVertConnectivity = mesh.getEltConnectivity();    //tableau de connectivite Element/Sommet
    const int *eltVertIndex = mesh.getEltConnectivityIndex();      //tableau d'index Element/Sommet


    std::vector<int> stockTabTmp;                                  //tableau temporaire de stockage de valeur
    std::vector<int> hachTab;                                      //table de hachage

    int nElt = mesh.getNElts();                                    //nombre d'elements du maillage
    int nVert = mesh.getNVertex();                                 //nombre de sommets du maillage
    int nVertElt;                                                  //nombre de sommet par element
    int sumVertKey;                                                //cle utilise pour la table de hachage (somme des indices des deux sommets)
    int sumVertKeyMax = 0;                                         //cle maximale (borne sup de la somme des indices)
    int swapTmp;                                                   //variable utilisee pour un swap
    bool doubleEdge;                                               //true si une arete est en double

    stockTabTmp.resize(2*nVert); 

    nEdge = 0;        


    /********   Algorithme *********/
    /*          ----------
       Calcul du tableau de connectivite Arete/Noeud
       ---------------------------------------------

      - Pour chaque element on calcule la somme des indices de chaque arete et on enregistre le 
        nombre d'arete par somme d'indice dans le tableau de stockage temporaire.
        Chaque arete est compte deux fois pour chaque sommet
      - On transforme le tableau de stockage en un tableau d'index par rappot a la somme des indices
      - Une fois le tableau d'index créé, on peut remplir le tableau de connectivite en evitant
        de rajouter deux fois la meme arete et on stocke dans la table de hachage le nombre 
        d'aretes rajoutees pour chaque somme d'indice
      - A partir de la table de hachage on supprime les cases vides du tableau de connectivite
        Edge/Vertex.
     */

    for(int iElt = 0; iElt < nElt ; iElt ++){

      nVertElt = eltVertIndex[iElt + 1] - eltVertIndex[iElt];

      for(int iVert = 0; iVert < nVertElt ; iVert ++){

        if(iVert == nVertElt -1)
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt]];
        else
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt] + iVert + 1];
        
        sumVertKeyMax = sumVertKey > sumVertKeyMax ? sumVertKey : sumVertKeyMax;
        
        stockTabTmp[sumVertKey] ++;
        
      }
    }


    stockTabTmp.resize(sumVertKeyMax + 1);

    for(int iEdge = 0; iEdge < sumVertKeyMax + 1; iEdge ++){
      swapTmp = stockTabTmp[iEdge]; 
      stockTabTmp[iEdge] = nEdge;
      nEdge += swapTmp ;
    }


    edgeVertexConnectivity.resize(2*nEdge);
    hachTab.resize(sumVertKeyMax + 2);

    for(int iElt = 0; iElt < nElt ; iElt ++){

      nVertElt = eltVertIndex[iElt + 1] - eltVertIndex[iElt];

      for(int iVert = 0; iVert < nVertElt ; iVert ++){
        doubleEdge = false;

        if(iVert == nVertElt -1){
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt]];
          
          for (int i = 0; i < 2*hachTab[sumVertKey] ; i ++){

            if(eltVertConnectivity[eltVertIndex[iElt] + iVert] == edgeVertexConnectivity[2*(stockTabTmp[sumVertKey]) + i]){
              doubleEdge = true;
              break;
            }

          }

          if(!doubleEdge){
            edgeVertexConnectivity[2*(stockTabTmp[sumVertKey] + hachTab[sumVertKey])] = eltVertConnectivity[eltVertIndex[iElt] + iVert];
            edgeVertexConnectivity[2*(stockTabTmp[sumVertKey] + hachTab[sumVertKey]) + 1] = eltVertConnectivity[eltVertIndex[iElt]];

            hachTab[sumVertKey] = hachTab[sumVertKey] + 1;                    
          }
              
        }
        else{
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt] + iVert + 1];

          for (int i = 0; i < 2*hachTab[sumVertKey] ; i ++){

            if(eltVertConnectivity[eltVertIndex[iElt] + iVert] == edgeVertexConnectivity[2*(stockTabTmp[sumVertKey]) + i]){
              doubleEdge = true;
              break;
            }
          }
          
          if(!doubleEdge){
            edgeVertexConnectivity[2*(stockTabTmp[sumVertKey] + hachTab[sumVertKey])] = eltVertConnectivity[eltVertIndex[iElt] + iVert];
            edgeVertexConnectivity[2*(stockTabTmp[sumVertKey] + hachTab[sumVertKey]) + 1] = eltVertConnectivity[eltVertIndex[iElt] + iVert + 1];
            hachTab[sumVertKey] ++;        
          }              
        }
      }
    }
    
    nEdge = 0;

    for(int iIndex = 0 ; iIndex <= sumVertKeyMax ; iIndex ++){

      for(int i = 0 ; i < hachTab[iIndex] ; i ++){
        edgeVertexConnectivity[2*(nEdge + i)] = edgeVertexConnectivity[2*(stockTabTmp[iIndex] + i)];
        edgeVertexConnectivity[2*(nEdge + i) + 1] = edgeVertexConnectivity[2*(stockTabTmp[iIndex] + i) + 1];
      }

      swapTmp = hachTab[iIndex];
      hachTab[iIndex] = nEdge;
      nEdge += swapTmp;
    }
    
    hachTab[sumVertKeyMax + 1] = nEdge;
    
    edgeVertexConnectivity.resize(2*nEdge);
    eltEdgeConnectivity.resize(eltVertIndex[nElt]);


    /********   Algorithme *********/
    /*          ----------
       Calcul du tableau de connectivite Element/Arete
       -----------------------------------------------

       - Pour chaque element on recupere la somme des indices de deux sommets consecutifs
       - A partir de cette somme on accede a la table de hachage a la bonne adresse du tableau 
         de connectivite calcule precedemment.
       - On parcourt chaque arete qui a la meme somme d'indice, et s'il y a un sommet en commun,
         on en deduit que l'on est sur la bonne arete et on rajoute dans le tableau de connectivite
         Element/Arete, l'indice de cette arete. On met un signe + si l'arete est dans le meme sens
         et un signe - si l'arete du tableau precedent et de l'element sont dans le sens contraire
    */

    for(int iElt = 0; iElt < nElt ; iElt ++){
      
      nVertElt = eltVertIndex[iElt + 1] - eltVertIndex[iElt];
      
      for(int iVert = 0; iVert < nVertElt ; iVert ++){
        if(iVert == nVertElt -1)
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt]];
        else
          sumVertKey = eltVertConnectivity[eltVertIndex[iElt] + iVert]
            + eltVertConnectivity[eltVertIndex[iElt] + iVert + 1];

        
        for (int i = 0; i < hachTab[sumVertKey + 1] - hachTab[sumVertKey] ; i ++){          
          if(eltVertConnectivity[eltVertIndex[iElt] + iVert] 
             == edgeVertexConnectivity[2*(hachTab[sumVertKey] + i)]){
            eltEdgeConnectivity[eltVertIndex[iElt] + iVert] = hachTab[sumVertKey] + i + 1;  //peut etre erreur!!!!!!!
            break;
          }
          else if(eltVertConnectivity[eltVertIndex[iElt] + iVert] 
                  == edgeVertexConnectivity[2*(hachTab[sumVertKey] + i) + 1]){
            eltEdgeConnectivity[eltVertIndex[iElt] + iVert] = - (hachTab[sumVertKey] + i + 1); //peut etre erreur!!!!!!!
            break;
          }
          
          
        }
      }
    }

    stockTabTmp.clear();
    hachTab.clear();


    /********   Algorithme *********/
    /*          ----------
       Calcul du tableau de connectivite Noeud/Arete
       -----------------------------------------------
      - Il suffit de parcourir une premiere fois le tableau Arete/Noeud et de compter pour chaque
        noeud le nombre d'arete.
      - On transforme ce tableau en tableau d'index.
      - On reparcourt le tableau Arete/Noeud pour remplir convenablement le tableau Noeud/Arete       
    */


    vertexEdgeIndex.resize(nVert + 1);

    for(int iEdge = 0; iEdge < nEdge ; iEdge++){
      vertexEdgeIndex[edgeVertexConnectivity[2*iEdge] - 1]++;
      vertexEdgeIndex[edgeVertexConnectivity[2*iEdge + 1]- 1] ++;
    }

    vertexEdgeConnectivity.resize(2*nEdge);//taille exacte
    
    int sumTmp;
    sumVertKey = 0;
    
    if(stockTabTmp.size() <= nVert )
      stockTabTmp.resize(nVert);
    

    for(int iVert = 0; iVert < nVert + 1; iVert ++){
      sumTmp = vertexEdgeIndex[iVert];
      vertexEdgeIndex[iVert] = sumVertKey;
      sumVertKey += sumTmp;
      stockTabTmp[iVert] = 0;
    }

    for (int iEdge = 0; iEdge < nEdge ; iEdge ++){
      vertexEdgeConnectivity[vertexEdgeIndex[edgeVertexConnectivity[2*iEdge] - 1]
                             + stockTabTmp[edgeVertexConnectivity[2*iEdge] - 1] ]
                             = (iEdge + 1);

      vertexEdgeConnectivity[vertexEdgeIndex[edgeVertexConnectivity[2*iEdge + 1] - 1]
                             + stockTabTmp[edgeVertexConnectivity[2*iEdge + 1] - 1] ]
                             = -(iEdge + 1);

      stockTabTmp[edgeVertexConnectivity[2*iEdge] - 1] =  stockTabTmp[edgeVertexConnectivity[2*iEdge] - 1] + 1;
      stockTabTmp[edgeVertexConnectivity[2*iEdge + 1] - 1] =  stockTabTmp[edgeVertexConnectivity[2*iEdge + 1] - 1] +1;
      
                             
    }

    stockTabTmp.clear();
  }


  void ConservativeMesh::computeMinDist(const std::vector<int>& edgeVertexConnectivitySM,
                                        const std::vector<int>& vertexEdgeConnectivitySM,
                                        const std::vector<int>& vertexEdgeIndexSM,
                                        const std::vector<int>& edgeVertexConnectivityTM,
                                        const std::vector<int>& vertexEdgeConnectivityTM,
                                        const std::vector<int>& vertexEdgeIndexTM,
                                        std::vector<double>& distMinVertex
                                        ){


    const double* vertexCoordSM = _sourceMesh.getVertexCoords();   //Coordonnees des sommets du maillage source
    const double* vertexCoordTM = _targetMesh.getVertexCoords();   //Coordonnees des sommets du maillage cible

    int nVertSM = _sourceMesh.getNVertex();                        //Nombre de sommets du maillage source
    int nVertTM = _targetMesh.getNVertex();                        //Nombre de sommets du maillage cible
    int vert1;                                                     //Indice du premier sommet de l'arete
    int vert2;                                                     //Indice du second sommet de l'arete
    int indEdge;                                                   //Indice de l'arete
    double normMinEdge;                                            //Norme minimale des aretes conectees a un sommet
    double normEdge;                                               //Norme d'une arete

    std::vector<double> edge(3);                                   //Coordonnees du vecteur de l'arete

    distMinVertex.resize(nVertTM + nVertSM);
    

    /********   Algorithme *********/
    /*          ----------
      Pour chaque maillage on parcourt tous les noeuds et a partir du tableau de 
      connectivite Sommet/Arete on calcule l'arete minimale connectee au sommet. La 
      distance minimale est calculée en multipliant l'arete minimale par la tolerance.

     */




    for(int iVert = 0 ; iVert < nVertSM ; iVert ++){
      indEdge = std::abs(vertexEdgeConnectivitySM[vertexEdgeIndexSM[iVert]]) - 1;
      vert1 = edgeVertexConnectivitySM[2*indEdge] - 1;
      vert2 = edgeVertexConnectivitySM[2*indEdge + 1] - 1;

      edge[0] = vertexCoordSM[3*vert2] - vertexCoordSM[3*vert1] ;
      edge[1] = vertexCoordSM[3*vert2 + 1] - vertexCoordSM[3*vert1 + 1] ;
      edge[2] = vertexCoordSM[3*vert2 + 2] - vertexCoordSM[3*vert1 + 2] ;

      normMinEdge = norm(&edge[0]);

      for(int iEdge = 1; iEdge < vertexEdgeIndexSM[iVert + 1] - vertexEdgeIndexSM[iVert] ; iEdge ++){
        indEdge = std::abs(vertexEdgeConnectivitySM[vertexEdgeIndexSM[iVert] + iEdge]) - 1;
        vert1 = edgeVertexConnectivitySM[2*indEdge] - 1;
        vert2 = edgeVertexConnectivitySM[2*indEdge + 1] - 1;
        
        edge[0] = vertexCoordSM[3*vert2] - vertexCoordSM[3*vert1] ;
        edge[1] = vertexCoordSM[3*vert2 + 1] - vertexCoordSM[3*vert1 + 1] ;
        edge[2] = vertexCoordSM[3*vert2 + 2] - vertexCoordSM[3*vert1 + 2] ;
        
        normEdge = norm(&edge[0]);
        
        if(normEdge < normMinEdge)
          normMinEdge = normEdge;

      } 
      distMinVertex[iVert] = _tolerance*normMinEdge;
    }


    for(int iVert = 0 ; iVert < nVertTM ; iVert ++){
      indEdge = std::abs(vertexEdgeConnectivityTM[vertexEdgeIndexTM[iVert]]) - 1;
      vert1 = edgeVertexConnectivityTM[2*indEdge] - 1;
      vert2 = edgeVertexConnectivityTM[2*indEdge + 1] - 1;

      edge[0] = vertexCoordTM[3*vert2] - vertexCoordTM[3*vert1] ;
      edge[1] = vertexCoordTM[3*vert2 + 1] - vertexCoordTM[3*vert1 + 1] ;
      edge[2] = vertexCoordTM[3*vert2 + 2] - vertexCoordTM[3*vert1 + 2] ;

      normMinEdge = norm(&edge[0]);

      for(int iEdge = 1; iEdge < vertexEdgeIndexTM[iVert + 1] - vertexEdgeIndexTM[iVert] ; iEdge ++){
        indEdge = std::abs(vertexEdgeConnectivityTM[vertexEdgeIndexTM[iVert] + iEdge]) - 1 ;
        vert1 = edgeVertexConnectivityTM[2*indEdge] - 1;
        vert2 = edgeVertexConnectivityTM[2*indEdge + 1] - 1;
        
        edge[0] = vertexCoordTM[3*vert2] - vertexCoordTM[3*vert1] ;
        edge[1] = vertexCoordTM[3*vert2 + 1] - vertexCoordTM[3*vert1 + 1] ;
        edge[2] = vertexCoordTM[3*vert2 + 2] - vertexCoordTM[3*vert1 + 2] ;
        
        normEdge = norm(&edge[0]);
        
        if(normEdge < normMinEdge)
          normMinEdge = normEdge;

      } 
      distMinVertex[iVert + nVertSM] = _tolerance*normMinEdge;

    }

    edge.clear();
  }


  bool ConservativeMesh::cutMesh(
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
                                 int& nVertex)
  {
    
    const int *eltVertIndexSM = _sourceMesh.getEltConnectivityIndex();           //Tableau d'index Element/Sommet du maillage source
    const int *eltVertConnectivitySM = _sourceMesh.getEltConnectivity();         //Tableau de connectivite Element/Sommet du maillage source

    const int *eltVertIndexTM = _targetMesh.getEltConnectivityIndex();           //Tableau d'index Element/Sommet du maillage cible
    const int *eltVertConnectivityTM = _targetMesh.getEltConnectivity();         //Tableau de connectivite Element/Sommet du maillage cible

    const double* vertexCoordsSM = _sourceMesh.getVertexCoords();                //Coordonnees des sommets du maillage source
    const double* vertexCoordsTM = _targetMesh.getVertexCoords();                //Coordonnees des sommets du maillage cible
    
    const double* p1EdgeSM;                                                      //Coordonnees du premier point de l'arete source
    const double* p2EdgeSM;                                                      //Coordonnees du second point de l'arete source
    const double* p1EdgeTM;                                                      //Coordonnees du premier point de l'arete cible
    const double* p2EdgeTM;                                                      //Coordonnees du second point de l'arete cible


    std::vector<int> testEdge;                                                   //Tableau contenant les aretes cibles a tester pour les intersections d'aretes
    std::vector<int> edgeVertexConnectivitySMTmp (edgeVertexConnectivitySM);     //Tableau de connectivite temporaire Arete/Sommet du maillage source
    std::vector<int> edgeVertexConnectivityTMTmp (edgeVertexConnectivityTM);     //Tableau de connectivite temporaire Arete/Sommet du maillage cible     
    std::vector<int> oldVertNewVert;                                             //Tableau Ancien sommet/Nouveau sommet avec d'abord les sommets presents dans 
                                                                                 //les deux maillages, puis ceux presents dans le maillage cible puis dans le maillage source
    std::vector<int> newEdgeOldEdgeSM;                                           //Tableau Nouvelle arete/Ancienne Arete appartenant au maillage source
    std::vector<int> newEdgeOldEdgeTM;                                           //Tableau Nouvelle arete/Ancienne Arete appartenant au maillage cible
    std::vector<int> stockTabTmp;                                                //Tableau de stockage temporaire
    std::vector<int> hachTab;                                                    //Table de hachage
    std::vector<int> classifyEdge;                                               //Tableau contenant les nouveaux indices des aretes apres les avoir trie avec la somme des indices des sommets
                                                                                 //classifyEdge[i] = j ->i ancien indice dans le tableau de connectivite temporaire
                                                                                 //                      j nouvel indice dans le tableau de connectivite du maillage intersecte
    std::vector<int> doubleEdge;                                                 //Tableau contenant le nombre d'aretes en double pour une somme d'indice donnee
    std::vector<double> distMinVertex;                                           //Tableau contenant les distances minimales avec au debut celles des noeuds du maillage source, puis cible, puis intersecte
    std::vector<double> locCoordInters;                                          //Tableau contenant les coordonnees locales d'une ou des intersections entre deux aretes

    
    int nVertSM = _sourceMesh.getNVertex();                                      //Nombre de sommets du maillage source
    int nVertTM = _targetMesh.getNVertex();                                      //Nombre de sommets du maillage cible
    int nIntersGlob = 0;                                                             //Nombre d'intersections totales
    int nIntersLoc;                                                              //Nombre d'intersections entre deux aretes
    int nEdgeInters = 0;                                                         //Nombre d'intersections Arete/Sommet ou Arete/Arete pour une arete cible avec toutes les aretes sources
    int nEdgeIntersTmp;                                                          //Nombre d'intersections Arete/Sommet ou Arete/Arete pour une arete cible avec une arete source
    int nCutEdgeSMTmp;                                                           //Nombre d'intersections Arete/Sommet ou Arete/Arete de toutes les aretes sources avec une arete cible
    int nCutEdgeTM = 0;                                                          //Nombre d'aretes du maillage cible apres le decoupage
    int nCutEdgeSM = nEdgeSM;                                                    //Nombre d'aretes du maillage source apres le decoupage
    int nCutEdgeSMOld;
    int indP1SM;                                                                 //Indice du premier sommet de l'arete source
    int indP2SM;                                                                 //Indice du second sommet de l'arete source
    int indP1TM;                                                                 //Indice du premier sommet de l'arete cible
    int indP2TM;                                                                 //Indice du second sommet de l'arete cible
    int nDoubleEdge = 0;                                                         //Nombre d'aretes en double
    int sumVertKey;                                                              //Cle pour la table de hachage (somme des indices)
    int sumVertKeyMax;                                                           //Somme de sommets maximale
    int sumTmp;                                                                  //Valeur temporaire


    double distMinP1SM;                                                          //Distance minimale du premier sommet de l'arete source
    double distMinP2SM;                                                          //Distance minimale du second sommet de l'arete source
    double distMinP1TM;                                                          //Distance minimale du premier sommet de l'arete cible
    double distMinP2TM;                                                          //Distance minimale du second sommet de l'arete cible
    double signe;                                                                //Signe utilise pour savoir si un sommet est le premier ou le second de l'arete
    bool doubleVertex;                                                           //True si le sommet a deja ete enregistree
    bool alreadyCreated;                                                         //True si l'arete a deja ete enregistree


    testEdge.resize(10);
    edgeVertexConnectivitySMTmp.resize(2*edgeVertexConnectivitySMTmp.size());
    edgeVertexConnectivityTMTmp.resize(2*(nEdgeSM + nEdgeTM));
    oldVertNewVert.resize(2*(nVertSM + nVertTM));
    newEdgeOldEdgeSM.resize(3*nEdgeSM);
    newEdgeOldEdgeTM.resize(3*nEdgeTM);
    distMinVertex.resize(2*(nVertSM + nVertTM));
    locCoordInters.resize(4);

    _coordsIM.resize(3*(nVertSM + nVertTM));


    /**********   Initialisation   **********/

    computeMinDist(edgeVertexConnectivitySM,
                   vertexEdgeConnectivitySM,
                   vertexEdgeIndexSM,
                   edgeVertexConnectivityTM,
                   vertexEdgeConnectivityTM,
                   vertexEdgeIndexTM,
                   distMinVertex);


    for (int iEdge = 0; iEdge < nEdgeSM ; iEdge ++)
	newEdgeOldEdgeSM[iEdge] = iEdge + 1;

    for (int iEdge = 0; iEdge < nEdgeTM ; iEdge ++)
	newEdgeOldEdgeTM[iEdge] = iEdge + 1;



    /**********   Intersections   **********/

    /********   Algorithme    ***********/
    /* 
      - On boucle sur les aretes cibles
        - On rajoute l'arete dans le tableau a tester
        - On boucle sur les aretes sources
          - On recupere les indices des sommets de l'arete source
          - Selon l'indice on recupere les coordonnees dans le tableau correspondant
          - On boucle sur les aretes a tester
            - On recupere les indices des sommets de l'arete cible
            - Selon l'indice on recupere les coordonnees dans le tableau correspondant
            - On calcule le nombre d'intersections via la fonction precedente
            - Si deux sommets ont le meme indice, on ignore cette intersection.
            - S'il y a deux intersections on a 18 types d'intersections :

              - A1 x--------x B1            deux sommets confondus
                A2            B2

              - A1 x--------x----x B2       un sommet confondu + arete cible coupee en deux
                A2          B1

              - A1 x--------x----x B1       un sommet confondu + arete source coupee en deux       
                A2          B2

              - A1 x--------x B1            deux sommets confondus
                B2            A2

              - A1 x--------x----x A2       un sommet confondu + arete cible coupee en deux
                B2          B1

              - A1 x--------x----x B1       un sommet confondu + arete source coupee en deux
                B2          A2

              - B2 x----x--------x B1       un sommet confondu + arete cible coupee en deux
                        A1         A2
                
              - A2 x----x--------x B1       un sommet confondu + arete cible coupee en deux
                        A1         B2

              - A2 x----x-----x----x B2     arete cible coupee en trois
                        A1    B1     

              - B2 x----x-----x----x A2     arete cible coupee en trois
                        A1    B1     

              - B2 x----x-----x----x B1     arete source coupee en deux + arete cible coupee en deux
                        A1    A2     

              - A2 x----x-----x----x B1     arete source coupee en deux + arete cible coupee en deux
                        A1    B2     

              - B1 x--------x----x A1       un sommet confondu + arete source coupee en deux
                A2          B2 

              - B1 x--------x----x A1       un sommet confondu + arete source coupee en deux
                B2          A2 

              - B2 x----x-----x----x A1     arete source coupee en deux + arete cible coupee en deux
                        B1    A2     

              - A2 x----x-----x----x A1     arete source coupee en deux + arete cible coupee en deux
                        B1    B2     

              - B1 x----x-----x----x A1     arete source coupee en trois
                        B2    A2     

              - B1 x----x-----x----x A1     arete source coupee en trois
                        A2    B2     

            - S'il y a une intersection, on a 9 cas possibles :
            
              -          B1                 un sommet confondu
                         x
                        /             
                       /
                      /
                     /
                 A1 x---------x B2
                 A2

              -          B1                 un sommet confondu
                         x
                        /
                       /
                      /
                     /
                 A1 x---------x A2
                 B2

              -          B1                 arete cible coupee en deux
                         x
                        /
                       /
                      /
                     /
            A2 x----x----x B2
                   A1


              -          A2                 un sommet confondu
                         B1
                         x
                        / \
                       /   \
                      /     \
                     /       \
                 A1 x         x B2
                 

              -          B2                 un sommet confondu
                         B1
                         x
                        / \
                       /   \
                      /     \
                     /       \
                 A1 x         x A2


              -          
                         B1
                 A2 x----x----x B2          arete cible coupee en deux
                        / 
                       /   
                      /     
                     /       
                 A1 x   


              -              B1
                            x
                           /
                          /
                      A2 x--------x B2      arete source coupee en deux
                        / 
                       /     
                      /       
                  A1 x   


              -              B1
                            x
                             \
                              \
                   A2 x--------x B2         arete source coupee en deux
                                \
                                 \
                                  \
                                   x A1


              -            A1 
                           x
                           |
                           |
                   A2x-----|------xB2       arete source et arete cible coupees en deux
                           |
                           |
                           X
                           B1


            - Si sommet confondu :
              -------------------
              
              - On cree un nouveau sommet resultant de la moyenne des deux sommets et 
                on enregistre ses coordonnees.
              - On cree le lien dans le tableau Ancien sommet/Nouveau sommet
              - On incremente le nombre d'intersection totale


            - Si arete source coupee en deux :
              -------------------------------

              - On cree une nouvelle arete dans le tableau de connectivite Edge/Sommet 
                temporaire du maillage source
              - On met a jour les deux aretes en inserant le sommet du maillage cible 
                qui coupe l'arete en deux
              - Pour le sommet du maillage cible on met a jour la distance minimale.
              - On incremente le nombre d'intersections pour le maillage source

            - Si arete cible coupee en deux :
              -------------------------------

              - On cree une nouvelle arete dans le tableau des aretes a tester du 
                maillage cible
              - On met a jour les deux aretes en inserant le sommet du maillage source 
                qui coupe l'arete en deux
              - Pour le sommet du maillage source on met a jour la distance minimale.
              - On incremente le nombre d'intersections pour le maillage cible
           
           - Si arete cible et arete source coupees en deux :
             ------------------------------------------------

             - On cree un nouveau point representant l'intersection entre les deux aretes
               et on enregistre les coordonnees.
             - On lui associe une distance minimale a partir des aretes qui lui sont
               associees.
             - On applique les deux algorithmes precedents que l'on utilise lorsqu'une
               arete est coupee en deux.

           - si les deux aretes ont des sommets avec des indices identiques, on ignore l'intersection

         - Une fois que toutes les intersections entre les aretes sources et l'arete 
            cible ont ete calculees, on ajoute au tableau de connectivite Edge/Sommet 
            temporaire du maillage cible les aretes du tableau des aretes a tester.
            
     */


   for(int iEdgeTM = 0 ; iEdgeTM < nEdgeTM ; iEdgeTM ++){      

      testEdge[0] = edgeVertexConnectivityTM[2*iEdgeTM] + nVertSM;
      testEdge[1] = edgeVertexConnectivityTM[2*iEdgeTM + 1] + nVertSM;

      nEdgeInters = 1;


      for(int iEdgeSM = 0 ; iEdgeSM < nCutEdgeSM ; iEdgeSM ++){                

        indP1SM = edgeVertexConnectivitySMTmp[2*iEdgeSM] - 1;

        while(oldVertNewVert[indP1SM] > 0)
          indP1SM = oldVertNewVert[indP1SM] - 1;

        if( (indP1SM >= 0) && (indP1SM < nVertSM))
          p1EdgeSM =& (vertexCoordsSM[3*indP1SM]);
        else if((indP1SM >= nVertSM ) && (indP1SM < nVertSM + nVertTM))
          p1EdgeSM =& (vertexCoordsTM[3*(indP1SM - nVertSM)]);            
        else
          p1EdgeSM =& (_coordsIM[3*(indP1SM - nVertSM - nVertTM)]);

        indP2SM = edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] - 1;


        while(oldVertNewVert[indP2SM] > 0)
          indP2SM = oldVertNewVert[indP2SM] - 1;

        if( (indP2SM >= 0) && (indP2SM < nVertSM))
          p2EdgeSM =& (vertexCoordsSM[3*indP2SM]);          
        else if((indP2SM >= nVertSM ) && (indP2SM < nVertSM + nVertTM))
          p2EdgeSM =& (vertexCoordsTM[3*(indP2SM - nVertSM)]);            
        else
          p2EdgeSM =& (_coordsIM[3*(indP2SM - nVertSM - nVertTM)]);
                
        if(testEdge.size() <= 2*(nEdgeInters + 1))
          testEdge.resize(4*nEdgeInters);

        nEdgeIntersTmp = 0;

        for(int iEdgeInters = 0 ; iEdgeInters < nEdgeInters; iEdgeInters ++){
          if(edgeVertexConnectivitySMTmp.size() <= 2*(nCutEdgeSM + 2))
            edgeVertexConnectivitySMTmp.resize(4*(nCutEdgeSM ));

          if(newEdgeOldEdgeSM.size() <= nCutEdgeSM + 1)
            newEdgeOldEdgeSM.resize(2*nCutEdgeSM);



          if(edgeVertexConnectivityTMTmp.size() <= 2*(nEdgeInters + nEdgeIntersTmp + 2))
            edgeVertexConnectivityTMTmp.resize(4*(nEdgeInters + nEdgeIntersTmp) );

          if(_coordsIM.size() <= 3*(nIntersGlob + 2))
            _coordsIM.resize(6*(nIntersGlob + 2));

          indP1TM = testEdge[2*iEdgeInters] - 1;
 
          while(oldVertNewVert[indP1TM ] > 0)
           indP1TM = oldVertNewVert[indP1TM] - 1;
         
          if( (indP1TM >= 0) && (indP1TM < nVertSM))
            p1EdgeTM =& (vertexCoordsSM[3*indP1TM]);          
          else if((indP1TM >= nVertSM ) && (indP1TM < nVertSM + nVertTM))
            p1EdgeTM =& (vertexCoordsTM[3*(indP1TM - nVertSM)]);            
          else
            p1EdgeTM =& (_coordsIM[3*(indP1TM - nVertSM - nVertTM)]);
          
          
          indP2TM = testEdge[2*iEdgeInters + 1] - 1;


          while(oldVertNewVert[indP2TM] > 0)
            indP2TM = oldVertNewVert[indP2TM] - 1;
          
          if( (indP2TM >= 0) && (indP2TM < nVertSM))
            p2EdgeTM =& (vertexCoordsSM[3*indP2TM]);          
          else if((indP2TM >= nVertSM ) && (indP2TM < nVertSM + nVertTM))
            p2EdgeTM =& (vertexCoordsTM[3*(indP2TM - nVertSM)]);            
          else
            p2EdgeTM =& (_coordsIM[3*(indP2TM - nVertSM - nVertTM)]);
          
          nIntersLoc = intersectionEdge(distMinVertex[indP1SM],
                                        distMinVertex[indP2SM],
                                        distMinVertex[indP1TM],
                                        distMinVertex[indP2TM],
                                        p1EdgeSM,
                                        p2EdgeSM,
                                        p1EdgeTM,
                                        p2EdgeTM,
                                        0.01,
                                        &(locCoordInters[0]));

              
          if(indP1SM == indP1TM || indP1SM == indP2TM){
            if(nIntersLoc == 2){
              locCoordInters[0] = locCoordInters[2];
              locCoordInters[1] = locCoordInters[3];
            }
            nIntersLoc --;
          }
          
          if(indP2SM == indP1TM || indP2SM == indP2TM){
            doubleVertex = locCoordInters[0] < 0 ;
            if(nIntersLoc == 2 &&  !doubleVertex){              
              locCoordInters[0] = locCoordInters[2];
              locCoordInters[1] = locCoordInters[3];
            }
            nIntersLoc --;
          }
          

          /***************************************************************************/
          /*          if(nIntersLoc > 0){
            
            bftc_printf(" indP1SM %d indP2SM %d  indP1TM %d indP2TM %d \n",
                            indP1SM + 1,indP2SM + 1,
                          indP1TM + 1 ,indP2TM + 1);
            bftc_printf("coordsSM  P1 %f  %f  %f \n",p1EdgeSM[0],p1EdgeSM[1],p1EdgeSM[2]);
            bftc_printf("coordsSM  P2 %f  %f  %f \n",p2EdgeSM[0],p2EdgeSM[1],p2EdgeSM[2]);
            bftc_printf("coordsTM  P1 %f  %f  %f \n",p1EdgeTM[0],p1EdgeTM[1],p1EdgeTM[2]);
            bftc_printf("coordsTM  P2 %f  %f  %f \n",p2EdgeTM[0],p2EdgeTM[1],p2EdgeTM[2]);
            bftc_printf("distMin SM1  %f  SM2  %f  TM1  %f  TM2  %f \n",
                            distMinVertex[indP1SM],
                            distMinVertex[indP2SM],
                            distMinVertex[indP1TM],
                            distMinVertex[indP2TM]);
            
            bftc_printf("nbr inters %d \n",nIntersLoc);
          
            if(nIntersLoc == 2)
              bftc_printf("intersection   s1  %f  t1  %f \n s2  %f  t2  %f \n \n",
                              locCoordInters[0],
                              locCoordInters[1],
                              locCoordInters[2],
                              locCoordInters[3]);
            else if(nIntersLoc == 1)
              bftc_printf("intersection   s1  %f  t1  %f \n \n",
                              locCoordInters[0],
                            locCoordInters[1]);
                            }*/
          
          /***************************************************************************/

          
          
          if(nIntersLoc == 2){
            
            if( oldVertNewVert.size() <= nEdgeSM + nEdgeTM + nIntersGlob)
              oldVertNewVert.resize(2*(nEdgeSM + nEdgeTM + nIntersGlob),-1);              
              
              if(locCoordInters[0] < 0){
                if(locCoordInters[1] < 0){
                  if(locCoordInters[2] < 0){
                    bftc_error(__FILE__,
                                   __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                   locCoordInters[0],
                                   locCoordInters[1],
                                   locCoordInters[2],
                                   locCoordInters[3]);
                    
                    return false;
                  }
                  else if (locCoordInters[2] > 1){
                    if(locCoordInters[3] < 0){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                      return false;
                    }

                    else if (locCoordInters[3] > 1){
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;
                      
                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP1TM]);


                      nIntersGlob++;

                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP2SM],distMinVertex[indP2TM]);

                      
                      nIntersGlob++;

                                            
                    }
                    else{
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP1TM]);


                      nIntersGlob++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                      
                      nEdgeIntersTmp ++;
                      
                    }
                  }
                  else{
                    if (locCoordInters[3] > 1){
                      
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP1TM]);


                      nIntersGlob++;

                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM)]  =
                        indP2TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                      
                      nCutEdgeSM ++;
                      
                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);

                      return false;
                    }                    
                  }
                }
                else if (locCoordInters[1] > 1){
                  if (locCoordInters[2] < 0){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                    return false;
                  }
                  else if (locCoordInters[2] > 1){                                      
                    if(locCoordInters[3] < 0){
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;
                      
                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP2TM]);


                      nIntersGlob++;

                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;
                      
                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP2SM],distMinVertex[indP1TM]);


                      nIntersGlob++;                                           
                    }
                    else if (locCoordInters[3] > 1){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                      return false;
                    }
                    else{
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP2TM]);


                      nIntersGlob++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                      
                      nEdgeIntersTmp ++;

                    }
                  }
                  else{
                    if(locCoordInters[3] < 0){
                      _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP2TM]);


                      nIntersGlob++;

                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP1TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                                            
                      nCutEdgeSM ++;

                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                   "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);

                      return false;
                    }
                  }
                }
                else{
                  if (locCoordInters[2] < 0){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                      
                      return false;
                  }
                  else if (locCoordInters[2] > 1){
                    if(locCoordInters[3] < 0){
                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP1SM],distMinVertex[indP2TM]);


                      nIntersGlob++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;
                      
                      nEdgeIntersTmp ++;


                    }
                    else if (locCoordInters[3] > 1){
                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP2SM],distMinVertex[indP2TM]);


                      nIntersGlob++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;
                      
                      nEdgeIntersTmp ++;

                    }
                    else if(locCoordInters[1] < locCoordInters[3]){                      
                      
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp + 1) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp + 1)] = indP2SM + 1;

                      nEdgeIntersTmp +=2;
                      
                    }
                    else if(locCoordInters[3] < locCoordInters[1]){
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp + 1) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp + 1)] = indP1SM + 1;
                               
                      nEdgeIntersTmp +=2;
                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                      
                      return false;
                    }
                  }
                  else{
                    if(locCoordInters[3] < 0){
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP1TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                      
                      nCutEdgeSM ++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;
                      
                      nEdgeIntersTmp ++;

                    }
                    else if (locCoordInters[3] > 1){
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP2TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                      
                      nCutEdgeSM ++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;

                      nEdgeIntersTmp ++;
                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);

                      return false;
                    }
                  }                  
                }
              }
              else if (locCoordInters[0] > 1){
                if(locCoordInters[1] < 0){
                  if(locCoordInters[2] < 0 || locCoordInters[2] > 1){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                    return false;
                  }
                  else{
                    if(locCoordInters[3] > 1){
                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p1EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p1EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p1EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;

                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP2SM],distMinVertex[indP1TM]);


                      nIntersGlob++;

                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP2TM + 1;

                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;

                      nCutEdgeSM ++;



                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                      
                      return false;
                    }
                  }
                }

                else if (locCoordInters[1] > 1){

                  if(locCoordInters[2] < 0 || locCoordInters [2] > 1){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                    return false;
                  }
                  else{
                    if(locCoordInters[3] < 0){
                      
                      _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p2EdgeTM[0])/2 ;              
                      _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p2EdgeTM[1])/2 ;
                      _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p2EdgeTM[2])/2 ;                
                      
                      oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;
                      oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;
                      
                      distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                        std::min(distMinVertex[indP2SM],distMinVertex[indP1TM]);
                      

                      nIntersGlob++;

                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP1TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;

                      nCutEdgeSM ++;

                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                      
                      return false;
                    }
                  }
                }
                else{
                  if (locCoordInters[2] < 0 || locCoordInters[2] > 1){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                    return false;
                  }
                  else{
                    if(locCoordInters[3] < 0){
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP1TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                      

                      nCutEdgeSM ++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                                           
                      nEdgeIntersTmp ++;
                      
                      
                    }
                    else if (locCoordInters[3] > 1){
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                      edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;
                      edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                        indP2TM + 1;
                      
                      newEdgeOldEdgeSM[nCutEdgeSM ] =
                        newEdgeOldEdgeSM[iEdgeSM] ;
                      
                      nCutEdgeSM ++;

                      testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                        = testEdge[2*iEdgeInters + 1] ;
                      testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                      testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                      
                      nEdgeIntersTmp ++;
                      
                    }
                    else{
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                      return false;
                    }
                  }
                }
              }
              else{
                if(locCoordInters[1] < 0){
                  if (locCoordInters[2] < 0 || locCoordInters[2] > 1){
                      bftc_error(__FILE__,
                                     __LINE__, 0,
                                     "Combinaison impossible coords locales %f %f %f %f \n",
                                     locCoordInters[0],
                                     locCoordInters[1],
                                     locCoordInters[2],
                                     locCoordInters[3]);
                    
                    return false;
                  }
                  
                  else{
                    if (locCoordInters[3] > 1){
                      if(locCoordInters[0] < locCoordInters[2] ){
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM 
                                                        + 1) + 1] =  
                          edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                        
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                        
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                          indP1TM + 1;

                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1]  =
                          indP2TM + 1;

                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM  + 1)]  =
                          indP2TM + 1;

                        
                        newEdgeOldEdgeSM[nCutEdgeSM ] =
                          newEdgeOldEdgeSM[iEdgeSM] ;

                        newEdgeOldEdgeSM[nCutEdgeSM  + 1] =
                          newEdgeOldEdgeSM[iEdgeSM] ;

                        
                        nCutEdgeSM +=2;
                      }
                      else if(locCoordInters[2] < locCoordInters[0] ){
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM 
                                                        + 1) + 1] =  
                          edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];                        
                        edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;                        
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  =
                          indP2TM + 1;
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1]  =
                          indP1TM + 1;
                        edgeVertexConnectivitySMTmp[2*(nCutEdgeSM  + 1)]  =
                          indP1TM + 1;
                        
                        newEdgeOldEdgeSM[nCutEdgeSM ] =
                          newEdgeOldEdgeSM[iEdgeSM] ;

                        newEdgeOldEdgeSM[nCutEdgeSM  + 1] =
                          newEdgeOldEdgeSM[iEdgeSM] ;

                        
                        nCutEdgeSM +=2;
                      }
                      else{
                        bftc_error(__FILE__,
                                       __LINE__, 0,
                                       "Combinaison impossible coords locales %f %f %f %f \n",
                                       locCoordInters[0],
                                       locCoordInters[1],
                                       locCoordInters[2],
                                       locCoordInters[3]);
                        
                        return false;
                      }                      
                    }
                    else{
                        bftc_error(__FILE__,
                                       __LINE__, 0,
                                       "Combinaison impossible coords locales %f %f %f %f \n",
                                       locCoordInters[0],
                                       locCoordInters[1],
                                       locCoordInters[2],
                                       locCoordInters[3]);
                                              
                      return false;
                    }                      
                  }
                }
                else{
                  bftc_error(__FILE__,
                                 __LINE__, 0,
                                 "Combinaison impossible coords locales %f %f %f %f \n",
                                 locCoordInters[0],
                                 locCoordInters[1],
                                 locCoordInters[2],
                                 locCoordInters[3]);
                  
                  return false;
                }                      
              }              
            }
            
            
            else if(nIntersLoc == 1){
              if( oldVertNewVert.size() <= nEdgeSM + nEdgeTM + nIntersGlob)
                oldVertNewVert.resize(2*oldVertNewVert.size(),-1);
              
              
              if(locCoordInters[0] < 0){
                if(locCoordInters[1] < 0){
                  _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p1EdgeTM[0])/2 ;              
                  _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p1EdgeTM[1])/2 ;
                  _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p1EdgeTM[2])/2 ;                
                  
                  oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;
                  oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;
                                    
                  distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                    std::min(distMinVertex[indP1TM],distMinVertex[indP1SM]);
                                    
                  nIntersGlob++;
                }
                else if(locCoordInters[1] > 1){
                  _coordsIM[3*nIntersGlob] = (p1EdgeSM[0] + p2EdgeTM[0])/2 ;              
                  _coordsIM[3*nIntersGlob + 1] = (p1EdgeSM[1] + p2EdgeTM[1])/2 ;
                  _coordsIM[3*nIntersGlob + 2] = (p1EdgeSM[2] + p2EdgeTM[2])/2 ;
                  
                  oldVertNewVert[indP1SM] = nIntersGlob + nVertSM + nVertTM + 1;
                  oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;

                  distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                    std::min(distMinVertex[indP1SM],distMinVertex[indP2TM]);
                 
                  nIntersGlob++;
                  
                }
                else {
                  
                  testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] =
                    testEdge[2*iEdgeInters + 1] ;
                  testEdge[2*iEdgeInters + 1] = indP1SM + 1;
                  testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP1SM + 1;
                                    
                  nEdgeIntersTmp ++;
                }
              }
              else if(locCoordInters[0] > 1){
                if(locCoordInters[1] < 0){
                  _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p1EdgeTM[0])/2 ;              
                  _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p1EdgeTM[1])/2 ;
                  _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p1EdgeTM[2])/2 ;
                  
                  oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;
                  oldVertNewVert[indP1TM] = nIntersGlob + nVertSM + nVertTM + 1;                
                  
                  distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                    std::min(distMinVertex[indP2SM],distMinVertex[indP1TM]);

                  nIntersGlob++;
                  
                }
                else if(locCoordInters[1] > 1){
                  _coordsIM[3*nIntersGlob] = (p2EdgeSM[0] + p2EdgeTM[0])/2 ;              
                  _coordsIM[3*nIntersGlob + 1] = (p2EdgeSM[1] + p2EdgeTM[1])/2 ;
                  _coordsIM[3*nIntersGlob + 2] = (p2EdgeSM[2] + p2EdgeTM[2])/2 ;
                  
                  oldVertNewVert[indP2SM] = nIntersGlob + nVertSM + nVertTM + 1;
                  oldVertNewVert[indP2TM] = nIntersGlob + nVertSM + nVertTM + 1;

                  distMinVertex[nIntersGlob + nVertSM + nVertTM] = 
                    std::min(distMinVertex[indP2SM],distMinVertex[indP2TM]);

                  nIntersGlob++;
                  
                }
                else {
                  testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] = testEdge[2*iEdgeInters + 1] ;
                  testEdge[2*iEdgeInters + 1] = indP2SM + 1;
                  testEdge[2*(nEdgeInters + nEdgeIntersTmp)] = indP2SM + 1;
                                   
                  nEdgeIntersTmp ++;                
                }
              }
              else if (locCoordInters[1] < 0 ){
                if((locCoordInters[0] > 0) && (locCoordInters[0] < 1)){
                  
                  edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] =  edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                  edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP1TM + 1;
                  edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  = indP1TM + 1;
                  
                  newEdgeOldEdgeSM[nCutEdgeSM ] = newEdgeOldEdgeSM[iEdgeSM] ;

                  nCutEdgeSM ++;
                }
              }
              else if (locCoordInters[1] > 1 ){
                if( (locCoordInters[0] > 0) && (locCoordInters[0] < 1)){
                  edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] 
                    =  edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];                  
                  edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] = indP2TM + 1;
                  edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  
                    = indP2TM + 1;
                  
                  newEdgeOldEdgeSM[nCutEdgeSM ] 
                    = newEdgeOldEdgeSM[iEdgeSM] ;

                  nCutEdgeSM ++;
                }
              }
              
              else{
                
                if(distMinVertex.size() <= nIntersGlob + nVertSM + nVertTM)
                  distMinVertex.resize( 2*distMinVertex.size());
                
                _coordsIM[3*nIntersGlob] = (1 - locCoordInters[0])*p1EdgeSM[0] 
                  + locCoordInters[0]* p2EdgeSM[0] ;              
                _coordsIM[3*nIntersGlob + 1] = (1 - locCoordInters[0])*p1EdgeSM[1] 
                  + locCoordInters[0]* p2EdgeSM[1] ;              
                _coordsIM[3*nIntersGlob + 2] = (1 - locCoordInters[0])*p1EdgeSM[2] 
                  + locCoordInters[0] * p2EdgeSM[2] ;
                
                
                oldVertNewVert[nVertTM + nVertSM + nIntersGlob] = 0;
                
                testEdge[2*(nEdgeInters + nEdgeIntersTmp) + 1] 
                  = testEdge[2*iEdgeInters + 1] ; 
                testEdge[2*iEdgeInters + 1] = nVertSM + nVertTM + nIntersGlob + 1;
                testEdge[2*(nEdgeInters + nEdgeIntersTmp)] 
                  = nVertSM + nVertTM + nIntersGlob + 1;
                
                edgeVertexConnectivitySMTmp[2*(nCutEdgeSM ) + 1] 
                  =  edgeVertexConnectivitySMTmp[2*iEdgeSM + 1];
                edgeVertexConnectivitySMTmp[2*iEdgeSM + 1] 
                  = nVertSM + nVertTM + nIntersGlob + 1;
                edgeVertexConnectivitySMTmp[2*(nCutEdgeSM )]  
                  = nVertSM + nVertTM + nIntersGlob + 1;
                
                newEdgeOldEdgeSM[nCutEdgeSM ] 
                  = newEdgeOldEdgeSM[iEdgeSM] ;
                
                
                
                distMinP1SM = norm(&(_coordsIM[3*nIntersGlob]),
                                                     p1EdgeSM);

                distMinP2SM = norm(&(_coordsIM[3*nIntersGlob]),
                                               p2EdgeSM);

                distMinP1TM = norm(&(_coordsIM[3*nIntersGlob]),
                                               p1EdgeTM);

                distMinP2TM = norm(&(_coordsIM[3*nIntersGlob]),
                                               p2EdgeTM);
                
                distMinVertex[nVertTM + nVertSM + nIntersGlob] =
                  _tolerance*
                  std::min(distMinP1SM,
                      std::min(distMinP2SM,
                          std::min(distMinP1TM,distMinP2TM)));


                nCutEdgeSM ++;
                nEdgeIntersTmp ++;
                nIntersGlob ++;              
                
              }
            }
        }

        nEdgeInters += nEdgeIntersTmp;              

        if(nEdgeIntersTmp > 0) 
          iEdgeSM = 0;      
      }
           
      if( newEdgeOldEdgeTM.size() <= nCutEdgeTM + nEdgeInters)
        newEdgeOldEdgeTM.resize(2*(nCutEdgeTM + nEdgeInters));

      if(edgeVertexConnectivityTMTmp.size() <=  2*(nCutEdgeTM + nEdgeInters + 1))
        edgeVertexConnectivityTMTmp.resize(2*edgeVertexConnectivityTMTmp.size());

      for (int iEdge = 0 ; iEdge < nEdgeInters ; iEdge ++){
        edgeVertexConnectivityTMTmp[2*(nCutEdgeTM + iEdge)] = testEdge[2*iEdge];
        edgeVertexConnectivityTMTmp[2*(nCutEdgeTM + iEdge) + 1] = testEdge[2*iEdge + 1];

        newEdgeOldEdgeTM[nCutEdgeTM + iEdge] = iEdgeTM + 1;
      }
      
      nCutEdgeTM += nEdgeInters;
    }

    nVertex = 0;

    
    /**********   Construction des tableaux   **********/

    /********  ALGO  Coordonnees    ***********/
    /* 
      - On copie dans un tableau temporaire les coordonnees trouvees lors d'intersection 
      Arete/Arete ou Sommet/Sommet.

      - On parcourt d'abord le tableau des coordonnees du maillage source ; si dans le 
      tableau de mise a jour des sommets, on a un indice nul, alors on sait que nous 
      avons un sommet du maillage intersecte.

      - On fait la meme chose pour le tableau des coordonnes du maillage cible et le
      tableau temporaire.
     */

    std::vector<double> coordsInters(_coordsIM);


    if(_coordsIM.size() <= 3*(nIntersGlob + nVertSM + nVertTM))
      _coordsIM.resize(3*(nIntersGlob + nVertSM + nVertTM));

    for(int iVert = 0 ; iVert < nVertSM ; iVert ++){
      if( oldVertNewVert[iVert] == 0) {
        _coordsIM[3*nVertex] = vertexCoordsSM[3*iVert];          
        _coordsIM[3*nVertex + 1] = vertexCoordsSM[3*iVert + 1];          
        _coordsIM[3*nVertex + 2] = vertexCoordsSM[3*iVert + 2];   
       
        oldVertNewVert[iVert] = - nVertex - 1;

        nVertex ++;
      }
    }

    for(int iVert = 0 ; iVert < nVertTM ; iVert ++){
      if( oldVertNewVert[iVert + nVertSM] == 0) {
        _coordsIM[3*nVertex] = vertexCoordsTM[3*iVert];          
        _coordsIM[3*nVertex + 1] = vertexCoordsTM[3*iVert + 1];          
        _coordsIM[3*nVertex + 2] = vertexCoordsTM[3*iVert + 2];          

        oldVertNewVert[iVert + nVertSM] = - nVertex - 1;

        nVertex ++;
      }
    }

    for(int iVert = 0 ; iVert < nIntersGlob; iVert ++){
      if( oldVertNewVert[iVert + nVertSM + nVertTM] == 0) {
        _coordsIM[3*nVertex] = coordsInters[3*iVert];          
        _coordsIM[3*nVertex + 1] = coordsInters[3*iVert + 1];          
        _coordsIM[3*nVertex + 2] = coordsInters[3*iVert + 2];          

        oldVertNewVert[iVert + nVertSM + nVertTM] = - nVertex - 1;

         nVertex ++;
      }
    }

    coordsInters.clear();
 


    /***********************************/

    /*    bftc_printf("\n coords IM \n");

    for(int i = 0 ; i < nVertex ; i++)
      bftc_printf(" %d   %f  %f  %f \n",i + 1,
                      _coordsIM[3*i],_coordsIM[3*i + 1],_coordsIM[3*i + 2]);

                      bftc_printf("------------------end coords--------------------- \n");*/

   /************************************/







    /********   ALGO  Tableau de connectivite Arete/Sommet    ***********/
   /*
     - On met tout d'abord a jour les indices de sommet dans les tableaux de connectivite
     Arete/Sommet temporaires du maillage source et cible.
     - On utilise ensuite une table de hachage afin de trier les aretes des deux tableaux
       de connectivite temporaire. On applique la meme methode que nous avions fait
       precedemment.
     - Une fois que la table de hachage est cree et transformee en un tableau d'index 
       selon les sommes des indices des aretes, on remplit le tableau de connectivite
       Arete/Sommet du maillage intersecte. Or nous aurons besoin de garder un lien entre
       le tableau de connectivite temporaire et le nouveau ; pour cela nous avons cree
       le tableau classifyEdge qui renvoie le nouvel indice de l'arete, avec d'abord
       les aretes du maillage cible, puis celles du maillage sources. 
     - Des plus on elimine les mauvaises aretes (aretes doubles et aretes avec un seul 
       point dus a l'algorithme) :
       - Si l'arete est en double, on "copie" les informations de l'arete deja enregistree
         et on fait attention au cas ou les deux aretes sont dans le sens contraire.
       - Si l'arete possede un seul point, on l'elimine en ne l'enregistrant pas dans le
         tableau et on initialise la valeur de son arete mere a 0.
     - La derniere partie consiste a eliminer les trous dans le nouveau tableau de 
       connectivite qui sont dus aux mauvaises aretes non enregistrees et a mettre a jour
       le tableau classifyEdge puisque les indices du nouveau tableau ont change.
    */


    vertexEdgeIndexIM.resize(nVertex + 1);
    edgeVertexConnectivityTMTmp.resize(2*(nCutEdgeSM + nCutEdgeTM + nIntersGlob));
    classifyEdge.resize(nCutEdgeSM + nCutEdgeTM);

    for(int iEdge = 0 ; iEdge < nCutEdgeTM ; iEdge ++){

      indP1TM = edgeVertexConnectivityTMTmp[2*iEdge];

      while(oldVertNewVert[indP1TM - 1] > 0)
        indP1TM = oldVertNewVert[indP1TM - 1];

      if(oldVertNewVert[indP1TM - 1] < 0)
        indP1TM = - oldVertNewVert[indP1TM - 1];

      edgeVertexConnectivityTMTmp[2*iEdge] = indP1TM;
            
      indP2TM = edgeVertexConnectivityTMTmp[2*iEdge + 1];

      while(oldVertNewVert[indP2TM - 1] > 0)
        indP2TM = oldVertNewVert[indP2TM - 1];

      if(oldVertNewVert[indP2TM - 1] < 0)
        indP2TM = -oldVertNewVert[indP2TM - 1];

      edgeVertexConnectivityTMTmp[2*iEdge + 1] = indP2TM;
      
    }


    for (int iEdge = 0 ; iEdge < nCutEdgeSM ; iEdge ++){

      indP1SM = edgeVertexConnectivitySMTmp[2*iEdge];

      while(oldVertNewVert[indP1SM - 1] > 0)
        indP1SM = oldVertNewVert[indP1SM - 1];

      if(oldVertNewVert[indP1SM - 1] < 0)
        indP1SM = - oldVertNewVert[indP1SM - 1];

      edgeVertexConnectivitySMTmp[2*iEdge] = indP1SM;

      indP2SM = edgeVertexConnectivitySMTmp[2*iEdge + 1];

      while(oldVertNewVert[indP2SM - 1] > 0)
        indP2SM = oldVertNewVert[indP2SM - 1];

      if(oldVertNewVert[indP2SM - 1] < 0)
        indP2SM = - oldVertNewVert[indP2SM - 1];

      edgeVertexConnectivitySMTmp[2*iEdge + 1] = indP2SM;      

    }

    sumVertKeyMax = 0;
    
    hachTab.resize(2*nVertex + 1);
    edgeVertexConnectivitySMTmp.resize(2*nCutEdgeSM);
    edgeVertexConnectivityIM.resize(2*(nCutEdgeSM + nCutEdgeTM));

    for(int iEdge = 0; iEdge < nCutEdgeTM ; iEdge ++){
      sumVertKey = edgeVertexConnectivityTMTmp[2*iEdge] + 
        edgeVertexConnectivityTMTmp[2*iEdge + 1];

      hachTab[sumVertKey - 1] ++;
      
      sumVertKeyMax = sumVertKey > sumVertKeyMax ? sumVertKey : sumVertKeyMax;
    }

    for(int iEdge = 0; iEdge < nCutEdgeSM ; iEdge ++){
      sumVertKey = edgeVertexConnectivitySMTmp[2*iEdge] + 
        edgeVertexConnectivitySMTmp[2*iEdge + 1];

      hachTab[sumVertKey - 1] ++;
      
      sumVertKeyMax = sumVertKey > sumVertKeyMax ? sumVertKey : sumVertKeyMax;
    }

    hachTab.resize(sumVertKeyMax + 2);

    nEdge = 0;
    
    for(int i = 0 ; i < sumVertKeyMax + 2; i ++){
      sumTmp = hachTab[i];
      hachTab[i] = nEdge;
      nEdge += sumTmp;      
    }

    stockTabTmp.resize(sumVertKeyMax + 2);
    doubleEdge.resize(sumVertKeyMax + 2);
    
    nEdge = 0;
    
    for(int iEdge = 0 ; iEdge < nCutEdgeTM ; iEdge ++){
      sumVertKey = edgeVertexConnectivityTMTmp[2*iEdge] + 
        edgeVertexConnectivityTMTmp[2*iEdge + 1];

      alreadyCreated = false;

      if(edgeVertexConnectivityTMTmp[2*iEdge] == 
         edgeVertexConnectivityTMTmp[2*iEdge + 1]){
          alreadyCreated = true;
          
          classifyEdge[iEdge] = hachTab[sumVertKey - 1]  + 1;

          newEdgeOldEdgeTM[iEdge] = 0;

          doubleEdge[sumVertKey] ++;

      }           
      else
        for(int i = 0 ; i < stockTabTmp[sumVertKey - 1] ; i++)
          if(edgeVertexConnectivityTMTmp[2*iEdge] == 
             edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i)] 
             || edgeVertexConnectivityTMTmp[2*iEdge] == 
             edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i) + 1]){
            
            alreadyCreated = true;
            
            classifyEdge[iEdge ] = hachTab[sumVertKey - 1] + i + 1;
            
            doubleEdge[sumVertKey] ++;

            if(edgeVertexConnectivityTMTmp[2*iEdge] == 
               edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i) + 1])
              newEdgeOldEdgeTM[iEdge] = -newEdgeOldEdgeTM[iEdge];
          }
      
      
      if(!alreadyCreated){
        edgeVertexConnectivityIM[ 2*(hachTab[sumVertKey - 1] 
                                     + stockTabTmp[sumVertKey - 1])] 
          = edgeVertexConnectivityTMTmp[2*iEdge];
        
        edgeVertexConnectivityIM[ 2*(hachTab[sumVertKey - 1] 
                                     + stockTabTmp[sumVertKey - 1]) + 1] 
          = edgeVertexConnectivityTMTmp[2*iEdge + 1];

        classifyEdge[iEdge] = hachTab[sumVertKey - 1] + stockTabTmp[sumVertKey - 1] + 1;
        
        stockTabTmp[sumVertKey - 1 ] ++;            
          
        nEdge ++;
      }        
    }

    
    for(int iEdge = 0 ; iEdge < nCutEdgeSM ; iEdge ++){
      sumVertKey = edgeVertexConnectivitySMTmp[2*iEdge] + 
        edgeVertexConnectivitySMTmp[2*iEdge + 1];

      alreadyCreated = false;

      if(edgeVertexConnectivitySMTmp[2*iEdge] == 
         edgeVertexConnectivitySMTmp[2*iEdge + 1]){

        alreadyCreated = true;
          
        classifyEdge[iEdge + nCutEdgeTM] = hachTab[sumVertKey - 1] +  1;

        newEdgeOldEdgeSM[iEdge] = 0;
            
        doubleEdge[sumVertKey] ++;
          

      }

      else
        for(int i = 0 ; i < stockTabTmp[sumVertKey - 1] ; i++)
          if(edgeVertexConnectivitySMTmp[2*iEdge] == 
             edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i)] 
             || edgeVertexConnectivitySMTmp[2*iEdge] == 
             edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i) + 1]){
            
            alreadyCreated = true;
            
            classifyEdge[iEdge + nCutEdgeTM] = hachTab[sumVertKey - 1] + i + 1;
            
            doubleEdge[sumVertKey] ++;
            
            if(edgeVertexConnectivitySMTmp[2*iEdge] == 
               edgeVertexConnectivityIM[2*(hachTab[sumVertKey - 1] + i) + 1])
              newEdgeOldEdgeSM[iEdge] = -newEdgeOldEdgeSM[iEdge];
            
          }
      

      if(!alreadyCreated){
        edgeVertexConnectivityIM[ 2*(hachTab[sumVertKey - 1] 
                                     + stockTabTmp[sumVertKey - 1])] 
          = edgeVertexConnectivitySMTmp[2*iEdge];
        
        edgeVertexConnectivityIM[ 2*(hachTab[sumVertKey - 1] 
                                     + stockTabTmp[sumVertKey - 1]) + 1] 
          = edgeVertexConnectivitySMTmp[2*iEdge + 1];
        
        classifyEdge[iEdge + nCutEdgeTM] = 
          hachTab[sumVertKey - 1] + stockTabTmp[sumVertKey - 1] + 1;
        
        stockTabTmp[sumVertKey - 1] ++;            
        
        nEdge++;
      }        
          
    }

    for(int iSum = 0 ; iSum < sumVertKeyMax + 1 ; iSum++){      
      for(int iEdge = 0 ; iEdge < stockTabTmp[iSum]; iEdge ++){

        edgeVertexConnectivityIM[ 2*(hachTab[iSum] 
                                     + iEdge - doubleEdge[iSum])] =
          edgeVertexConnectivityIM[ 2*(hachTab[iSum] 
                                       + iEdge)];

        edgeVertexConnectivityIM[ 2*(hachTab[iSum] 
                                     + iEdge - doubleEdge[iSum]) + 1] =
          edgeVertexConnectivityIM[ 2*(hachTab[iSum] 
                                       + iEdge) + 1];

      }
      doubleEdge[iSum + 1] += doubleEdge[iSum];
    }


    edgeVertexConnectivityTMTmp.clear();   
    edgeVertexConnectivitySMTmp.clear();   
    

    for (int iEdge = 0 ; iEdge < nCutEdgeTM ; iEdge ++){
      sumVertKey = edgeVertexConnectivityTMTmp[2*iEdge] + 
        edgeVertexConnectivityTMTmp[2*iEdge + 1];
      
      classifyEdge[iEdge] -= doubleEdge[sumVertKey - 1]; 
    }


    for (int iEdge = 0 ; iEdge < nCutEdgeSM ; iEdge ++){
      sumVertKey = edgeVertexConnectivitySMTmp[2*iEdge] + 
        edgeVertexConnectivitySMTmp[2*iEdge + 1];

      classifyEdge[iEdge + nCutEdgeTM] -= doubleEdge[sumVertKey - 1]; 
    }


    stockTabTmp.clear();


    /*****************************************************************/

    /*    bftc_printf("edgeVertexConnectivity \n");
    for(int i = 0; i < nEdge ; i++)
      bftc_printf("%d  %d  %d \n",
             i + 1, 
             edgeVertexConnectivityIM[2*i],
             edgeVertexConnectivityIM[2*i + 1]);
             bftc_printf("--------end edgeVertexConnectivity --------\n");*/


    /*****************************************************************/



    /********   ALGO  Tableau de connectivite Arete/Sommet    ***********/
   /*
     - On cherche a construire le tableau Ancienne arete / Nouvelle Arete qui est 
       necessaire pour la construction des elements du maillage intersecte.
     - On calcule d'abord le tableau d'index a partir du tableau 
       Nouvelle arete/Ancienne arete en ne prenant pas en compte les indices nuls (arete
       a un sommet).
     - A partir de ce meme tableau on en deduit le tableau voulu en faisant attention
       au signe de l'arete mere et fille qui correspondent a leur sens.
    */
          
    oldEdgeNewEdgeIndexIM.resize(nEdgeSM + nEdgeTM + 1,0);
    edgeMeshTag.resize(nEdge);

    for (int iEdge = 0; iEdge < nCutEdgeTM ; iEdge ++)
      if(std::abs(newEdgeOldEdgeTM[iEdge]) > 0)
        oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeTM[iEdge]) - 1] ++;

    for (int iEdge = 0; iEdge < nCutEdgeSM ; iEdge ++)
      if(std::abs(newEdgeOldEdgeSM[iEdge]) > 0)
        oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeSM[iEdge]) + nEdgeTM - 1 ] ++;


    int sum = 0;
    int swapTmp;

    for(int iEdge = 0; iEdge < nEdgeSM + nEdgeTM + 1; iEdge ++ ){
      swapTmp = oldEdgeNewEdgeIndexIM[iEdge];
      oldEdgeNewEdgeIndexIM[iEdge] = sum;
      sum += swapTmp ;      
    }
    

    stockTabTmp.resize(nEdgeSM + nEdgeTM,0);

    oldEdgeNewEdgeConnectivityIM.resize(2*(nCutEdgeTM + nCutEdgeSM));
    

    for (int iEdge = 0; iEdge < nCutEdgeTM ; iEdge ++)
      if(std::abs(newEdgeOldEdgeTM[iEdge]) > 0){
        signe = newEdgeOldEdgeTM[iEdge] > 0 ? 1 : -1;
        
        oldEdgeNewEdgeConnectivityIM[ oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeTM[iEdge]) - 1] +
                                      stockTabTmp[std::abs(newEdgeOldEdgeTM[iEdge]) - 1]] 
          = signe * classifyEdge[iEdge];

        if(oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeTM[iEdge])] - oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeTM[iEdge]) - 1] > 1)
          edgeMeshTag[std::abs(classifyEdge[iEdge]) - 1] = 3;
        else
          edgeMeshTag[std::abs(classifyEdge[iEdge]) - 1] = 2;

        stockTabTmp[std::abs(newEdgeOldEdgeTM[iEdge]) - 1] ++;
      }

    for (int iEdge = 0; iEdge < nCutEdgeSM ; iEdge ++)
      if(std::abs(newEdgeOldEdgeSM[iEdge]) > 0){

        signe = newEdgeOldEdgeSM[iEdge] > 0 ? 1 : -1;
        
        oldEdgeNewEdgeConnectivityIM[ oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeSM[iEdge]) + 
                                                            nEdgeTM - 1] + 
                                    stockTabTmp[std::abs(newEdgeOldEdgeSM[iEdge]) 
                                                     + nEdgeTM - 1]] 
          = signe*classifyEdge[iEdge  + nCutEdgeTM];


        if( oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeSM[iEdge]) + 
                                                            nEdgeTM] 
            - oldEdgeNewEdgeIndexIM[std::abs(newEdgeOldEdgeSM[iEdge]) + 
                                                            nEdgeTM - 1] > 1 
            || edgeMeshTag[std::abs(classifyEdge[iEdge + nCutEdgeTM]) - 1] == 2)

          edgeMeshTag[std::abs(classifyEdge[iEdge + nCutEdgeTM]) - 1] = 3;
        else
          edgeMeshTag[std::abs(classifyEdge[iEdge + nCutEdgeTM]) - 1] = 1;

        
        stockTabTmp[std::abs(newEdgeOldEdgeSM[iEdge])  + nEdgeTM - 1] ++;
      }
    
    stockTabTmp.clear();
    

   /********   ALGO  Tableau de connectivite Sommet/Arete    ***********/
   /*
     - On construit le tableau d'index comme dans les cas precedents a l'aide du tableau
       de connectivite Arete/Sommet.
     - On calcule ensuite le tableau de connectivite Sommet/Arete avec un indice d'arete 
       positif pour le premier sommet et negatif pour le second.

   */

    for(int iEdge = 0 ; iEdge < nEdge ; iEdge ++){
      vertexEdgeIndexIM[edgeVertexConnectivityIM[2*iEdge] - 1] ++;
      vertexEdgeIndexIM[edgeVertexConnectivityIM[2*iEdge + 1] - 1] ++;      
    }

    sum = 0;
    
    for(int iVert = 0; iVert < nVertex + 1 ; iVert ++ ){
      swapTmp = vertexEdgeIndexIM[iVert];
      vertexEdgeIndexIM[iVert] = sum;
      sum += swapTmp ;      
    }
    
    stockTabTmp.resize(nVertex);

   
    vertexEdgeConnectivityIM.resize(vertexEdgeIndexIM[nVertex]);

    for (int iEdge = 0; iEdge < nEdge ; iEdge ++){
      
      vertexEdgeConnectivityIM[ vertexEdgeIndexIM[edgeVertexConnectivityIM[2*iEdge] - 1]
                                + stockTabTmp[edgeVertexConnectivityIM[2*iEdge] - 1]]
        = (iEdge + 1);

      stockTabTmp[edgeVertexConnectivityIM[2*iEdge] - 1] ++;
                                 
      vertexEdgeConnectivityIM[ vertexEdgeIndexIM[edgeVertexConnectivityIM[2*iEdge + 1] 
                                                  - 1]
                                + stockTabTmp[edgeVertexConnectivityIM[2*iEdge + 1] 
                                                   - 1]]
        = -(iEdge + 1);

      stockTabTmp[edgeVertexConnectivityIM[2*iEdge + 1] - 1]++;

    }
    



    /**************   Libération des donnees  ************************/

    testEdge.clear();
    locCoordInters.clear();
    edgeVertexConnectivitySMTmp.clear();
    oldVertNewVert.clear();
    newEdgeOldEdgeSM.clear();
    newEdgeOldEdgeTM.clear();
    distMinVertex.clear();
    stockTabTmp.clear();
    hachTab.clear();

    return true;

  }




  
  bool ConservativeMesh::buildIntersectionMesh(
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
                                               std::vector<int>& _eltVertConnectivityIM,
                                               std::vector<int>& _eltVertIndexIM,
                                               int& nElts) {
    


    const int *eltVertIndexSM = _sourceMesh.getEltConnectivityIndex();           //Tableau d'index Element/Sommet du maillage source
    const int *eltVertConnectivitySM = _sourceMesh.getEltConnectivity();         //Tableau de connectivite Element/Sommet du maillage source

    const int *eltVertIndexTM = _targetMesh.getEltConnectivityIndex();           //Tableau d'index Element/Sommet du maillage cible
    const int *eltVertConnectivityTM = _targetMesh.getEltConnectivity();         //Tableau de connectivite Element/Sommet du maillage cible

    std::vector<double> normalFaceSM(3);                                         //Tableau de coordonnees des normales d'une face du maillage source
    std::vector<double> normalFaceTM(3);                                         //Tableau de coordonnees des normales d'une face du maillage cible
    std::vector<int> tagEdge;                                                    //Tableau de marquage d'aretes. Indice pair dans le meme sens, 
                                                                                 //indice impair dans le sens contraire
    std::vector<int> stackEdge;                                                  //Tableau de stockage de donnees
    std::vector<bool> borderEdge(nEdge,false);                                   //True si l'arete i est au bord de l'element pere
    std::vector<bool> edgeElementTag(nEdge,false);
    std::vector<int> doubleEdge(nEdge);                                          //Indice des aretes en double.
    std::vector<int> quadrangleElt;                                              //Tableau de connectivite Element quadrangle/Sommets 
    std::vector<int> newEltOldEltQuad;                                           //Tableau de liaison entre element quadrangle fils et element pere
    std::vector<int> polygonEltIndex;                                            //Tableau d'index Element polygonal/Sommets
    std::vector<int> polygonElt;                                                 //Tableau de connectivite Element polygonal/Sommets
    std::vector<int> newEltOldEltPol;                                            //Tableau de liaison entre element polygonal fils et element pere
    std::vector<int> indElt;                                                     //Tableau de liaison entre indice general et indice du tableau de connectivite Element/Sommets
    std::vector<int> newElt(10);                                                 //Tableau contenant les indices de l'element en creation
    std::vector<double> vectCrossProduct(3);                                         //Tableau de coordonnees du produit vectoriel
    std::vector<double> vect1(3);                                                //Tableau de coordonnees du vecteur 1
    std::vector<double> vect2(3);                                                //Tableau de coordonnees du vecteur 2
    std::vector<double> vect1Proj(3);                                            //Tableau de coordonnees de la projection du vecteur 1
    std::vector<double> vect2Proj(3);                                            //Tableau de coordonnees de la projection du vecteur 2


    int nEltSM = _sourceMesh.getNElts();                                         //Nombre d'elements du maillage source
    int nEltTM = _targetMesh.getNElts();                                         //Nombre d'elements du maillage cible
    int nEdgeElt;                                                                //Nombre d'aretes dans l'element
    int nStackEdge;                                                              //Nombre d'arete dans le tableau de stockage
    int iStackEdge;                                                              //Position de l'indice dans le tableau de stockage
    int oldEdge;                                                                 //Indice de l'arete mere
    int nNewEdge;                                                                //Nombre d'aretes filles pour une arete mere
    int prevVert;                                                                //Indice du sommet precedent
    int intVert;                                                                 //Indice du sommet actuel
    int nextVert;                                                                //Indice du sommet suivant 
    int nextVertTmp;                                                             //Indice du sommet suivant potentiel
    int firstVertElt;                                                            //Premier sommet de l'element en creation
    int signe;                                                                   //Signe
    int newEdge;                                                                 //Indice de l'arete
    int nextEdge;                                                                //Indice de l'arete suivante
    int nextEdgeTmp;
    int nIntersEdge;                                                             //Nombre d'aretes connectees a un sommet
    int direction;                                                               //Direction de l'arete actuelle
    int directionLeft;                                                           //Direction de l'arete la plus a gauche par rapport a l'arete actuelle
    int directionBorder;                                                         //Direction de l'arete au bord de l'element pere
    int nVertNewElt;                                                             //Nombre de sommets dans le nouvel element
    int nTrianElt = 0;                                                           //Nombre d'elements triangulaires
    int nQuadElt = 0;                                                            //Nombre d'elements quadrangles
    int nPolElt = 0;                                                             //Nombre d'elements polygonaux 
    int indTag;                                                                  //Valeur du tableau tagEdge

    double cosLeft;                                                              //Valeur du cosinus entre l'arete actuelle et l'arete la plus a gauche
    double cosBorder;                                                            //Valeur du cosinus entre l'arete actuelle et l'arete au bord de l'element
    double normeNormaleFace;                                                     //Norme de la normale a l'element pere
    double cosV1V2;                                                              //Valeur du cosinus entre le vecteur 1 et 2
    double cosNormalNext;                                                        //Valeur du cosinus entre l'arete actuelle et l'arete suivante
    double cosNormalNextTmp;                                                     //Valeur du cosinus entre l'arete actuelle et l'arete suivante potentielle
    double psV1Normal;                                                           //Valeur du produit scalaire entre le vecteur 1 et la normale a l'element pere
    double psV2Normal;                                                           //Valeur du produit scalaire entre le vecteur 2 et la normale a l'element pere

    bool alreadyCreated;                                                         //True si l'arete existe deja
    bool isBorderEdge;                                                           //True s'il existe une arete bord connectee a l'arete actuelle
    bool normalTest;                                                             //True si le test sur la normale est vraie
    bool areteConfondue;                                                         //True s'il y a une arete confondue


    /**********   Initialisation   **********/

    _newEltOldElt.resize(6*(nEltTM + nEltSM),-1);
    newEltOldEltQuad.resize(4*(nEltTM + nEltSM),-1);
    newEltOldEltPol.resize(6*(nEltTM + nEltSM),-1);
    tagEdge.resize(2*nEdge,-1);
    stackEdge.resize(40);
    _eltVertConnectivityIM.resize(6*(nEltTM + nEltSM));
    _eltVertIndexIM.resize((nEltTM + nEltSM));
    indElt.resize(2*((nEltTM + nEltSM)));
    quadrangleElt.resize(4*(nEltTM + nEltSM));
    polygonElt.resize(6*(nEltTM + nEltSM));
    polygonEltIndex.resize(nEltTM + nEltSM);
   

    nElts = 0;



    /********  ALGO  Coordonnees    ***********/
    /* - On boucle sur les elements peres du maillage cible
         - On calcule la normale pour chacun de l'element
         - On boucle sur chaque arete de l'element
           - On recupere la decomposition de l'arete
           - On boucle sur chaque arete fille, que l'on rajoute au tableau de stockage.

      ATTENTION : - Au signe des aretes meres et filles.
      ----------  - De verifier de ne pas rajouter l'arete fille dans un sens et dans l'autre
                    (lie a l'algorithme de decoupage)

         - Tant que le tableau de stockage n'est pas vide :
           - On recupere un indice d'arete non nul
           - Si cette arete n'a pas deja ete marque (valeur negative) on cree un nouvel element
             - On enregistre le premier sommet de l'element
             - Tant que le second sommet de l'arete traitee est different de ce premier sommet :
               - On recupere le nombre d'aretes connectees au second sommet
               - On projette l'arete actuelle sur le plan orthogonale a la normale de l'element pere
               - On boucle sur chaque arete suivante potentielle pour verifier s'il y a une arete bord de l'element pere. Si oui :
                 - On projette l'arete bord sur le plan orthogonale a la normale de l'element pere
                 - On calcule le cosinus et la direction (sens direct ou indirect) avec l'arete actuelle
                 - On garde en memoire l'arete bord la plus a gauche 
                   (cosinus le plus petit avec un sinus > 0 ou cosinus le plus grand avec un sinus < 0)
               - On boucle sur chaque arete suivante potentielle pour choisir l'arete la plus a gauche
                 - On ignore l'arete actuelle comme candidate
                 - Afin d'eviter des cas pathologiques, on effectue un test entre la normale a l'element
                   pere et la normale au plan forme par l'arete actuelle et l'arete suivante potentielle.
                   On souhaite que ces deux aretes soient un minimum colineaire.
                 - On garde en memoire l'arete adjacente la plus a gauche 
                   (cosinus le plus petit avec un sinus > 0 ou cosinus le plus grand avec un sinus < 0)
                 - Si l'arete suivante potentielle est a gauche de l'arete de bord (s'il y en a une), on
                   ajoute l'arete potentielle dans le sens contraire (pour ne pas ajouter une seconde fois l'arete elue)
                   a la pile des aretes a tester.
               - On marque l'arete suivante choisie
               - On met a jour les sommets de la nouvelle arete actuelle
            - Selon le nombre de sommets du nouvel element cree on remplit les tableaux concernes.
              - Le tableau de connectivite Element/Sommets
              - Le tableau de liaison NouvelElement/Ancien Element
              - Le tableau de liaison Numerotation Globale/Numerotation Locale 

           - On passe a l'arete suivante dans le tableau de stockage



       - On applique la meme chose sur les elements peres du maillage source en verifiant
         que si l'arete est deja marque pour un element du maillage cible, il faut mettre
         a jour le tableau Nouvel element/Ancien Element correspondant en remplissant
         la partie pour les elements peres source.

       - On remplit definitivement les tableaux suivants :
         - Tableau d'index Element/Sommets du maillage intersecte
         - Tableau de connectivite Element/Sommets du maillage intersecte
         - Tableau Nouvel element/Ancien element (indice pair element pere source
                                                  indice impair element pere cible)

     */
    
    double tolerance_loc = 1e-19*_tolerance;


    for(int iElt = 0; iElt < nEltTM ; iElt ++){

      normeNormaleFace = norm(&_targetMesh.getNormalFace()[3*iElt]);
      normalFaceTM[0] = _targetMesh.getNormalFace()[3*iElt]/normeNormaleFace;
      normalFaceTM[1] = _targetMesh.getNormalFace()[3*iElt + 1]/normeNormaleFace;
      normalFaceTM[2] = _targetMesh.getNormalFace()[3*iElt + 2]/normeNormaleFace;
      
      nEdgeElt = eltVertIndexTM[iElt + 1] - eltVertIndexTM[iElt];
      
      nStackEdge = 0;
      iStackEdge = 0;
      
      for(int iEdge = 0 ; iEdge < nEdgeElt ; iEdge ++){
        
        oldEdge = eltEdgeConnectivityTM[eltVertIndexTM[iElt] + iEdge];
        signe = oldEdge > 0 ? 1 : -1;
        oldEdge = std::abs(oldEdge) - 1;
        
        nNewEdge = oldEdgeNewEdgeIndexIM[oldEdge + 1] - oldEdgeNewEdgeIndexIM[oldEdge];
        
        for(int iCutEdge = 0 ; iCutEdge < nNewEdge; iCutEdge ++){
          newEdge = signe*
            oldEdgeNewEdgeConnectivityIM[oldEdgeNewEdgeIndexIM[oldEdge] + iCutEdge];
          
          borderEdge[std::abs(newEdge) - 1] = true;
          
          if(doubleEdge[std::abs(newEdge) - 1] == 0){
            if(stackEdge.size() <= nStackEdge)
              stackEdge.resize(2*stackEdge.size());
            
            stackEdge[nStackEdge] = newEdge;
            
            doubleEdge[std::abs(newEdge) - 1] 
              = nStackEdge + 1;            
            
            nStackEdge ++;
          }
          else{
            stackEdge[doubleEdge[std::abs(newEdge) - 1] - 1] 
              = 0;
            
            doubleEdge[std::abs(newEdge) - 1] = 0;            
          }
        }
      }
      
      for(int i = 0 ; i < nStackEdge ; i++)        
        if(stackEdge[i] != 0)
          doubleEdge[std::abs(stackEdge[i]) - 1] = 0;
      
      
      /******************************************************************/
      
      
      while(iStackEdge < nStackEdge){
        alreadyCreated = false;
        
        newEdge = stackEdge[iStackEdge];          
        while(newEdge == 0){
          iStackEdge ++;
          newEdge = stackEdge[iStackEdge];          
        }
        
        if(newEdge > 0){
          
          firstVertElt = edgeVertexConnectivityIM[2*(newEdge - 1)] - 1;
          intVert = edgeVertexConnectivityIM[2*(newEdge - 1) + 1] - 1;
          indTag = tagEdge[2*(newEdge - 1)];
          
          if(indTag < 0)
            tagEdge[2*(newEdge - 1)] = nElts;
          else
            alreadyCreated = true;
        }
        else{
          firstVertElt = edgeVertexConnectivityIM[-2*(newEdge + 1) + 1] - 1;
          intVert = edgeVertexConnectivityIM[-2*(newEdge + 1) ] - 1;
          indTag = tagEdge[-2*(newEdge + 1) + 1];
          if(indTag < 0)
            tagEdge[-2*(newEdge + 1) + 1] = nElts;
          else
            alreadyCreated = true;
        }
        
        
        /****************************************************************/
        
        
        
        if(!alreadyCreated){
          
          if(_eltVertIndexIM.size() <= nElts + 1){
            _eltVertIndexIM.resize(2*_eltVertIndexIM.size());
            indElt.resize(2*indElt.size());
            _eltVertConnectivityIM.resize(2*_eltVertConnectivityIM.size());
          }
          
          if(_newEltOldElt.size() <= 2*(nElts + 1))
            _newEltOldElt.resize(4*(nElts + 1));
          
          prevVert = firstVertElt;
          
          nVertNewElt = 1;
          newElt[0] = firstVertElt + 1;
          
          
          nIntersEdge 
            = vertexEdgeIndexIM[intVert + 1] - vertexEdgeIndexIM[intVert];
          
          vect1[0] = _coordsIM[3*intVert] - _coordsIM[3*prevVert];
          vect1[1] = _coordsIM[3*intVert + 1] - _coordsIM[3*prevVert + 1];
          vect1[2] = _coordsIM[3*intVert + 2] - _coordsIM[3*prevVert + 2];
          
          vect1Proj[0] = vect1[0] - psV1Normal * normalFaceTM[0];
          vect1Proj[1] = vect1[1] - psV1Normal * normalFaceTM[1];
          vect1Proj[2] = vect1[2] - psV1Normal * normalFaceTM[2];
          
          areteConfondue = false;

          for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){
            
            nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                   + iEdge];
            
            if(nextEdgeTmp > 0)
              nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
            else
              nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;
            
            if(intVert != nextVertTmp && nextVertTmp != prevVert ){
              
              vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
              vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
              vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];
              
              
              vect2Proj[0] = vect2[0] - psV2Normal * normalFaceTM[0];
              vect2Proj[1] = vect2[1] - psV2Normal * normalFaceTM[1];
              vect2Proj[2] = vect2[2] - psV2Normal * normalFaceTM[2];                
              
              cosV1V2 = cosinus(&vect1Proj[0], &vect2Proj[0], tolerance_loc);
                            
              if(fabs(1 + cosV1V2) < 1e-15 && 
                 norm(&vect2Proj[0]) < norm(&vect1Proj[0]))
           
                areteConfondue = true;
                
              
            }
            
          }
          
          if(!areteConfondue){
            
            while(firstVertElt != intVert){
              
              if (edgeElementTag[std::abs(newEdge) - 1 ]){
                
                bftc_error(__FILE__,
                               __LINE__, 0,
                               "Probleme de reconstruction (arete utilisee plusieurs fois) \n");
                
                return false;
              }
              else
                edgeElementTag[std::abs(newEdge) - 1 ] = true;
              
              isBorderEdge = false;
              
              nIntersEdge 
                = vertexEdgeIndexIM[intVert + 1] - vertexEdgeIndexIM[intVert];
              
              vect1[0] = _coordsIM[3*intVert] - _coordsIM[3*prevVert];
              vect1[1] = _coordsIM[3*intVert + 1] - _coordsIM[3*prevVert + 1];
              vect1[2] = _coordsIM[3*intVert + 2] - _coordsIM[3*prevVert + 2];
              
              psV1Normal = dotProduct(&vect1[0],&normalFaceTM[0]);
              
              vect1Proj[0] = vect1[0] - psV1Normal * normalFaceTM[0];
              vect1Proj[1] = vect1[1] - psV1Normal * normalFaceTM[1];
              vect1Proj[2] = vect1[2] - psV1Normal * normalFaceTM[2];
              
              
              for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){
                
                if(borderEdge[std::abs(vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                           + iEdge]) - 1]){
                  
                  nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                         + iEdge];
                  
                  if(nextEdgeTmp > 0)
                    nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
                  else
                    nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;
                  
                  vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
                  vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
                  vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];
                  
                  psV2Normal = dotProduct(&vect2[0],&normalFaceTM[0]);
                  
                  vect2Proj[0] = vect2[0] - psV2Normal * normalFaceTM[0];
                  vect2Proj[1] = vect2[1] - psV2Normal * normalFaceTM[1];
                  vect2Proj[2] = vect2[2] - psV2Normal * normalFaceTM[2];                
                  
                  cosV1V2 = cosinus(&vect1Proj[0], &vect2Proj[0], tolerance_loc);
                  normalizedCrossProduct(&vect1Proj[0],
                                         &vect2Proj[0],
                                         tolerance_loc,
                                         &vectCrossProduct[0]);                  
                  
                  if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance)
                    direction = cosV1V2 > 0 ? 1 : -1;
                  else
                    direction = signDotProduct(&vectCrossProduct[0],&normalFaceTM[0]);                
                  
                  if(!isBorderEdge){
                    cosBorder = cosV1V2;
                    directionBorder = direction;                
                    
                  }
                  else{
                    if(directionBorder < 0){
                      if(direction > 0){
                        cosBorder = cosV1V2;
                        directionBorder = direction;
                      }
                      else{
                        if (cosV1V2> cosBorder){
                          cosBorder = cosV1V2;
                          directionBorder = direction;
                        }
                      }
                    }
                    else{
                      if (direction > 0 && cosV1V2 < cosBorder){
                        cosBorder = cosV1V2;
                        directionBorder = direction;
                      }
                    }
                  }
                  
                  isBorderEdge = true;
                  
                }
              }
              
              
              if(stackEdge.size() <= nStackEdge + nIntersEdge)
                stackEdge.resize(2*(nStackEdge + nIntersEdge));            
              
              nextEdge = 0;
              
              
              for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){
                
                nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                       + iEdge];
                
                if(nextEdgeTmp > 0)
                  nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
                else
                  nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;
                
                if(intVert != nextVertTmp && nextVertTmp != prevVert){
                  
                  vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
                  vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
                  vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];
                  
                  normalizedCrossProduct(&vect1[0],&vect2[0],tolerance_loc,&vectCrossProduct[0]);
                  
                  normalTest = false;
                  
                  if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance){
                    if(cosinus(&vect1[0],&vect2[0],tolerance_loc) > 0)
                      normalTest = true;
                  }
                  else
                    normalTest = (fabs(cosinus(&vectCrossProduct[0],
                                              &normalFaceTM[0], tolerance_loc)) > (1-_tolerance));
                  
                  if( normalTest ||
                      (fabs(dotProduct(&vect1[0],&normalFaceTM[0])) < 0.001*_tolerance && 
                       fabs(dotProduct(&vect2[0],&normalFaceTM[0])) < 0.001*_tolerance)){
                    
                    psV2Normal = dotProduct(&vect2[0],&normalFaceTM[0]);                                 
                    
                    vect2Proj[0] = vect2[0] - psV2Normal * normalFaceTM[0];
                    vect2Proj[1] = vect2[1] - psV2Normal * normalFaceTM[1];
                    vect2Proj[2] = vect2[2] - psV2Normal * normalFaceTM[2];                
                    
                    cosV1V2 = cosinus(&vect1Proj[0],&vect2Proj[0], tolerance_loc);
                    normalizedCrossProduct(&vect1Proj[0],&vect2Proj[0],tolerance_loc,&vectCrossProduct[0]);
                    
                    if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance)
                      direction = cosV1V2 > 0 ? 1 : -1;
                    else
                      direction = signDotProduct(&vectCrossProduct[0],&normalFaceTM[0]);                
                    
                    cosNormalNextTmp = cosinus(&vect2[0],&normalFaceTM[0],tolerance_loc);
                    
                    if(nextEdge == 0){
                      cosLeft = cosV1V2;
                      directionLeft = direction;
                      nextVert = nextVertTmp;
                      nextEdge = nextEdgeTmp;
                      cosNormalNext = cosNormalNextTmp;
                    }
                    else if( directionLeft - direction == 0 
                             && fabs(cosLeft - cosV1V2) < 1e-15){
                      
                      if(fabs(cosNormalNextTmp) < fabs(cosNormalNext) || 
                         
                         (fabs(fabs(cosNormalNextTmp) - fabs(cosNormalNext)) < 1e-15*_tolerance
                          && norm(&vect2Proj[0]) 
                          < norm(&_coordsIM[3*intVert],&_coordsIM[3*nextVert]))){
                        
                        cosLeft = cosV1V2;
                        directionLeft = direction;
                        nextVert = nextVertTmp;
                        nextEdge = nextEdgeTmp;
                        cosNormalNext = cosNormalNextTmp;
                      }
                    }
                    else{
                      if(directionLeft < 0){
                        if(direction > 0){
                          cosLeft = cosV1V2;
                          directionLeft = direction ;
                          nextVert = nextVertTmp;
                          nextEdge = nextEdgeTmp;
                          cosNormalNext = cosNormalNextTmp;
                        }
                        else{                      
                          if ( cosV1V2 > cosLeft){
                            cosLeft = cosV1V2;
                            directionLeft = direction;
                            nextVert = nextVertTmp;
                            nextEdge = nextEdgeTmp;
                            cosNormalNext = cosNormalNextTmp;
                          }
                        }
                      }                  
                      else{
                        if (direction > 0 && cosV1V2 < cosLeft){
                          cosLeft = cosV1V2;
                          directionLeft = direction ;
                          nextVert = nextVertTmp;
                          nextEdge = nextEdgeTmp;
                          cosNormalNext = cosNormalNextTmp;
                        }
                      }                    
                    }
                    
                    if(isBorderEdge){
                      if(directionBorder < 0){
                        if(direction > 0){
                          stackEdge[nStackEdge] = -(nextEdgeTmp);
                          nStackEdge ++;
                        }
                        else{
                          if (cosV1V2> cosBorder){
                            stackEdge[nStackEdge] = -(nextEdgeTmp);
                            nStackEdge ++;
                          }
                        }
                      }
                      else{
                        if (direction > 0 && cosV1V2 < cosBorder){
                          stackEdge[nStackEdge] = -(nextEdgeTmp);
                          nStackEdge ++;
                        }
                      }
                    }
                    else {
                      stackEdge[nStackEdge] = -(nextEdgeTmp);
                      nStackEdge ++;                           
                    }
                  }
                }
              }
              if(nextEdge == 0){
                bftc_error(__FILE__,
                               __LINE__, 0,
                               "Impossible de trouver une arete a gauche \n Verifiez si les sommets de l'element n° %d du maillage cible sont coplanaires \n",iElt + 1);
                
                return false;
              }
              
              if(edgeMeshTag[std::abs(nextEdge) - 1] == 2 && !(borderEdge[std::abs(nextEdge) - 1])){
                bftc_error(__FILE__,
                               __LINE__, 0,
                               "Reconstruction d'un element avec une arete d'un autre element du meme maillage \n");
              }
              
              if(nextEdge > 0)
                tagEdge[2*(nextEdge - 1)] = nElts;
              else
                tagEdge[-2*(nextEdge + 1) + 1] = nElts;              
              
              if(newElt.size() <= nVertNewElt)
                newElt.resize(2*nVertNewElt);
              
              newElt[nVertNewElt] = intVert + 1;
              
              prevVert = intVert;
              intVert = nextVert;
              newEdge = nextEdge;
              
              nVertNewElt ++;                   
            }
            
            
            
            /****************************************************************/
            
            
            if(nVertNewElt == 3){
              
              if(_eltVertConnectivityIM.size() <= 3*(nTrianElt + 1))
                _eltVertConnectivityIM.resize(2*_eltVertConnectivityIM.size());
              
              for(int iVert = 0; iVert < 3 ; iVert ++)                
                _eltVertConnectivityIM[3*nTrianElt + iVert] = newElt[iVert];
              
              if(_newEltOldElt.size() <= 2*(nTrianElt + 1))
                _newEltOldElt.resize(2*_newEltOldElt.size());
              
              _newEltOldElt[2*nTrianElt] = iElt + 1;
              
              indElt[2*nElts] = 3;
              indElt[2*nElts + 1] = nTrianElt;
              
              nTrianElt ++;
            }
            
            else if(nVertNewElt == 4){
              
              if(quadrangleElt.size() <= 4*(nQuadElt + 1))
                quadrangleElt.resize(2*quadrangleElt.size());
              
              for(int iVert = 0; iVert < 4 ; iVert ++)                
                quadrangleElt[4*nQuadElt + iVert] = newElt[iVert];
              
              if(newEltOldEltQuad.size() <= 2*(nQuadElt + 1))
                newEltOldEltQuad.resize(2*newEltOldEltQuad.size(),-1);
              
              newEltOldEltQuad[2*nQuadElt] = iElt + 1;
              
              indElt[2*nElts] = 4;
              indElt[2*nElts + 1] = nQuadElt;
              
              nQuadElt ++;
            }
            else if (nVertNewElt > 4){
              
              if(polygonElt.size() <= polygonEltIndex[nPolElt] + nVertNewElt)
                polygonElt.resize(2*polygonElt.size());
              
              
              if(polygonEltIndex.size() <= nPolElt + 1)
                polygonEltIndex.resize(2*polygonEltIndex.size());
              
              polygonEltIndex[nPolElt + 1] = polygonEltIndex[nPolElt] + nVertNewElt;
              
              for(int iVert = 0; iVert < nVertNewElt ; iVert ++)
                polygonElt[polygonEltIndex[nPolElt] + iVert] = newElt[iVert];
              
              if(newEltOldEltPol.size() <= 2*(nPolElt + 1))
                newEltOldEltPol.resize(2*newEltOldEltPol.size(),-1);
              
              newEltOldEltPol[2*nPolElt] = iElt + 1;
              
              indElt[2*nElts] = nVertNewElt;
              indElt[2*nElts + 1] = nPolElt;
              
              nPolElt ++;
            }
            else{
              bftc_error(__FILE__,
                             __LINE__, 0,"error element compose de %d \n ",nVertNewElt);
              return false;
            }
            
            edgeElementTag.assign(nEdge,false);

            nElts ++;
          }
        }
        
        iStackEdge ++;
        
      }
    }
    
    int nEltsTM = nElts;
    
    borderEdge.clear();
    borderEdge.resize(nEdge,false);    
    
    printf("Edge TM done \n");
    
    /****************************************************************/
    
    for(int iElt = 0; iElt < nEltSM ; iElt ++){
      
      normeNormaleFace = norm(&_sourceMesh.getNormalFace()[3*iElt]);
      normalFaceSM[0] = _sourceMesh.getNormalFace()[3*iElt]/normeNormaleFace;
      normalFaceSM[1] = _sourceMesh.getNormalFace()[3*iElt + 1]/normeNormaleFace;
      normalFaceSM[2] = _sourceMesh.getNormalFace()[3*iElt + 2]/normeNormaleFace;
      
      nEdgeElt = eltVertIndexSM[iElt + 1] - eltVertIndexSM[iElt];
      
      nStackEdge = 0;
      iStackEdge = 0;
      
      for(int iEdge = 0 ; iEdge < nEdgeElt ; iEdge ++){
        
        oldEdge = eltEdgeConnectivitySM[eltVertIndexSM[iElt] + iEdge];
        
        signe = oldEdge > 0 ? 1 : -1;
        
        oldEdge = std::abs(oldEdge) - 1 + nEdgeTM;
        
        nNewEdge = oldEdgeNewEdgeIndexIM[oldEdge + 1] - oldEdgeNewEdgeIndexIM[oldEdge];
        
        for(int iCutEdge = 0 ; iCutEdge < nNewEdge; iCutEdge ++){
          
          newEdge = signe*
            oldEdgeNewEdgeConnectivityIM[oldEdgeNewEdgeIndexIM[oldEdge] + iCutEdge];
          
          borderEdge[std::abs(newEdge) - 1] = true;
          
          if(doubleEdge[std::abs(newEdge) - 1] == 0){
            if(stackEdge.size() <= nStackEdge)
              stackEdge.resize(2*stackEdge.size());
            stackEdge[nStackEdge] = newEdge;
            
            doubleEdge[std::abs(newEdge) - 1] 
              = nStackEdge + 1;            
            
            nStackEdge ++;
          }
          else{
            stackEdge[doubleEdge[std::abs(newEdge) - 1] - 1] 
              = 0;
            doubleEdge[std::abs(newEdge) - 1] = 0;
            
          }
        }
      }
      
      for(int i = 0 ; i < nStackEdge ; i++)
        if(stackEdge[i] != 0)
          doubleEdge[std::abs(stackEdge[i]) - 1] = 0;
      
      
      /****************************************************************/
      
      
      while(iStackEdge < nStackEdge){
        alreadyCreated = false;
        
        newEdge = stackEdge[iStackEdge];
        
        while(newEdge == 0){
          iStackEdge ++;
          newEdge = stackEdge[iStackEdge];          
        }
        
        if(newEdge > 0){
          firstVertElt = edgeVertexConnectivityIM[2*(newEdge - 1)] - 1;
          intVert = edgeVertexConnectivityIM[2*(newEdge - 1) + 1] - 1;
          indTag = tagEdge[2*(newEdge - 1)];
          
          if(indTag < 0)
            tagEdge[2*(newEdge - 1)] = nElts;
          else {
            if(indTag < nEltsTM){
              if(indElt[2*indTag] == 3 )
                _newEltOldElt[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
              else if(indElt[2*indTag] == 4 )
                newEltOldEltQuad[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
              else if (indElt[2*indTag] > 4)
                newEltOldEltPol[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
            }
            alreadyCreated = true;
            
          }
        }
        else{
          firstVertElt = edgeVertexConnectivityIM[-2*(newEdge + 1) + 1] - 1;
          intVert = edgeVertexConnectivityIM[-2*(newEdge + 1) ] - 1;
          indTag = tagEdge[-2*(newEdge + 1) + 1];
          
          if(indTag < 0)
            tagEdge[-2*(newEdge + 1) + 1] = nElts;
          else{
            if (indTag < nEltsTM ){
              if(indElt[2*indTag] == 3 )
                _newEltOldElt[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
              else if(indElt[2*indTag] == 4 )
                newEltOldEltQuad[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
              else if (indElt[2*indTag] > 4){
                newEltOldEltPol[2*indElt[2*indTag + 1] + 1] = 
                  iElt + 1;
                
              }
            }
            
            alreadyCreated = true;
            
          }
        }
        
        
        /****************************************************************/
        

        if(!alreadyCreated){

          if(_eltVertIndexIM.size() <= nElts + 1){
            _eltVertIndexIM.resize(2*_eltVertIndexIM.size());
            indElt.resize(2*indElt.size(),-1);
            _eltVertConnectivityIM.resize(2*_eltVertConnectivityIM.size());
          }
          
          if(_newEltOldElt.size() <= 2*(nElts + 1))
            _newEltOldElt.resize(4*(nElts + 1));
          
          prevVert = firstVertElt;
          
          nVertNewElt = 1;
          newElt[0] = firstVertElt + 1;
          

          nIntersEdge 
            = vertexEdgeIndexIM[intVert + 1] - vertexEdgeIndexIM[intVert];
          
          vect1[0] = _coordsIM[3*intVert] - _coordsIM[3*prevVert];
          vect1[1] = _coordsIM[3*intVert + 1] - _coordsIM[3*prevVert + 1];
          vect1[2] = _coordsIM[3*intVert + 2] - _coordsIM[3*prevVert + 2];
          
          vect1Proj[0] = vect1[0] - psV1Normal * normalFaceTM[0];
          vect1Proj[1] = vect1[1] - psV1Normal * normalFaceTM[1];
          vect1Proj[2] = vect1[2] - psV1Normal * normalFaceTM[2];
          
          areteConfondue = false;
          
          for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){
            
            nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                   + iEdge];
            
            if(nextEdgeTmp > 0)
              nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
            else
              nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;
            
            if(intVert != nextVertTmp && nextVertTmp != prevVert){
              
              vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
              vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
              vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];
              
              
              vect2Proj[0] = vect2[0] - psV2Normal * normalFaceTM[0];
              vect2Proj[1] = vect2[1] - psV2Normal * normalFaceTM[1];
              vect2Proj[2] = vect2[2] - psV2Normal * normalFaceTM[2];                
              
              cosV1V2 = cosinus(&vect1Proj[0],&vect2Proj[0],tolerance_loc);
              
              
              if(fabs(1 + cosV1V2) < 1e-15 && 
                 norm(&vect2Proj[0]) < norm(&vect1Proj[0])){
                areteConfondue = true;
                
              }
            }
            
          }
          
          if(!areteConfondue){
            
            while(firstVertElt != intVert){

              if (edgeElementTag[std::abs(newEdge) - 1 ]){
                
                bftc_error(__FILE__,
                               __LINE__, 0,
                               "Probleme de reconstruction (arete utilisee plusieurs fois) \n");
                
                return false;
              }
              else
                edgeElementTag[std::abs(newEdge) - 1 ] = true;
              
              isBorderEdge = false;
            
              nIntersEdge 
                = vertexEdgeIndexIM[intVert + 1] - vertexEdgeIndexIM[intVert];
              
              vect1[0] = _coordsIM[3*intVert] - _coordsIM[3*prevVert];
              vect1[1] = _coordsIM[3*intVert + 1] - _coordsIM[3*prevVert + 1];
              vect1[2] = _coordsIM[3*intVert + 2] - _coordsIM[3*prevVert + 2];
              
              psV1Normal = dotProduct(&vect1[0],&normalFaceSM[0]);
              
              vect1Proj[0] = vect1[0] - psV1Normal * normalFaceSM[0];
              vect1Proj[1] = vect1[1] - psV1Normal * normalFaceSM[1];
              vect1Proj[2] = vect1[2] - psV1Normal * normalFaceSM[2];
              
              for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){
                
                if(borderEdge[std::abs(vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                           + iEdge]) - 1]){
                
                  nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                         + iEdge];
                  
                  if(nextEdgeTmp > 0)
                    nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
                  else
                    nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;
                  
                       
                  vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
                  vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
                  vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];
                  
                  psV2Normal = dotProduct(&vect2[0],&normalFaceSM[0]);
                  
                  vect2Proj[0] = vect2[0] - psV2Normal * normalFaceSM[0];
                  vect2Proj[1] = vect2[1] - psV2Normal * normalFaceSM[1];
                  vect2Proj[2] = vect2[2] - psV2Normal * normalFaceSM[2];
                  
                  cosV1V2 = cosinus(&vect1[0],&vect2[0],tolerance_loc);
                  crossProduct(&vect1Proj[0],&vect2Proj[0],&vectCrossProduct[0]);
                  
                  if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance)
                    direction = cosV1V2 > 0 ? 1 : -1;
                  else
                    direction = signDotProduct(&vectCrossProduct[0],&normalFaceSM[0]);                
                  
                  
                  
                  
                if(!isBorderEdge){
                  cosBorder = cosV1V2;
                  directionBorder = direction;                
                }
                else{
                  if(directionBorder < 0){
                    if(direction > 0){
                      cosBorder = cosV1V2;
                      directionBorder = direction;
                    }
                    else{
                      if (cosV1V2> cosBorder){
                        cosBorder = cosV1V2;
                        directionBorder = direction;
                      }
                    }
                  }
                  else{
                    if (direction > 0 && cosV1V2 < cosBorder){
                      cosBorder = cosV1V2;
                      directionBorder = direction;
                    }
                  }
                }
                
                isBorderEdge = true;
                
              }
            }
            
            if(stackEdge.size() <= nStackEdge + nIntersEdge)
              stackEdge.resize(2*(nStackEdge + nIntersEdge));            
           
            
            nextEdge = 0;

            for(int iEdge = 0; iEdge < nIntersEdge ; iEdge ++){

              nextEdgeTmp = vertexEdgeConnectivityIM[vertexEdgeIndexIM[intVert]
                                                 + iEdge];
                            
              if(nextEdgeTmp > 0)
                nextVertTmp = edgeVertexConnectivityIM[2*(nextEdgeTmp - 1) + 1] - 1;
              else
                nextVertTmp = edgeVertexConnectivityIM[-2*(nextEdgeTmp + 1)] - 1;


              if(intVert != nextVertTmp && nextVertTmp != prevVert){
                                         
                vect2[0] = _coordsIM[3*nextVertTmp] - _coordsIM[3*intVert];
                vect2[1] = _coordsIM[3*nextVertTmp + 1] - _coordsIM[3*intVert + 1];
                vect2[2] = _coordsIM[3*nextVertTmp + 2] - _coordsIM[3*intVert + 2];

                normalizedCrossProduct(&vect1[0],&vect2[0],tolerance_loc,&vectCrossProduct[0]);
                                                   
                normalTest = false;

                if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance)
                  normalTest = true;
                else
                  normalTest = (fabs(cosinus(&vectCrossProduct[0],&normalFaceTM[0],tolerance_loc)) 
                                > 0.7*(1-_tolerance));
                
                if( normalTest || 
                    (fabs(dotProduct(&vect1[0],&normalFaceSM[0])) < 0.001*_tolerance && 
                     fabs(dotProduct(&vect2[0],&normalFaceSM[0])) < 0.001*_tolerance) ){
                  

                    psV2Normal = dotProduct(&vect2[0],&normalFaceSM[0]);
                  
                  vect2Proj[0] = vect2[0] - psV2Normal * normalFaceSM[0];
                  vect2Proj[1] = vect2[1] - psV2Normal * normalFaceSM[1];
                  vect2Proj[2] = vect2[2] - psV2Normal * normalFaceSM[2];


                
                  cosV1V2 = cosinus(&vect1Proj[0],&vect2Proj[0], tolerance_loc);

                  normalizedCrossProduct(&vect1Proj[0],&vect2Proj[0],tolerance_loc,&vectCrossProduct[0]);
                  
                  if(norm(&vectCrossProduct[0]) < 1e-15*_tolerance)
                    direction = cosV1V2 > 0 ? 1 : -1;
                  else
                    direction = signDotProduct(&vectCrossProduct[0],&normalFaceSM[0]);                
                   
                  if(nextEdge == 0){
                    cosLeft = cosV1V2;
                    directionLeft = direction;
                    nextVert = nextVertTmp;
                    nextEdge = nextEdgeTmp;
                    cosNormalNext = cosNormalNextTmp;
                  }
                  else if( std::abs(directionLeft - direction) == 0
                           && fabs(cosLeft - cosV1V2) < 1e-15){
                    if(fabs(cosNormalNextTmp) < fabs(cosNormalNext)
                       ||                        
                       (fabs(fabs(cosNormalNextTmp) - fabs(cosNormalNext)) < 1e-15
                        && norm(&vect2Proj[0]) 
                        < norm(&_coordsIM[3*intVert],&_coordsIM[3*nextVert]))){

                      cosLeft = cosV1V2;
                      directionLeft = direction;
                      nextVert = nextVertTmp;
                      nextEdge = nextEdgeTmp;
                      cosNormalNext = cosNormalNextTmp;
                    }
                  }                  
                  else{
                    if(directionLeft < 0){
                      if(direction > 0){
                        cosLeft = cosV1V2;
                        directionLeft = direction ;
                        nextVert = nextVertTmp;
                        nextEdge = nextEdgeTmp;
                        cosNormalNext = cosNormalNextTmp;
                      }
                      else{                      
                        if ( cosV1V2 > cosLeft){
                          cosLeft = cosV1V2;
                          directionLeft = direction;
                          nextVert = nextVertTmp;
                          nextEdge = nextEdgeTmp;
                          cosNormalNext = cosNormalNextTmp;
                        }
                      }
                    }
                    else{
                      if (direction > 0 && cosV1V2 < cosLeft){
                        cosLeft = cosV1V2; 
                        directionLeft = direction ;
                        nextVert = nextVertTmp;
                        nextEdge = nextEdgeTmp;
                        cosNormalNext = cosNormalNextTmp;
                      }
                    }
                  }
                                
                if(isBorderEdge){
                  if(directionBorder < 0){
                    if(direction > 0){
                      
                      stackEdge[nStackEdge] = -(nextEdgeTmp);                        
                      
                      nStackEdge ++;
                    }
                    else{
                      if (cosV1V2 > cosBorder){
                        stackEdge[nStackEdge] = -(nextEdgeTmp);
                        
                        nStackEdge ++;
                      }
                    }
                  }
                  else{
                    if (direction > 0 && cosV1V2 < cosBorder){
                      stackEdge[nStackEdge] = -(nextEdgeTmp);
                      
                      nStackEdge ++;
                    }
                  }
                }
                else {
                  stackEdge[nStackEdge] = -(nextEdgeTmp);
                  
                  nStackEdge ++;                
                }
                }  
              }           
            }

            if(nextEdge == 0){
              
              return false;
            }
            
            if(edgeMeshTag[std::abs(nextEdge) - 1] == 1 && !(borderEdge[std::abs(nextEdge) - 1])){
              bftc_error(__FILE__,
                             __LINE__, 0,
                             "Reconstruction d'un element avec une arete d'un autre element du meme maillage \n");
            }
            
            if(nextEdge > 0)
              tagEdge[2*(nextEdge - 1)] = nElts;
            else
              tagEdge[-2*(nextEdge + 1) + 1] = nElts;              
            
            if(newElt.size() <= nVertNewElt)
              newElt.resize(2*nVertNewElt);
            
            newElt[nVertNewElt] = intVert + 1;

            prevVert = intVert;
            intVert = nextVert;
            newEdge = nextEdge;

            nVertNewElt ++;

                   
          }
          

          /****************************************************************/
          

          if(nVertNewElt == 3){
            
            if(_eltVertConnectivityIM.size() <= 3*(nTrianElt + 1))
              _eltVertConnectivityIM.resize(2*_eltVertConnectivityIM.size());
            
            for(int iVert = 0; iVert < 3 ; iVert ++)                
              _eltVertConnectivityIM[3*nTrianElt + iVert] = newElt[iVert];


            if(_newEltOldElt.size() <= 2*(nTrianElt + 1))
              _newEltOldElt.resize(2*_newEltOldElt.size());
            
            _newEltOldElt[2*nTrianElt + 1] = iElt + 1;

            nTrianElt ++;
          }
          
          else if(nVertNewElt == 4){
            
            if(quadrangleElt.size() <= 4*(nQuadElt + 1))
              quadrangleElt.resize(2*quadrangleElt.size());
            
            for(int iVert = 0; iVert < 4 ; iVert ++)                
              quadrangleElt[4*nQuadElt + iVert] = newElt[iVert];

            if(newEltOldEltQuad.size() <= 2*(nQuadElt + 1))
              newEltOldEltQuad.resize(2*newEltOldEltQuad.size(),-1);
            
            newEltOldEltQuad[2*nQuadElt + 1] = iElt + 1;


            nQuadElt ++;
          }
          else if (nVertNewElt > 4){

            if(polygonElt.size() <= polygonEltIndex[nPolElt] + nVertNewElt)
              polygonElt.resize(2*polygonElt.size());
            
            
            if(polygonEltIndex.size() <= nPolElt + 1)
              polygonEltIndex.resize(2*polygonEltIndex.size());
            
            polygonEltIndex[nPolElt + 1] = polygonEltIndex[nPolElt] + nVertNewElt;
            
            for(int iVert = 0; iVert < nVertNewElt ; iVert ++)                
              polygonElt[polygonEltIndex[nPolElt] + iVert] = newElt[iVert];


            if(newEltOldEltPol.size() <= 2*(nPolElt + 1))
              newEltOldEltPol.resize(2*newEltOldEltPol.size(),-1);
            
            newEltOldEltPol[2*nPolElt + 1] = iElt + 1;

            nPolElt ++;
          }
          else{
             bftc_error(__FILE__,
                            __LINE__, 0,"error element compose de %d \n ",nVertNewElt);

             return false;
          }

          edgeElementTag.assign(nEdge,false);

          nElts ++;
          }
        }

        iStackEdge ++;
      }
    }


    printf("Edge SM done \n");

    _eltVertIndexIM.resize(nElts + 1);
    _eltVertConnectivityIM.resize(3*nTrianElt + 4*nQuadElt + polygonEltIndex[nPolElt]);
    _eltVertIndexIM[0] = 0;
    _newEltOldElt.resize(2*nElts);

    /****************************************************************/


    for(int iTri = 0 ; iTri < nTrianElt ; iTri++)
      _eltVertIndexIM[iTri + 1] = _eltVertIndexIM[iTri] + 3;

    for(int iQuad = 0 ; iQuad < nQuadElt ; iQuad++){
      _eltVertIndexIM[iQuad + nTrianElt + 1] = _eltVertIndexIM[iQuad + nTrianElt] + 4;

      _newEltOldElt[2*(nTrianElt + iQuad)] = newEltOldEltQuad[2*iQuad];
      _newEltOldElt[2*(nTrianElt + iQuad) + 1] = newEltOldEltQuad[2*iQuad + 1];
      }

    for(int iQuad = 0 ; iQuad < nQuadElt ; iQuad++)    
      for(int iVert = 0 ; iVert < 4 ; iVert++)
        _eltVertConnectivityIM[_eltVertIndexIM[iQuad + nTrianElt] + iVert]
          = quadrangleElt[4*iQuad + iVert];
    
    for(int iPol = 0 ; iPol < nPolElt ; iPol++){
      _eltVertIndexIM[iPol + nTrianElt + nQuadElt + 1] 
        = _eltVertIndexIM[nTrianElt + nQuadElt] + polygonEltIndex[iPol + 1];

      _newEltOldElt[2*(nTrianElt + nQuadElt + iPol)] = newEltOldEltPol[2*iPol];
      _newEltOldElt[2*(nTrianElt + nQuadElt + iPol) + 1] = newEltOldEltPol[2*iPol + 1];
    }
    
    for(int iPol = 0 ; iPol < nPolElt ; iPol++)
      for(int iVert = 0 ;
          iVert <  polygonEltIndex[iPol + 1] - polygonEltIndex[iPol] ;
          iVert++){       
        /* printf("size  %d  ind  %d \n",(int) polygonElt.size(), 
           polygonEltIndex[iPol] + iVert);*/

        _eltVertConnectivityIM[_eltVertIndexIM[iPol + nTrianElt + nQuadElt] + iVert]
          = polygonElt[polygonEltIndex[iPol] + iVert];
      }

    /*******************************************************************/
    /*    bftc_printf("triangles %d  quadrangles %d polygones %d  \n",
                    nTrianElt,nQuadElt,nPolElt);


    bftc_printf("_eltVertConnectivityIM  \n");
    for(int i = 0 ; i < nElts ; i++){
      bftc_printf("element n° %d  ",i + 1);
      for(int j = 0; j < _eltVertIndexIM[i + 1] - _eltVertIndexIM[i] ; j++)
        bftc_printf("%d ",_eltVertConnectivityIM[_eltVertIndexIM[i] + j]);
      bftc_printf("\n");
    }
    bftc_printf("--------end------- \n");*/
    /******************************************************************  /

    /**************   Libération des donnees  ************************/

    tagEdge.clear();
    stackEdge.clear();
    vect1.clear();
    vect2.clear();
    newElt.clear();
    borderEdge.clear();
    polygonElt.clear();
    quadrangleElt.clear();
    indElt.clear();
    doubleEdge.clear();
    edgeElementTag.clear();

    return true;
  }
}

