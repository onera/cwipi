---
title: 'CWIPI: Coupling With Interpolation Parallel Interface'
tags:
  - parallel computing
  - coupling
  - multi-physics
  - MPI

authors:
  - name: Eric Quémerais
    orcid: 0009-0007-2018-6358
    affiliation: "1"

  - name: Bastien Andrieu
    orcid: 0009-0000-9937-0244
    affiliation: "1"

  - name: Karmijn Hoogveld
    orcid: 0009-0007-8251-1152
    affiliation: "1"

affiliations:
 - index: 1
   name: DMPE, ONERA, Université Paris-Saclay, 92320 Châtillon, France

date: 27 January 2026

bibliography: paper.bib

---

# Summary

CWIPI is an open-source library dedicated to the parallel coupling of independent simulation codes for multi-physics research, a field also covered by tools such as MUI [@Tang2015] and preCICE [@Chourdakis2022]. It is aimed at anyone wishing to conduct studies and research in this field. It provides MPI-based communication and fully distributed data exchange on non-coincident meshes, including on-the-fly spatial interpolation. CWIPI stands out for its support of a wide variety of spatial interpolation families and a broad range of geometric element types. This extensive spatial capability allows users to leverage solver-specific numerical methods, including finite volumes, finite elements, discontinuous Galerkin, and spectral methods. Advanced users can further combine CWIPI's geometric algorithms with solvers' internal interpolations through callbacks, improving accuracy across highly heterogeneous solvers.

CWIPI supports a wide range of interface geometries, from linear, surface, and volumetric elements to higher-order elements with curved edges and faces. Performance-critical algorithms have been specifically optimized for handling coupling cases involving dynamic or adaptive meshes. Recent developments in the ParaDiGM library [@Quemerais2026] now provide dynamic load and memory balancing for large-scale simulations. Ongoing work includes GPU acceleration of geometric algorithms [@Cazalbou2024] and the development of temporal interpolation methods [@Simon2025].

# Statement of Need

Multi-physics problems can be approached either through monolithic methods — where all physical phenomena are solved simultaneously within a single solver — or through partitioned methods, where independent solvers exchange data across shared interfaces. Partitioned approaches are particularly attractive in high-performance computing (HPC) contexts, as they allow the reuse of existing specialized codes while maintaining their native parallel architecture.

Robust and efficient coupling tools are essential to enable data exchange between solvers operating on distinct meshes, geometries, and discretization schemes. CWIPI was initially developed in 2009 to provide an open-source solution complementary to MPCCI [@Joppich2006], whose black-box nature limited flexibility and transparency for multi-physics research.

To ensure a straightforward and efficient deployment of coupled simulations, an HPC coupling tool must abstract the complexity of parallelism by interacting solely with data local to each MPI rank. Furthermore, it must deliver high-performance spatial interpolation algorithms whose execution time remains negligible compared to a single solver iteration, thereby avoiding any penalty on the overall simulation wall-clock time. Alongside these, it must provide temporal coupling schemes—or at least the essential building blocks to implement them. While meeting these requirements is sufficient for standard use cases, the ability to implement specific spatial interpolation methods is critical to achieve a global accuracy order close to the solvers' internal precision. This necessity was demonstrated by [@Leger2012] for an isophysical configuration based on the linearized Euler equations: by exposing CWIPI's geometric algorithms to the solvers' internal interpolation capabilities via callbacks, and by employing the same Runge-Kutta scheme with coupling performed at the sub-iteration level, a global accuracy order of 3.5 was preserved when coupling a 4th-order finite difference code with a 4th-order discontinuous Galerkin solver.

# State of the Field

Existing coupling frameworks differ significantly in scope, geometry support, and interpolation capabilities. MUI [@Tang2015] primarily supports point clouds and focuses on general-purpose temporal and spatial interpolation. preCICE [@Chourdakis2022] supports point clouds or meshes composed of planar triangular or quadrilateral elements, provides a broad ecosystem of adapters for established simulation codes, and integrates black-box temporal coupling algorithms. In contrast to these tools, OpenPALM [@Duchaine2015] features a dedicated supervision process that orchestrates temporal coupling algorithms, which are defined beforehand by the user through the PrePALM graphical user interface. In this architecture, OpenPALM operates as a higher-level coupler, relying entirely on CWIPI as its underlying spatial interpolation and low-level communication engine. Although a comparative study of coupling tools within an HPC framework [@Rubin2022] has been conducted, it is not exhaustive; it mentions other tools besides those cited above but only discusses CWIPI through its utilization via OpenPALM.

Regarding temporal coupling algorithms, and unlike black-box solutions such as preCICE, CWIPI does not offer any turn-key solutions; instead, it solely provides the building blocks to define the orchestration of exchanges. Consequently, the temporal coupling algorithm is defined by the user in a decentralized manner.

Regarding spatial interpolation of fields, a concept designated as "mapping" or "sampling" depending on the framework, CWIPI, due to its HPC nature, has oriented itself exclusively towards spatial interpolation methods based on local neighboring data. This contrasts notably with preCICE which, in addition to local interpolation methods, early on proposed a set of global RBF-type interpolation methods, involving all nodes of the interface mesh.

CWIPI offers greater versatility than the solutions mentioned in terms of geometric representation of interfaces. While preCICE's "mapping" or MUI's "sampling" rely mainly on point cloud-based interpolation, CWIPI supports a broad spectrum of interface representations. This includes linear, surface and volumetric elements, polygons and polyhedra, as well as higher-order elements with curved edges and faces. To exploit these geometries, CWIPI provides three families of interpolation methods based on an underlying geometric algorithm: one based on point clouds (similar to preCICE and MUI), one based on locating target points within the source mesh, and one based on the intersection of source and target meshes. For each of these, a default method is provided, supplemented by a callback interface allowing users to implement the interpolation method best suited to their specific context. If MUI also allows users to write a custom callback, this is limited to point clouds; conversely, CWIPI extends this capability to the three geometric families. Within a CWIPI callback, parallelism management is entirely transparent: the user works exclusively with local data, without having to manage remote communications, except in the specific case where they invoke internal functions of their own code requiring their intracommunicator. It should be noted that the family based on locating target points was the first to be implemented, which explains its privileged use in the vast majority of CWIPI applications.

Strategies for launching and communication differ significantly across frameworks. MUI and CWIPI favor an approach based on MPI MPMD (Multiple Program, Multiple Data), whereas preCICE opts for a decoupled approach where each code is launched via its own independent MPI command (MPI Port), with the connection between processes established a posteriori. Although CWIPI offers a TCP/IP communication mode similar to preCICE, it is rarely used, as the requirements of HPC users almost exclusively necessitate the use of MPI.

CWIPI stands out for its high flexibility in distributing codes across MPI ranks. Unlike approaches that mandate disjoint communicators, CWIPI allows a code to be placed on any rank. Consequently, the intersection of the MPI intra-communicators of different codes is not necessarily empty. In the limiting case where codes share all their ranks, their intra-communicators become duplications of the same communicator. This architecture offers total freedom: the same coupling case can be handled with disjoint intra-communicators or fully overlapping ones, without any assumptions regarding the location of the coupling interface. In practice, although this capability is available in C and Fortran, it will mainly be used in Python, since the flexibility of this language makes it easy to load the module for a given coupled code depending on the MPI rank number, and thus to distribute the codes as desired. This capability therefore makes it possible to adjust this distribution in order to improve load balancing across the whole coupled simulation. CWIPI does not provide a capability for searching for an optimal distribution, but this could constitute a future research direction.

# Software Design

Since version 1.0, CWIPI is built on top of the ParaDiGM library [@Quemerais2026], originally created specifically to meet the needs of CWIPI. ParaDiGM provides the low-level geometric algorithms and parallel data structures. More broadly, ParaDiGM was designed to mutualize parallel geometric algorithms that share common requirements across other contexts, such as dynamic mesh adaptation and, more generally, any geometric co-processing tasks performed by solvers during runtime.

CWIPI provides three APIs — C/C++, Fortran, and Python/NumPy — making it compatible with a wide range of scientific codes. The integration of CWIPI into an existing solver is designed to be minimally intrusive: it does not require modifying the code structure, and the API relies solely on simple arrays, keeping the coupling layer lightweight and straightforward to adopt.

To support various discretization schemes, CWIPI's three interpolation families provide distinct default approaches and data structures for custom implementations. The first family, based on searching for the k-nearest source nodes for each target point, operates on nodes, cell centers, or user-defined points for both source and target locations. Its default method relies on a least squares approach, and its callback interface exposes, for each target point, the local data containing the k-nearest source points, their distances, and the associated source field values. The second family focuses on the location of target points within the source mesh. It maps source fields defined at nodes or cell centers onto target fields located at nodes, cell centers, or user-defined points. By default, it uses the Lagrange basis functions of the source elements; for custom interpolations, the callback interface provides the data local to each source element, specifically the target points it contains (or their projections) and the source element fields. The third family is based on the intersection of source and target meshes, where both fields are strictly located at cell centers. Its default method performs a weighting by the volumes of the intersecting polyhedra, and the callback interface supplies each target element with its exact intersection volumes and the fields defined on the intersecting source polyhedra.

At each coupling step, the execution sequence typically proceeds as follows: first, the geometric algorithm corresponding to the selected interpolation family is invoked; second, the physical fields are exchanged across the communication graph established by the preceding geometric step; and finally, the fields are interpolated onto the target mesh. However, if the chosen family is based on the location of target points within the source mesh, the field interpolation is performed directly by the source code, which inverts the order of the last two steps. Furthermore, in cases where the mesh interface remains fixed throughout the simulation, the geometric algorithm is invoked only once during the initialization phase, and its results are subsequently stored in memory for all remaining coupling steps.

Within this coupling communicator, the geometric algorithms linked to the interpolation families generate optimized point-to-point communication graphs between the processes of the coupled codes. Each geometric algorithm produces its own communication graph. This graph also depends on the location of the fields' degrees of freedom (nodes, cell centers, etc.). Since the location of the degrees of freedom and the choice of interpolation family are properties specific to each field, CWIPI is capable of managing several distinct communication graphs between the codes for a single coupling instance.

A major design challenge lies in the **performance of geometric algorithms**, especially when the coupling interface moves during the simulation or when solvers employ dynamic mesh adaptation. In the 0.x branch of CWIPI, algorithms relied on the FVM library from Code_Saturne [@Archambeau2004], which later evolved into PLE, a lightweight coupling library [@Fournier2020]. New algorithms [@Andrieu2026] developed within ParaDiGM provide **dynamic load and memory balancing** at each coupling step, which is critical because coupling problems are inherently unbalanced across MPI processes. These developments bring significant performance gains for large-scale and complex simulations.


## Future Directions

Two main research directions are currently pursued. The first concerns the **GPU porting** of the most computationally expensive geometric algorithms, for which new algorithms have recently been designed [@Cazalbou2024] and will be integrated in an upcoming release. The second focuses on the development of **temporal interpolation methods** [@Simon2025; @Francois2023], complementary to those provided by existing frameworks such as preCICE, aiming to improve accuracy and flexibility for dynamic interfaces and heterogeneous solvers.

# Research Impact Statement

From a scientific perspective, CWIPI is primarily used within the aeronautics and aerospace community from which it originated. The library is integrated into several codes developed at ONERA, including elsA [@Cambier2013] for computational fluid dynamics, CEDRE [@Refloch2011] for multiphysics simulations, Z-set [@Garaud2019] for structural and materials mechanics, MoDeTheC [@Dellinger2024] for thermal degradation of materials, and SPACE [@Delorme2005] for acoustic propagation. Through these tools, CWIPI enables advanced multiphysics simulations in aerothermodynamics [@Errera2019], aeroacoustics [@Langenais2019; @MoratillaVega2022], magnetohydrodynamics [@Rocamora2025], electric arc propagation, and fire safety [@Dellinger2023].

CWIPI is also employed as both an internal and external coupling library in leading combustion simulation environments such as YALES2 [@Moureau2011] and AVBP [@Schonfeld1999]. A notable HPC application [@Dombard2018] combining CWIPI with AVBP on a complex geometry demonstrated the robustness and scalability of the coupling approach on massively parallel architectures. Its capabilities were highlighted as early as 2017 [@Duchaine2017] and confirmed by an invitation to the ExCALIBUR workshop [@Andrieu2021]. Recent work in aeroelasticity with YALES2 [@Fabbri2023] further demonstrates its versatility. CWIPI's integration into the software stack of national HPC facilities such as the CCRT confirms its recognition at a national level.

Beyond its original ecosystem, the open-source community has adopted CWIPI through GitHub-based interfacing projects with OpenFOAM [@Weller1998; @MoratillaVega2022_2; @nd922025] and Nektar++ [@Cantwell2015], enabling applications in aeroacoustics and data assimilation [@Valero2026]. Finally, CWIPI serves as the interpolation engine of the OpenPALM coupler [@Duchaine2015], illustrating its central role in the broader multiphysics and HPC coupling ecosystem.

# AI Usage Disclosure

Generative AI tools were used to assist with the writing of this paper: translation, formatting, and improving the fluency of the text. No AI was used in the development of the CWIPI software. All scientific content was written and validated by the authors.

# Acknowledgements

The authors would like to thank the following people for their contributions to the development, testing, and dissemination of the CWIPI library within the community:

Romain Casta¹, Nicolas Dellinger², Florent Duchaine¹, Bastien Frisulli², Jean-Didier Garaud², Thomas Hennion², Stéphanie Lala², Xavier Lamboley¹², Bruno Maugars², Thierry Morel¹, Christophe Peyret²

¹ CERFACS — ² ONERA

The authors also wish to thank BPI France, the Directorate General for Civil Aviation (DGAC), and the General Scientific Directorate of ONERA for their financial support.

# Author Contributions

The contributions to this software are listed according to the CRediT taxonomy:

- **Eric Quémerais**: Conceptualization, Methodology, Software, Validation, Writing – Original Review & Editing, Project Administration, Funding Acquisition, Supervision
- **Bastien Andrieu**: Software, Validation, Writing – Original Draft
- **Karmijn Hoogveld**: Software, Validation, Writing – Original Draft

# References
