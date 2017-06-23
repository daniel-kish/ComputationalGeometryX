# ComputationalGeometryX

## About
This project is a part of my master's thesis (defense took place in June, 2017). The thesis is dedicated to a research of the so-called Delaunay refinement algorithms.

## Rationale
The Delaunay refinement meshing techniques are distinguished among all other approaches by its mathematical strictness and great amount of proofs and formal analysis. Resulting meshes can be used in **Finite Element** or **Finite Volume Analysis**. I have been studying in CAE field for six years now and the subject of computational geometry and meshing for that matter was always a mystery to me. Mostly this is due to unawareness of even simple meshing techniques among our tutors and professors. I have chosen such a subject simply to fully educate myself.

## Contents
Here I've implemented some of the most successful 2D Delaunay refinement algorithms: Ruppert's alg., the second Chew's alg. and a variant of the latter with Alper Ungor's offcenters.
Of course no Delaunay meshing approach can be implemented without a Delaunay triangulation algorithms. Here I've implemented certain foundational things like:
* Quad-Edge data structure
* Divide-and-conquer Delaunay triangulation algo
* Incremental Delaunay triangulation
* Constrained Delaunay triangulation
* Some other *computational geometry* trifles
