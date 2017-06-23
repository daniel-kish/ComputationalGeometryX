# ComputationalGeometryX

## About
This project is a part of my master's thesis (defense took place in June, 2017). The thesis is dedicated to a research of the so-called Delaunay refinement algorithms.

## Rationale
The Delaunay refinement meshing techniques are distinguished among all other approaches by its mathematical strictness and great amount of proofs and formal analysis. Resulting meshes can be used in **Finite Element** or **Finite Volume Analysis**. 
I have been studying in CAD/CAE field for six years now and the subject of computational geometry (and meshing for that matter) was always a mystery to me. Mostly this is due to unawareness of even simple meshing techniques among our tutors and professors. I have chosen such a subject simply to fully educate myself.

## Contents
Here I've implemented some of the most successful 2D Delaunay refinement algorithms: **Ruppert's alg.**, the **second Chew's alg.** and a variant of the latter with **Alper Ungor's offcenters**.
Of course no Delaunay meshing approach can be implemented without a Delaunay triangulation algorithms. Here I've implemented certain foundational things like:
* Quad-Edge data structure
* Divide-and-conquer Delaunay triangulation algo
* Incremental Delaunay triangulation
* Constrained Delaunay triangulation
* Some other computational geometry trifles

Aforementioned Delaunay meshing algorithms were programmed on the top of that CG algos and data structures. The main feature of interest in the programmatic part of the project is probably Jonathan Shewchuk's *robust adaptive arithmetic* which has helped hugely in achieving numerical stability of the code.

## Results
The most interesting meshes are presented below for your visual amusement :^)

Here's what we start with - a Delaunay triangulation of the points constrained with a number of edges:

![init](/1del_init.png)

And here comes the refined version of it:

![result](/2del_res.png)

Here's how A. Ungors technique helps with mesh size optimization:

![ruppert](/3ruppert_5.png)

![alper](/4alper_5.png)

