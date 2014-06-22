Geometry Radio-Frequency Library by Chris Tralie
===========

This library started out as a geometry library to support some work with computational electromagnetics for radar, but it has turned into a much more general purpose 3D geometry library.  Features include

Features
--------------
* Support for 3D primitives: Vectors, Points, Rays, Planes, 
* Support for 3D primitive transformations
* Support for 3D polygon meshes, including geometry methods (PCA, plane slicing) and some topology methods (triangle subdivision, basic no-frills filling).  Can load and save .off or .obj files with color
* Basic 3D mesh viewer with a polar camera using PyOpenGL (meshView.py)

Algorithms Implemented
--------------
* Iterative closest points
* Laplacian Mesh Editing
* Image sources calculation for an arbitrary polygon mesh


Algorithms in Development (sort of working but buggy...hopefully items will slowly mature)
--------------
* Polygon beam tracing
* 3D Planar Reflective Symmetry Transform
* Fast marching for geodesic distances
* Generalized Multidimensional Scaling
