Geometry Radio-Frequency Library by Chris Tralie
===========

This library started out as a geometry library to support some work with computational electromagnetics for radar, but it has turned into a much more general purpose 3D geometry library.  Please visit the wiki for some screenshots and videos of this software in action

Features
--------------
* Support for 3D primitives and primitive transformations: Vectors, Points, Rays, Planes, etc
* Support for 3D polygon meshes, including geometry methods (PCA, slice by plane) and some topology methods (triangle subdivision, basic no-frills hole filling).  Can load and save .off or .obj files with color
* Basic 3D mesh viewer with a polar camera using PyOpenGL (meshView.py)

Algorithms Implemented
--------------
* Iterative closest points
* Laplacian Mesh Editing
* Image sources calculation for an arbitrary polygon mesh


Algorithms in Development
--------------
This stuff is sort of working but buggy...hopefully items will slowly mature

* Polygon beam tracing (numerical precision problems)
* 3D Planar Reflective Symmetry Transform (sampling problems)
* Fast marching for geodesic distances (need to handle obtuse triangles)
* Generalized Multidimensional Scaling (boundary case problems cause it to get stuck)
