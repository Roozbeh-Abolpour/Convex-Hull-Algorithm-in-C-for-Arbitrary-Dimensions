# Convex Hull Algorithm in C for Arbitrary Dimensions
This project implements a safe and efficient convex hull construction algorithm in pure C, based on the Quickhull method. It supports point clouds in arbitrary dimensions.

## Features
- Supports n-dimensional point sets (2D, 3D, and higher).
- Based on the Quickhull algorithm (Barber et al., 1996).
- Custom lightweight matrix library (MatLib).
- Incremental convex hull construction.
- Memory-safe and well-structured C code.

## Examples

### 2D Convex Hull
Input points:
[2D Input](examples/Fig (1).jpg)
Computed convex hull:
![2D Hull](examples/Fig (2).jpg)

### 3D Convex Hull
Input points:
![3D Input](examples/Fig (3).jpg)
Computed convex hull:
![3D Hull](examples/Fig (4).jpg)

### Code Example: Compute Convex Hull from a File

The file (examples/example_from_file.c) shows how to:
- Read a matrix of points from file1.txt
- Compute the convex hull
- Write the facets and neighbor indices to output files
- 
Inside the examples/ :
-file1.txt: Contains the input points matrix.  
  Example (2D points).
  
When you run the example program (example_from_file.c), it generates:
-file2.txt: Convex hull facets (each facet as a set of point indices)
-file3.txt: Neighbor relationships between facets
-file4.txt: Outpoints information (points lying outside each facet)

