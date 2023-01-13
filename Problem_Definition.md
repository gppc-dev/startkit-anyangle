# Problem statement: any-angle setting
Your task is to find a path on a given graph, and the solution is evaluated based on optimality, time performance and space cost.

## From Grid to Euclidean Plane
Input grid format is specified at [MovingAI](https://movingai.com/benchmarks/formats.html).

The grid is mapped to the Euclidean plane, details below:

  1. Euclidean plane `(x,y)` coordinates increases to the right (east) and downards (south) respectively.

  2. Each grid cell maps to the Euclidean plane with their `(x,y)` coordinate coorisponding to the Euclidean plane.

  3. The grid cell `(x,y)` translates to Euclidean plane square `(x,y)` to `(x+1,y+1)`, i.e. a grid coordinate is the top-left corner.

  4. The blocked grid cell squares translates to non-traversable space on the Euclidean plane.

### Example Map
Grid:

    ...#.
    #.#..
    #..#.
    ##...

Euclidean palne:

<p align="center">
<img src="./figs/grid_plane.png" height="200">
</p>

## Agent

1. agent at `(x, y)` is a point on the Euclidean plane

2. agent start and target positions are the top-left corner of non-blocked grid cell squares. For example, all blue positions in the following plane are valid start and target positions, and all red positions are invalid.
  <p align="center">
    <img src="figs/euclidean_start_target.png" height="200" width="200"> 
  </p>

## Path
For a query `(s, t)`, a valid path is a sequence of **Integer Coordinates** `p=(s,v1,...vn,t)`, any two adjacent coordinates `a` and `b` on the path must be a **valid path segment** (see details in the next section).


**When start and target are the same position, the path must be empty, the length must be `0`.**

## Valid Path Segments

For any given two adjacent coordinates `a` and `b` in a path, where the agent moves from position `a` to position `b`, the following constraints must be followed in general:
  
The straight line path segment defined by `a` and `b`, must not traverse through:
1. Any non-traversable area, excluding the boundaries and corners of the area.

2. Any boundary of a non-traversable area that is also the boundary of an adjacent non-traversable area.

3. Any corner of a non-traversable area that is also the corner of a diagonal adjacent non-traversable area (No Double Corner Cutting).

In the following example, path segment `(a, b)` is invalid if `b` is located on red circles, and is valid if `b` is located on blue cicles:
  <p align="center">
    <img src="figs/invalid_segments.png" height="200" width="200">
  </p>


### When `a` is a start position on a double conner

A special case is that `a` is a start position located on a double corner. In this case, the outgoing direction must between `EAST` and  `SOUTH`.

In the following example, path segment `(a, b)` is invalid if `b` is located on red circles, and is valid if `b` is located on blue positions:
  <p align="center">
    <img src="figs/invalid_segments_start.png" height="200" width="200">
  </p>

### When `b` is a target position on a double conner

Another special case is that `b` is a target position located on a double corner. In this case, the incoming direction must come between `EAST` and  `SOUTH`.

In the following example, path segment `(a, b)` is invalid if `a` is located on red circles, and is valid if `a` is located on blue positions:
  <p align="center">
    <img src="figs/invalid_segments_target.png" height="200" width="200">
  </p>


