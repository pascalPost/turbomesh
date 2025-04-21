# Introduction

## Mathematical Background

## Single Block Smoothing

## Multi Block Smoothing

Different possibilities to construct the system matrix.

- All mesh DOF in the matrix.
- Only Solved DOF in the matrix -> more complexity in essembly

# Details

# System matrix

We need to generate a system matrix.

This is a sparse matrix. There are different sparse matrix solvers and formats.

We use umfpack from SuiteSparse.

// UMFPACK docs: https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK/Doc

The packages uses a Compressed Column format, but also allows to create the matrix in a Compressed Row format.
A row based asembly seems to be more natural for us.

Loop over all blocks and add all block internal points.

For each point, we want to add data for each neighbor that is also in the matrix (to make the matrix as implicit as possible -> less iterations needed).

After each block internal point, the block boundary points are processed.


## Derive a standard stencil

## Show what happens with the stencil on the boundary

## Connections

TODO: discussion about different options

- a single DOF (or matrix entry) for each connection point
- a DOF for each point in the mesh (that means 2 or more entries for each point pair connected by connections)
- Ghost Cells with additional DOF

## Periodic connections

Number of internal points per block (equal to the rows in the system matrix)

| block number | block name | rows in the matrix |
|--|--|--|
| 0 | ss | 0..13568 |
| 1 | ps | 13568..18656 |
| 2 | in | ..24584 |
| 3 | out | ..29164 |
| 4 | down | 29164..32831 |
| 5 | up |  ..37077 |
| 6 | upstream | 37077..40060 |
| 7 | downstream | 40060..45228 |

We look at the first connection

```zig
try mesh.connections.append(.{ .data = .{
    .{ .block = down_id, .side = boundary.Side.j_min, .start = self.num_cells.down_j, .end = 0 },
    .{ .block = upstream_id, .side = boundary.Side.j_max, .start = 0, .end = self.num_cells.down_j },
} });
```

```
              |
              x
              |
              |
              |
              |
upstream      |             down -> i
              |                  |
              |                  V
              |                  j
              x
```

|connection point | point index block 0 (down) | point index block 1 (upstream) | matrix index |
|--|--|--|--|
| 0 | (0, 20) = [20] = {} | (20, 0) = [3180] | - |
| 1 | (0, 19) = [19] = {} | (20, 1) = [3181] | 45228 |
| 2 | | | 25229 |
| 3 |
| 4 |
| 5 |
| 6 |
| 7 |
| 8 |
| 9 |
| 10 |
| 11 |
| 12 |
| 13 |
| 14 |
| 15 |
| 16 |
| 17 |
| 18 |
| 19 | (0, 1) = [1] = {} | (20, 19) = [3199] |
| 20 | (0, 0) = [0] = {} | (20, 20) = [3200]

column entries for row 45228:

natural order:
29182
29181
39939
39940
45228
45229

but the column entries must always given in ascending order. Thus the 6 non-zero columns for row 45228 are:

29181
29182
39939
39940
45228
45229

For the next row 45229:
29180
29181
29182
39939
39940
39941
45228
45229
45230


## Overlapping identical and periodic points

For all other overlapping points, identity to the smoothed point is enforced

$$
\v{x}_i - \v{x}_j = 0
$$
or for periodic points
$$
\v{x}_i - \v{x}_j = \v{p}_j
$$

## Laplacian smoothing
$$
\gdef\v#1{\vec{#1}}
$$

As soon as we have multiple block corner points that overlap, we make use of Laplacian smoothing for these corner points,
as we generally cannot construct a 9 point stencil for these corner points.

Multi-block structured grids can generally not be merged into a single block structured grid.

This is the case since for points in which 3 blocks touch, it is not clear in which direction the structured block lines
evolve.


In the current implementation we apply the Laplacian smoothing for the overlapping corner point with the lowest global index.
For all other overlapping points, identity to the smoothed point is enforced as described before.

We use Laplacian smoothing to compute the new position of point $\v{x}$
$$
\v{x} = \frac{1}{N} \sum_{i} \v{x}_i
$$
by summing the positions of the $N$ adjacent points with coordinates $\v{x}_i$.

For periodic points, we can take the periodicity of each adjacent point $\v{p}_j$ into account
$$
\v{x}_i = \frac{1}{N} \sum_{j} \left( \v{x}_j + \v{p}_j \right) \, ,
$$
where $\v{p}_j = \v{0}$ for non-periodic points.
This can be transformed to
$$
-N \v{x}_i + \sum_{j} \v{x}_j = \sum_{j} \v{p}_j \, ,
$$
which yields the following entries for point $i$ in the system matrix $a_{ii} = -N$ as well as
$a_{ij} = 1$ for all accounted adjacent points and for the RHS $b_i = \sum_{j} \v{p}_j$.

In the current implementation we do not account for all adjacent points, but just for the diagnal points.
This has the big advantage that we do not have to include an algorithm that checks for overlapping
adjacent points. A short testing has shown that this approach yields basically the same results compared
to accounting for all adjacent points. However, at some point we might want to extend the treatment.
Another check that is not performed yet is weather a laplacian smoothing shall be conducted. Right now,
this is done for every corner that an overlapping point can be found for; but there might be situation in
non turbomachinery specific configrations where we have overlapping corners that are not block internal and
where a smoothing in this way cannot be done since the laplacian smoothed point needs to be sourrounded by
adjacent points to give resonable results. For other situation, a laplacian smoothing might only make sense
on a vector line or something like this.

This is the situation that the current implementation assumes:
```
x | x
  o-----
x | x
```
This might be something that could be found in other circumstances:
```
x
o----
x
```

Algorithm:

- we first collect all connection endpints in a contigouse array: $\v{x}_{00}, \v{x}_{01}, \v{x}_{10}, \v{x}_{11},...$
- for each endpoint we search the following points for identical indices. If we have a match, we need to check if there is already an overlapping point containg this id is identified. If not, we create a new point that will now have up to 4 overlapping point ids (we check that we do not collect ids multiple times). For exising points, we add try to the id of the connected endpoint to the match (will be the point just before for uneven and the point after for even indices in the entpoint array). We also know the periodicity, from the connection id (endpoint id divided by 4).
- Having found all overlapping points and their periodicity, we can identify the global index of the point to be smoothed and the indices that we need to enforce identity with.
- Collect the adjacent (diagonal) block internal point for each overlapping point (can be 1 point for corners and 2 for points on edges).
- Add the point itself to the stencil and sort in ascending order for the CSC format.
- Compute the RHS based on the collected periodicity info (for each periodically overlapping point, the adjacent block internal point is also periodic and we need to account for that on the RHS).

Example:
Global point 31324

```
*-*-o-*-*     block #n: (n), o: overlapping points, ~~~: periodicity, X: adjacent points to include in stencil, *: block points
    |
* X * X *
    |
(6) | (5)    overlapping points on this side (o) : (6; 20, 158) [43417] & (5; 194, 23) [40078]

~~~~~~~~~

(6) | (4)    overlapping points on this side (o) : (6; 20, 0) [43259] & (4; 0, 20) [31324]
    |
* X * X *
    |
*-*-o-*-*
```

Here, the overlapping points are: 31324, 40078, 43259, 43417.

The stencil is: (4; 0, 20) 31324 (the point to be smoothed), (4; 1, 19) 31344, (5; 193, 22) 40053, (6; 19, 1) 43101, (6; 19, 157) 43257

And the RHS must be: $2 * \v{p}$ since 40053 and 43257 are $\v{p}$ periodic w.r.t. 31324.

For 31324, the coefficients in the system matrix and on the RHS are computed during the initialization.
All other points enforce identity with 31324, with $\v{b}_{40078} = \v{b}_{43417} = \v{p}$ and $\v{b}_{43259} = \v{0}$.
