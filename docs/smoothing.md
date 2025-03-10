# Introduction

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
