const std = @import("std");
const discrete = @import("discrete.zig");
const types = @import("types.zig");
const umfpack = @import("umfpack.zig");
const boundary = @import("boundary.zig");

// This module provides a block-structured elliptic grid generation algorithm.
// It takes a set of boundary points and iteratively adjusts interior points to create
// a smooth grid that conforms to the boundaries.
//
// The algorithm:
// 1. Constructs a sparse, banded matrix representing the discrete Laplace equations
// 2. Uses a compressed row storage format for efficiency with:
//    - Non-zero matrix values
//    - Row/column indices of non-zeros
//    - Cumulative count of non-zeros per row
// 3. Iteratively solves the system:
//    - Assembles matrix coefficients based on current grid point positions
//    - Uses UMFPACK to solve the linear system for new interior point positions
//    - Updates grid points and checks convergence
//
// The matrix structure assembles all internal block points first
// - Corner points (4 neighbors w/o connections)
// - Edge points (6 neighbors w/o connections)
// - Interior points (9 neighbors)
// Then connection data is added as additional rows.
//
// Parameters:
// - points: Current grid point positions
// - iterations: Number of smoothing iterations to perform
// - s,t: Stretching parameters for grid control
//
// This is an example of how to get info for the compressed row format:
// {
//     const i = 18;
//     const j = 1;
//     std.debug.print("Debug ({}, {}) : {}\n", .{ i, j, points.index(.{ i, j }) });
//     std.debug.print("Point: {}\n", .{points.data[points.index(.{ i, j })]});
//     const col_size = points.size[1] - 2;
//     const matrix_idx = (j - 1) + col_size * (i - 1);
//     std.debug.print("Matrix Row: {}\n", .{matrix_idx});
//     std.debug.print("Non Zero Entries: {} bis {}\n", .{ lhs_p[matrix_idx], lhs_p[matrix_idx + 1] });
//     const start: usize = @intCast(lhs_p[matrix_idx]);
//     const end: usize = @intCast(lhs_p[matrix_idx + 1]);
//     for (start..end) |idx| {
//         std.debug.print("  {} : {} : {}\n", .{ idx, lhs_i[idx], lhs_values[idx] });
//     }
//     std.debug.print("RHS: ({}, {})\n", .{ rhs_x[matrix_idx], rhs_y[matrix_idx] });
// }

/// The 9 point stencil data for the point (i,j). The values are stored in an array where the index in
/// the array can be found in the contained index enum (internally going from 0 to 8). For easy access,
/// a get function is provided that allows access based on the index enum.
const StencilData = struct {
    data: [9]f64,

    /// this enum provides the index access for the data
    const index = enum {
        i_j,
        ip1_j,
        im1_j,
        i_jp1,
        i_jm1,
        ip1_jp1,
        ip1_jm1,
        im1_jp1,
        im1_jm1,
    };

    /// returns the contained data based on the index enum, e.g `get(.i_j)`
    fn get(self: StencilData, i: index) f64 {
        return self.data[@intFromEnum(i)];
    }

    fn init(im1_j: types.Vec2d, ip1_j: types.Vec2d, i_jm1: types.Vec2d, i_jp1: types.Vec2d, s: f64, t: f64) StencilData {
        const x_xi = 0.5 * (ip1_j.data[0] - im1_j.data[0]);
        const x_eta = 0.5 * (i_jp1.data[0] - i_jm1.data[0]);
        const y_xi = 0.5 * (ip1_j.data[1] - im1_j.data[1]);
        const y_eta = 0.5 * (i_jp1.data[1] - i_jm1.data[1]);

        const p = x_eta * x_eta + y_eta * y_eta;
        const q = x_xi * x_eta + y_xi * y_eta;
        const r = x_xi * x_xi + y_xi * y_xi;

        return .{
            .data = .{
                -2.0 * p - 2.0 * r, // get(.i_j)
                p + 0.5 * s, //ip1_j
                p - 0.5 * s, //im1_j
                r + 0.5 * t, //i_jp1
                r - 0.5 * t, //i_jm1
                -0.5 * q, //ip1_jp1
                0.5 * q, //ip1_jm1
                0.5 * q, //im1_jp1
                -0.5 * q, //im1_jm1
            },
        };
    }
};

const RowCompressedMatrixSystem2d = struct {
    mesh: *discrete.Mesh,
    allocator: std.mem.Allocator,
    buffer_int: []c_int,
    buffer_float: []f64,

    // cum sum of non-zero entries in columns (size DOF + 1); first value must be zero
    lhs_p: []c_int,

    // row indicies with non-zero values (size non-zero entries)
    lhs_i: []c_int,

    lhs_values: []f64,

    rhs_x: []f64,
    rhs_y: []f64,

    // solution vectors
    x_new: []f64,
    y_new: []f64,

    fn init(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh) !RowCompressedMatrixSystem2d {
        // check that the connection points connect
        const abs_tol = 1e-15;
        for (mesh_data.connections.items, 0..) |connection, connection_idx| {
            const range_0 = connection.data[0];
            const range_1 = connection.data[1];

            var it_0 = range_0.iterate(mesh_data);
            var it_1 = range_1.iterate(mesh_data);

            var point_idx: usize = 0;
            while (true) {
                const p_0 = it_0.next() orelse {
                    std.debug.assert(it_1.next() == null);
                    break;
                };
                const p_1 = it_1.next().?;

                const x_0 = mesh_data.blocks.items[range_0.block].points.data[p_0];
                const x_1 = mesh_data.blocks.items[range_1.block].points.data[p_1];

                if (!types.eqlApprox(x_0, x_1, abs_tol)) {
                    std.debug.panic("non matching points for connection {} point {}:\n\t{}\n\t{}\n", .{ connection_idx, point_idx, x_0, x_1 });
                }

                point_idx += 1;
            }
        }

        // // create boundary point buffer
        // // boundary points | side [0] neighbors | side [1] neighbors
        // var buffer_len: usize = 0;
        // for (mesh_data.connections.items) |connection| {
        //     const len = connection.len();
        //     buffer_len += len;
        // }

        // var boundary_point_buffer = try allocator.alloc(types.Vec2d, buffer_len * 3);
        // defer allocator.free(boundary_point_buffer);

        // var idx: usize = 0;
        // for (mesh_data.connections.items) |connection| {
        //     // boundary points (taken from side 0)
        //     {
        //         const range_0 = connection.data[0];
        //         const point_data = mesh_data.blocks.items[range_0.block].points.data;
        //         var it = range_0.iterate(mesh_data);
        //         while (it.next()) |point_idx| {
        //             boundary_point_buffer[idx] = point_data[point_idx];
        //             idx += 1;
        //         }
        //     }

        //     // boundary point neighbors from side 0 and 1
        //     for (0..2) |side_idx| {
        //         const range = connection.data[side_idx];
        //         const points = mesh_data.blocks.items[range.block].points;

        //         const shift: isize = switch (range.side) {
        //             .i_min => @intCast(points.size[1]),
        //             .i_max => -@as(isize, @intCast(points.size[1])),
        //             .j_min => 1,
        //             .j_max => -1,
        //         };

        //         var it = range.iterate(mesh_data);
        //         while (it.next()) |point_idx| {
        //             const i: isize = @as(isize, @intCast(point_idx)) + shift;
        //             boundary_point_buffer[idx] = points.data[@intCast(i)];
        //         }
        //     }
        // }

        // // collect boundary point connections
        // var boundary_point_connections = try boundary.PointData(std.BoundedArray(usize, 4)).init(allocator, mesh_data);
        // defer boundary_point_connections.deinit();

        // for (mesh_data.connections.items, 0..) |connection, connection_idx| {
        //     for (0..2) |side_idx| {
        //         var it = boundary_point_connections.iterateRange(connection.data[side_idx]);
        //         while (it.nextPtr()) |point_connection_data| {
        //             try point_connection_data.append(connection_idx);
        //         }
        //     }
        // }

        // // tag boundary points based on connections
        // var boundary_point_kind = try boundary.PointData(BoundaryPointProp).init(allocator, mesh_data);
        // defer boundary_point_kind.deinit();

        // for (boundary_point_connections.buffer, boundary_point_kind.buffer) |connections, *kind| {
        //     switch (connections.len) {
        //         0 => kind.* = .fix,
        //         1 => kind.* = .interface,
        //         else => kind.* = .junction,
        //     }
        // }

        // first we add block internal points,
        // then we add block boundary points to be solved based on conenction data

        var dof: usize = 0;
        var non_zero_entries_capacity: usize = 0;

        for (mesh_data.blocks.items) |b| {
            const points = b.points;
            const block_dof = (points.size[0] - 2) * (points.size[1] - 2);
            dof += block_dof;

            // this is a max limit of non zero entries per point
            // the true number will be less dependent on how many
            // neighboring points need to bes solved
            const block_non_zero_entries = dof * 9;
            non_zero_entries_capacity += block_non_zero_entries;
        }

        for (mesh_data.connections.items) |connection| {
            // for now we assume that only connection internal points will be solve
            const internal_points = connection.lenInternal();
            dof += internal_points;
            non_zero_entries_capacity += internal_points * 9;
        }

        const buffer_int = try allocator.alloc(c_int, (dof + 1) + non_zero_entries_capacity);

        // left hand side (lhs) in compressed row form
        var lhs_p = buffer_int[0 .. dof + 1]; // cum sum of non-zero entries in columns
        lhs_p[0] = 0; // first value must be zero
        const lhs_i = buffer_int[dof + 1 ..]; // row indicies with non-zero values
        std.debug.assert(lhs_i.len == non_zero_entries_capacity);

        const buffer_float = try allocator.alloc(f64, 4 * dof + non_zero_entries_capacity);

        const x_new = buffer_float[0..dof];
        const y_new = buffer_float[dof .. 2 * dof];

        const rhs_x = buffer_float[2 * dof .. 3 * dof];
        const rhs_y = buffer_float[3 * dof .. 4 * dof];

        const lhs_values = buffer_float[4 * dof ..];

        try RowCompressedMatrixSystem2d.nonZeroMatrixEntries(allocator, mesh_data, lhs_p, lhs_i, dof, non_zero_entries_capacity);

        return .{
            .mesh = mesh_data,
            .allocator = allocator,
            .buffer_int = buffer_int,
            .buffer_float = buffer_float,
            .lhs_p = lhs_p,
            .lhs_i = lhs_i,
            .lhs_values = lhs_values,
            .rhs_x = rhs_x,
            .rhs_y = rhs_y,
            .x_new = x_new,
            .y_new = y_new,
        };
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        defer self.allocator.free(self.buffer_int);
        defer self.allocator.free(self.buffer_float);
    }

    fn nonZeroMatrixEntries(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh, lhs_p: []c_int, lhs_i: []c_int, dof: usize, non_zero_entries: usize) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_i[0..]));
        var row_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries);

        var matrix_idx: c_int = -1;

        for (mesh_data.blocks.items) |b| {
            const point_size_j: c_int = @intCast(b.points.size[1]);
            const col_size = point_size_j - 2;

            // corner A[0, 0] = (1, 1)
            {
                matrix_idx += 1; // point_idx = point_size_j + 1;
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // TODO make entries and RHS dependent on the block boundaries

            // edge A[0, j] = (1, j)
            for (2..b.points.size[1] - 2) |_| {
                matrix_idx += 1; // points.size[1] + j
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // corner A[0, m] = (1, J-2)
            {
                matrix_idx += 1; // points.size[1] + points.size[1] - 2
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // loop over i
            for (2..b.points.size[0] - 2) |_| {

                // edge A[i, 0] = (i, 0)
                {
                    matrix_idx += 1;
                    try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                    try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                    try row_entries.append(matrix_idx); // A[i, j]
                    try row_entries.append(matrix_idx + 1); // A[i, j+1]
                    try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                    try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                    try row_count.append(@intCast(row_entries.items.len));
                }

                // internal points
                for (2..b.points.size[1] - 2) |_| {
                    // loop over j
                    matrix_idx += 1;
                    try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                    try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                    try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                    try row_entries.append(matrix_idx - 1); // A[i, j-1]
                    try row_entries.append(matrix_idx); // A[i, j]
                    try row_entries.append(matrix_idx + 1); // A[i, j+1]
                    try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                    try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                    try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                    try row_count.append(@intCast(row_entries.items.len));
                }

                // edge A[i, m] = (i, J-2)
                {
                    matrix_idx += 1;
                    try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                    try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                    try row_entries.append(matrix_idx - 1); // A[i, j-1]
                    try row_entries.append(matrix_idx); // A[i, j]
                    try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                    try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                    try row_count.append(@intCast(row_entries.items.len));
                }
            }

            // corner A[n, 0] = (I-2, 1)
            {
                matrix_idx += 1;
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // edge A[n, j] = (I-2, j)
            for (2..b.points.size[1] - 2) |_| {
                matrix_idx += 1; // points.size[1] + j
                try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // corner A[n, m] = (I-1, J-1)
            {
                matrix_idx += 1; // points.size[1] + points.size[1] - 2
                try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_count.append(@intCast(row_entries.items.len));
            }
        }

        // create a map that gives the first entry in the system matrix for a given block.
        // This entry id corresponds to the row index of the first internal point of the matrix.
        // All following internal points are ordered in column major format (j-index first).
        // This map allows to compute the matrix index for a given internal point
        //
        // TODO move this to an own structure that allows the mapping of internal points to matrix indices
        var block_2_matrix_start_idx_map = try allocator.alloc(usize, mesh_data.blocks.items.len);
        defer allocator.free(block_2_matrix_start_idx_map);

        {
            var internal_dof: usize = 0;
            for (mesh_data.blocks.items, 0..) |block, block_idx| {
                block_2_matrix_start_idx_map[block_idx] = internal_dof;
                const points = block.points;
                const block_dof = (points.size[0] - 2) * (points.size[1] - 2);
                internal_dof += block_dof;
            }
        }

        for (mesh_data.connections.items) |connection| {

            // we assume that all connection internal points must be solved

            // we rely on block 0 matrix entries to be ahead of the block 1 matrix entries,
            // otherwise we run into an invalid matrix error
            // We could handle this here dynamically in the future:
            // - move index 0 to a variable lower_idx_block
            // - move index 1 to a variable higher_idx_block
            std.debug.assert(connection.data[0].block < connection.data[1].block);

            // we need at least the two extreme points handled in the current implementation;
            // it should be simple to add handling of the edge cases.
            std.debug.assert(connection.lenInternal() > 3);

            // ignore first and last point (assuming these points to be fixed)
            const internal_ranges = connection.internalRanges();

            var it_0 = RangeNeighborMatrixIndexIterator.init(internal_ranges[0], mesh_data, block_2_matrix_start_idx_map);
            var it_1 = RangeNeighborMatrixIndexIterator.init(internal_ranges[1], mesh_data, block_2_matrix_start_idx_map);

            // NOTE: we need to make sure that the matrix entries are in ascending order
            // as required by the compressed row format
            const matrix_increment_0: c_int = @intCast(it_0.increment);
            const matrix_increment_1: c_int = @intCast(it_1.increment);

            // first internal point has first point fixed.
            {
                //           |
                //       connection
                //           |
                //           V
                //
                // [fixed]                [fixed]          [fixed]
                // (matrix_idx_0)         (matrix_idx)     (matrix_idx_1)
                // (matrix_idx_0 + inc_0) (matrix_idx + 1) (matrix_idx_1 + inc_1)

                const matrix_idx_0: c_int = @intCast(it_0.next().?);
                const matrix_idx_1: c_int = @intCast(it_1.next().?);

                // connection point's  internal points neighbor matrix entries
                const entries_for_block_0 = [2]c_int{ matrix_idx_0, matrix_idx_0 + matrix_increment_0 };
                const entries_for_block_1 = [2]c_int{ matrix_idx_1, matrix_idx_1 + matrix_increment_1 };
                try row_entries.append(@min(entries_for_block_0[0], entries_for_block_0[1]));
                try row_entries.append(@max(entries_for_block_0[0], entries_for_block_0[1]));
                try row_entries.append(@min(entries_for_block_1[0], entries_for_block_1[1]));
                try row_entries.append(@max(entries_for_block_1[0], entries_for_block_1[1]));

                matrix_idx += 1;

                // new entries for the connection (this point and connection neighbor)
                try row_entries.append(matrix_idx);
                try row_entries.append(matrix_idx + 1);

                try row_count.append(@intCast(row_entries.items.len));
            }

            const matrix_increment_0_abs: c_int = @intCast(@abs(matrix_increment_0));
            const matrix_increment_1_abs: c_int = @intCast(@abs(matrix_increment_1));

            // loop edge
            for (0..it_0.count - 1) |_| {
                //           |
                //       connection
                //           |
                //           V
                //
                // (matrix_idx_0 - inc_0) (matrix_idx - 1) (matrix_idx_1 - inc_1)
                // (matrix_idx_0)           (matrix_idx)   (matrix_idx_1)
                // (matrix_idx_0 + inc_0) (matrix_idx + 1) (matrix_idx_1 + inc_1)

                const matrix_idx_0: c_int = @intCast(it_0.next().?);
                const matrix_idx_1: c_int = @intCast(it_1.next().?);

                try row_entries.append(matrix_idx_0 - matrix_increment_0_abs);
                try row_entries.append(matrix_idx_0);
                try row_entries.append(matrix_idx_0 + matrix_increment_0_abs);

                try row_entries.append(matrix_idx_1 - matrix_increment_1_abs);
                try row_entries.append(matrix_idx_1);
                try row_entries.append(matrix_idx_1 + matrix_increment_1_abs);

                matrix_idx += 1;
                try row_entries.append(matrix_idx - 1);
                try row_entries.append(matrix_idx);
                try row_entries.append(matrix_idx + 1);

                try row_count.append(@intCast(row_entries.items.len));
            }

            // second to last point has last point fixed.
            {
                //           |
                //       connection
                //           |
                //           V
                //
                // (matrix_idx_0 - inc_0) (matrix_idx - 1) (matrix_idx_1 - inc_1)
                // (matrix_idx_0)         (matrix_idx)     (matrix_idx_1)
                // [fixed]                [fixed]          [fixed]

                const matrix_idx_0: c_int = @intCast(it_0.next().?);
                const matrix_idx_1: c_int = @intCast(it_1.next().?);

                // connection point's internal points neighbor matrix entries
                const entries_for_block_0 = [2]c_int{ matrix_idx_0, matrix_idx_0 - matrix_increment_0 };
                const entries_for_block_1 = [2]c_int{ matrix_idx_1, matrix_idx_1 - matrix_increment_1 };
                try row_entries.append(@min(entries_for_block_0[0], entries_for_block_0[1]));
                try row_entries.append(@max(entries_for_block_0[0], entries_for_block_0[1]));
                try row_entries.append(@min(entries_for_block_1[0], entries_for_block_1[1]));
                try row_entries.append(@max(entries_for_block_1[0], entries_for_block_1[1]));

                matrix_idx += 1;

                // new entries for the connection (this point and connection neighbor)
                try row_entries.append(matrix_idx - 1);
                try row_entries.append(matrix_idx);

                try row_count.append(@intCast(row_entries.items.len));
            }
        }
    }

    fn fillAndSolve(self: *RowCompressedMatrixSystem2d) !void {
        var lhs_values = self.lhs_values;
        var rhs_x = self.rhs_x;
        var rhs_y = self.rhs_y;

        var matrix_idx: usize = 0;
        var rhs_idx: usize = 0;

        // laplace conditions (zero intialization)
        const s = 0.0;
        const t = 0.0;

        for (self.mesh.blocks.items) |b| {
            const points = b.points;

            const point_size_j = points.size[1];
            var points_idx = point_size_j + 1;

            // corner A[0, 0] = (1, 1)
            {
                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.i_j); // A[i, j]
                lhs_values[matrix_idx + 1] = values.get(.i_jp1); // A[i, j+1]
                lhs_values[matrix_idx + 2] = values.get(.ip1_j); // A[i+1, j]
                lhs_values[matrix_idx + 3] = values.get(.ip1_jp1); // A[i+1, j+1]

                matrix_idx += 4;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j-1] A[i+1, j-1]

                rhs_x[rhs_idx + 0] = -(values.get(.im1_jm1) * im1_jm1.data[0] +
                    values.get(.im1_j) * im1_j.data[0] +
                    values.get(.im1_jp1) * im1_jp1.data[0] +
                    values.get(.i_jm1) * i_jm1.data[0] +
                    values.get(.ip1_jm1) * ip1_jm1.data[0]);

                rhs_y[rhs_idx + 0] = -(values.get(.im1_jm1) * im1_jm1.data[1] +
                    values.get(.im1_j) * im1_j.data[1] +
                    values.get(.im1_jp1) * im1_jp1.data[1] +
                    values.get(.i_jm1) * i_jm1.data[1] +
                    values.get(.ip1_jm1) * ip1_jm1.data[1]);

                rhs_idx += 1;
            }

            // edge A[0, j] = (1, j)
            for (2..points.size[1] - 2) |_| {
                points_idx += 1; // points.size[1] + j

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_j = points.data[points_idx + points.size[1]];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.i_jm1); // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.get(.i_j); // A[i, j]
                lhs_values[matrix_idx + 2] = values.get(.i_jp1); // A[i, j+1]
                lhs_values[matrix_idx + 3] = values.get(.ip1_jm1); // A[i+1, j-1]
                lhs_values[matrix_idx + 4] = values.get(.ip1_j); // A[i+1, j]
                lhs_values[matrix_idx + 5] = values.get(.ip1_jp1); // A[i+1, j+1]

                matrix_idx += 6;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[0] +
                    values.get(.im1_j) * im1_j.data[0] +
                    values.get(.im1_jp1) * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[1] +
                    values.get(.im1_j) * im1_j.data[1] +
                    values.get(.im1_jp1) * im1_jp1.data[1]);

                rhs_idx += 1;
            }

            // corner A[0, m] = (1, J-2)
            {
                points_idx += 1; // points.size[1] + points.size[1] - 2

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.i_jm1); // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.get(.i_j); // A[i, j]
                lhs_values[matrix_idx + 2] = values.get(.ip1_jm1); // A[i+1, j-1]
                lhs_values[matrix_idx + 3] = values.get(.ip1_j); // A[i+1, j]

                matrix_idx += 4;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j+1] A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[0] +
                    values.get(.im1_j) * im1_j.data[0] +
                    values.get(.im1_jp1) * im1_jp1.data[0] +
                    values.get(.i_jp1) * i_jp1.data[0] +
                    values.get(.ip1_jp1) * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[1] +
                    values.get(.im1_j) * im1_j.data[1] +
                    values.get(.im1_jp1) * im1_jp1.data[1] +
                    values.get(.i_jp1) * i_jp1.data[1] +
                    values.get(.ip1_jp1) * ip1_jp1.data[1]);

                rhs_idx += 1;
            }

            // loop over i
            for (2..points.size[0] - 2) |_| {
                points_idx += 2;

                // edge A[i, 0] = (i, 0)
                {
                    points_idx += 1;

                    const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                    const im1_j = points.data[points_idx - points.size[1]];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                    const ip1_j = points.data[points_idx + points.size[1]];

                    const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.get(.im1_j); // A[i-1, j]
                    lhs_values[matrix_idx + 1] = values.get(.im1_jp1); // A[i-1, j+1]
                    lhs_values[matrix_idx + 2] = values.get(.i_j); // A[i, j]
                    lhs_values[matrix_idx + 3] = values.get(.i_jp1); // A[i, j+1]
                    lhs_values[matrix_idx + 4] = values.get(.ip1_j); // A[i+1, j]
                    lhs_values[matrix_idx + 5] = values.get(.ip1_jp1); // A[i+1, j+1]

                    matrix_idx += 6;

                    // RHS: A[i-1, j-1], A[i, j-1], A[i+1, j-1]

                    rhs_x[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[0] +
                        values.get(.i_jm1) * i_jm1.data[0] +
                        values.get(.ip1_jm1) * ip1_jm1.data[0]);

                    rhs_y[rhs_idx] = -(values.get(.im1_jm1) * im1_jm1.data[1] +
                        values.get(.i_jm1) * i_jm1.data[1] +
                        values.get(.ip1_jm1) * ip1_jm1.data[1]);

                    rhs_idx += 1;
                }

                // internal points
                for (2..points.size[1] - 2) |_| {
                    // loop over j
                    points_idx += 1;

                    const im1_j = points.data[points_idx - points.size[1]];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_j = points.data[points_idx + points.size[1]];

                    const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.get(.im1_jm1); // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.get(.im1_j); // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.get(.im1_jp1); // A[i-1, j+1]
                    lhs_values[matrix_idx + 3] = values.get(.i_jm1); // A[i, j-1]
                    lhs_values[matrix_idx + 4] = values.get(.i_j); // A[i, j]
                    lhs_values[matrix_idx + 5] = values.get(.i_jp1); // A[i, j+1]
                    lhs_values[matrix_idx + 6] = values.get(.ip1_jm1); // A[i+1, j-1]
                    lhs_values[matrix_idx + 7] = values.get(.ip1_j); // A[i+1, j]
                    lhs_values[matrix_idx + 8] = values.get(.ip1_jp1); // A[i+1, j+1]

                    matrix_idx += 9;

                    rhs_x[rhs_idx] = 0;
                    rhs_y[rhs_idx] = 0;

                    rhs_idx += 1;
                }

                // edge A[i, m] = (i, J-2)
                {
                    points_idx += 1;

                    const im1_j = points.data[points_idx - points.size[1]];
                    const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_j = points.data[points_idx + points.size[1]];
                    const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                    const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.get(.im1_jm1); // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.get(.im1_j); // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.get(.i_jm1); // A[i, j-1]
                    lhs_values[matrix_idx + 3] = values.get(.i_j); // A[i, j]
                    lhs_values[matrix_idx + 4] = values.get(.ip1_jm1); // A[i+1, j-1]
                    lhs_values[matrix_idx + 5] = values.get(.ip1_j); // A[i+1, j]

                    matrix_idx += 6;

                    // RHS: A[i-1, j+1], A[i, j+1], A[i+1, j+1]

                    rhs_x[rhs_idx] = -(values.get(.im1_jp1) * im1_jp1.data[0] +
                        values.get(.i_jp1) * i_jp1.data[0] +
                        values.get(.ip1_jp1) * ip1_jp1.data[0]);

                    rhs_y[rhs_idx] = -(values.get(.im1_jp1) * im1_jp1.data[1] +
                        values.get(.i_jp1) * i_jp1.data[1] +
                        values.get(.ip1_jp1) * ip1_jp1.data[1]);

                    rhs_idx += 1;
                }
            }

            points_idx += 2;

            // corner A[n, 0] = (I-2, 1)
            {
                points_idx += 1;

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.im1_j); // A[i-1, j]
                lhs_values[matrix_idx + 1] = values.get(.im1_jp1); // A[i-1, j+1]
                lhs_values[matrix_idx + 2] = values.get(.i_j); // A[i, j]
                lhs_values[matrix_idx + 3] = values.get(.i_jp1); // A[i, j+1]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j-1], A[i-1, j-1]

                rhs_x[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[0] +
                    values.get(.ip1_j) * ip1_j.data[0] +
                    values.get(.ip1_jp1) * ip1_jp1.data[0] +
                    values.get(.i_jm1) * i_jm1.data[0] +
                    values.get(.im1_jm1) * im1_jm1.data[0]);

                rhs_y[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[1] +
                    values.get(.ip1_j) * ip1_j.data[1] +
                    values.get(.ip1_jp1) * ip1_jp1.data[1] +
                    values.get(.i_jm1) * i_jm1.data[1] +
                    values.get(.im1_jm1) * im1_jm1.data[1]);

                rhs_idx += 1;
            }

            // edge A[n, j] = (I-2, j)
            for (2..points.size[1] - 2) |_| {
                points_idx += 1; // points.size[1] + j

                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.im1_jm1); // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.get(.im1_j); // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.get(.im1_jp1); // A[i-1, j+1]
                lhs_values[matrix_idx + 3] = values.get(.i_jm1); // A[i, j-1]
                lhs_values[matrix_idx + 4] = values.get(.i_j); // A[i, j]
                lhs_values[matrix_idx + 5] = values.get(.i_jp1); // A[i, j+1]

                matrix_idx += 6;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[0] +
                    values.get(.ip1_j) * ip1_j.data[0] +
                    values.get(.ip1_jp1) * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[1] +
                    values.get(.ip1_j) * ip1_j.data[1] +
                    values.get(.ip1_jp1) * ip1_jp1.data[1]);

                rhs_idx += 1;
            }

            // corner A[n, m] = (I-2, J-2)
            {
                points_idx += 1; // points.size[1] + points.size[1] - 2

                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.get(.im1_jm1); // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.get(.im1_j); // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.get(.i_jm1); // A[i, j-1]
                lhs_values[matrix_idx + 3] = values.get(.i_j); // A[i, j]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j+1], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[0] +
                    values.get(.ip1_j) * ip1_j.data[0] +
                    values.get(.ip1_jp1) * ip1_jp1.data[0] +
                    values.get(.i_jp1) * i_jp1.data[0] +
                    values.get(.im1_jp1) * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.get(.ip1_jm1) * ip1_jm1.data[1] +
                    values.get(.ip1_j) * ip1_j.data[1] +
                    values.get(.ip1_jp1) * ip1_jp1.data[1] +
                    values.get(.i_jp1) * i_jp1.data[1] +
                    values.get(.im1_jp1) * im1_jp1.data[1]);

                rhs_idx += 1;
            }
        }

        for (self.mesh.connections.items) |connection| {
            const point_data = [2][]types.Vec2d{
                self.mesh.blocks.items[connection.data[0].block].points.data,
                self.mesh.blocks.items[connection.data[1].block].points.data,
            };

            var it = RangeFillMatrixIterator.init(connection, self.mesh);

            // ignore first point (assuming it is fixed)
            _ = it.next();

            //           |
            //       connection
            //           |
            //           V
            //
            // A(i-1,j-1) A(i-1,j) A(i-1,j+1)
            // A(i,j-1)   A(i,j)   A(i,j+1)
            // A(i+1,j-1) A(i+1,j) A(i+1,j+1)
            //
            // or
            //
            // A(matrix_idx_0 - inc_0) A(matrix_idx - 1) A(matrix_idx_1 - inc_1)
            // A(matrix_idx_0)         A(matrix_idx)     A(matrix_idx_1)
            // A(matrix_idx_0 + inc_0) A(matrix_idx + 1) A(matrix_idx_1 + inc_1)

            // first internal point has first point fixed.
            {
                const boundary_idx = it.next().?;
                std.debug.assert(types.eqlApprox(point_data[0][boundary_idx[0]], point_data[1][boundary_idx[1]], 1e-12));
                const boundary_idx_0: isize = @intCast(boundary_idx[0]);
                const boundary_idx_1: isize = @intCast(boundary_idx[1]);

                const idx_im1_jm1: usize = @intCast(boundary_idx_0 - it.in_connection_direction_shift[0] + it.first_internal_point_shift[0]);
                const idx_im1_j: usize = @intCast(boundary_idx_0 - it.in_connection_direction_shift[0]);

                std.debug.print("{}\n", .{it});
                std.debug.print("point idx: {} {}\n", .{ idx_im1_jm1, idx_im1_j });

                const im1_jm1 = point_data[0][idx_im1_jm1];
                const im1_j = point_data[0][idx_im1_j];
                const i_jm1 = point_data[0][@intCast(boundary_idx_0 + it.first_internal_point_shift[0])];
                const ip1_j = point_data[0][@intCast(boundary_idx_0 + it.in_connection_direction_shift[0])];
                const i_jp1 = point_data[1][@intCast(boundary_idx_1 + it.first_internal_point_shift[1])];
                const im1_jp1 = point_data[1][@intCast(boundary_idx_1 - it.in_connection_direction_shift[1] + it.first_internal_point_shift[1])];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                if (it.in_connection_direction_shift[0] > 0) {
                    lhs_values[matrix_idx + 0] = values.get(.i_jm1); // A[matrix_idx, matrix_idx_0]
                    lhs_values[matrix_idx + 1] = values.get(.ip1_jm1); // A[matrix_idx, matrix_idx_0 + matrix_increment_0]
                } else {
                    lhs_values[matrix_idx + 0] = values.get(.ip1_jm1); // A[matrix_idx, matrix_idx_0 + matrix_increment_0]
                    lhs_values[matrix_idx + 1] = values.get(.i_jm1); // A[matrix_idx, matrix_idx_0]
                }

                if (it.in_connection_direction_shift[1] > 0) {
                    lhs_values[matrix_idx + 2] = values.get(.i_jp1); // A[matrix_idx, matrix_idx_1]
                    lhs_values[matrix_idx + 3] = values.get(.ip1_jp1); // A[matrix_idx, matrix_idx_1 + matrix_increment_1]
                } else {
                    lhs_values[matrix_idx + 2] = values.get(.ip1_jp1); // A[matrix_idx, matrix_idx_1 + matrix_increment_1]
                    lhs_values[matrix_idx + 3] = values.get(.i_jp1); // A[matrix_idx, matrix_idx_1]
                }

                lhs_values[matrix_idx + 4] = values.get(.i_j); // A[matrix_idx, matrix_idx]
                lhs_values[matrix_idx + 5] = values.get(.ip1_j); // A[matrix_idx, matrix_idx + 1]

                matrix_idx += 6;

                // RHS: fixed first point of the connection including its neighbors (i-1,j-1) (i-1,j) (i-1,j+1)

                rhs_x[rhs_idx + 0] = -(values.get(.im1_jm1) * im1_jm1.data[0] + values.get(.im1_j) * im1_j.data[0] + values.get(.im1_jp1) * im1_jp1.data[0]);
                rhs_y[rhs_idx + 0] = -(values.get(.im1_jm1) * im1_jm1.data[1] + values.get(.im1_j) * im1_j.data[1] + values.get(.im1_jp1) * im1_jp1.data[1]);

                rhs_idx += 1;
            }

            // loop edge
            for (0..it.count - 1) |_| {
                const boundary_idx = it.next().?;
                std.debug.assert(types.eqlApprox(point_data[0][boundary_idx[0]], point_data[1][boundary_idx[1]], 1e-12));
                const boundary_idx_0: isize = @intCast(boundary_idx[0]);
                const boundary_idx_1: isize = @intCast(boundary_idx[1]);

                const im1_j = point_data[0][@intCast(boundary_idx_0 - it.in_connection_direction_shift[0])];
                const i_jm1 = point_data[0][@intCast(boundary_idx_0 + it.first_internal_point_shift[0])];
                const ip1_j = point_data[0][@intCast(boundary_idx_0 + it.in_connection_direction_shift[0])];
                const i_jp1 = point_data[1][@intCast(boundary_idx_1 + it.first_internal_point_shift[1])];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 2] = values.get(.im1_jm1); // A[matrix_idx, matrix_idx_0 - matrix_increment_0]
                lhs_values[matrix_idx + 1] = values.get(.i_jm1); // A[matrix_idx, matrix_idx_0]
                lhs_values[matrix_idx + 0] = values.get(.ip1_jm1); // A[matrix_idx, matrix_idx_0 + matrix_increment_0]

                lhs_values[matrix_idx + 3] = values.get(.im1_jp1); // A[matrix_idx, matrix_idx_1 - matrix_increment_1]
                lhs_values[matrix_idx + 4] = values.get(.i_jp1); // A[matrix_idx, matrix_idx_1]
                lhs_values[matrix_idx + 5] = values.get(.ip1_jp1); // A[matrix_idx, matrix_idx_1 + matrix_increment_1]

                lhs_values[matrix_idx + 6] = values.get(.im1_j); // A[matrix_idx, matrix_idx - 1]
                lhs_values[matrix_idx + 7] = values.get(.i_j); // A[matrix_idx, matrix_idx]
                lhs_values[matrix_idx + 8] = values.get(.ip1_j); // A[matrix_idx, matrix_idx + 1]

                matrix_idx += 9;

                rhs_x[rhs_idx] = 0;
                rhs_y[rhs_idx] = 0;

                rhs_idx += 1;
            }

            // second to last point has last point fixed.
            {
                const boundary_idx = it.next().?;
                std.debug.assert(types.eqlApprox(point_data[0][boundary_idx[0]], point_data[1][boundary_idx[1]], 1e-12));
                const boundary_idx_0: isize = @intCast(boundary_idx[0]);
                const boundary_idx_1: isize = @intCast(boundary_idx[1]);

                const im1_j = point_data[0][@intCast(boundary_idx_0 - it.in_connection_direction_shift[0])];
                const i_jm1 = point_data[0][@intCast(boundary_idx_0 + it.first_internal_point_shift[0])];
                const ip1_j = point_data[0][@intCast(boundary_idx_0 + it.in_connection_direction_shift[0])];
                const ip1_jm1 = point_data[0][@intCast(boundary_idx_0 + it.in_connection_direction_shift[0] + it.first_internal_point_shift[0])];
                const i_jp1 = point_data[1][@intCast(boundary_idx_1 + it.first_internal_point_shift[1])];
                const ip1_jp1 = point_data[1][@intCast(boundary_idx_1 + it.in_connection_direction_shift[1] + it.first_internal_point_shift[1])];

                const values = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 1] = values.get(.im1_jm1); // A[matrix_idx, matrix_idx_0 - matrix_increment_0]
                lhs_values[matrix_idx + 0] = values.get(.i_jm1); // A[matrix_idx, matrix_idx_0]
                lhs_values[matrix_idx + 2] = values.get(.im1_jp1); // A[matrix_idx, matrix_idx_1 - matrix_increment_1]
                lhs_values[matrix_idx + 3] = values.get(.i_jp1); // A[matrix_idx, matrix_idx_1]
                lhs_values[matrix_idx + 4] = values.get(.im1_j); // A[matrix_idx, matrix_idx - 1]
                lhs_values[matrix_idx + 5] = values.get(.i_j); // A[matrix_idx, matrix_idx]

                matrix_idx += 6;

                // RHS: fixed first point of the connection including its neighbors (i-1,j-1) (i-1,j) (i-1,j+1)

                rhs_x[rhs_idx + 0] = -(values.get(.ip1_jm1) * ip1_jm1.data[0] + values.get(.ip1_j) * ip1_j.data[0] + values.get(.ip1_jp1) * ip1_jp1.data[0]);
                rhs_y[rhs_idx + 0] = -(values.get(.ip1_jm1) * ip1_jm1.data[1] + values.get(.ip1_j) * ip1_j.data[1] + values.get(.ip1_jp1) * ip1_jp1.data[1]);

                rhs_idx += 1;
            }

            // std.debug.assert(it.count == 1);
        }

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, lhs_values, rhs_x, self.x_new[0..], rhs_y, self.y_new[0..]);
    }
};

const BoundaryPointProp = enum {
    fix,
    interface,
    junction,

    periodic,

    // sliding
};

pub fn mesh(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh, iterations: usize) !void {
    var system = try RowCompressedMatrixSystem2d.init(allocator, mesh_data);
    defer system.deinit();

    // iterate and fill matrix values
    for (0..iterations) |n| {
        std.debug.print("  iteration: {}\n", .{n});

        try system.fillAndSolve();

        var x_norm_sqr: f64 = 0.0;
        var y_norm_sqr: f64 = 0.0;

        // TODO inefficient - enhance
        {
            var matrix_idx: usize = 0;
            for (mesh_data.blocks.items) |b| {
                const points = b.points;
                for (1..points.size[0] - 1) |i| {
                    for (1..points.size[1] - 1) |j| {
                        const p = points.data[points.index(.{ i, j })];
                        const dx = p.data[0] - system.x_new[matrix_idx];
                        const dy = p.data[1] - system.y_new[matrix_idx];

                        x_norm_sqr += dx * dx;
                        y_norm_sqr += dy * dy;

                        matrix_idx += 1;
                    }
                }
            }

            for (mesh_data.connections.items) |connection| {
                // we only need to use one block here as the connected boundary points are identical
                const point_data = mesh_data.blocks.items[connection.data[0].block].points.data;

                var it = connection.data[0].iterate(mesh_data);

                // ignore first point (assuming it is fixed)
                _ = it.next();

                for (0..it.count) |_| {
                    const point_idx = it.next().?;

                    const dx = point_data[point_idx].data[0] - system.x_new[matrix_idx];
                    const dy = point_data[point_idx].data[1] - system.y_new[matrix_idx];

                    x_norm_sqr += dx * dx;
                    y_norm_sqr += dy * dy;

                    matrix_idx += 1;
                }

                // ignore last point (assuming it is fixed)
                std.debug.assert(it.count == 0 and it.idx != null);
            }
        }

        const norm = (x_norm_sqr + y_norm_sqr) * (x_norm_sqr + y_norm_sqr);
        std.debug.print("\tresidual: {}\n", .{norm});

        // copy into coordinate field
        {
            var matrix_idx: usize = 0;
            for (mesh_data.blocks.items) |b| {
                const points = b.points;
                for (1..points.size[0] - 1) |i| {
                    for (1..points.size[1] - 1) |j| {

                        // TODO inefficient - enhance
                        const points_idx = points.index(.{ i, j });

                        points.data[points_idx] = types.Vec2d.init(system.x_new[matrix_idx], system.y_new[matrix_idx]);

                        matrix_idx += 1;
                    }
                }
            }

            for (mesh_data.connections.items) |connection| {
                // we iterate over both connected blocks individually
                const matrix_idx_connection_start = matrix_idx;
                for (0..2) |side_idx| {
                    const point_data = mesh_data.blocks.items[connection.data[side_idx].block].points.data;
                    var it = connection.data[side_idx].iterate(mesh_data);
                    matrix_idx = matrix_idx_connection_start;

                    // ignore first point (assuming it is fixed)
                    _ = it.next();

                    for (0..it.count) |_| {
                        const point_idx = it.next().?;
                        point_data[point_idx] = types.Vec2d.init(system.x_new[matrix_idx], system.y_new[matrix_idx]);
                        matrix_idx += 1;
                    }

                    // ignore last point (assuming it is fixed)
                    std.debug.assert(it.count == 0 and it.idx != null);
                }
            }
        }
    }
}

const RangeNeighborMatrixIndexIterator = struct {
    count: usize,
    increment: isize,
    position: usize,

    fn next(self: *@This()) ?usize {
        if (self.count == 0) return null;
        self.count -= 1;
        const matrix_index = self.position;
        self.position = @intCast(@as(isize, @intCast(self.position)) + self.increment);
        return matrix_index;
    }

    fn init(
        range: boundary.Range,
        mesh_data: *const discrete.Mesh,
        block_2_matrix_start_idx_map: []const usize,
    ) RangeNeighborMatrixIndexIterator {
        const block_size_i = mesh_data.blocks.items[range.block].points.size[0];
        const block_size_j = mesh_data.blocks.items[range.block].points.size[1];

        var first_internal_point_index: types.Index2d = undefined;
        var increment: isize = undefined;
        switch (range.side) {
            .i_min => {
                first_internal_point_index = .{ range.start, 1 };
                increment = @intCast(block_size_j - 2);
            },
            .i_max => {
                first_internal_point_index = .{ range.start, block_size_j - 2 };
                increment = @intCast(block_size_j - 2);
            },
            .j_min => {
                first_internal_point_index = .{ 1, range.start };
                increment = 1;
            },
            .j_max => {
                first_internal_point_index = .{ block_size_i - 2, range.start };
                increment = 1;
            },
        }

        const matrix_idx_start = block_2_matrix_start_idx_map[range.block];
        const matrix_index = matrix_idx_start +
            (first_internal_point_index[0] - 1) * (block_size_j - 2) + (first_internal_point_index[1] - 1);

        var count: usize = 1;
        if (range.start <= range.end) {
            count += range.end - range.start;
        } else {
            increment = -increment;
            count += range.start - range.end;
        }

        return .{
            .count = count,
            .increment = increment,
            .position = matrix_index,
        };
    }
};

const RangeFillMatrixIterator = struct {
    count: usize,
    first_internal_point_shift: [2]isize,
    in_connection_direction_shift: [2]isize,
    position: [2]usize,

    fn next(self: *@This()) ?[2]usize {
        if (self.count == 0) return null;
        const position = self.position;

        self.count -= 1;

        self.position = .{
            @intCast(@as(isize, @intCast(self.position[0])) + self.in_connection_direction_shift[0]),
            @intCast(@as(isize, @intCast(self.position[1])) + self.in_connection_direction_shift[1]),
        };

        return position;
    }

    fn init(
        connection: boundary.Connection,
        mesh_data: *const discrete.Mesh,
    ) RangeFillMatrixIterator {
        var data: RangeFillMatrixIterator = undefined;

        for (0..2) |side_idx| {
            const block_idx = connection.data[side_idx].block;
            const points = mesh_data.blocks.items[block_idx].points;
            const start = connection.data[side_idx].start;
            const end = connection.data[side_idx].end;
            switch (connection.data[side_idx].side) {
                .i_min => {
                    data.first_internal_point_shift[side_idx] = 1;
                    data.in_connection_direction_shift[side_idx] = @intCast(points.size[1]);
                    data.position[side_idx] = points.index(.{ start, 0 });
                },
                .i_max => {
                    data.first_internal_point_shift[side_idx] = -1;
                    data.in_connection_direction_shift[side_idx] = @intCast(points.size[1]);
                    data.position[side_idx] = points.index(.{ start, points.size[1] - 1 });
                },
                .j_min => {
                    data.first_internal_point_shift[side_idx] = @intCast(points.size[1]);
                    data.in_connection_direction_shift[side_idx] = 1;
                    data.position[side_idx] = points.index(.{ 0, start });
                },
                .j_max => {
                    data.first_internal_point_shift[side_idx] = -@as(isize, @intCast(points.size[1]));
                    data.in_connection_direction_shift[side_idx] = 1;
                    data.position[side_idx] = points.index(.{ points.size[0] - 1, start });
                },
            }
            if (start > end) {
                data.in_connection_direction_shift[side_idx] = -data.in_connection_direction_shift[side_idx];
                data.count = start - end;
            } else {
                data.count = end - start;
            }
        }

        return data;
    }
};
