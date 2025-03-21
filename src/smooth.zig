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
// TODO: rework w.r.t. new handling of block connections.

pub fn mesh(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh, iterations: usize) !void {
    var system = try RowCompressedMatrixSystem2d.init(allocator, mesh_data);
    defer system.deinit();

    // iterate and fill matrix values
    for (0..iterations) |n| {
        std.debug.print("  iteration: {}\n", .{n});

        try system.fillAndSolve();

        var x_norm_sqr: f64 = 0.0;
        var y_norm_sqr: f64 = 0.0;

        // TODO: inefficient - enhance
        {
            var row_block_start_idx: usize = 0;
            for (mesh_data.blocks.items) |block| {
                for (1..block.points.size[0] - 1) |i| {
                    for (1..block.points.size[1] - 1) |j| {
                        const point_idx = block.points.index(.{ i, j });
                        const row_idx = row_block_start_idx + point_idx;
                        const p = block.points.data[point_idx];
                        const dx = p.data[0] - system.x_new[row_idx];
                        const dy = p.data[1] - system.y_new[row_idx];

                        x_norm_sqr += dx * dx;
                        y_norm_sqr += dy * dy;
                    }
                }
                row_block_start_idx += block.points.size[0] * block.points.size[1];
            }
        }

        const norm = (x_norm_sqr + y_norm_sqr) * (x_norm_sqr + y_norm_sqr);
        std.debug.print("\tresidual: {}\n", .{norm});

        // copy into coordinate field
        // TODO: inefficient - enhance
        {
            var row_block_start_idx: usize = 0;
            for (mesh_data.blocks.items) |block| {
                for (1..block.points.size[0] - 1) |i| {
                    for (1..block.points.size[1] - 1) |j| {
                        const point_idx = block.points.index(.{ i, j });
                        const row_idx = row_block_start_idx + point_idx;
                        block.points.data[point_idx] = types.Vec2d.init(system.x_new[row_idx], system.y_new[row_idx]);
                    }
                }
                row_block_start_idx += block.points.size[0] * block.points.size[1];
            }
        }
    }
}

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

/// loops all connection and checks that the data contains idential point locations.
/// TODO: remove panic with proper error handling.
fn connectionDataCheck(mesh_data: *const discrete.Mesh) void {
    const abs_tol = 1e-15;
    for (mesh_data.connections.items, 0..) |connection, connection_idx| {
        const range_0 = connection.ranges[0];
        const range_1 = connection.ranges[1];

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
}

const RowCompressedMatrixSystem2d = struct {
    mesh: *discrete.Mesh,
    allocator: std.mem.Allocator,
    buffer_int: []c_int,
    buffer_float: []f64,

    /// cum sum of non-zero entries in columns (size DOF + 1); first value must be zero
    lhs_p: []c_int,

    /// row indicies with non-zero values (size non-zero entries)
    lhs_i: []c_int,

    /// sparse matrix coefficients w.r.t. lhs_i
    lhs_values: []f64,

    rhs_x: []f64,
    rhs_y: []f64,

    // solution vectors
    x_new: []f64,
    y_new: []f64,

    // added for handling multi-block smoothing
    // NOTE: we could save an allocation if changed to c_int and merged with the other allocation.
    //
    /// row index range start for each block (global point index)
    row_idx_range_start_for_each_block: []usize,

    fn init(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh) !RowCompressedMatrixSystem2d {
        connectionDataCheck(mesh_data);

        var dof: usize = 0;
        var non_zero_entries_capacity: usize = 0;

        var row_idx_range_start_for_each_block = try allocator.alloc(usize, mesh_data.blocks.items.len);
        errdefer allocator.free(row_idx_range_start_for_each_block);

        for (mesh_data.blocks.items, 0..) |block, block_idx| {
            const points = block.points;
            const block_dof = points.size[0] * points.size[1];
            dof += block_dof;

            // this is a max limit of non zero entries per point
            // the true number will be less dependent on how many
            // neighboring points need to bes solved
            const block_non_zero_entries = dof * 9;
            non_zero_entries_capacity += block_non_zero_entries;

            row_idx_range_start_for_each_block[block_idx] = dof;
        }

        const buffer_int = try allocator.alloc(c_int, (dof + 1) + non_zero_entries_capacity);
        errdefer allocator.free(buffer_int);

        // left hand side (lhs) in compressed row form
        var lhs_p = buffer_int[0 .. dof + 1]; // cum sum of non-zero entries in columns
        lhs_p[0] = 0; // first value must be zero
        const lhs_i = buffer_int[dof + 1 ..]; // row indicies with non-zero values
        std.debug.assert(lhs_i.len == non_zero_entries_capacity);

        const buffer_float = try allocator.alloc(f64, 4 * dof + non_zero_entries_capacity);
        errdefer allocator.free(buffer_float);

        const x_new = buffer_float[0..dof];
        const y_new = buffer_float[dof .. 2 * dof];

        const rhs_x = buffer_float[2 * dof .. 3 * dof];
        const rhs_y = buffer_float[3 * dof .. 4 * dof];

        const lhs_values = buffer_float[4 * dof ..];

        var system = RowCompressedMatrixSystem2d{
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
            .row_idx_range_start_for_each_block = row_idx_range_start_for_each_block,
        };

        // Handling of the inter-block connections:
        // For now, we will add rows for each connected point to enforce identity of the connected points.
        // Each connected point adds his part of stencil data to the system matrix.
        // This renders the stencil info most easy (of the considered option) at the expense of additional
        // matrix entries.
        const index_converter = try IndexConverter.init(mesh_data, allocator);
        defer index_converter.deinit();

        const connected_points = try BlockBoundaryPointConnections.init(allocator, mesh_data, index_converter);
        defer connected_points.deinit();

        try system.nonZeroMatrixEntries(dof, non_zero_entries_capacity, connected_points);

        return system;
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        defer self.allocator.free(self.buffer_int);
        defer self.allocator.free(self.buffer_float);
        defer self.allocator.free(self.row_idx_range_start_for_each_block);
    }

    /// returns the non zero entries range for the given row (the row index range start for each block can be retrieved from an array in the struct)
    fn nonZeroEntriesRangeStart(self: RowCompressedMatrixSystem2d, row_idx: c_int) usize {
        return @intCast(self.lhs_p[@intCast(row_idx + 1)]);
    }

    fn addOrAllocNonZeroMatrixEntriesForBoundaryPoint(
        non_zero_entries: *std.ArrayList(c_int),
        row_count: *std.ArrayList(c_int),
        boundary_point_idx: *usize,
        connected_points: BlockBoundaryPointConnections,
        row_idx: *c_int, // equivalent to the global point index
    ) !void {
        // for all points on the boundary we need to check how this point is to be added to the system matrix:
        // - fixed point: only a single entry with the current coordinates on the RHS
        // - smoothed: 9 entries for the full stencil
        // - laplacian: TBD (probably also 9 entries) TODO: fill this info.

        switch (connected_points.getKind(boundary_point_idx.*, @intCast(row_idx.*))) {
            .fix => {
                try non_zero_entries.append(row_idx.*); // A[i, j]
            },
            .smooth => {
                // smooth the block boundary point based on the same 9 point stancil that is used for block
                // internal points. The actual data is set in a connection based loop after this undefined
                // allocation.
                _ = try non_zero_entries.addManyAsSlice(9);
            },
            .connect => {
                // enforce the equality for the connected points with two entries.
                try non_zero_entries.append(row_idx.*);
                try non_zero_entries.append(connected_points.data.buffer[boundary_point_idx.*].get(0));
            },
            .junction => {
                unreachable;
                // TODO: connect all points to the index with lowest value.
                // add laplace smoothing for the lowest index: will have as many non zero entries for this row
                // equal to the number of connected points.
            },
        }
        boundary_point_idx.* += 1;

        try row_count.append(@intCast(non_zero_entries.items.len));
        row_idx.* += 1;
    }

    fn addOrAllocNonZeroMatrixEntries(
        self: @This(),
        non_zero_entries: *std.ArrayList(c_int),
        row_count: *std.ArrayList(c_int),
        connected_points: BlockBoundaryPointConnections,
    ) !void {
        var boundary_point_idx: usize = 0;
        var row_idx: c_int = 0;
        for (self.mesh.blocks.items) |block| {
            const col_size: c_int = @intCast(block.points.size[1]);

            // we loop first over j_min then the the middle (including i_min and i_max parts) and then j_max
            //
            //       |  |----------|  |
            //       |  |          |  |
            //       |  |          |  |
            // j_min |  |  middle  |  | j_max
            //       |  |          |  |
            //       |  |          |  |
            //       |  |----------|  |

            // edge j_min
            for (0..block.points.size[1]) |_| {
                try addOrAllocNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, connected_points, &row_idx);
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                try addOrAllocNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, connected_points, &row_idx);

                // internal points with full 9 point stencil
                for (1..block.points.size[1] - 1) |_| {
                    try non_zero_entries.append(row_idx - col_size - 1); // A[i-1, j-1]
                    try non_zero_entries.append(row_idx - col_size); // A[i-1, j]
                    try non_zero_entries.append(row_idx - col_size + 1); // A[i-1, j+1]
                    try non_zero_entries.append(row_idx - 1); // A[i, j-1]
                    try non_zero_entries.append(row_idx); // A[i, j]
                    try non_zero_entries.append(row_idx + 1); // A[i, j+1]
                    try non_zero_entries.append(row_idx + col_size - 1); // A[i+1, j-1]
                    try non_zero_entries.append(row_idx + col_size); // A[i+1, j]
                    try non_zero_entries.append(row_idx + col_size + 1); // A[i+1, j+1]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row_idx += 1;
                }

                // edge i_max
                try addOrAllocNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, connected_points, &row_idx);
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                try addOrAllocNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, connected_points, &row_idx);
            }
        }
    }

    fn fillNonZeroMatrixEntriesForConnections(
        self: @This(),
        non_zero_entries: *std.ArrayList(c_int),
    ) void {
        for (self.mesh.connections.items) |connection| {

            // NOTE: we rely on block 0 matrix entries to be ahead of the block 1 matrix entries,
            // otherwise we run into an invalid matrix error. We could handle this here dynamically in the future:
            // - move index 0 to a variable lower_idx_block
            // - move index 1 to a variable higher_idx_block
            std.debug.assert(connection.ranges[0].block < connection.ranges[1].block);

            // we need at least the two extreme points handled in the current implementation;
            // it should be simple to add handling of the edge cases.
            std.debug.assert(connection.lenInternal() > 3);

            // NOTE: we need to make sure that the matrix entries are in ascending order
            // as required by the compressed row format
            var it = RangeFillMatrixIterator.init(connection, self.mesh);

            // TODO: remove this hardcoded handling
            // we set point 0 to fixed and point 1 to connected to 0
            {
                _ = it.next().?;
                // const connected_boundary_points = it.next().?;

                // const boundary_idx_0: c_int = @intCast(connected_boundary_points[0]);
                // const boundary_idx_1: c_int = @intCast(connected_boundary_points[1]);
                //
                // const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
                // non_zero_entries.items[row_non_zero_entries_start_idx] = boundary_idx_0;
                // non_zero_entries.items[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
            }

            for (0..it.count - 1) |_| {
                // TODO: introduce a seperate type for the different index types.
                const connected_boundary_points_local_idx = it.next().?;

                // local to global index

                const boundary_idx_0: c_int = @intCast(connected_boundary_points_local_idx[0]);
                const boundary_idx_1: c_int = @intCast(connected_boundary_points_local_idx[1]);

                // smooth 1st point
                {
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_0);

                    // TODO: here we just set the data. The allocation must have been done before!
                    // TODO: compute this dynamically. For now this is hardcoded for the first connection.

                    non_zero_entries.items[row_non_zero_entries_start_idx + 0] = boundary_idx_0 + it.first_internal_point_shift[0] - it.in_connection_direction_shift[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 1] = boundary_idx_0 + it.first_internal_point_shift[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 2] = boundary_idx_0 + it.first_internal_point_shift[0] + it.in_connection_direction_shift[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 3] = boundary_idx_0 - it.in_connection_direction_shift[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 4] = boundary_idx_0;
                    non_zero_entries.items[row_non_zero_entries_start_idx + 5] = boundary_idx_0 + it.in_connection_direction_shift[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 6] = boundary_idx_1 + it.first_internal_point_shift[1] + it.in_connection_direction_shift[1];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 7] = boundary_idx_1 + it.first_internal_point_shift[1];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 8] = boundary_idx_1 + it.first_internal_point_shift[1] - it.in_connection_direction_shift[1];

                    // check that we have all entries ascending as required by the row compressed format.
                    std.debug.assert(blk: {
                        for (0..8) |i| {
                            if (non_zero_entries.items[row_non_zero_entries_start_idx + i] > non_zero_entries.items[row_non_zero_entries_start_idx + i + 1]) break :blk true;
                        }
                        break :blk false;
                    });
                }

                // // enforce 2nd point to conform to 1st
                // {
                //     const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
                //     non_zero_entries.items[row_non_zero_entries_start_idx] = boundary_idx_0;
                //     non_zero_entries.items[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
                // }
            }

            // // TODO: remove this hardcoded handling
            // // we set point 0 to fixed and point 1 to connected to 0
            // {
            //     const connected_boundary_points = it.next().?;
            //
            //     const boundary_idx_0: c_int = @intCast(connected_boundary_points[0]);
            //     const boundary_idx_1: c_int = @intCast(connected_boundary_points[1]);
            //
            //     const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
            //     non_zero_entries.items[row_non_zero_entries_start_idx] = boundary_idx_0;
            //     non_zero_entries.items[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
            // }
        }
    }

    fn nonZeroMatrixEntries(
        self: *RowCompressedMatrixSystem2d,
        dof: usize,
        non_zero_entries_capacity: usize,
        connected_points: BlockBoundaryPointConnections,
    ) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_i[0..]));
        var non_zero_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

        try self.addOrAllocNonZeroMatrixEntries(&non_zero_entries, &row_count, connected_points);
        self.fillNonZeroMatrixEntriesForConnections(&non_zero_entries);
    }

    fn fillBlockInternalPointData(lhs: []f64, rhs_x: []f64, rhs_y: []f64, s: f64, t: f64, mesh_data: *const discrete.Mesh) void {
        var non_zero_entry_idx: usize = 0;
        var row_idx: usize = 0;

        for (mesh_data.blocks.items) |block| {
            var point_idx: usize = 0;

            // edge j_min
            for (0..block.points.size[1]) |_| {
                lhs[non_zero_entry_idx] = 1;
                rhs_x[row_idx] = block.points.data[point_idx].data[0];
                rhs_y[row_idx] = block.points.data[point_idx].data[1];

                non_zero_entry_idx += 1;
                point_idx += 1;
                row_idx += 1;
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                {
                    lhs[non_zero_entry_idx] = 1;
                    rhs_x[row_idx] = block.points.data[point_idx].data[0];
                    rhs_y[row_idx] = block.points.data[point_idx].data[1];

                    non_zero_entry_idx += 1;
                    point_idx += 1;
                    row_idx += 1;
                }

                // internal points with full 9 point stencil
                for (1..block.points.size[1] - 1) |_| {
                    const im1_j = block.points.data[point_idx - block.points.size[1]];
                    const i_jm1 = block.points.data[point_idx - 1];
                    const i_jp1 = block.points.data[point_idx + 1];
                    const ip1_j = block.points.data[point_idx + block.points.size[1]];

                    const stencil = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs[non_zero_entry_idx + 0] = stencil.get(.im1_jm1); // A[i-1, j-1]
                    lhs[non_zero_entry_idx + 1] = stencil.get(.im1_j); // A[i-1, j]
                    lhs[non_zero_entry_idx + 2] = stencil.get(.im1_jp1); // A[i-1, j+1]
                    lhs[non_zero_entry_idx + 3] = stencil.get(.i_jm1); // A[i, j-1]
                    lhs[non_zero_entry_idx + 4] = stencil.get(.i_j); // A[i, j]
                    lhs[non_zero_entry_idx + 5] = stencil.get(.i_jp1); // A[i, j+1]
                    lhs[non_zero_entry_idx + 6] = stencil.get(.ip1_jm1); // A[i+1, j-1]
                    lhs[non_zero_entry_idx + 7] = stencil.get(.ip1_j); // A[i+1, j]
                    lhs[non_zero_entry_idx + 8] = stencil.get(.ip1_jp1); // A[i+1, j+1]

                    rhs_x[row_idx] = 0;
                    rhs_y[row_idx] = 0;

                    non_zero_entry_idx += 9;
                    point_idx += 1;
                    row_idx += 1;
                }

                // edge i_max
                {
                    lhs[non_zero_entry_idx] = 1;
                    rhs_x[row_idx] = block.points.data[point_idx].data[0];
                    rhs_y[row_idx] = block.points.data[point_idx].data[1];

                    non_zero_entry_idx += 1;
                    point_idx += 1;
                    row_idx += 1;
                }
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                lhs[non_zero_entry_idx] = 1;
                rhs_x[row_idx] = block.points.data[point_idx].data[0];
                rhs_y[row_idx] = block.points.data[point_idx].data[1];

                non_zero_entry_idx += 1;
                point_idx += 1;
                row_idx += 1;
            }
        }
    }

    fn fillBlockConnectionData(self: @This(), lhs: []f64, rhs_x: []f64, rhs_y: []f64, s: f64, t: f64) void {
        for (self.mesh.connections.items) |connection| {
            const point_data = [2][]types.Vec2d{
                self.mesh.blocks.items[connection.ranges[0].block].points.data,
                self.mesh.blocks.items[connection.ranges[1].block].points.data,
            };

            var it = RangeFillMatrixIterator.init(connection, self.mesh);

            // TODO: check how to solve each point.
            // TODO: solve range 0 and connect range 1 to connected range 0 point
            // TODO: make sure that the corner points are conencted!

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

            // const matrix_idx_inc_intern_0: [3]usize = if (it.in_connection_direction_shift[0] > 0) .{ 0, 1, 2 } else .{ 2, 1, 0 };
            // const matrix_idx_inc_intern_1: [3]usize = if (it.in_connection_direction_shift[1] > 0) .{ 3, 4, 5 } else .{ 5, 4, 3 };

            while (it.next()) |connected_boundary_points| {
                std.debug.assert(types.eqlApprox(point_data[0][connected_boundary_points[0]], point_data[1][connected_boundary_points[1]], 1e-12));
                const boundary_idx_0: c_int = @intCast(connected_boundary_points[0]);
                const boundary_idx_1: c_int = @intCast(connected_boundary_points[1]);

                // we add the stencil data to the 1st point and force equality to it's solution for the 2nd point
                {
                    const row_idx: usize = @intCast(boundary_idx_0);

                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_0);

                    const im1_j = point_data[0][@intCast(boundary_idx_0 - it.in_connection_direction_shift[0])];
                    const i_jm1 = point_data[0][@intCast(boundary_idx_0 + it.first_internal_point_shift[0])];
                    const ip1_j = point_data[0][@intCast(boundary_idx_0 + it.in_connection_direction_shift[0])];
                    const i_jp1 = point_data[1][@intCast(boundary_idx_1 + it.first_internal_point_shift[1])];

                    const stencil = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    // add 9 point stencil data in the right order (associated points must be ascending)

                    // TODO: compute the offset dynamically (loop at matrix_idx_inc_intern as inspiration.); right now the indices are hard coded for the first connection.

                    lhs[row_non_zero_entries_start_idx + 2] = stencil.get(.im1_jm1); // A[matrix_idx, matrix_idx_0 - matrix_increment_0]
                    lhs[row_non_zero_entries_start_idx + 1] = stencil.get(.i_jm1); // A[matrix_idx, matrix_idx_0]
                    lhs[row_non_zero_entries_start_idx + 0] = stencil.get(.ip1_jm1); // A[matrix_idx, matrix_idx_0 + matrix_increment_0]

                    lhs[row_non_zero_entries_start_idx + 5] = stencil.get(.im1_j); // A[matrix_idx, matrix_idx - 1]
                    lhs[row_non_zero_entries_start_idx + 4] = stencil.get(.i_j); // A[matrix_idx, matrix_idx]
                    lhs[row_non_zero_entries_start_idx + 3] = stencil.get(.ip1_j); // A[matrix_idx, matrix_idx + 1]

                    lhs[row_non_zero_entries_start_idx + 8] = stencil.get(.im1_jp1); // A[matrix_idx, matrix_idx_1 - matrix_increment_1]
                    lhs[row_non_zero_entries_start_idx + 7] = stencil.get(.i_jp1); // A[matrix_idx, matrix_idx_1]
                    lhs[row_non_zero_entries_start_idx + 6] = stencil.get(.ip1_jp1); // A[matrix_idx, matrix_idx_1 + matrix_increment_1]

                    rhs_x[row_idx] = 0;
                    rhs_y[row_idx] = 0;
                }

                // enforce equality of the connected 2nd point with the 1st
                // NOTE: it might be more efficient to do this in a 2nd loop to reduce cache misses.
                {
                    const row_idx: usize = @intCast(boundary_idx_1);
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);

                    lhs[row_non_zero_entries_start_idx] = 1.0;
                    lhs[row_non_zero_entries_start_idx + 1] = -1.0;
                    rhs_x[row_idx] = 0;
                    rhs_y[row_idx] = 0;
                }
            }
        }
    }

    fn fillAndSolve(self: *RowCompressedMatrixSystem2d) !void {
        var lhs = self.lhs_values;

        var rhs_x = self.rhs_x;
        var rhs_y = self.rhs_y;

        // laplace conditions (zero intialization)
        const s = 0.0;
        const t = 0.0;

        // TODO: adjust naming: where are the fixed points set!? The second function handles only connected boundary points !?
        fillBlockInternalPointData(lhs[0..], rhs_x[0..], rhs_y[0..], s, t, self.mesh);
        self.fillBlockConnectionData(lhs[0..], rhs_x[0..], rhs_y[0..], s, t);

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, lhs, rhs_x, self.x_new[0..], rhs_y, self.y_new[0..]);
    }
};

const BlockBoundaryPointKind = enum {
    fix, // not smoothed

    // interface points are either smoothed like interior points or connected to a smoothed point
    smooth,
    connect,

    junction, // smoothed using Laplacian smoothing
};

/// contains for each block boundary point all connected points in a flat array. The number of connections allows
/// to categorize each boundary point w.r.t. how to smooth it.
const BlockBoundaryPointConnections = struct {
    // NOTE: we use a c_int here only to match what is needed by the matrix lib; usize would be the natural choice.
    data: boundary.PointData(std.BoundedArray(c_int, 4)),

    fn init(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh, index_converter: IndexConverter) !BlockBoundaryPointConnections {
        var connected_points = BlockBoundaryPointConnections{
            .data = try .init(allocator, mesh_data),
        };
        errdefer connected_points.deinit();

        for (connected_points.data.buffer[0..]) |*boundary_point| boundary_point.* = try .init(0);

        for (mesh_data.connections.items) |connection| {
            var it = connection.iterate(mesh_data);
            while (it.next()) |connected_points_local_indices| {
                // TODO: is there a nicer way to doing this?

                // local to global point indices
                const global_idx_0 = index_converter.globalIndex(connection.ranges[0].block, connected_points_local_indices[0]);
                const global_idx_1 = index_converter.globalIndex(connection.ranges[1].block, connected_points_local_indices[1]);

                const local_idx_0 = IndexConverter.index(connection.ranges[0].block, connected_points_local_indices[0], mesh_data);
                const local_idx_1 = IndexConverter.index(connection.ranges[1].block, connected_points_local_indices[1], mesh_data);

                const buffer_idx_0 = try connected_points.data.bufferIndex(local_idx_0, mesh_data.blocks.items[connection.ranges[0].block].points.size);

                // const buffer_idx_0 = try connected_points.data.bufferIndex(.{ .block = connection.ranges[0].block, .idx = local_idx_0 }, mesh_data.blocks.items[connection.ranges[0].block].points.size);
                const buffer_idx_1 = try connected_points.data.bufferIndex(.{ .block = connection.ranges[1].block, .idx = local_idx_1 }, mesh_data.blocks.items[connection.ranges[1].block].points.size);

                connected_points.data.buffer[buffer_idx_0].append(@intCast(global_idx_1));
                connected_points.data.buffer[buffer_idx_1].append(@intCast(global_idx_0));
            }
        }

        // TODO: either check that all boundary points are defined (would work with a tagged union) or handle end points of connections seperately (if
        // these are not junction points, set the 1st to fixed and connect the 2nd).

        // TODO: remove this hard coding of the end points of the 1st connection.
        {
            const size = mesh_data.blocks.items[4].points.size;
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ 4, .{ 0, 0 } }, size)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ 4, .{ 0, 20 } }, size)].clear();
        }

        return connected_points;
    }

    fn deinit(self: BlockBoundaryPointConnections) void {
        self.data.deinit();
    }

    /// returns the kind (i.e. how to treat the point w.r.t. smoothing) for the given boundary point index.
    fn getKind(self: BlockBoundaryPointConnections, bounday_point_index: usize, global_point_index: usize) BlockBoundaryPointKind {
        // const point_data = self.connected_points.get(index);
        const point_data = self.data.buffer[bounday_point_index];

        return switch (point_data.len) {
            0 => .fix,
            1 => if (point_data.get(0) > global_point_index) .smooth else .connect, // smooth lower index point and connect higher index point
            else => .junction,
        };
    }
};

const RangeFillMatrixIterator = struct {
    count: usize,
    first_internal_point_shift: [2]c_int,
    in_connection_direction_shift: [2]c_int,
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
            const block_idx = connection.ranges[side_idx].block;
            const points = mesh_data.blocks.items[block_idx].points;
            const start = connection.ranges[side_idx].start;
            const end = connection.ranges[side_idx].end;
            switch (connection.ranges[side_idx].side) {
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
                    data.first_internal_point_shift[side_idx] = -@as(c_int, @intCast(points.size[1]));
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

const GlobalIndex = usize;
const LocalIndex = usize;

/// Convert index types necessary for smoothing.
const IndexConverter = struct {
    allocator: std.mem.Allocator,
    global_point_index_range_start: []GlobalIndex,

    fn init(mesh_data: *const discrete.Mesh, allocator: std.mem.Allocator) !IndexConverter {
        var global_point_index_range_start = try allocator.alloc(GlobalIndex, mesh_data.blocks.items.len);

        var total_num_points: usize = 0;
        for (mesh_data.blocks.items, 0..) |block, block_idx| {
            global_point_index_range_start[block_idx] = total_num_points;
            total_num_points += block.points.size[0] * block.points.size[1];
        }

        return .{
            .allocator = allocator,
            .global_point_index_range_start = global_point_index_range_start,
        };
    }

    fn deinit(self: IndexConverter) void {
        self.allocator.free(self.global_point_index_range_start);
    }

    fn globalIndex(self: IndexConverter, block: usize, local_idx: LocalIndex) GlobalIndex {
        return self.global_point_index_range_start[block] + local_idx;
    }

    fn localIndex(self: IndexConverter, global_idx: GlobalIndex) struct {
        block: usize,
        local_idx: LocalIndex,
    } {
        var block_idx: usize = self.blocks.len - 1;
        while (global_idx < self.global_point_index_range_start[block_idx]) : (block_idx -= 1) {}
        const local_idx = global_idx - self.global_point_index_range_start[block_idx];
        return .{ .block = block_idx, .local_idx = local_idx };
    }

    fn index(block_idx: usize, local_idx: LocalIndex, mesh_data: *const discrete.Mesh) struct {
        block: usize,
        idx: types.Index2d,
    } {
        const block = mesh_data.blocks.items[block_idx];
        const point_idx = blk: {
            const point_i = @divTrunc(local_idx, block.size[1]);
            const point_j = local_idx - point_i * block.size[1];
            break :blk types.Index2d{ point_i, point_j };
        };
        return .{
            .block = block,
            .idx = point_idx,
        };
    }
};
