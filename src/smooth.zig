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
                for (0..block.points.size[0]) |i| {
                    for (0..block.points.size[1]) |j| {
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
                for (0..block.points.size[0]) |i| {
                    for (0..block.points.size[1]) |j| {
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

        if (connection.periodicity) |periodicity| {
            var point_idx: usize = 0;
            while (true) {
                const p_0 = it_0.next() orelse {
                    std.debug.assert(it_1.next() == null);
                    break;
                };
                const p_1 = it_1.next().?;

                const x_0 = types.add(mesh_data.blocks.items[range_0.block].points.data[p_0], periodicity);
                const x_1 = mesh_data.blocks.items[range_1.block].points.data[p_1];

                if (!types.eqlApprox(x_0, x_1, abs_tol)) {
                    std.debug.panic("non matching points for connection {} point {}:\n\t{}\n\t{}\n", .{ connection_idx, point_idx, x_0, x_1 });
                }

                point_idx += 1;
            }
        } else {
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

    // TODO: replace this with a boundary buffer of just the kind!
    connected_points: BlockBoundaryPointConnections,

    fn init(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh) !RowCompressedMatrixSystem2d {
        connectionDataCheck(mesh_data);

        var dof: usize = 0;
        var non_zero_entries_capacity: usize = 0;

        var row_idx_range_start_for_each_block = try allocator.alloc(usize, mesh_data.blocks.items.len);
        errdefer allocator.free(row_idx_range_start_for_each_block);

        for (mesh_data.blocks.items, 0..) |block, block_idx| {
            row_idx_range_start_for_each_block[block_idx] = dof;

            const points = block.points;
            const block_dof = points.size[0] * points.size[1];
            dof += block_dof;

            // this is a max limit of non zero entries per point
            // the true number will be less dependent on how many
            // neighboring points need to bes solved
            const block_non_zero_entries = dof * 9;
            non_zero_entries_capacity += block_non_zero_entries;
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

        // Handling of the inter-block connections:
        // For now, we will add rows for each connected point to enforce identity of the connected points.
        // Each connected point adds his part of stencil data to the system matrix.
        // This renders the stencil info most easy (of the considered option) at the expense of additional
        // matrix entries.
        const index_converter = try IndexConverter.init(mesh_data, allocator);
        defer index_converter.deinit();

        const connected_points = try BlockBoundaryPointConnections.init(allocator, mesh_data, index_converter);

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
            .connected_points = connected_points,
        };

        try system.nonZeroMatrixEntries(dof, non_zero_entries_capacity, connected_points, index_converter);

        return system;
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        self.allocator.free(self.buffer_int);
        self.allocator.free(self.buffer_float);
        self.allocator.free(self.row_idx_range_start_for_each_block);
        self.connected_points.deinit();
    }

    /// returns the non zero entries range for the given row (the row index range start for each block can be retrieved from an array in the struct)
    fn nonZeroEntriesRangeStart(self: RowCompressedMatrixSystem2d, row_idx: c_int) usize {
        return @intCast(self.lhs_p[@intCast(row_idx)]);
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
                // NOTE: the ordering must be ascending. Since we smooth the smaller global index, it must come first.
                try non_zero_entries.append(connected_points.data.buffer[boundary_point_idx.*].get(0));
                try non_zero_entries.append(row_idx.*);
            },
            .junction => {
                // TODO: remove hard coding.
                std.debug.assert(row_idx.* == 26221);
                _ = try non_zero_entries.appendSlice(&[_]c_int{
                    26104,
                    26105,
                    26220,
                    26221,
                    31305,
                    31326,
                    43119,
                    43120,
                    43121,
                });
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
        index_converter: IndexConverter,
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

            // it.first_internal_point_shift[0] > 0

            // NOTE: this is needed to add the point indices in ascending order depending on the direction of the connection.
            const direction_modifier: [2]c_int = .{
                if (it.in_connection_direction_shift[0] > 0) 1.0 else -1.0,
                if (it.in_connection_direction_shift[1] > 0) 1.0 else -1.0,
            };

            // NOTE: this is needed to have the right ordering of the points taken from the 0th block of the connection.
            const shift: [2]c_int = if (it.first_internal_point_shift[0] > 0) .{ 0, it.first_internal_point_shift[0] } else .{ it.first_internal_point_shift[0], 0 };

            // TODO: remove this hardcoded handling
            // we set point 0 to fixed and point 1 to connected to 0
            _ = it.next().?;

            // smooth 1st point
            for (0..it.count) |_| {
                const connected_points_local_idx = it.next().?;
                const connected_points_global_idx = [2]c_int{
                    @intCast(index_converter.globalIndex(connection.ranges[0].block, connected_points_local_idx[0]).value),
                    @intCast(index_converter.globalIndex(connection.ranges[1].block, connected_points_local_idx[1]).value),
                };

                const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(connected_points_global_idx[0]);

                // NOTE: here we set the data based on the connections; the allocation happend in the previous function based on a point loop.

                // TODO: do this w/o switch; e.g., via the indices set once for a connection.
                switch (connection.ranges[0].side) {
                    .i_min, .i_max => {
                        non_zero_entries.items[row_non_zero_entries_start_idx + 0] = connected_points_global_idx[0] + shift[0] - it.in_connection_direction_shift[0] * direction_modifier[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 2] = connected_points_global_idx[0] + shift[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 4] = connected_points_global_idx[0] + shift[0] + it.in_connection_direction_shift[0] * direction_modifier[0];

                        non_zero_entries.items[row_non_zero_entries_start_idx + 1] = connected_points_global_idx[0] + shift[1] - it.in_connection_direction_shift[0] * direction_modifier[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 3] = connected_points_global_idx[0] + shift[1];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 5] = connected_points_global_idx[0] + shift[1] + it.in_connection_direction_shift[0] * direction_modifier[0];
                    },
                    .j_min, .j_max => {
                        non_zero_entries.items[row_non_zero_entries_start_idx + 0] = connected_points_global_idx[0] + shift[0] - it.in_connection_direction_shift[0] * direction_modifier[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 1] = connected_points_global_idx[0] + shift[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 2] = connected_points_global_idx[0] + shift[0] + it.in_connection_direction_shift[0] * direction_modifier[0];

                        non_zero_entries.items[row_non_zero_entries_start_idx + 3] = connected_points_global_idx[0] + shift[1] - it.in_connection_direction_shift[0] * direction_modifier[0];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 4] = connected_points_global_idx[0] + shift[1];
                        non_zero_entries.items[row_non_zero_entries_start_idx + 5] = connected_points_global_idx[0] + shift[1] + it.in_connection_direction_shift[0] * direction_modifier[0];
                    },
                }

                // NOTE: the entries into the 2nd block always come last since the block index is higher and thus the global indices are higher.
                non_zero_entries.items[row_non_zero_entries_start_idx + 6] = connected_points_global_idx[1] + it.first_internal_point_shift[1] - it.in_connection_direction_shift[1] * direction_modifier[1];
                non_zero_entries.items[row_non_zero_entries_start_idx + 7] = connected_points_global_idx[1] + it.first_internal_point_shift[1];
                non_zero_entries.items[row_non_zero_entries_start_idx + 8] = connected_points_global_idx[1] + it.first_internal_point_shift[1] + it.in_connection_direction_shift[1] * direction_modifier[1];

                // check that we have all entries ascending as required by the row compressed format.
                std.debug.assert(blk: {
                    for (0..8) |i| {
                        if (non_zero_entries.items[row_non_zero_entries_start_idx + i] >= non_zero_entries.items[row_non_zero_entries_start_idx + i + 1]) {
                            std.debug.print("wrong ordering of connection stencil data for point {d}: {any}\n", .{ connected_points_global_idx[0], non_zero_entries.items[row_non_zero_entries_start_idx .. row_non_zero_entries_start_idx + 9] });
                            break :blk true;
                        }
                    }
                    break :blk true;
                });
            }
        }
    }

    fn nonZeroMatrixEntries(
        self: *RowCompressedMatrixSystem2d,
        dof: usize,
        non_zero_entries_capacity: usize,
        connected_points: BlockBoundaryPointConnections,
        index_converter: IndexConverter,
    ) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_i[0..]));
        var non_zero_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

        try self.addOrAllocNonZeroMatrixEntries(&non_zero_entries, &row_count, connected_points);
        self.fillNonZeroMatrixEntriesForConnections(&non_zero_entries, index_converter);
    }

    fn fillBlockBoundaryPointData(
        self: @This(),
        lhs: []f64,
        rhs_x: []f64,
        rhs_y: []f64,
        block: discrete.Block2d,
        boundary_point_idx: *usize,
        row_idx: *usize,
        point_idx: *usize,
        non_zero_entry_idx: *usize,
    ) void {
        // TODO: remove setting the RHS in this loop. It only needs to be set once before looping.

        switch (self.connected_points.getKind(boundary_point_idx.*, row_idx.*)) {
            .fix => {
                // TODO: we can also set RHS to 1 and only set the LHS.
                lhs[non_zero_entry_idx.*] = 1;
                non_zero_entry_idx.* += 1;

                rhs_x[row_idx.*] = block.points.data[point_idx.*].data[0];
                rhs_y[row_idx.*] = block.points.data[point_idx.*].data[1];
            },
            .smooth => {
                // NOTE: here, we only skip the necessary stencil data since the values are set next in a connection based loop.
                non_zero_entry_idx.* += 9;

                rhs_x[row_idx.*] = 0;
                rhs_y[row_idx.*] = 0;
            },
            .connect => {
                // NOTE: we would only need to fill this data one time.
                // TODO: do this in a one time loop up front.
                lhs[non_zero_entry_idx.*] = 1;
                lhs[non_zero_entry_idx.* + 1] = -1;
                non_zero_entry_idx.* += 2;

                // TODO: this needs to be set depending on weather it is a simple or periodic connection. Right now it only handles
                // simple connections.
                rhs_x[row_idx.*] = 0;
                rhs_y[row_idx.*] = 0;
            },
            .junction => {
                lhs[non_zero_entry_idx.*] = 1; // 26104
                lhs[non_zero_entry_idx.* + 1] = 1; // 26105
                lhs[non_zero_entry_idx.* + 2] = 1; // 26220
                lhs[non_zero_entry_idx.* + 3] = -8; // 26221
                lhs[non_zero_entry_idx.* + 4] = 1; // 31305
                lhs[non_zero_entry_idx.* + 5] = 1; // 31326
                lhs[non_zero_entry_idx.* + 6] = 1; // 43119
                lhs[non_zero_entry_idx.* + 7] = 1; // 43120
                lhs[non_zero_entry_idx.* + 8] = 1; // 43121
                non_zero_entry_idx.* += 9;

                rhs_x[row_idx.*] = 0;
                rhs_y[row_idx.*] = 0;
            },
        }

        row_idx.* += 1;
        boundary_point_idx.* += 1;
        point_idx.* += 1;
    }

    fn fillBlockInternalPointData(
        self: @This(),
        lhs: []f64,
        rhs_x: []f64,
        rhs_y: []f64,
        s: f64,
        t: f64,
    ) void {
        var non_zero_entry_idx: usize = 0;
        var row_idx: usize = 0;
        var boundary_point_idx: usize = 0;

        for (self.mesh.blocks.items) |block| {
            var point_idx: usize = 0;

            // edge j_min
            for (0..block.points.size[1]) |_| {
                self.fillBlockBoundaryPointData(lhs, rhs_x, rhs_y, block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx);
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                self.fillBlockBoundaryPointData(lhs, rhs_x, rhs_y, block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx);

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
                self.fillBlockBoundaryPointData(lhs, rhs_x, rhs_y, block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx);
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                self.fillBlockBoundaryPointData(lhs, rhs_x, rhs_y, block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx);
            }
        }
    }

    fn fillBlockConnectionData(
        self: RowCompressedMatrixSystem2d,
        lhs: []f64,
        s: f64,
        t: f64,
    ) void {
        for (self.mesh.connections.items) |connection| {
            const point_data = [2][]types.Vec2d{
                self.mesh.blocks.items[connection.ranges[0].block].points.data,
                self.mesh.blocks.items[connection.ranges[1].block].points.data,
            };

            var it = RangeFillMatrixIterator.init(connection, self.mesh);

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

            // NOTE: this is needed to add the point indices in ascending order depending on the direction of the connection (see the non zero matrix entries definition).
            const direction_modifier: [2]isize = .{
                if (it.in_connection_direction_shift[0] > 0) 1.0 else -1.0,
                if (it.in_connection_direction_shift[1] > 0) 1.0 else -1.0,
            };

            // NOTE: this is needed to have the right ordering of the points taken from the 0th block of the connection (entries 0 - 5):
            // the ordering is first the points inside block 0 ([i-1,j-1],[i,j-1],[i+1,j-2]) and then the points on the connection
            // ([i-1,j],[i,j],[i+1,j]).
            const point_stencil_idx: [6]usize = switch (connection.ranges[0].side) {
                .i_min => .{
                    @intCast(3 - 2 * direction_modifier[0]),
                    3,
                    @intCast(3 + 2 * direction_modifier[0]),
                    @intCast(2 - 2 * direction_modifier[0]),
                    2,
                    @intCast(2 + 2 * direction_modifier[0]),
                },
                .i_max => .{
                    @intCast(2 - 2 * direction_modifier[0]),
                    2,
                    @intCast(2 + 2 * direction_modifier[0]),
                    @intCast(3 - 2 * direction_modifier[0]),
                    3,
                    @intCast(3 + 2 * direction_modifier[0]),
                },
                .j_min => .{
                    @intCast(4 - direction_modifier[0]),
                    4,
                    @intCast(4 + direction_modifier[0]),
                    @intCast(1 - direction_modifier[0]),
                    1,
                    @intCast(1 + direction_modifier[0]),
                },
                .j_max => .{
                    @intCast(1 - direction_modifier[0]),
                    1,
                    @intCast(1 + direction_modifier[0]),
                    @intCast(4 - direction_modifier[0]),
                    4,
                    @intCast(4 + direction_modifier[0]),
                },
            };

            // TODO: remove this hard coding for the first point
            _ = it.next();

            for (0..it.count) |_| {
                const connected_points = it.next().?;
                const point_idx: [2]c_int = .{ @intCast(connected_points[0].value), @intCast(connected_points[1].value) };

                // TODO: this could be enhanced by computing this just once for the connection (plus +1 or -1 for the next connection point)
                const global_idx_0 = @as(c_int, @intCast(self.row_idx_range_start_for_each_block[connection.ranges[0].block])) + point_idx[0];
                const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(global_idx_0);

                const im1_j = point_data[0][@intCast(point_idx[0] - it.in_connection_direction_shift[0])];
                const i_jm1 = point_data[0][@intCast(point_idx[0] + it.first_internal_point_shift[0])];
                const ip1_j = point_data[0][@intCast(point_idx[0] + it.in_connection_direction_shift[0])];
                const i_jp1 = point_data[1][@intCast(point_idx[1] + it.first_internal_point_shift[1])];

                const stencil = StencilData.init(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                // points inside block 0
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[0]] = stencil.get(.im1_jm1);
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[1]] = stencil.get(.i_jm1);
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[2]] = stencil.get(.ip1_jm1);

                // points on boundary connection (taken from block 0)
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[3]] = stencil.get(.im1_j);
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[4]] = stencil.get(.i_j);
                lhs[row_non_zero_entries_start_idx + point_stencil_idx[5]] = stencil.get(.ip1_j);

                // points inside block 1
                lhs[row_non_zero_entries_start_idx + @as(usize, @intCast(7 - direction_modifier[1]))] = stencil.get(.im1_jp1);
                lhs[row_non_zero_entries_start_idx + 7] = stencil.get(.i_jp1);
                lhs[row_non_zero_entries_start_idx + @as(usize, @intCast(7 + direction_modifier[1]))] = stencil.get(.ip1_jp1);

                // NOTE: RHS is already set in the point based loop.
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
        self.fillBlockInternalPointData(lhs[0..], rhs_x[0..], rhs_y[0..], s, t);
        self.fillBlockConnectionData(lhs[0..], s, t);

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

// TODO: remove above kind enum and rename this type.
const Tags = enum {
    fixed,

    connection_end,

    smoothed,

    connected,

    periodic,

    laplacian,
};

/// contains for each block boundary point all connected points in a flat array. The number of connections allows
/// to categorize each boundary point w.r.t. how to smooth it.
const BlockBoundaryPointConnections = struct {
    // NOTE: we use a c_int here only to match what is needed by the matrix lib; usize would be the natural choice.
    // TODO: rename to connected points.
    data: boundary.PointData(std.BoundedArray(c_int, 4)),

    // TODO: we need to somehow tag the different points
    tags: boundary.PointData(Tags),

    fn init(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh, index_converter: IndexConverter) !BlockBoundaryPointConnections {
        var connected_points = BlockBoundaryPointConnections{
            .data = try .init(allocator, mesh_data),
            .tags = try .init(allocator, mesh_data),
        };
        errdefer connected_points.deinit();

        for (connected_points.data.buffer[0..]) |*boundary_point| boundary_point.* = try .init(0);

        for (mesh_data.connections.items) |connection| {
            var it = connection.iterate(mesh_data);
            while (it.next()) |connected_points_local_indices| {
                // TODO: is there a nicer way to doing this?

                // local to global point indices
                const global_idx_0 = index_converter.globalIndex(connection.ranges[0].block, .{ .value = connected_points_local_indices[0] });
                const global_idx_1 = index_converter.globalIndex(connection.ranges[1].block, .{ .value = connected_points_local_indices[1] });

                const local_idx_0 = IndexConverter.index(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = connected_points_local_indices[0] } }, mesh_data);
                const local_idx_1 = IndexConverter.index(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = connected_points_local_indices[1] } }, mesh_data);

                const buffer_idx_0 = try connected_points.data.bufferIndex(local_idx_0, mesh_data);
                const buffer_idx_1 = try connected_points.data.bufferIndex(local_idx_1, mesh_data);

                try connected_points.data.buffer[buffer_idx_0].append(@intCast(global_idx_1.value));
                try connected_points.data.buffer[buffer_idx_1].append(@intCast(global_idx_0.value));
            }
        }

        for (mesh_data.connections.items) |connection| {
            var it = connection.iterate(mesh_data);

            // first end point
            {
                const local_flat_indices = it.next();

                // TODO: move all of this index transformation into a seperate function!

                const global_idx_0 = index_converter.globalIndex(connection.ranges[0].block, .{ .value = connected_points_local_indices[0] });
                const global_idx_1 = index_converter.globalIndex(connection.ranges[1].block, .{ .value = connected_points_local_indices[1] });

                const local_idx_0 = IndexConverter.index(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = connected_points_local_indices[0] } }, mesh_data);
                const local_idx_1 = IndexConverter.index(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = connected_points_local_indices[1] } }, mesh_data);

                const buffer_idx_0 = try connected_points.data.bufferIndex(local_idx_0, mesh_data);
                const buffer_idx_1 = try connected_points.data.bufferIndex(local_idx_1, mesh_data);

                try connected_points.tags.buffer[buffer_idx_0] = .connection_end;
                try connected_points.tags.buffer[buffer_idx_1] = .connection_end;
            }

            // connection internal points
            const internal_point_tag = if (connection.periodicity) |_| .periodic else .connected;
            for (0..it.count) |_| {
                const local_flat_indices = it.next();

                const global_idx_0 = index_converter.globalIndex(connection.ranges[0].block, .{ .value = connected_points_local_indices[0] });
                const global_idx_1 = index_converter.globalIndex(connection.ranges[1].block, .{ .value = connected_points_local_indices[1] });

                const local_idx_0 = IndexConverter.index(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = connected_points_local_indices[0] } }, mesh_data);
                const local_idx_1 = IndexConverter.index(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = connected_points_local_indices[1] } }, mesh_data);

                const buffer_idx_0 = try connected_points.data.bufferIndex(local_idx_0, mesh_data);
                const buffer_idx_1 = try connected_points.data.bufferIndex(local_idx_1, mesh_data);

                try connected_points.tags.buffer[buffer_idx_0] = .smoothed;
                try connected_points.tags.buffer[buffer_idx_1] = internal_point_tag;
            }

            // second end point
            {
                const local_flat_indices = it.next();

                const global_idx_0 = index_converter.globalIndex(connection.ranges[0].block, .{ .value = connected_points_local_indices[0] });
                const global_idx_1 = index_converter.globalIndex(connection.ranges[1].block, .{ .value = connected_points_local_indices[1] });

                const local_idx_0 = IndexConverter.index(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = connected_points_local_indices[0] } }, mesh_data);
                const local_idx_1 = IndexConverter.index(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = connected_points_local_indices[1] } }, mesh_data);

                const buffer_idx_0 = try connected_points.data.bufferIndex(local_idx_0, mesh_data);
                const buffer_idx_1 = try connected_points.data.bufferIndex(local_idx_1, mesh_data);

                try connected_points.tags.buffer[buffer_idx_0] = .connection_end;
                try connected_points.tags.buffer[buffer_idx_1] = .connection_end;
            }
        }

        // TODO: collect all connecting points
        // TODO: remove duplicats.

        // TODO: either check that all boundary points are defined (would work with a tagged union) or handle end points of connections seperately (if
        // these are not junction points, set the 1st to fixed and connect the 2nd).

        // TODO: remove this hard coding of the end points of the 1st connection.
        {
            // try connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 4, .point = .{ 0, 0 } }, mesh_data)].append(26221);

            // connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 4, .point = .{ 0, 0 } }, size)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 4, .point = .{ 0, 20 } }, mesh_data)].clear();
            // connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 4, .point = .{ 115, 0 } }, size)].clear();
            // connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 2, .point = .{ 53, 115 } }, mesh_data)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 2, .point = .{ 53, 0 } }, mesh_data)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 2, .point = .{ 0, 115 } }, mesh_data)].clear();

            // TODO: adjust points for periodic connection!
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 5, .point = .{ 0, 23 } }, mesh_data)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 5, .point = .{ 194, 23 } }, mesh_data)].clear();
            connected_points.data.buffer[try connected_points.data.bufferIndex(.{ .block = 4, .point = .{ 194, 20 } }, mesh_data)].clear();
        }

        // sort junction point data in ascending order
        for (connected_points.data.buffer, 0..) |*point_data, idx| {
            if (point_data.len > 1) {
                // add yourself to have a consistent ordering for all points
                const global_idx = connected_points.data.pointIndex(idx, mesh_data);
                try point_data.append(@intCast(global_idx));

                std.mem.sort(c_int, point_data.slice(), {}, std.sort.asc(c_int));
            }
        }

        return connected_points;
    }

    fn deinit(self: BlockBoundaryPointConnections) void {
        self.data.deinit();
    }

    /// returns the kind (i.e. how to treat the point w.r.t. smoothing) for the given boundary point index.
    fn getKind(self: BlockBoundaryPointConnections, bounday_point_index: usize, global_point_index: usize) BlockBoundaryPointKind {
        const point_data = self.data.buffer[bounday_point_index];

        return switch (point_data.len) {
            0 => .fix,
            // TODO: remove this test out of here: return a struct that has a func to test this.
            1 => if (point_data.get(0) > global_point_index) .smooth else .connect, // smooth lower index point and connect higher index point
            else => {
                for (point_data.slice()[1..]) |p| std.debug.assert(point_data.get(0) < p);
                return if (point_data.get(0) == global_point_index) .junction else .connect; // treat first point as junction and the others as connected to the first point
            },
        };
    }
};

const RangeFillMatrixIterator = struct {
    count: usize,
    first_internal_point_shift: [2]c_int,
    in_connection_direction_shift: [2]c_int,
    position: [2]LocalIndex,

    fn next(self: *RangeFillMatrixIterator) ?[2]LocalIndex {
        if (self.count == 0) return null;
        const position = self.position;

        self.count -= 1;

        self.position = .{
            .{ .value = @intCast(@as(isize, @intCast(self.position[0].value)) + self.in_connection_direction_shift[0]) },
            .{ .value = @intCast(@as(isize, @intCast(self.position[1].value)) + self.in_connection_direction_shift[1]) },
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
                    data.position[side_idx] = .{ .value = points.index(.{ start, 0 }) };
                },
                .i_max => {
                    data.first_internal_point_shift[side_idx] = -1;
                    data.in_connection_direction_shift[side_idx] = @intCast(points.size[1]);
                    data.position[side_idx] = .{ .value = points.index(.{ start, points.size[1] - 1 }) };
                },
                .j_min => {
                    data.first_internal_point_shift[side_idx] = @intCast(points.size[1]);
                    data.in_connection_direction_shift[side_idx] = 1;
                    data.position[side_idx] = .{ .value = points.index(.{ 0, start }) };
                },
                .j_max => {
                    data.first_internal_point_shift[side_idx] = -@as(c_int, @intCast(points.size[1]));
                    data.in_connection_direction_shift[side_idx] = 1;
                    data.position[side_idx] = .{ .value = points.index(.{ points.size[0] - 1, start }) };
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

// NOTE: the different index types are introduce to help express the required index type in a type save manner.
// I do not know if this is wearth it TBH, but perhaps it will be helpful in the future to not only rely in sane
// variable naming...
// When it comes to implementation, it seems resonable to do it this way as the tuple type index notatin can leed
// to confusion with the array unpacking indexing.

// DOCS: explain the different types of indices with examples. Explain why type safety is added here like this.

/// A (flat) global unique point identifier.
const GlobalIndex = struct { value: usize };

/// A (flat) block local unique point identifier: index into a block array.
const LocalIndex = struct { value: usize };

const BlockAndLocalIndex = struct { block: usize, local_idx: LocalIndex };

/// Convert index types necessary for smoothing.
const IndexConverter = struct {
    allocator: std.mem.Allocator,
    global_point_index_range_start: []GlobalIndex,

    fn init(mesh_data: *const discrete.Mesh, allocator: std.mem.Allocator) !IndexConverter {
        var global_point_index_range_start = try allocator.alloc(GlobalIndex, mesh_data.blocks.items.len);

        var total_num_points: usize = 0;
        for (mesh_data.blocks.items, 0..) |block, block_idx| {
            global_point_index_range_start[block_idx] = .{ .value = total_num_points };
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
        return .{ .value = self.global_point_index_range_start[block].value + local_idx.value };
    }

    fn localIndex(self: IndexConverter, global_idx: GlobalIndex) BlockAndLocalIndex {
        var block_idx: usize = self.global_point_index_range_start.len - 1;
        while (global_idx.value < self.global_point_index_range_start[block_idx].value) : (block_idx -= 1) {}
        const local_idx = global_idx.value - self.global_point_index_range_start[block_idx].value;
        return .{ .block = block_idx, .local_idx = .{ .value = local_idx } };
    }

    fn index(block_and_local_idx: BlockAndLocalIndex, mesh_data: *const discrete.Mesh) types.MeshIndex2d {
        const block = mesh_data.blocks.items[block_and_local_idx.block];
        const point_idx = blk: {
            const point_i = @divTrunc(block_and_local_idx.local_idx.value, block.points.size[1]);
            const point_j = block_and_local_idx.local_idx.value - point_i * block.points.size[1];
            break :blk types.Index2d{ point_i, point_j };
        };
        return .{ .block = block_and_local_idx.block, .point = point_idx };
    }
};
