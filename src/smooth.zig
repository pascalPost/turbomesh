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

        for (mesh_data.blocks.items, 0..1) |block, block_idx| {
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

        // Handling of the inter-block connections:
        // For now, we will add rows for each connected point to enforce identity of the connected points.
        // Each connected point adds his part of stencil data to the system matrix.
        // This renders the stencil info most easy (of the considered option) at the expense of additional
        // matrix entries.

        const block_boundary_points_kind = categorizeBlockBoundaryPoints(allocator, mesh_data);
        defer block_boundary_points_kind.deinit();

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

        try system.nonZeroMatrixEntries(dof, non_zero_entries_capacity);

        return system;
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        defer self.allocator.free(self.buffer_int);
        defer self.allocator.free(self.buffer_float);
        defer self.allocator.free(self.block_point_index_range_start);
    }

    /// returns the non zero entries range for the given row (the row index range start for each block can be retrieved from an array in the struct)
    fn nonZeroEntriesRangeStart(self: RowCompressedMatrixSystem2d, row: usize) usize {
        return self.lhs_p[row + 1];
    }

    fn nonZeroMatrixEntries(self: *RowCompressedMatrixSystem2d, dof: usize, non_zero_entries_capacity: usize) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_i[0..]));
        var non_zero_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

        // TODO: move this to own function
        var row_idx: c_int = 0;
        for (self.mesh.blocks.items, 0..) |block, block_idx| {
            const col_size: c_int = @intCast(block.points.size[1]);

            //
            //
            // TODO: this info needs to go before to be able to compute the DOF

            // for all points on the boundary we need to check how this point is to be added to the system matrix:
            // - fixed point: only a single entry with the current coordinates on the RHS
            // - smoothed: 9 entries for the full stencil
            // - laplacian: TBD (probably also 9 entries) TODO: fill this info.

            _ = block_idx;
            var block_boundary_point_idx: usize = 0;
            const boundary_points_kind: []BlockBoundaryPointKind = undefined;

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
                // NOTE: here we assume a fixed point...
                // The actual number of non zero entries for this point's row depends on weather it will be smoothed or not.
                // That means for each boundary point we need to check the point status.

                // We need to access the block boundary properties
                const block_boundary_point_kind = boundary_points_kind[block_boundary_point_idx];
                block_boundary_point_idx += 1;

                switch (block_boundary_point_kind) {
                    .fix => {
                        try non_zero_entries.append(row_idx); // A[i, j]
                    },
                    .interface => {
                        // TODO: we also need to know which is the point that needs to be solved for!
                        // TODO: here we need all 9 point stencil matrix indices (that are equal to the global point indices)

                        // Option: split the contirbutions - but we would need to add another DOF to enforce the points to be identical
                        // Option: try to collect the info that is needed to assemble the full coefficients for one point and enforce
                        // identity for the second point

                        try non_zero_entries.append(row_idx - col_size - 1); // A[i-1, j-1]
                        try non_zero_entries.append(row_idx - col_size); // A[i-1, j]
                        try non_zero_entries.append(row_idx - col_size + 1); // A[i-1, j+1]
                        try non_zero_entries.append(row_idx - 1); // A[i, j-1]
                        try non_zero_entries.append(row_idx); // A[i, j]
                        try non_zero_entries.append(row_idx + 1); // A[i, j+1]
                        try non_zero_entries.append(row_idx + col_size - 1); // A[i+1, j-1]
                        try non_zero_entries.append(row_idx + col_size); // A[i+1, j]
                        try non_zero_entries.append(row_idx + col_size + 1); // A[i+1, j+1]

                        // TODO: if point is to be smoothed add 9 undefiened entries. The values will be computed later based on the connection data.
                        // TODO: if point is to be connected add 2 undefiened entries. The values will be computed later based on the connection data.
                    },
                    .junction => {
                        unreachable;
                    },
                }

                try row_count.append(@intCast(non_zero_entries.items.len));
                row_idx += 1;
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                {
                    try non_zero_entries.append(row_idx); // A[i, j]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row_idx += 1;
                }

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
                {
                    try non_zero_entries.append(row_idx); // A[i, j]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row_idx += 1;
                }
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                try non_zero_entries.append(row_idx); // A[i, j]
                try row_count.append(@intCast(non_zero_entries.items.len));
                row_idx += 1;
            }
        }

        // NOTE: it might be better to process the connections and split it into
        // ranges to be smoothed and points to be handled by Laplacian smoothing.

        // TODO: move this to own function.
        {
            for (self.mesh.connections.items) |connection| {

                // NOTE: we rely on block 0 matrix entries to be ahead of the block 1 matrix entries,
                // otherwise we run into an invalid matrix error. We could handle this here dynamically in the future:
                // - move index 0 to a variable lower_idx_block
                // - move index 1 to a variable higher_idx_block
                std.debug.assert(connection.data[0].block < connection.data[1].block);

                // we need at least the two extreme points handled in the current implementation;
                // it should be simple to add handling of the edge cases.
                std.debug.assert(connection.lenInternal() > 3);

                // NOTE: we need to make sure that the matrix entries are in ascending order
                // as required by the compressed row format
                var it = RangeFillMatrixIterator.init(connection, self.mesh);

                // TODO: remove this hardcoded handling of the first point
                {
                    const connected_boundary_points = it.next().?;

                    const boundary_idx_0: isize = @intCast(connected_boundary_points[0]);
                    const boundary_idx_1: isize = @intCast(connected_boundary_points[1]);

                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
                    non_zero_entries[row_non_zero_entries_start_idx] = boundary_idx_0;
                    non_zero_entries[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
                }

                for (0..it.count - 1) |_| {
                    const connected_boundary_points = it.next().?;

                    const boundary_idx_0: isize = @intCast(connected_boundary_points[0]);
                    const boundary_idx_1: isize = @intCast(connected_boundary_points[1]);

                    // smooth 1st point
                    {
                        const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_0);

                        // TODO: here we just set the data. The allocation must have been done before!
                        // TODO: compute this dynamically. For now this is hardcoded for the first connection.

                        non_zero_entries[row_non_zero_entries_start_idx + 0] = boundary_idx_0 + it.first_internal_point_shift[0] - it.in_connection_direction_shift[0];
                        non_zero_entries[row_non_zero_entries_start_idx + 1] = boundary_idx_0 + it.first_internal_point_shift[0];
                        non_zero_entries[row_non_zero_entries_start_idx + 2] = boundary_idx_0 + it.first_internal_point_shift[0] + it.in_connection_direction_shift[0];
                        non_zero_entries[row_non_zero_entries_start_idx + 3] = boundary_idx_0 - it.in_connection_direction_shift[0];
                        non_zero_entries[row_non_zero_entries_start_idx + 4] = boundary_idx_0;
                        non_zero_entries[row_non_zero_entries_start_idx + 5] = boundary_idx_0 + it.in_connection_direction_shift[0];
                        non_zero_entries[row_non_zero_entries_start_idx + 6] = boundary_idx_1 + it.first_internal_point_shift[1] + it.in_connection_direction_shift[1];
                        non_zero_entries[row_non_zero_entries_start_idx + 7] = boundary_idx_1 + it.first_internal_point_shift[1];
                        non_zero_entries[row_non_zero_entries_start_idx + 8] = boundary_idx_1 + it.first_internal_point_shift[1] - it.in_connection_direction_shift[1];

                        // check that we have all entries ascending as required by the row compressed format.
                        std.debug.assert(blk: {
                            for (0..8) |i| {
                                if (non_zero_entries[row_non_zero_entries_start_idx + i] > non_zero_entries[row_non_zero_entries_start_idx + i + 1]) break :blk true;
                            }
                            break :blk false;
                        });
                    }

                    // enforce 2nd point to conform to 1st
                    {
                        const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
                        non_zero_entries[row_non_zero_entries_start_idx] = boundary_idx_0;
                        non_zero_entries[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
                    }
                }

                // TODO: remove this hardcoded handling of the first point
                {
                    const connected_boundary_points = it.next().?;

                    const boundary_idx_0: isize = @intCast(connected_boundary_points[0]);
                    const boundary_idx_1: isize = @intCast(connected_boundary_points[1]);

                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(boundary_idx_1);
                    non_zero_entries[row_non_zero_entries_start_idx] = boundary_idx_0;
                    non_zero_entries[row_non_zero_entries_start_idx + 1] = boundary_idx_1;
                }
            }
        }
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

    fn fillBlockConnectionData(mesh_data: *const discrete.Mesh) void {
        for (mesh_data.connections.items) |connection| {
            const point_data = [2][]types.Vec2d{
                mesh_data.blocks.items[connection.data[0].block].points.data,
                mesh_data.blocks.items[connection.data[1].block].points.data,
            };

            var it = RangeFillMatrixIterator.init(connection, mesh_data);

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
                const boundary_idx_0: isize = @intCast(connected_boundary_points[0]);
                const boundary_idx_1: isize = @intCast(connected_boundary_points[1]);

                // we add the stencil data to the 1st point and force equality to it's solution for the 2nd point
                {
                    const row_idx = boundary_idx_0;

                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(row_idx);

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
                    const row_idx = boundary_idx_1;
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(row_idx);

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
        fillBlockConnectionData();

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, lhs, rhs_x, self.x_new[0..], rhs_y, self.y_new[0..]);
    }
};

const BlockBoundaryPointKind = enum {
    fix, // not smoothed

    // interface points are either smoothed or connected to a smoothed point
    smooth,
    connect,

    interface, // smoothed just like interior points

    junction, // smoothed using Laplacian smoothing
};

/// provides for every block boundary point the point type info (fixed, ...) based on the number of connections found for each block boundary point
const BlockBoundaryPointsInfo = struct {
    // /// contains all connections a point is contained in (in place buffer of given size)
    // point_connections: boundary.PointData(std.BoundedArray(usize, 4)),

    /// contains for each point a buffer with all connected points
    connected_points: boundary.PointData(std.BoundedArray(usize, 4)),

    fn init(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh) BlockBoundaryPointsInfo {
        // var point_connections = try boundary.PointData(std.BoundedArray(usize, 4)).init(allocator, mesh_data);
        var info = BlockBoundaryPointsInfo{ .connected_points = try .init(allocator, mesh_data) };

        for (mesh_data.connections.items) |connection| {
            var it_0 = info.connected_points.iterateRange(connection.data[0]);
            var it_1 = info.connected_points.iterateRange(connection.data[1]);

            _ = it_0;
            _ = it_1;

            // while (true) {
            //
            // }

            // while(it.next()) |connected_points| {
            //
            // }
        }

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

        // TODO: remove this: hard coded fixed
        boundary_point_props.set(4, .{ 0, 0 }, .fix);
        boundary_point_props.set(4, .{ 0, 0 }, .fix);

        // info.point_connections.set(4, .{0, 0})

        return info;
    }

    fn deinit(self: BlockBoundaryPointsInfo) void {
        self.point_connections.deinit();
    }

    fn getKind(
        self: BlockBoundaryPointsInfo,
        block_idx: usize,
        side: boundary.Side,
        point_idx: usize,
    ) BlockBoundaryPointKind {
        self.connections_count.get(block_idx, side, point_idx);
    }
};

fn categorizeBlockBoundaryPoints(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh) boundary.PointData(BlockBoundaryPointKind) {

    // (1) collect for all block boundary points the number of connections the point is contained in.
    // (2) Define how the point is to be put into the system matrix based on the number of connections.

    var boundary_point_connections = try boundary.PointData(std.BoundedArray(usize, 4)).init(allocator, mesh_data);
    defer boundary_point_connections.deinit();

    for (mesh_data.connections.items, 0..) |connection, connection_idx| {
        for (0..2) |side_idx| {
            var it = boundary_point_connections.iterateRange(connection.data[side_idx]);
            while (it.nextPtr()) |point_connection_data| {
                try point_connection_data.append(connection_idx);
            }
        }
    }

    // tag boundary points based on connections
    var boundary_point_props = try boundary.PointData(BlockBoundaryPointProp).init(allocator, mesh_data);
    for (boundary_point_connections.buffer, boundary_point_props.buffer) |connections, *prop| {
        switch (connections.len) {
            0 => prop.* = .fix,
            1 => {
                const connection_idx = connections[0];
                const connection = mesh_data.connections.items[connection_idx];
                // set one side to solve
                // set the other side to connect
            },
            else => prop.* = .junction,
        }
    }

    return boundary_point_kind;

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

        // TODO: fix this! We need the global point idx and the start of the non zero entires of this point!
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
