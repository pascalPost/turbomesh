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

        var dof: usize = 0;
        var non_zero_entries_capacity: usize = 0;

        for (mesh_data.blocks.items) |b| {
            const points = b.points;
            const block_dof = points.size[0] * points.size[1];
            dof += block_dof;

            // this is a max limit of non zero entries per point
            // the true number will be less dependent on how many
            // neighboring points need to bes solved
            const block_non_zero_entries = dof * 9;
            non_zero_entries_capacity += block_non_zero_entries;
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

        try RowCompressedMatrixSystem2d.nonZeroMatrixEntries(mesh_data, lhs_p, lhs_i, dof, non_zero_entries_capacity);

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

    fn nonZeroMatrixEntries(mesh_data: *const discrete.Mesh, lhs_p: []c_int, lhs_i: []c_int, dof: usize, non_zero_entries_capacity: usize) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_i[0..]));
        var non_zero_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

        var row: c_int = 0;

        for (mesh_data.blocks.items) |block| {
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
                try non_zero_entries.append(row); // A[i, j]
                try row_count.append(@intCast(non_zero_entries.items.len));
                row += 1;
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                {
                    try non_zero_entries.append(row); // A[i, j]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row += 1;
                }

                // internal points with full 9 point stencil
                for (1..block.points.size[1] - 1) |_| {
                    try non_zero_entries.append(row - col_size - 1); // A[i-1, j-1]
                    try non_zero_entries.append(row - col_size); // A[i-1, j]
                    try non_zero_entries.append(row - col_size + 1); // A[i-1, j+1]
                    try non_zero_entries.append(row - 1); // A[i, j-1]
                    try non_zero_entries.append(row); // A[i, j]
                    try non_zero_entries.append(row + 1); // A[i, j+1]
                    try non_zero_entries.append(row + col_size - 1); // A[i+1, j-1]
                    try non_zero_entries.append(row + col_size); // A[i+1, j]
                    try non_zero_entries.append(row + col_size + 1); // A[i+1, j+1]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row += 1;
                }

                // edge i_max
                {
                    try non_zero_entries.append(row); // A[i, j]
                    try row_count.append(@intCast(non_zero_entries.items.len));
                    row += 1;
                }
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                try non_zero_entries.append(row); // A[i, j]
                try row_count.append(@intCast(non_zero_entries.items.len));
                row += 1;
            }
        }
    }

    fn fillAndSolve(self: *RowCompressedMatrixSystem2d) !void {
        var lhs = self.lhs_values;

        var rhs_x = self.rhs_x;
        var rhs_y = self.rhs_y;

        var non_zero_entry_idx: usize = 0;
        var row_idx: usize = 0;

        // laplace conditions (zero intialization)
        const s = 0.0;
        const t = 0.0;

        for (self.mesh.blocks.items) |block| {
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

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, lhs, rhs_x, self.x_new[0..], rhs_y, self.y_new[0..]);
    }
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
        // TODO inefficient - enhance
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
