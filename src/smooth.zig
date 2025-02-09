const std = @import("std");
const discrete = @import("discrete.zig");
const types = @import("types.zig");
const umfpack = @import("umfpack.zig");

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
// The matrix structure handles special cases for:
// - Corner points (4 neighbors)
// - Edge points (6 neighbors)
// - Interior points (9 neighbors)
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

fn computeMatrixValues(im1_j: types.Vec2d, ip1_j: types.Vec2d, i_jm1: types.Vec2d, i_jp1: types.Vec2d, s: f64, t: f64) struct {
    a_i_j: f64,
    a_ip1_j: f64,
    a_im1_j: f64,
    a_i_jp1: f64,
    a_i_jm1: f64,
    a_ip1_jp1: f64,
    a_ip1_jm1: f64,
    a_im1_jp1: f64,
    a_im1_jm1: f64,
} {
    const x_xi = 0.5 * (ip1_j.data[0] - im1_j.data[0]);
    const x_eta = 0.5 * (i_jp1.data[0] - i_jm1.data[0]);
    const y_xi = 0.5 * (ip1_j.data[1] - im1_j.data[1]);
    const y_eta = 0.5 * (i_jp1.data[1] - i_jm1.data[1]);

    const p = x_eta * x_eta + y_eta * y_eta;
    const q = x_xi * x_eta + y_xi * y_eta;
    const r = x_xi * x_xi + y_xi * y_xi;

    return .{
        .a_i_j = -2.0 * p - 2.0 * r,
        .a_ip1_j = p + 0.5 * s,
        .a_im1_j = p - 0.5 * s,
        .a_i_jp1 = r + 0.5 * t,
        .a_i_jm1 = r - 0.5 * t,
        .a_ip1_jp1 = -0.5 * q,
        .a_ip1_jm1 = 0.5 * q,
        .a_im1_jp1 = 0.5 * q,
        .a_im1_jm1 = -0.5 * q,
    };
}

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
        var dof: usize = 0;
        var non_zero_entries: usize = 0;

        for (mesh_data.blocks.items) |b| {
            const points = b.points;
            const block_dof = (points.size[0] - 2) * (points.size[1] - 2);
            dof += block_dof;

            const block_non_zero_entries =
                // corners
                4 * 4 +
                // edges
                ((points.size[0] - 4) + (points.size[1] - 4)) * 2 * 6 +
                // internal
                (points.size[0] - 4) * (points.size[1] - 4) * 9;
            non_zero_entries += block_non_zero_entries;
        }

        const buffer_int = try allocator.alloc(c_int, (dof + 1) + non_zero_entries);

        // left hand side (lhs) in compressed row form
        var lhs_p = buffer_int[0 .. dof + 1]; // cum sum of non-zero entries in columns
        lhs_p[0] = 0; // first value must be zero
        const lhs_i = buffer_int[dof + 1 ..]; // row indicies with non-zero values
        std.debug.assert(lhs_i.len == non_zero_entries);

        const buffer_float = try allocator.alloc(f64, 4 * dof + non_zero_entries);

        const x_new = buffer_float[0..dof];
        const y_new = buffer_float[dof .. 2 * dof];

        const rhs_x = buffer_float[2 * dof .. 3 * dof];
        const rhs_y = buffer_float[3 * dof .. 4 * dof];

        const lhs_values = buffer_float[4 * dof ..];

        try RowCompressedMatrixSystem2d.nonZeroMatrixEntries(mesh_data, lhs_p, lhs_i, dof, non_zero_entries);

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

    fn nonZeroMatrixEntries(mesh_data: *const discrete.Mesh, lhs_p: []c_int, lhs_i: []c_int, dof: usize, non_zero_entries: usize) !void {
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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 1] = values.a_i_jp1; // A[i, j+1]
                lhs_values[matrix_idx + 2] = values.a_ip1_j; // A[i+1, j]
                lhs_values[matrix_idx + 3] = values.a_ip1_jp1; // A[i+1, j+1]

                matrix_idx += 4;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j-1] A[i+1, j-1]

                rhs_x[rhs_idx + 0] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0] +
                    values.a_i_jm1 * i_jm1.data[0] +
                    values.a_ip1_jm1 * ip1_jm1.data[0]);

                rhs_y[rhs_idx + 0] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1] +
                    values.a_i_jm1 * i_jm1.data[1] +
                    values.a_ip1_jm1 * ip1_jm1.data[1]);

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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 2] = values.a_i_jp1; // A[i, j+1]
                lhs_values[matrix_idx + 3] = values.a_ip1_jm1; // A[i+1, j-1]
                lhs_values[matrix_idx + 4] = values.a_ip1_j; // A[i+1, j]
                lhs_values[matrix_idx + 5] = values.a_ip1_jp1; // A[i+1, j+1]

                matrix_idx += 6;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1]);

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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 2] = values.a_ip1_jm1; // A[i+1, j-1]
                lhs_values[matrix_idx + 3] = values.a_ip1_j; // A[i+1, j]

                matrix_idx += 4;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j+1] A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0] +
                    values.a_i_jp1 * i_jp1.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1] +
                    values.a_i_jp1 * i_jp1.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1]);

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

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 1] = values.a_im1_jp1; // A[i-1, j+1]
                    lhs_values[matrix_idx + 2] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 3] = values.a_i_jp1; // A[i, j+1]
                    lhs_values[matrix_idx + 4] = values.a_ip1_j; // A[i+1, j]
                    lhs_values[matrix_idx + 5] = values.a_ip1_jp1; // A[i+1, j+1]

                    matrix_idx += 6;

                    // RHS: A[i-1, j-1], A[i, j-1], A[i+1, j-1]

                    rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                        values.a_i_jm1 * i_jm1.data[0] +
                        values.a_ip1_jm1 * ip1_jm1.data[0]);

                    rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                        values.a_i_jm1 * i_jm1.data[1] +
                        values.a_ip1_jm1 * ip1_jm1.data[1]);

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

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.a_im1_jp1; // A[i-1, j+1]
                    lhs_values[matrix_idx + 3] = values.a_i_jm1; // A[i, j-1]
                    lhs_values[matrix_idx + 4] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 5] = values.a_i_jp1; // A[i, j+1]
                    lhs_values[matrix_idx + 6] = values.a_ip1_jm1; // A[i+1, j-1]
                    lhs_values[matrix_idx + 7] = values.a_ip1_j; // A[i+1, j]
                    lhs_values[matrix_idx + 8] = values.a_ip1_jp1; // A[i+1, j+1]

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

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.a_i_jm1; // A[i, j-1]
                    lhs_values[matrix_idx + 3] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 4] = values.a_ip1_jm1; // A[i+1, j-1]
                    lhs_values[matrix_idx + 5] = values.a_ip1_j; // A[i+1, j]

                    matrix_idx += 6;

                    // RHS: A[i-1, j+1], A[i, j+1], A[i+1, j+1]

                    rhs_x[rhs_idx] = -(values.a_im1_jp1 * im1_jp1.data[0] +
                        values.a_i_jp1 * i_jp1.data[0] +
                        values.a_ip1_jp1 * ip1_jp1.data[0]);

                    rhs_y[rhs_idx] = -(values.a_im1_jp1 * im1_jp1.data[1] +
                        values.a_i_jp1 * i_jp1.data[1] +
                        values.a_ip1_jp1 * ip1_jp1.data[1]);

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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 1] = values.a_im1_jp1; // A[i-1, j+1]
                lhs_values[matrix_idx + 2] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 3] = values.a_i_jp1; // A[i, j+1]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j-1], A[i-1, j-1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0] +
                    values.a_i_jm1 * i_jm1.data[0] +
                    values.a_im1_jm1 * im1_jm1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1] +
                    values.a_i_jm1 * i_jm1.data[1] +
                    values.a_im1_jm1 * im1_jm1.data[1]);

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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.a_im1_jp1; // A[i-1, j+1]
                lhs_values[matrix_idx + 3] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 4] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 5] = values.a_i_jp1; // A[i, j+1]

                matrix_idx += 6;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1]);

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

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 3] = values.a_i_j; // A[i, j]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j+1], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0] +
                    values.a_i_jp1 * i_jp1.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1] +
                    values.a_i_jp1 * i_jp1.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1]);

                rhs_idx += 1;
            }
        }

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, lhs_values, rhs_x, self.x_new[0..], rhs_y, self.y_new[0..]);
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
        }
    }
}
