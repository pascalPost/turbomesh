// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const builtin = @import("builtin");
const config = @import("config");
const discrete = @import("../discrete.zig");
const types = @import("../types.zig");
const boundary = @import("../boundary.zig");
const wall_control_function = @import("wall_control_function.zig");
const solver = @import("solver.zig");

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
//    - Uses the selected sparse solver backend to update interior point positions
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

const log = std.log.scoped(.smoothing);

fn writeFile(name: []const u8, field: anytype, buffer: []u8) !void {
    var file = try std.fs.cwd().createFile(name, .{});
    defer file.close();

    var file_writer = file.writer(buffer);
    var writer = &file_writer.interface;

    for (field) |item| {
        try writer.print("{}\n", .{item});
    }

    try writer.flush();
}

pub fn mesh(
    allocator: std.mem.Allocator,
    mesh_data: *discrete.Mesh,
    iterations: usize,
    solver_option: solver.Option,
    control_function_algorithm: wall_control_function.Algorithm,
) !void {
    const time_start =
        if (comptime builtin.cpu.arch != .wasm32)
            try std.time.Instant.now()
        else
            void;

    var system = try RowCompressedMatrixSystem2d.init(allocator, mesh_data, control_function_algorithm);
    defer system.deinit();

    // debug output
    // if (builtin.cpu.arch != .wasm32) {
    //     var buffer: [1024]u8 = undefined;
    //     try writeFile("Ap.txt", system.lhs_p, &buffer);
    //     try writeFile("Ai.txt", system.lhs_i, &buffer);
    //     try writeFile("Ax.txt", system.lhs_values, &buffer);
    //     try writeFile("rhs_x.txt", system.rhs_x, &buffer);
    //     try writeFile("rhs_y.txt", system.rhs_y, &buffer);
    // }

    var s = try solver.Solver.init(solver_option, system);
    defer s.deinit();

    // iterate and fill matrix values
    for (0..iterations) |n| {
        log.info("iteration: {}", .{n});

        // fill common parts of the matrix
        try system.fill(n);

        try s.solve();

        var x_norm_sqr: f64 = 0.0;
        var y_norm_sqr: f64 = 0.0;

        {
            var row_block_start_id: usize = 0;
            for (mesh_data.blocks.items) |block| {
                var point_id: usize = 0;
                for (0..block.points.size[0]) |_| {
                    for (0..block.points.size[1]) |_| {
                        const row_idx = row_block_start_id + point_id;
                        const p = block.points.data[point_id];
                        const dx = p.data[0] - system.x_new[row_idx];
                        const dy = p.data[1] - system.y_new[row_idx];

                        x_norm_sqr += dx * dx;
                        y_norm_sqr += dy * dy;

                        point_id += 1;
                    }
                }
                row_block_start_id += block.points.size[0] * block.points.size[1];
            }
        }

        const norm = (x_norm_sqr + y_norm_sqr) * (x_norm_sqr + y_norm_sqr);
        log.info("\tresidual: {}", .{norm});

        // copy into coordinate field
        {
            var row_block_start_id: usize = 0;
            for (mesh_data.blocks.items) |block| {
                var point_id: usize = 0;
                for (0..block.points.size[0]) |_| {
                    for (0..block.points.size[1]) |_| {
                        const row_idx = row_block_start_id + point_id;
                        block.points.data[point_id] = types.Vec2d.init(system.x_new[row_idx], system.y_new[row_idx]);
                        point_id += 1;
                    }
                }
                row_block_start_id += block.points.size[0] * block.points.size[1];
            }
        }
    }

    if (comptime builtin.cpu.arch != .wasm32) {
        const time_end = try std.time.Instant.now();
        const time_delta: f32 = @floatFromInt(time_end.since(time_start));
        log.info("elapsed time for smoothing: {d:.2} s\n", .{time_delta / std.time.ns_per_s});
    }

    // TODO: remove asap (e.g. with an option)
    if (config.use_cgns) {
        try system.write("smooth.cgns");
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

    fn init(im1_j: types.Vec2d, ip1_j: types.Vec2d, i_jm1: types.Vec2d, i_jp1: types.Vec2d, P: f64, Q: f64) StencilData {
        const x_xi = 0.5 * (ip1_j.data[0] - im1_j.data[0]);
        const x_eta = 0.5 * (i_jp1.data[0] - i_jm1.data[0]);
        const y_xi = 0.5 * (ip1_j.data[1] - im1_j.data[1]);
        const y_eta = 0.5 * (i_jp1.data[1] - i_jm1.data[1]);

        const g22 = x_eta * x_eta + y_eta * y_eta;
        const g12 = x_xi * x_eta + y_xi * y_eta;
        const g11 = x_xi * x_xi + y_xi * y_xi;

        return .{
            .data = .{
                -2.0 * g22 - 2.0 * g11, // i_j
                g22 * (1 + 0.5 * P), //ip1_j
                g22 * (1 - 0.5 * P), //im1_j
                g11 * (1 + 0.5 * Q), //i_jp1
                g11 * (1 - 0.5 * Q), //i_jm1
                -0.5 * g12, //ip1_jp1
                0.5 * g12, //ip1_jm1
                0.5 * g12, //im1_jp1
                -0.5 * g12, //im1_jm1
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
                    if (it_1.next() != null) {
                        std.debug.panic("connection {} point {}", .{ connection_idx, point_idx });
                    }
                    break;
                };
                const p_1 = it_1.next() orelse {
                    std.debug.panic("connection {} point {}", .{ connection_idx, point_idx });
                };

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
                const p_1 = it_1.next() orelse {
                    std.debug.panic("non matching connection data for connection {}.", .{connection_idx});
                };

                const x_0 = mesh_data.blocks.items[range_0.block].points.data[p_0];
                const x_1 = mesh_data.blocks.items[range_1.block].points.data[p_1];

                if (!types.eqlApprox(x_0, x_1, abs_tol)) {
                    std.debug.panic("non matching points for connection {} point {}:\n\t{}\n\t{}\n", .{ connection_idx, point_idx, x_0, x_1 });
                }

                point_idx += 1;
            }
        }
    }

    // TODO: check that endpoints match! That allows to ease the connection handling!
}

pub const RowCompressedMatrixSystem2d = struct {
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

    boundary_points: BlockBoundaryPoints,
    control_function: wall_control_function.ControlFunction,
    index_converter: IndexConverter,

    fn init(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh, control_function_algorithm: wall_control_function.Algorithm) !RowCompressedMatrixSystem2d {
        connectionDataCheck(mesh_data);

        // TODO: check that endpoints match! That allows to ease the connection handling!

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

        const index_converter = try IndexConverter.init(mesh_data, allocator);
        errdefer index_converter.deinit();

        const boundary_points = try BlockBoundaryPoints.init(allocator, index_converter, mesh_data);
        errdefer boundary_points.deinit();

        const control_function = try wall_control_function.ControlFunction.init(allocator, dof, mesh_data.*, control_function_algorithm);
        errdefer control_function.deinit();

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
            .boundary_points = boundary_points,
            .control_function = control_function,
            .index_converter = index_converter,
        };

        try system.initNonZeroMatrixEntries(dof, non_zero_entries_capacity);
        system.initBoundaryData();

        return system;
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        self.allocator.free(self.buffer_int);
        self.allocator.free(self.buffer_float);
        self.allocator.free(self.row_idx_range_start_for_each_block);
        self.boundary_points.deinit();
        self.control_function.deinit();
        self.index_converter.deinit();
    }

    fn write(self: RowCompressedMatrixSystem2d, filename: [:0]const u8) !void {
        if (config.use_cgns) {
            const cgns = @import("../cgns.zig");
            // buffer
            var size: usize = 0;
            for (self.mesh.blocks.items) |b| {
                size = @max(size, b.points.data.len);
            }
            var buffer = try self.allocator.alloc(types.Float, size);
            defer self.allocator.free(buffer);

            // block data
            var block_points = try self.allocator.alloc(types.Mat2d, self.mesh.blocks.items.len);
            defer self.allocator.free(block_points);
            for (self.mesh.blocks.items, block_points[0..]) |b, *data| data.* = b.points;

            try cgns.write(filename, self.mesh.names.items, block_points, buffer[0..], self.control_function.data);
        } else return error.OutputFormatNotEnabled;
    }

    /// returns the non zero entries range for the given row (the row index range start for each block can be retrieved from an array in the struct)
    fn nonZeroEntriesRangeStart(self: RowCompressedMatrixSystem2d, row_idx: c_int) usize {
        return @intCast(self.lhs_p[@intCast(row_idx)]);
    }

    fn initNonZeroMatrixEntriesForBoundaryPoint(
        non_zero_entries: *std.array_list.Managed(c_int),
        row_count: *std.array_list.Managed(c_int),
        boundary_point_id: *usize,
        boundary_points: BlockBoundaryPoints,
        row_idx: *c_int, // equivalent to the global point index
        laplacian_count: *usize,
    ) !void {
        switch (boundary_points.kind.buffer[boundary_point_id.*]) {
            .fixed => {
                // single entry with the current coordinates on the RHS
                try non_zero_entries.append(row_idx.*); // A[i, j]
            },
            .smoothed => {
                // smooth the block boundary point based on the same 9 point stancil that is used for block
                // internal points. The actual data is set in a connection based loop after this undefined
                // allocation.
                _ = try non_zero_entries.addManyAsSlice(9);
            },
            .connected => {
                // enforce the equality for the connected points with two entries.
                _ = try non_zero_entries.addManyAsSlice(2);
            },
            .laplacian_smoothed => {
                const stencil_ids = boundary_points.laplacian_points.items[laplacian_count.*].stencil_ids.slice();
                try non_zero_entries.appendSlice(stencil_ids);
                laplacian_count.* += 1;
            },
            .sliding_circ => {
                // we need 2 entries for the sliding circular boundary points [Neumann BC] (for a fixed axial direction, i.e. Dirichlet BC, we only need 1)
                _ = try non_zero_entries.addManyAsSlice(2);
            },
        }
        boundary_point_id.* += 1;

        try row_count.append(@intCast(non_zero_entries.items.len));
        row_idx.* += 1;
    }

    fn initNonZeroMatrixEntriesPointBased(
        self: @This(),
        non_zero_entries: *std.array_list.Managed(c_int),
        row_count: *std.array_list.Managed(c_int),
    ) !void {
        var boundary_point_idx: usize = 0;
        var row_idx: c_int = 0;
        var laplacian_count: usize = 0; // needed for accessing the right laplacian point info
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
                try initNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, self.boundary_points, &row_idx, &laplacian_count);
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {

                // edge i_min
                try initNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, self.boundary_points, &row_idx, &laplacian_count);

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
                try initNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, self.boundary_points, &row_idx, &laplacian_count);
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                try initNonZeroMatrixEntriesForBoundaryPoint(non_zero_entries, row_count, &boundary_point_idx, self.boundary_points, &row_idx, &laplacian_count);
            }
        }
    }

    fn computeConnectionStencilPositions(connection: boundary.Connection, it: RangeFillMatrixIterator) [9]usize {
        // NOTE: this is needed to add the point indices in ascending order depending on the direction of the connection (see the non zero matrix entries definition).

        // handle connections to within a block
        if (connection.ranges[0].block == connection.ranges[1].block) {
            if (connection.ranges[0].side == .i_min and connection.ranges[1].side == .i_max) {
                if (it.in_connection_direction_shift[0] > 0) {

                    // x----x----x  i_max
                    //
                    // x    x    x
                    // 2    5    8
                    //
                    // ~~~~~~~~~~~~
                    //
                    //
                    // x    x    x
                    // 1    4    7
                    //
                    // x----x----x -> connection i_min
                    // 0    3    6

                    std.debug.assert(it.in_connection_direction_shift[1] > 0);
                    return .{
                        1, 4, 7, // internal points on i_min side
                        0, 3, 6, // points on connection
                        2, 5, 8, // internal points on i_max side
                    };
                } else {

                    // transpose of the points above as the connection runs in the opposite direction

                    std.debug.assert(it.in_connection_direction_shift[1] < 0);
                    return .{
                        7, 4, 1, // internal points on i_min side
                        6, 3, 0, // points on connection
                        8, 5, 2, // internal points on i_max side
                    };
                }
            } else {
                unreachable;
            }
        }

        std.debug.assert(connection.ranges[0].block < connection.ranges[1].block);

        const direction_modifier: [2]c_int = .{
            if (it.in_connection_direction_shift[0] > 0) 1.0 else -1.0,
            if (it.in_connection_direction_shift[1] > 0) 1.0 else -1.0,
        };

        // NOTE: this is needed to have the right ordering of the points taken from the 0th block of the connection (entries 0 - 5):
        // the ordering is first the points inside block 0 ([i-1,j-1],[i,j-1],[i+1,j-1]) and then the points on the connection
        // ([i-1,j],[i,j],[i+1,j]).

        var positions: [9]usize = undefined;

        switch (connection.ranges[0].side) {
            .i_min => {
                positions[0] = @intCast(3 - 2 * direction_modifier[0]);
                positions[1] = 3;
                positions[2] = @intCast(3 + 2 * direction_modifier[0]);
                positions[3] = @intCast(2 - 2 * direction_modifier[0]);
                positions[4] = 2;
                positions[5] = @intCast(2 + 2 * direction_modifier[0]);
            },
            .i_max => {
                positions[0] = @intCast(2 - 2 * direction_modifier[0]);
                positions[1] = 2;
                positions[2] = @intCast(2 + 2 * direction_modifier[0]);
                positions[3] = @intCast(3 - 2 * direction_modifier[0]);
                positions[4] = 3;
                positions[5] = @intCast(3 + 2 * direction_modifier[0]);
            },
            .j_min => {
                positions[0] = @intCast(4 - direction_modifier[0]);
                positions[1] = 4;
                positions[2] = @intCast(4 + direction_modifier[0]);
                positions[3] = @intCast(1 - direction_modifier[0]);
                positions[4] = 1;
                positions[5] = @intCast(1 + direction_modifier[0]);
            },
            .j_max => {
                positions[0] = @intCast(1 - direction_modifier[0]);
                positions[1] = 1;
                positions[2] = @intCast(1 + direction_modifier[0]);
                positions[3] = @intCast(4 - direction_modifier[0]);
                positions[4] = 4;
                positions[5] = @intCast(4 + direction_modifier[0]);
            },
        }

        // positions in block 1
        positions[6] = @intCast(7 - direction_modifier[1]);
        positions[7] = 7;
        positions[8] = @intCast(7 + direction_modifier[1]);

        return positions;
    }

    fn initNonZeroMatrixEntriesConnectionBased(
        self: @This(),
        non_zero_entries: *std.array_list.Managed(c_int),
    ) !void {
        for (self.mesh.connections.items) |connection| {
            // NOTE: we rely on block 0 matrix entries to be ahead of the block 1 matrix entries,
            // otherwise we run into an invalid matrix error. We could handle this here dynamically in the future:
            // - move index 0 to a variable lower_idx_block
            // - move index 1 to a variable higher_idx_block
            std.debug.assert(connection.ranges[0].block <= connection.ranges[1].block);

            // we need at least the two extreme points handled in the current implementation;
            // it should be simple to add handling of the edge cases.
            std.debug.assert(connection.lenInternal() > 3);

            // NOTE: we need to make sure that the matrix entries are in ascending order
            // as required by the compressed row format
            var it = RangeFillMatrixIterator.init(connection, self.mesh);

            const point_stencil_idx = computeConnectionStencilPositions(connection, it);

            try self.initNonZeroMatrixForConnectionEndpoint(it.next().?, connection, non_zero_entries);

            // middle of the connection
            for (0..it.count - 1) |_| {
                const connected_points_local_idx = it.next().?;
                const connected_points_global_idx = [2]c_int{
                    @intCast(self.index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = connected_points_local_idx[0] }).value),
                    @intCast(self.index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = connected_points_local_idx[1] }).value),
                };

                // connect 2nd point to 1st
                {
                    std.debug.assert(connected_points_global_idx[0] < connected_points_global_idx[1]);
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(connected_points_global_idx[1]);
                    non_zero_entries.items[row_non_zero_entries_start_idx] = connected_points_global_idx[0];
                    non_zero_entries.items[row_non_zero_entries_start_idx + 1] = connected_points_global_idx[1];
                }

                // smooth 1st point
                {
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(connected_points_global_idx[0]);

                    // NOTE: here we set the data based on the connections; the allocation happend in the previous function based on a point loop.

                    // points inside block 0
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[0]] = connected_points_global_idx[0] - it.in_connection_direction_shift[0] + it.first_internal_point_shift[0]; // i-1,j-1
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[1]] = connected_points_global_idx[0] + it.first_internal_point_shift[0]; // i,j-1
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[2]] = connected_points_global_idx[0] + it.in_connection_direction_shift[0] + it.first_internal_point_shift[0]; //i+1,j-1

                    // points on boundary connection (taken from block 0)
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[3]] = connected_points_global_idx[0] - it.in_connection_direction_shift[0]; // i-1,j
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[4]] = connected_points_global_idx[0]; // i,j
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[5]] = connected_points_global_idx[0] + it.in_connection_direction_shift[0]; // i+1,j

                    // points inside block 1
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[6]] = connected_points_global_idx[1] - it.in_connection_direction_shift[1] + it.first_internal_point_shift[1]; // i-1,j+1
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[7]] = connected_points_global_idx[1] + it.first_internal_point_shift[1]; // i,j+1
                    non_zero_entries.items[row_non_zero_entries_start_idx + point_stencil_idx[8]] = connected_points_global_idx[1] + it.in_connection_direction_shift[1] + it.first_internal_point_shift[1]; // i+1,j+1

                    // check that we have all entries ascending as required by the row compressed format.
                    std.debug.assert(blk: {
                        for (0..8) |i| {
                            if (non_zero_entries.items[row_non_zero_entries_start_idx + i] >= non_zero_entries.items[row_non_zero_entries_start_idx + i + 1]) {
                                // std.debug.print("wrong ordering of connection stencil data for point {d}: {any}\n", .{ connected_points_global_idx[0], non_zero_entries.items[row_non_zero_entries_start_idx .. row_non_zero_entries_start_idx + 9] });
                                break :blk false;
                            }
                        }
                        break :blk true;
                    });
                }
            }

            try self.initNonZeroMatrixForConnectionEndpoint(it.position, connection, non_zero_entries);
        }
    }

    fn initNonZeroMatrixForConnectionEndpoint(
        self: @This(),
        local_ids: [2]LocalIndex,
        connection: boundary.Connection,
        non_zero_entries: *std.array_list.Managed(c_int),
    ) !void {
        const local_id_2d = self.index_converter.index2d(.{ .block = connection.ranges[0].block, .local_idx = local_ids[0] });
        const buffer_id = try self.boundary_points.index_converter.bufferIndex(local_id_2d);

        switch (self.boundary_points.kind.buffer[buffer_id]) {
            .fixed, .sliding_circ => {
                // connect the 2nd point to the 1st
                const connected_points_global_idx = [2]c_int{
                    @intCast(self.index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = local_ids[0] }).value),
                    @intCast(self.index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = local_ids[1] }).value),
                };

                std.debug.assert(connected_points_global_idx[0] < connected_points_global_idx[1]);
                const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(connected_points_global_idx[1]);
                non_zero_entries.items[row_non_zero_entries_start_idx] = connected_points_global_idx[0];
                non_zero_entries.items[row_non_zero_entries_start_idx + 1] = connected_points_global_idx[1];
            },
            .laplacian_smoothed => {}, // already set
            .connected => {},
            else => unreachable,
        }
    }

    fn initNonZeroMatrixEntries(
        self: *RowCompressedMatrixSystem2d,
        dof: usize,
        non_zero_entries_capacity: usize,
    ) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_p[1..]));
        var row_count = try std.array_list.Managed(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_i[0..]));
        var non_zero_entries = try std.array_list.Managed(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

        try self.initNonZeroMatrixEntriesPointBased(&non_zero_entries, &row_count);

        // init laplacian points

        for (self.boundary_points.laplacian_points.items) |laplacian_point| {
            const smoothed_id = laplacian_point.globalId();

            // set all other points to connected to 1st
            for (laplacian_point.overlapping_points.slice()[1..]) |overlapping_point| {
                const non_zero_entries_range_start = self.nonZeroEntriesRangeStart(@intCast(overlapping_point.global_id));
                non_zero_entries.items[non_zero_entries_range_start] = @intCast(smoothed_id);
                non_zero_entries.items[non_zero_entries_range_start + 1] = @intCast(overlapping_point.global_id);
            }
        }

        try self.initNonZeroMatrixEntriesConnectionBased(&non_zero_entries);

        for (self.mesh.boundary_conditions.items) |bc| {
            switch (bc.kind) {
                .inlet, .outlet => {
                    const first_internal_point_shift = bc.range.firstInternalPointShift(self.mesh);

                    var it = bc.range.iterate(self.mesh);
                    while (it.next()) |local_id| {
                        const block_and_local_id = BlockAndLocalIndex{ .block = bc.range.block, .local_idx = .{ .value = local_id } };
                        const local_id_2d = self.index_converter.index2d(block_and_local_id);
                        const boundary_id = try self.boundary_points.index_converter.bufferIndex(local_id_2d);
                        if (self.boundary_points.kind.buffer[boundary_id] != .sliding_circ) continue;

                        const global_id = self.index_converter.globalIndex(block_and_local_id).value;
                        const non_zero_entries_range_start = self.nonZeroEntriesRangeStart(@intCast(global_id));

                        if (first_internal_point_shift > 0) {
                            non_zero_entries.items[non_zero_entries_range_start] = @intCast(global_id);
                            non_zero_entries.items[non_zero_entries_range_start + 1] = @as(c_int, @intCast(global_id)) + first_internal_point_shift;
                        } else {
                            non_zero_entries.items[non_zero_entries_range_start] = @as(c_int, @intCast(global_id)) + first_internal_point_shift;
                            non_zero_entries.items[non_zero_entries_range_start + 1] = @intCast(global_id);
                        }
                    }
                },
                else => unreachable,
            }
        }
    }

    fn initBoundaryPointData(
        self: @This(),
        block: discrete.Block2d,
        boundary_point_id: *usize,
        row_id: *usize,
        point_id: *usize,
        non_zero_entry_id: *usize,
        laplacian_count: *usize,
    ) void {
        switch (self.boundary_points.kind.buffer[boundary_point_id.*]) {
            .fixed => {
                self.lhs_values[non_zero_entry_id.*] = 1;
                non_zero_entry_id.* += 1;

                self.rhs_x[row_id.*] = block.points.data[point_id.*].data[0];
                self.rhs_y[row_id.*] = block.points.data[point_id.*].data[1];
            },
            .smoothed => {
                // NOTE: here, we only skip the necessary stencil data since the values are set next in a connection based loop.
                non_zero_entry_id.* += 9;

                self.rhs_x[row_id.*] = 0;
                self.rhs_y[row_id.*] = 0;
            },
            .connected => {
                self.lhs_values[non_zero_entry_id.*] = 1;
                self.lhs_values[non_zero_entry_id.* + 1] = -1;
                non_zero_entry_id.* += 2;

                // NOTE: for periodic connection this is not zero and set in a connection based loop.
                self.rhs_x[row_id.*] = 0;
                self.rhs_y[row_id.*] = 0;
            },
            .laplacian_smoothed => {
                const laplacian_point = self.boundary_points.laplacian_points.items[laplacian_count.*];

                const laplacian_point_id = laplacian_point.globalId();

                var point_position_in_stencil: usize = 0;
                for (laplacian_point.stencil_ids.slice()) |global_id| {
                    if (global_id == laplacian_point_id) break;
                    point_position_in_stencil += 1;
                }

                // account for all stencil entries (the point itself is adjusted next)
                for (0..laplacian_point.stencil_ids.len) |i| self.lhs_values[non_zero_entry_id.* + i] = 1;

                // account for the point itself
                self.lhs_values[non_zero_entry_id.* + point_position_in_stencil] = -@as(f64, @floatFromInt(laplacian_point.stencil_ids.len)) + 1;

                non_zero_entry_id.* += laplacian_point.stencil_ids.len;

                self.rhs_x[row_id.*] = 0;
                self.rhs_y[row_id.*] = 0;

                laplacian_count.* += 1;
            },
            .sliding_circ => {
                // For y we need 2 entries, for x we only need one entry.
                // We have to update the coefficients for the x and y directions separately.

                non_zero_entry_id.* += 2;

                // NOTE the matrix needs to be changed for x and for y; this is done in the iteration loop.

                // for x:
                // self.lhs_values[non_zero_entry_id.*] = 1;
                // self.lhs_values[non_zero_entry_id.* + 1] = 0;
                //
                // for y:
                // self.lhs_values[non_zero_entry_id.*] = 1;
                // self.lhs_values[non_zero_entry_id.* + 1] = -1;

                const global_id = row_id.*;
                const local = self.index_converter.localIndex(.{ .value = global_id });
                const x = self.mesh.blocks.items[local.block].points.data[local.local_idx.value].data[0];

                self.rhs_x[row_id.*] = x;
                self.rhs_y[row_id.*] = 0.0;
            },
        }

        row_id.* += 1;
        boundary_point_id.* += 1;
        point_id.* += 1;
    }

    fn initBoundaryData(self: @This()) void {
        var non_zero_entry_idx: usize = 0;
        var row_idx: usize = 0;
        var boundary_point_idx: usize = 0;
        var laplacian_count: usize = 0;

        for (self.mesh.blocks.items) |block| {
            var point_idx: usize = 0;

            // edge j_min
            for (0..block.points.size[1]) |_| {
                self.initBoundaryPointData(block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx, &laplacian_count);
            }

            for (1..block.points.size[0] - 1) |_| {
                // edge i_min
                self.initBoundaryPointData(block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx, &laplacian_count);

                // middle
                // TODO: remove this loop and just add the right count in one step!!!!
                for (1..block.points.size[1] - 1) |_| {
                    non_zero_entry_idx += 9;
                    point_idx += 1;
                    row_idx += 1;
                }

                // edge i_max
                self.initBoundaryPointData(block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx, &laplacian_count);
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                self.initBoundaryPointData(block, &boundary_point_idx, &row_idx, &point_idx, &non_zero_entry_idx, &laplacian_count);
            }
        }

        // account for periodicity
        for (self.mesh.connections.items) |connection| {
            if (connection.periodicity) |periodicity| {
                var it = RangeFillMatrixIterator.init(connection, self.mesh);
                while (it.next()) |connected_points| {
                    const global_idx_1 = self.row_idx_range_start_for_each_block[connection.ranges[1].block] + connected_points[1].value;

                    // adjust RHS for periodicity for the connected point
                    self.rhs_x[global_idx_1] = -periodicity.data[0];
                    self.rhs_y[global_idx_1] = -periodicity.data[1];
                }
            }
        }

        for (self.boundary_points.laplacian_points.items) |laplacian_point| {
            self.rhs_x[laplacian_point.globalId()] = laplacian_point.rhs.data[0];
            self.rhs_y[laplacian_point.globalId()] = laplacian_point.rhs.data[1];
        }
    }

    fn fillBlockInternalPointData(self: RowCompressedMatrixSystem2d) void {
        var row_idx: usize = 0;

        // const control_function = must be either constant 0 w/0 data access or get value by row_id and field access

        for (self.mesh.blocks.items) |block| {
            var point_idx: usize = 0;

            // edge j_min
            for (0..block.points.size[1]) |_| {
                // TODO: enhance performance!
                point_idx += 1;
                row_idx += 1;
            }

            // middle
            for (1..block.points.size[0] - 1) |_| {
                // edge i_min
                row_idx += 1;
                point_idx += 1;

                // internal points with full 9 point stencil
                for (1..block.points.size[1] - 1) |_| {
                    const im1_j = block.points.data[point_idx - block.points.size[1]];
                    const i_jm1 = block.points.data[point_idx - 1];
                    const i_jp1 = block.points.data[point_idx + 1];
                    const ip1_j = block.points.data[point_idx + block.points.size[1]];

                    const cf = self.control_function.data[row_idx];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        cf.data[0],
                        cf.data[1],
                    );

                    var non_zero_entry_idx = self.nonZeroEntriesRangeStart(@intCast(row_idx));

                    self.lhs_values[non_zero_entry_idx + 0] = stencil.get(.im1_jm1); // A[i-1, j-1]
                    self.lhs_values[non_zero_entry_idx + 1] = stencil.get(.im1_j); // A[i-1, j]
                    self.lhs_values[non_zero_entry_idx + 2] = stencil.get(.im1_jp1); // A[i-1, j+1]
                    self.lhs_values[non_zero_entry_idx + 3] = stencil.get(.i_jm1); // A[i, j-1]
                    self.lhs_values[non_zero_entry_idx + 4] = stencil.get(.i_j); // A[i, j]
                    self.lhs_values[non_zero_entry_idx + 5] = stencil.get(.i_jp1); // A[i, j+1]
                    self.lhs_values[non_zero_entry_idx + 6] = stencil.get(.ip1_jm1); // A[i+1, j-1]
                    self.lhs_values[non_zero_entry_idx + 7] = stencil.get(.ip1_j); // A[i+1, j]
                    self.lhs_values[non_zero_entry_idx + 8] = stencil.get(.ip1_jp1); // A[i+1, j+1]

                    self.rhs_x[row_idx] = 0;
                    self.rhs_y[row_idx] = 0;

                    non_zero_entry_idx += 9;
                    point_idx += 1;
                    row_idx += 1;
                }

                // edge i_max
                row_idx += 1;
                point_idx += 1;
            }

            // edge j_max
            for (0..block.points.size[1]) |_| {
                // TODO: enhance performance!
                row_idx += 1;
            }
        }
    }

    fn fillBlockConnectionData(self: RowCompressedMatrixSystem2d) void {
        for (self.mesh.connections.items) |connection| {
            const point_data = [2][]types.Vec2d{
                self.mesh.blocks.items[connection.ranges[0].block].points.data,
                self.mesh.blocks.items[connection.ranges[1].block].points.data,
            };

            var it = RangeFillMatrixIterator.init(connection, self.mesh);
            it.limitToRangeInternalPoints();

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

            const point_stencil_idx = computeConnectionStencilPositions(connection, it);

            if (connection.periodicity) |periodicity| {
                while (it.next()) |connected_points| {
                    const point_idx: [2]c_int = .{ @intCast(connected_points[0].value), @intCast(connected_points[1].value) };

                    // TODO: this could be enhanced by computing this just once for the connection (plus +1 or -1 for the next connection point)
                    const global_idx_0 = @as(c_int, @intCast(self.row_idx_range_start_for_each_block[connection.ranges[0].block])) + point_idx[0];
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(global_idx_0);

                    const im1_j = point_data[0][@intCast(point_idx[0] - it.in_connection_direction_shift[0])];
                    const i_jm1 = point_data[0][@intCast(point_idx[0] + it.first_internal_point_shift[0])];
                    const ip1_j = point_data[0][@intCast(point_idx[0] + it.in_connection_direction_shift[0])];
                    const i_jp1 = types.add(point_data[1][@intCast(point_idx[1] + it.first_internal_point_shift[1])], types.neg(periodicity));

                    const cf = self.control_function.data[@intCast(global_idx_0)];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        cf.data[0],
                        cf.data[1],
                    );

                    // points inside block 0
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[0]] = stencil.get(.im1_jm1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[1]] = stencil.get(.i_jm1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[2]] = stencil.get(.ip1_jm1);

                    // points on boundary connection (taken from block 0)
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[3]] = stencil.get(.im1_j);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[4]] = stencil.get(.i_j);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[5]] = stencil.get(.ip1_j);

                    // points inside block 1
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[6]] = stencil.get(.im1_jp1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[7]] = stencil.get(.i_jp1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[8]] = stencil.get(.ip1_jp1);

                    // adjust RHS for periodicity
                    self.rhs_x[@intCast(global_idx_0)] = periodicity.data[0] * (stencil.get(.im1_jp1) + stencil.get(.i_jp1) + stencil.get(.ip1_jp1));
                    self.rhs_y[@intCast(global_idx_0)] = periodicity.data[1] * (stencil.get(.im1_jp1) + stencil.get(.i_jp1) + stencil.get(.ip1_jp1));
                }
            } else {
                while (it.next()) |connected_points| {
                    const point_idx: [2]c_int = .{ @intCast(connected_points[0].value), @intCast(connected_points[1].value) };

                    // TODO: this could be enhanced by computing this just once for the connection (plus +1 or -1 for the next connection point)
                    const global_idx_0 = @as(c_int, @intCast(self.row_idx_range_start_for_each_block[connection.ranges[0].block])) + point_idx[0];
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(global_idx_0);

                    const im1_j = point_data[0][@intCast(point_idx[0] - it.in_connection_direction_shift[0])];
                    const i_jm1 = point_data[0][@intCast(point_idx[0] + it.first_internal_point_shift[0])];
                    const ip1_j = point_data[0][@intCast(point_idx[0] + it.in_connection_direction_shift[0])];
                    const i_jp1 = point_data[1][@intCast(point_idx[1] + it.first_internal_point_shift[1])];

                    const cf = self.control_function.data[@intCast(global_idx_0)];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        cf.data[1],
                        cf.data[0],
                    );

                    // points inside block 0
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[0]] = stencil.get(.im1_jm1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[1]] = stencil.get(.i_jm1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[2]] = stencil.get(.ip1_jm1);

                    // points on boundary connection (taken from block 0)
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[3]] = stencil.get(.im1_j);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[4]] = stencil.get(.i_j);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[5]] = stencil.get(.ip1_j);

                    // points inside block 1
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[6]] = stencil.get(.im1_jp1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[7]] = stencil.get(.i_jp1);
                    self.lhs_values[row_non_zero_entries_start_idx + point_stencil_idx[8]] = stencil.get(.ip1_jp1);

                    // NOTE: RHS is already set in the point based loop.
                }
            }
        }
    }

    fn fill(self: *RowCompressedMatrixSystem2d, iteration: usize) !void {
        if (iteration > 0) {
            self.control_function.update(self.mesh.*);
        }
        self.fillBlockInternalPointData();
        self.fillBlockConnectionData();
    }

    pub fn fillXSpecific(self: *RowCompressedMatrixSystem2d) !void {
        for (self.mesh.boundary_conditions.items) |bc| {
            switch (bc.kind) {
                .inlet, .outlet => {
                    const first_internal_point_shift = bc.range.firstInternalPointShift(self.mesh);

                    var it = bc.range.iterate(self.mesh);
                    while (it.next()) |local_id| {
                        const block_and_local_id = BlockAndLocalIndex{ .block = bc.range.block, .local_idx = .{ .value = local_id } };
                        const local_id_2d = self.index_converter.index2d(block_and_local_id);
                        const boundary_id = try self.boundary_points.index_converter.bufferIndex(local_id_2d);
                        if (self.boundary_points.kind.buffer[boundary_id] != .sliding_circ) continue;

                        const global_id = self.index_converter.globalIndex(block_and_local_id).value;
                        const non_zero_range_start = self.nonZeroEntriesRangeStart(@intCast(global_id));

                        if (first_internal_point_shift > 0) {
                            self.lhs_values[non_zero_range_start] = 1.0;
                            self.lhs_values[non_zero_range_start + 1] = 0.0;
                        } else {
                            self.lhs_values[non_zero_range_start] = 0.0;
                            self.lhs_values[non_zero_range_start + 1] = 1.0;
                        }
                    }
                },
                else => unreachable,
            }
        }
    }

    pub fn fillYSpecific(self: *RowCompressedMatrixSystem2d) !void {
        for (self.mesh.boundary_conditions.items) |bc| {
            switch (bc.kind) {
                .inlet, .outlet => {
                    var it = bc.range.iterate(self.mesh);
                    while (it.next()) |local_id| {
                        const block_and_local_id = BlockAndLocalIndex{ .block = bc.range.block, .local_idx = .{ .value = local_id } };
                        const local_id_2d = self.index_converter.index2d(block_and_local_id);
                        const boundary_id = try self.boundary_points.index_converter.bufferIndex(local_id_2d);
                        if (self.boundary_points.kind.buffer[boundary_id] != .sliding_circ) continue;

                        const global_id = self.index_converter.globalIndex(block_and_local_id).value;
                        const non_zero_range_start = self.nonZeroEntriesRangeStart(@intCast(global_id));
                        self.lhs_values[non_zero_range_start] = 1.0;
                        self.lhs_values[non_zero_range_start + 1] = -1.0;
                    }
                },
                else => unreachable,
            }
        }
    }
};

const BlockBoundaryPointKind = enum {
    fixed,
    smoothed,
    connected, // perhaps save the global id of the point that it is connected to?
    laplacian_smoothed, // id to access the neighbor info
    sliding_circ,
};

fn BoundedArray(comptime T: type, comptime buffer_capacity: usize) type {
    return struct {
        const Self = @This();

        buffer: [buffer_capacity]T = undefined,
        len: usize = 0,

        fn slice(self: anytype) switch (@TypeOf(&self.buffer)) {
            *[buffer_capacity]T => []T,
            *const [buffer_capacity]T => []const T,
            else => unreachable,
        } {
            return self.buffer[0..self.len];
        }

        fn append(self: *Self, item: T) error{Overflow}!void {
            try self.ensureUnusedCapacity(1);
            self.buffer[self.len] = item;
            self.len += 1;
        }

        fn appendSlice(self: *Self, items: []const T) error{Overflow}!void {
            try self.ensureUnusedCapacity(items.len);
            const old_len = self.len;
            self.len += items.len;
            @memcpy(self.buffer[old_len..self.len], items[0..]);
        }

        fn ensureUnusedCapacity(self: Self, additional_count: usize) error{Overflow}!void {
            if (self.len + additional_count > buffer_capacity) return error.Overflow;
        }
    };
}

/// contains for each block boundary point all connected points in a flat array. The number of connections allows
/// to categorize each boundary point w.r.t. how to smooth it.
const BlockBoundaryPoints = struct {
    kind: boundary.PointData(BlockBoundaryPointKind),
    index_converter: boundary.PointDataBufferIndexConverter,
    laplacian_points: std.array_list.Managed(LaplacianPoint),

    // TODO: consider enhancing the periodicity handling of laplacian points by e.g. moving AoS to SoA

    const LaplacianPoint = struct {
        /// collection of the overlapping points; the first point is the global id of the point.
        overlapping_points: BoundedArray(OverlappingPoint, 4),

        /// stencil ids that are used in the system matrix; this includes the point itself.
        stencil_ids: BoundedArray(c_int, 6),

        /// RHS that potentially accounts for periodicity (count of periodic points * periodicity); defaults to 0;
        rhs: types.Vec2d,

        fn globalId(self: LaplacianPoint) usize {
            return @intCast(self.overlapping_points.slice()[0].global_id);
        }
    };

    fn init(allocator: std.mem.Allocator, index_converter: IndexConverter, mesh_data: *const discrete.Mesh) !BlockBoundaryPoints {
        var boundary_points = BlockBoundaryPoints{
            .kind = try .init(allocator, mesh_data),
            .index_converter = try .init(allocator, mesh_data),
            .laplacian_points = try initLaplacianPoints(allocator, index_converter, mesh_data),
        };
        errdefer boundary_points.deinit();

        // set kind (how to smooth the point)
        @memset(boundary_points.kind.buffer, .fixed);

        // add laplacian points
        for (boundary_points.laplacian_points.items) |laplacian_point| {
            {
                // set lowest point to laplacian smoothed
                const global_id = laplacian_point.overlapping_points.slice()[0].global_id;
                const local_id = index_converter.localIndex(.{ .value = global_id });
                const local_id_2d = index_converter.index2d(local_id);
                const buffer_id = try boundary_points.index_converter.bufferIndex(local_id_2d);
                boundary_points.kind.buffer[buffer_id] = .laplacian_smoothed;
            }

            // all other point are connected to the laplacian smoothed point
            for (laplacian_point.overlapping_points.slice()[1..]) |overlapping_point| {
                const local_id = index_converter.localIndex(.{ .value = overlapping_point.global_id });
                const local_id_2d = index_converter.index2d(local_id);
                const buffer_id = try boundary_points.index_converter.bufferIndex(local_id_2d);
                boundary_points.kind.buffer[buffer_id] = .connected;
            }
        }

        for (mesh_data.boundary_conditions.items) |boundary_condition| {
            switch (boundary_condition.kind) {
                .inlet, .outlet => {
                    var range_iterator = boundary_condition.range.iterate(mesh_data);
                    while (range_iterator.next()) |local_id| {
                        const local_id_2d = index_converter.index2d(.{ .block = boundary_condition.range.block, .local_idx = .{ .value = local_id } });
                        const buffer_id = try boundary_points.index_converter.bufferIndex(local_id_2d);
                        boundary_points.kind.buffer[buffer_id] = .sliding_circ;
                    }
                },
                .wall => {},
            }
        }

        // set other connections
        for (mesh_data.connections.items) |connection| {
            var connected_points = connection.iterate(mesh_data);

            // NOTE: endpoints are either laplacian smoothed fixed or sliding.

            // check 1st endpoint
            {
                const local_ids = connected_points.next().?;
                const buffer_ids = try computeBufferIds(.{
                    .{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0] } },
                    .{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1] } },
                }, index_converter, boundary_points.index_converter);

                // set 2nd point to connect if the 1st point is fixed or sliding (laplacian points are already set!)
                switch (boundary_points.kind.buffer[buffer_ids[0]]) {
                    .fixed, .sliding_circ => {
                        boundary_points.kind.buffer[buffer_ids[1]] = .connected;
                    },
                    else => {},
                }
            }

            // middle of the connection is smoothed (1st) and connected (2nd to 1st)
            for (0..connected_points.data[0].count) |_| {
                const local_ids = connected_points.next().?;
                const buffer_ids = try computeBufferIds(.{
                    .{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0] } },
                    .{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1] } },
                }, index_converter, boundary_points.index_converter);
                boundary_points.kind.buffer[buffer_ids[0]] = .smoothed;
                boundary_points.kind.buffer[buffer_ids[1]] = .connected;
            }

            // check 2nd endpoint
            {
                const local_ids = connected_points.next().?;
                const buffer_ids = try computeBufferIds(.{
                    .{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0] } },
                    .{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1] } },
                }, index_converter, boundary_points.index_converter);

                // set 2nd point to connect if the 1st point is fixed or sliding (laplacian points are already set!)
                switch (boundary_points.kind.buffer[buffer_ids[0]]) {
                    .fixed, .sliding_circ => {
                        boundary_points.kind.buffer[buffer_ids[1]] = .connected;
                    },
                    else => {},
                }
            }
        }

        return boundary_points;
    }

    const OverlappingPoint = struct {
        global_id: usize,
        periodicity: types.Vec2d,
    };

    /// returns the necessary laplacian smoothed point data.
    fn initLaplacianPoints(
        allocator: std.mem.Allocator,
        index_converter: IndexConverter,
        mesh_data: *const discrete.Mesh,
    ) !std.array_list.Managed(LaplacianPoint) {
        // collect connection endpoint ids into a flat array (start_00, start_01, end_00, end_01, start_10, start_11, end_10, end_11, ...)
        var endpoint_ids = try allocator.alloc(usize, mesh_data.connections.items.len * 4);
        defer allocator.free(endpoint_ids);

        for (mesh_data.connections.items, 0..) |connection, connection_id| {
            const local_ids = .{ connection.ranges[0].endpoints(mesh_data), connection.ranges[1].endpoints(mesh_data) };
            const range_start = connection_id * 4;
            endpoint_ids[range_start] = index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0][0] } }).value;
            endpoint_ids[range_start + 1] = index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1][0] } }).value;
            endpoint_ids[range_start + 2] = index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0][1] } }).value;
            endpoint_ids[range_start + 3] = index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1][1] } }).value;
        }

        // std.debug.print("endpoint_ids: {any}\n", .{endpoint_ids});

        var laplacian_points = try std.array_list.Managed(LaplacianPoint).initCapacity(allocator, mesh_data.connections.items.len * 2);
        errdefer laplacian_points.deinit();

        // linear search for identical ids: if there are identical ids, these are junction points!
        for (0..endpoint_ids.len - 1) |endpoint_id| {
            const endpoint = endpoint_ids[endpoint_id];
            for (endpoint_id + 1..endpoint_ids.len) |endpoint_id_to_check| {
                const endpoint_other = endpoint_ids[endpoint_id_to_check];
                if (endpoint == endpoint_other) {

                    // see, if we need to merge it with an already existing laplacian point.
                    var exisitng_point_found = false;
                    for (laplacian_points.items) |*laplacian_point| {
                        for (laplacian_point.overlapping_points.slice()) |point| {
                            if (point.global_id == endpoint) {
                                exisitng_point_found = true;

                                // add the ID that is connected to the match
                                const endpoint_id_to_add = if (endpoint_id_to_check % 2 == 0) endpoint_id_to_check + 1 else endpoint_id_to_check - 1;

                                // get the periodicity
                                const connection_id = @divTrunc(endpoint_id_to_add, 4);
                                const periodicity = if (mesh_data.connections.items[connection_id].periodicity) |periodicity| periodicity else types.Vec2d.init(0, 0);

                                try appendIfUnique(&laplacian_point.overlapping_points, endpoint_ids[endpoint_id_to_add], periodicity);
                            }
                        }
                    }

                    if (!exisitng_point_found) {
                        const points = .{ @divTrunc(endpoint_id, 2), @divTrunc(endpoint_id_to_check, 2) };
                        std.debug.assert(points[0] != points[1]);

                        var laplacian_point = BoundedArray(OverlappingPoint, 4){};

                        // add first point
                        {
                            const connection_id = points[0] / 2;
                            const periodicity = if (mesh_data.connections.items[connection_id].periodicity) |periodicity| periodicity else types.Vec2d.init(0, 0);
                            try laplacian_point.appendSlice(&.{
                                .{ .global_id = endpoint_ids[points[0] * 2], .periodicity = .init(0, 0) },
                                .{ .global_id = endpoint_ids[points[0] * 2 + 1], .periodicity = periodicity },
                            });

                            std.debug.assert(laplacian_point.slice()[0].global_id != laplacian_point.slice()[1].global_id);

                            std.debug.assert(blk: {
                                const local_ids = .{
                                    index_converter.localIndex(.{ .value = laplacian_point.slice()[0].global_id }),
                                    index_converter.localIndex(.{ .value = laplacian_point.slice()[1].global_id }),
                                };

                                const local_ids_2d = .{ index_converter.index2d(local_ids[0]), index_converter.index2d(local_ids[1]) };

                                const coords = .{
                                    mesh_data.blocks.items[local_ids[0].block].points.getIndex(local_ids_2d[0].point),
                                    mesh_data.blocks.items[local_ids[1].block].points.getIndex(local_ids_2d[1].point),
                                };

                                break :blk types.eqlApprox(coords[0], coords[1], 1e-12);
                            });
                        }

                        // add 2nd point
                        {
                            const connection_id = points[1] / 2;
                            const periodicity = if (mesh_data.connections.items[connection_id].periodicity) |periodicity| periodicity else types.Vec2d.init(0, 0);
                            try appendIfUnique(&laplacian_point, endpoint_ids[points[1] * 2], periodicity);
                            try appendIfUnique(&laplacian_point, endpoint_ids[points[1] * 2 + 1], periodicity);
                        }

                        try laplacian_points.append(.{
                            .overlapping_points = laplacian_point,
                            .stencil_ids = .{},
                            .rhs = .init(0, 0),
                        });
                    }
                }
            }
        }

        // in place sort the laplacian points
        for (laplacian_points.items) |*laplacian_point| {
            std.mem.sort(OverlappingPoint, laplacian_point.overlapping_points.slice(), {}, comptime struct {
                pub fn inner(_: void, a: OverlappingPoint, b: OverlappingPoint) bool {
                    return a.global_id < b.global_id;
                }
            }.inner);
        }

        // sort laplacian points in order of the lowest id
        std.mem.sort(LaplacianPoint, laplacian_points.items, {}, struct {
            fn inner(_: void, a: LaplacianPoint, b: LaplacianPoint) bool {
                return a.overlapping_points.slice()[0].global_id < b.overlapping_points.slice()[0].global_id;
            }
        }.inner);

        // set the stencil ids
        for (laplacian_points.items) |*laplacian_point| {
            try laplacian_point.stencil_ids.append(@intCast(laplacian_point.globalId())); // add laplacian point to stencil

            // add neighboring points
            for (laplacian_point.overlapping_points.slice()) |overlapping_point| {
                const local_id = index_converter.localIndex(.{ .value = overlapping_point.global_id });
                const local_2d = index_converter.index2d(local_id);

                // check if it is a corner or an internal point
                const size = mesh_data.blocks.items[local_2d.block].points.size;

                const internal_points: [2]?types.Index2d = blk: {
                    if (local_2d.point[0] == 0) {
                        // j_min
                        if (local_2d.point[1] == 0) {
                            break :blk .{ .{ 1, 1 }, null };
                        } else if (local_2d.point[1] == size[1] - 1) {
                            break :blk .{ .{ 1, size[1] - 2 }, null };
                        } else {
                            break :blk .{ .{ 1, local_2d.point[1] - 1 }, .{ 1, local_2d.point[1] + 1 } };
                        }
                    } else if (local_2d.point[0] == size[0] - 1) {
                        // j_max
                        if (local_2d.point[1] == 0) {
                            break :blk .{ .{ size[0] - 2, 1 }, null };
                        } else if (local_2d.point[1] == size[1] - 1) {
                            break :blk .{ .{ size[0] - 2, size[1] - 2 }, null };
                        } else {
                            break :blk .{ .{ size[0] - 2, local_2d.point[1] - 1 }, .{ size[0] - 2, local_2d.point[1] + 1 } };
                        }
                    } else {
                        std.debug.assert(local_2d.point[1] == 0 or local_2d.point[1] == size[1] - 1);
                        if (local_2d.point[1] == 0) {
                            break :blk .{ .{ local_2d.point[0] - 1, 1 }, .{ local_2d.point[0] + 1, 1 } };
                        } else if (local_2d.point[1] == size[1] - 1) {
                            break :blk .{ .{ local_2d.point[0] - 1, local_2d.point[1] - 1 }, .{ local_2d.point[0] + 1, local_2d.point[1] - 1 } };
                        } else {
                            unreachable; // it is not possible to land here!
                        }
                    }
                };

                for (internal_points) |point| {
                    if (point) |p| {
                        const global_id = index_converter.globalIndex(index_converter.index(.{ .block = local_id.block, .point = p }));
                        try laplacian_point.stencil_ids.append(@intCast(global_id.value));
                        laplacian_point.rhs.add(overlapping_point.periodicity);
                    }
                }
            }

            // sort the ids
            std.mem.sort(c_int, laplacian_point.stencil_ids.slice(), {}, comptime std.sort.asc(c_int));
        }

        return laplacian_points;
    }

    fn appendIfUnique(overlapping_points: *BoundedArray(OverlappingPoint, 4), id: usize, periodicity: types.Vec2d) !void {
        for (overlapping_points.slice()) |overlapping_point| {
            if (overlapping_point.global_id == id) return;
        }

        try overlapping_points.append(.{ .global_id = id, .periodicity = periodicity });
    }

    fn deinit(self: BlockBoundaryPoints) void {
        self.kind.deinit();
        self.index_converter.deinit();
        self.laplacian_points.deinit();
    }
};

pub const RangeFillMatrixIterator = struct {
    count: usize,
    first_internal_point_shift: [2]c_int,
    in_connection_direction_shift: [2]c_int,
    position: [2]LocalIndex,

    pub fn next(self: *RangeFillMatrixIterator) ?[2]LocalIndex {
        if (self.count == 0) return null;
        const position = self.position;

        self.count -= 1;

        self.position = .{
            .{ .value = @intCast(@as(isize, @intCast(self.position[0].value)) + self.in_connection_direction_shift[0]) },
            .{ .value = @intCast(@as(isize, @intCast(self.position[1].value)) + self.in_connection_direction_shift[1]) },
        };

        return position;
    }

    fn limitToRangeInternalPoints(self: *RangeFillMatrixIterator) void {
        _ = self.next();
        self.count -= 1;
    }

    pub fn init(
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
                data.count = start - end + 1;
            } else {
                data.count = end - start + 1;
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
    mesh: *const discrete.Mesh,

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
            .mesh = mesh_data,
        };
    }

    fn deinit(self: IndexConverter) void {
        self.allocator.free(self.global_point_index_range_start);
    }

    fn globalIndex(self: IndexConverter, block_local_idx: BlockAndLocalIndex) GlobalIndex {
        return .{ .value = self.global_point_index_range_start[block_local_idx.block].value + block_local_idx.local_idx.value };
    }

    fn localIndex(self: IndexConverter, global_idx: GlobalIndex) BlockAndLocalIndex {
        var block_idx: usize = self.global_point_index_range_start.len - 1;
        while (global_idx.value < self.global_point_index_range_start[block_idx].value) : (block_idx -= 1) {}
        const local_idx = global_idx.value - self.global_point_index_range_start[block_idx].value;
        return .{ .block = block_idx, .local_idx = .{ .value = local_idx } };
    }

    fn index2d(self: IndexConverter, block_and_local_idx: BlockAndLocalIndex) types.MeshIndex2d {
        const block = self.mesh.blocks.items[block_and_local_idx.block];
        const point_idx = blk: {
            const point_i = @divTrunc(block_and_local_idx.local_idx.value, block.points.size[1]);
            const point_j = block_and_local_idx.local_idx.value - point_i * block.points.size[1];
            break :blk types.Index2d{ point_i, point_j };
        };
        return .{ .block = block_and_local_idx.block, .point = point_idx };
    }

    fn index(self: IndexConverter, point: types.MeshIndex2d) BlockAndLocalIndex {
        const size = self.mesh.blocks.items[point.block].points.size;
        return .{ .block = point.block, .local_idx = .{ .value = point.point[0] * size[1] + point.point[1] } };
    }
};

fn computeBufferIds(local_ids: [2]BlockAndLocalIndex, point_index_converter: IndexConverter, boundary_index_converter: boundary.PointDataBufferIndexConverter) ![2]usize {
    const local2d = .{ point_index_converter.index2d(local_ids[0]), point_index_converter.index2d(local_ids[1]) };
    return .{ try boundary_index_converter.bufferIndex(local2d[0]), try boundary_index_converter.bufferIndex(local2d[1]) };
}
