const std = @import("std");
const discrete = @import("discrete.zig");
const types = @import("types.zig");
const umfpack = @import("umfpack.zig");
const boundary = @import("boundary.zig");
const tfi = @import("tfi.zig");

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
                        const delta_x = system.x_new[row_idx] - block.points.data[point_idx].data[0];
                        const delta_y = system.y_new[row_idx] - block.points.data[point_idx].data[1];
                        block.points.data[point_idx].data[0] += 1.0 * delta_x;
                        block.points.data[point_idx].data[1] += 1.0 * delta_y;
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

    boundary_points: BlockBoundaryPoints,

    control_function: ControlFunction,

    fn init(allocator: std.mem.Allocator, mesh_data: *discrete.Mesh) !RowCompressedMatrixSystem2d {
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
        defer index_converter.deinit();

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
            .boundary_points = try .init(allocator, index_converter, mesh_data),
            .control_function = try .init(allocator, dof, mesh_data),
        };

        try system.initNonZeroMatrixEntries(dof, non_zero_entries_capacity, index_converter);
        system.initBoundaryData();

        return system;
    }

    fn deinit(self: RowCompressedMatrixSystem2d) void {
        self.allocator.free(self.buffer_int);
        self.allocator.free(self.buffer_float);
        self.allocator.free(self.row_idx_range_start_for_each_block);
        self.boundary_points.deinit();
        self.control_function.deinit();
    }

    // fn initControlFunctionKhamaysehEtAl(self: *RowCompressedMatrixSystem2d) !void {
    //     const buf_size = types.Index2d{
    //         @max(self.mesh.blocks.items[0].points.size[0], self.mesh.blocks.items[1].points.size[0]),
    //         @max(self.mesh.blocks.items[0].points.size[1], self.mesh.blocks.items[1].points.size[1]),
    //     };
    //     var edge_buf = try self.allocator.alloc(types.Vec2d, 2 * buf_size[0] + 2 * buf_size[1]);
    //     defer self.allocator.free(edge_buf);
    //
    //     var block_range_start: usize = 0;
    //     for (self.mesh.blocks.items[0..2]) |block| {
    //         const size = block.points.size;
    //
    //         @memset(edge_buf[0 .. 2 * size[0] + 2 * size[1]], .{ .data = .{ 0, 0 } });
    //
    //         const edge_i_min = edge_buf[0..size[0]];
    //         const edge_i_max = edge_buf[size[0] .. 2 * size[0]];
    //         const edge_j_min = edge_buf[2 * size[0] .. 2 * size[0] + size[1]];
    //         const edge_j_max = edge_buf[2 * size[0] + size[1] .. 2 * size[0] + 2 * size[1]];
    //
    //         // NOTE: hard coded for i_min right now!
    //
    //         // edge i_min
    //         {
    //             var local_id: usize = size[1];
    //             for (1..block.points.size[0] - 1) |edge_id| {
    //                 const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
    //                 const x_i_0, const y_i_0 = block.points.data[local_id].data;
    //                 const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
    //                 const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;
    //
    //                 // central differences
    //                 const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
    //                 const y_xi = 0.5 * (y_ip1_0 - y_im1_0);
    //                 const x_xi2 = x_ip1_0 - 2.0 * x_i_0 + x_im1_0;
    //                 const y_xi2 = y_ip1_0 - 2.0 * y_i_0 + y_im1_0;
    //                 const g11 = x_xi * x_xi + y_xi * y_xi;
    //
    //                 // ghost point computation
    //
    //                 // one sided approx normal to wall
    //                 const x_eta_0 = x_i_1 - x_i_0;
    //                 const y_eta_0 = y_i_1 - y_i_0;
    //
    //                 // eq. 6.13
    //                 const factor = (-y_xi * x_eta_0 + x_xi * y_eta_0) / g11;
    //                 const x_eta_fixed = factor * -y_xi;
    //                 const y_eta_fixed = factor * x_xi;
    //
    //                 const x_ghost = x_i_0 - x_eta_fixed;
    //                 const y_ghost = y_i_0 - y_eta_fixed;
    //                 // TODO: save ghost point positions
    //
    //                 const g22 = x_eta_fixed * x_eta_fixed + y_eta_fixed * y_eta_fixed;
    //
    //                 const x_eta2 = x_ghost - 2.0 * x_i_0 + x_i_1;
    //                 const y_eta2 = y_ghost - 2.0 * y_i_0 + y_i_1;
    //
    //                 const x_eta = x_i_1 - x_i_0;
    //                 const y_eta = y_i_1 - y_i_0;
    //
    //                 // eq. 6.10
    //                 const P = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
    //                 const Q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;
    //
    //                 edge_i_min[edge_id] = .init(P, Q);
    //
    //                 local_id += size[1];
    //             }
    //         }
    //
    //         const num_points = block.points.size[0] * block.points.size[1];
    //         try tfi.linear2d(self.control_function[block_range_start .. block_range_start + num_points], edge_i_min, edge_i_max, edge_j_min, edge_j_max);
    //         block_range_start += num_points;
    //     }
    // }

    fn initControlFunctionThomasAndMiddlecoff(self: *RowCompressedMatrixSystem2d) !void {

        // TODO: we assume the first 2 blocks to be the O blocks
        // TODO: for now, we just compute on the j edges

        // TODO: introduce global buffer to reduce number of allocations.
        const buf_size = types.Index2d{
            @max(self.mesh.blocks.items[0].points.size[0], self.mesh.blocks.items[1].points.size[0]),
            @max(self.mesh.blocks.items[0].points.size[1], self.mesh.blocks.items[1].points.size[1]),
        };
        var edge_buf = try self.allocator.alloc(types.Vec2d, 2 * buf_size[0] + 2 * buf_size[1]);
        defer self.allocator.free(edge_buf);

        var block_range_start: usize = 0;
        for (self.mesh.blocks.items[0..2]) |block| {
            const size = block.points.size;

            @memset(edge_buf[0 .. 2 * size[0] + 2 * size[1]], .{ .data = .{ 0, 0 } });

            const edge_i_min = edge_buf[0..size[0]];
            const edge_i_max = edge_buf[size[0] .. 2 * size[0]];
            const edge_j_min = edge_buf[2 * size[0] .. 2 * size[0] + size[1]];
            const edge_j_max = edge_buf[2 * size[0] + size[1] .. 2 * size[0] + 2 * size[1]];

            // edge i_min
            {
                var local_id: usize = size[1];
                for (1..block.points.size[0] - 1) |edge_id| {
                    const im1_j = block.points.data[local_id - size[1]];
                    const i_j = block.points.data[local_id];
                    const ip1_j = block.points.data[local_id + size[1]];

                    std.debug.assert(local_id - size[1] == block.points.index(.{ edge_id - 1, 0 }));
                    std.debug.assert(local_id == block.points.index(.{ edge_id, 0 }));
                    std.debug.assert(local_id + size[1] == block.points.index(.{ edge_id + 1, 0 }));

                    // central differences
                    const x_xi = types.scale(0.5, types.sub(ip1_j, im1_j));
                    const x_xi2 = types.sub(types.add(ip1_j, im1_j), types.scale(2.0, i_j));

                    // eq. 11
                    const phi = -(x_xi.data[0] * x_xi2.data[0] + x_xi.data[1] * x_xi2.data[1]) / (x_xi.data[0] * x_xi.data[0] + x_xi.data[1] * x_xi.data[1]);
                    const psi = 0.0;

                    edge_i_min[edge_id] = .init(phi, psi);

                    local_id += size[1];
                }
            }

            // edge i_max
            {
                var local_id: usize = 2 * size[1] - 1;
                for (1..block.points.size[0] - 1) |edge_id| {
                    const im1_j = block.points.data[local_id - size[1]];
                    const i_j = block.points.data[local_id];
                    const ip1_j = block.points.data[local_id + size[1]];

                    std.debug.assert(local_id - size[1] == block.points.index(.{ edge_id - 1, size[1] - 1 }));
                    std.debug.assert(local_id == block.points.index(.{ edge_id, size[1] - 1 }));
                    std.debug.assert(local_id + size[1] == block.points.index(.{ edge_id + 1, size[1] - 1 }));

                    // central differences
                    const x_xi = types.scale(0.5, types.sub(ip1_j, im1_j));
                    const x_xi2 = types.sub(types.add(ip1_j, im1_j), types.scale(2.0, i_j));

                    // eq. 11
                    const phi = -(x_xi.data[0] * x_xi2.data[0] + x_xi.data[1] * x_xi2.data[1]) / (x_xi.data[0] * x_xi.data[0] + x_xi.data[1] * x_xi.data[1]);
                    const psi = 0.0;

                    edge_i_max[edge_id] = .init(phi, psi);

                    local_id += size[1];
                }
            }

            // edge j_min
            {
                var local_id: usize = 1;
                for (1..block.points.size[1] - 1) |edge_id| {
                    const i_jm1 = block.points.data[local_id - 1];
                    const i_j = block.points.data[local_id];
                    const i_jp1 = block.points.data[local_id + 1];

                    std.debug.assert(local_id - 1 == block.points.index(.{ 0, edge_id - 1 }));
                    std.debug.assert(local_id == block.points.index(.{ 0, edge_id }));
                    std.debug.assert(local_id + 1 == block.points.index(.{ 0, edge_id + 1 }));

                    // central differences
                    const x_eta = types.scale(0.5, types.sub(i_jp1, i_jm1));
                    const x_eta2 = types.sub(types.add(i_jp1, i_jm1), types.scale(2.0, i_j));

                    // eq. 11
                    const phi = 0.0;
                    const psi = -(x_eta.data[0] * x_eta2.data[0] + x_eta.data[1] * x_eta2.data[1]) / (x_eta.data[0] * x_eta.data[0] + x_eta.data[1] * x_eta.data[1]);

                    edge_j_min[edge_id] = .init(phi, psi);

                    local_id += 1;
                }
            }

            // edge j_max
            {
                var local_id = (block.points.size[0] - 1) * block.points.size[1] + 1;
                for (1..block.points.size[1] - 1) |edge_id| {
                    const i_jm1 = block.points.data[local_id - 1];
                    const i_j = block.points.data[local_id];
                    const i_jp1 = block.points.data[local_id + 1];

                    std.debug.assert(local_id - 1 == block.points.index(.{ size[0] - 1, edge_id - 1 }));
                    std.debug.assert(local_id == block.points.index(.{ size[0] - 1, edge_id }));
                    std.debug.assert(local_id + 1 == block.points.index(.{ size[0] - 1, edge_id + 1 }));

                    // central differences
                    const x_eta = types.scale(0.5, types.sub(i_jp1, i_jm1));
                    const x_eta2 = types.sub(types.add(i_jp1, i_jm1), types.scale(2.0, i_j));

                    // eq. 11
                    const phi = 0.0;
                    const psi = -(x_eta.data[0] * x_eta2.data[0] + x_eta.data[1] * x_eta2.data[1]) / (x_eta.data[0] * x_eta.data[0] + x_eta.data[1] * x_eta.data[1]);

                    edge_j_max[edge_id] = .init(phi, psi);

                    local_id += 1;
                }
            }

            const num_points = block.points.size[0] * block.points.size[1];
            try tfi.linear2d(self.control_function[block_range_start .. block_range_start + num_points], edge_i_min, edge_i_max, edge_j_min, edge_j_max);
            block_range_start += num_points;
        }
    }

    /// returns the non zero entries range for the given row (the row index range start for each block can be retrieved from an array in the struct)
    fn nonZeroEntriesRangeStart(self: RowCompressedMatrixSystem2d, row_idx: c_int) usize {
        return @intCast(self.lhs_p[@intCast(row_idx)]);
    }

    fn initNonZeroMatrixEntriesForBoundaryPoint(
        non_zero_entries: *std.ArrayList(c_int),
        row_count: *std.ArrayList(c_int),
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
        }
        boundary_point_id.* += 1;

        try row_count.append(@intCast(non_zero_entries.items.len));
        row_idx.* += 1;
    }

    fn initNonZeroMatrixEntriesPointBased(
        self: @This(),
        non_zero_entries: *std.ArrayList(c_int),
        row_count: *std.ArrayList(c_int),
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
        non_zero_entries: *std.ArrayList(c_int),
        index_converter: IndexConverter,
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

            try self.initNonZeroMatrixForConnectionEndpoint(it.next().?, connection, non_zero_entries, index_converter);

            // middle of the connection
            for (0..it.count) |_| {
                const connected_points_local_idx = it.next().?;
                const connected_points_global_idx = [2]c_int{
                    @intCast(index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = connected_points_local_idx[0] }).value),
                    @intCast(index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = connected_points_local_idx[1] }).value),
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
                                std.debug.print("wrong ordering of connection stencil data for point {d}: {any}\n", .{ connected_points_global_idx[0], non_zero_entries.items[row_non_zero_entries_start_idx .. row_non_zero_entries_start_idx + 9] });
                                break :blk false;
                            }
                        }
                        break :blk true;
                    });
                }
            }

            try self.initNonZeroMatrixForConnectionEndpoint(it.position, connection, non_zero_entries, index_converter);
        }
    }

    fn initNonZeroMatrixForConnectionEndpoint(
        self: @This(),
        local_ids: [2]LocalIndex,
        connection: boundary.Connection,
        non_zero_entries: *std.ArrayList(c_int),
        index_converter: IndexConverter,
    ) !void {
        const local_id_2d = index_converter.index2d(.{ .block = connection.ranges[0].block, .local_idx = local_ids[0] });
        const buffer_id = try self.boundary_points.index_converter.bufferIndex(local_id_2d);

        switch (self.boundary_points.kind.buffer[buffer_id]) {
            .fixed => {
                // connect the 2nd point to the 1st
                const connected_points_global_idx = [2]c_int{
                    @intCast(index_converter.globalIndex(.{ .block = connection.ranges[0].block, .local_idx = local_ids[0] }).value),
                    @intCast(index_converter.globalIndex(.{ .block = connection.ranges[1].block, .local_idx = local_ids[1] }).value),
                };

                std.debug.assert(connected_points_global_idx[0] < connected_points_global_idx[1]);
                const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(connected_points_global_idx[1]);
                non_zero_entries.items[row_non_zero_entries_start_idx] = connected_points_global_idx[0];
                non_zero_entries.items[row_non_zero_entries_start_idx + 1] = connected_points_global_idx[1];
            },
            .laplacian_smoothed => {}, // already set
            else => undefined,
        }
    }

    fn initNonZeroMatrixEntries(
        self: *RowCompressedMatrixSystem2d,
        dof: usize,
        non_zero_entries_capacity: usize,
        index_converter: IndexConverter,
    ) !void {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(self.lhs_i[0..]));
        var non_zero_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

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

        try self.initNonZeroMatrixEntriesConnectionBased(&non_zero_entries, index_converter);
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

                // TODO: change this to a while loop!

                for (0..it.count) |_| {
                    const connected_points = it.next().?;
                    const global_idx_1 = self.row_idx_range_start_for_each_block[connection.ranges[1].block] + connected_points[1].value;

                    // adjust RHS for periodicity for the connected point
                    self.rhs_x[global_idx_1] = -periodicity.data[0];
                    self.rhs_y[global_idx_1] = -periodicity.data[1];
                }

                // TODO: fix the RangeFillMatrixIterator to also account for the 2nd endpoint!
                // Perhaps also create an internal iterator version and something that just gives the endpoints.

                // account for 2nd endpoint
                {
                    const connected_points = it.position;
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

                    // TODO: how can we effectively handle 0 w/o data access?
                    const control_function = self.control_function.data[row_idx];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        control_function.data[0],
                        control_function.data[1],
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

            _ = it.next();

            if (connection.periodicity) |periodicity| {
                for (0..it.count) |_| {
                    const connected_points = it.next().?;
                    const point_idx: [2]c_int = .{ @intCast(connected_points[0].value), @intCast(connected_points[1].value) };

                    // TODO: this could be enhanced by computing this just once for the connection (plus +1 or -1 for the next connection point)
                    const global_idx_0 = @as(c_int, @intCast(self.row_idx_range_start_for_each_block[connection.ranges[0].block])) + point_idx[0];
                    const row_non_zero_entries_start_idx = self.nonZeroEntriesRangeStart(global_idx_0);

                    const im1_j = point_data[0][@intCast(point_idx[0] - it.in_connection_direction_shift[0])];
                    const i_jm1 = point_data[0][@intCast(point_idx[0] + it.first_internal_point_shift[0])];
                    const ip1_j = point_data[0][@intCast(point_idx[0] + it.in_connection_direction_shift[0])];
                    const i_jp1 = types.add(point_data[1][@intCast(point_idx[1] + it.first_internal_point_shift[1])], .{ .data = .{ 0, -0.08836 } });

                    const control_function = self.control_function.data[@intCast(global_idx_0)];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        control_function.data[0],
                        control_function.data[1],
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

                    const control_function = self.control_function.data[@intCast(global_idx_0)];
                    const stencil = StencilData.init(
                        im1_j,
                        ip1_j,
                        i_jm1,
                        i_jp1,
                        control_function.data[0],
                        control_function.data[1],
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

    fn fillAndSolve(self: *RowCompressedMatrixSystem2d) !void {
        // try self.control_function.algorithm.update(self.control_function.data);
        self.fillBlockInternalPointData();
        self.fillBlockConnectionData();

        const dof = self.rhs_x.len;
        try umfpack.solve2(@intCast(dof), @intCast(dof), self.lhs_p, self.lhs_i, self.lhs_values, self.rhs_x, self.x_new, self.rhs_y, self.y_new);
    }
};

const BlockBoundaryPointKind = enum {
    fixed,
    smoothed,
    connected, // perhaps save the global id of the point that it is connected to?
    laplacian_smoothed, // id to access the neighbor info
};

/// contains for each block boundary point all connected points in a flat array. The number of connections allows
/// to categorize each boundary point w.r.t. how to smooth it.
const BlockBoundaryPoints = struct {
    kind: boundary.PointData(BlockBoundaryPointKind),
    index_converter: boundary.PointDataBufferIndexConverter,
    laplacian_points: std.ArrayList(LaplacianPoint),

    // TODO: consider enhancing the periodicity handling of laplacian points by e.g. moving AoS to SoA

    const LaplacianPoint = struct {
        /// collection of the overlapping points; the first point is the global id of the point.
        overlapping_points: std.BoundedArray(OverlappingPoint, 4),

        /// stencil ids that are used in the system matrix; this includes the point itself.
        stencil_ids: std.BoundedArray(c_int, 6),

        /// RHS that potentially accounts for periodicity (count of periodic points * periodicity); defaults to 0;
        rhs: types.Vec2d,

        fn globalId(self: LaplacianPoint) usize {
            return @intCast(self.overlapping_points.get(0).global_id);
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
                const global_id = laplacian_point.overlapping_points.get(0).global_id;
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

        // set other connections
        for (mesh_data.connections.items) |connection| {
            var connected_points = connection.iterate(mesh_data);

            // NOTE: endpoints are either laplacian smoothed or fixed.

            // check 1st endpoint
            {
                const local_ids = connected_points.next().?;
                const buffer_ids = try computeBufferIds(.{
                    .{ .block = connection.ranges[0].block, .local_idx = .{ .value = local_ids[0] } },
                    .{ .block = connection.ranges[1].block, .local_idx = .{ .value = local_ids[1] } },
                }, index_converter, boundary_points.index_converter);

                // set 2nd point to connect if the 1st point is fixed (laplacian points are already set!)
                if (boundary_points.kind.buffer[buffer_ids[0]] == .fixed) {
                    boundary_points.kind.buffer[buffer_ids[1]] = .connected;
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

                // set 2nd point to connect if the 1st point is fixed (laplacian points are already set!)
                if (boundary_points.kind.buffer[buffer_ids[0]] == .fixed) {
                    boundary_points.kind.buffer[buffer_ids[1]] = .connected;
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
    ) !std.ArrayList(LaplacianPoint) {
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

        var laplacian_points = try std.ArrayList(LaplacianPoint).initCapacity(allocator, mesh_data.connections.items.len * 2);
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

                        var laplacia_point = try std.BoundedArray(OverlappingPoint, 4).init(0);

                        // add first point
                        {
                            const connection_id = points[0] / 2;
                            const periodicity = if (mesh_data.connections.items[connection_id].periodicity) |periodicity| periodicity else types.Vec2d.init(0, 0);
                            try laplacia_point.appendSlice(&.{
                                .{ .global_id = endpoint_ids[points[0] * 2], .periodicity = .init(0, 0) },
                                .{ .global_id = endpoint_ids[points[0] * 2 + 1], .periodicity = periodicity },
                            });

                            std.debug.assert(laplacia_point.get(0).global_id != laplacia_point.get(1).global_id);

                            std.debug.assert(blk: {
                                const local_ids = .{
                                    index_converter.localIndex(.{ .value = laplacia_point.get(0).global_id }),
                                    index_converter.localIndex(.{ .value = laplacia_point.get(1).global_id }),
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
                            try appendIfUnique(&laplacia_point, endpoint_ids[points[1] * 2], periodicity);
                            try appendIfUnique(&laplacia_point, endpoint_ids[points[1] * 2 + 1], periodicity);
                        }

                        try laplacian_points.append(.{ .overlapping_points = laplacia_point, .stencil_ids = try .init(0), .rhs = .init(0, 0) });
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
                return a.overlapping_points.get(0).global_id < b.overlapping_points.get(0).global_id;
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

        // for (laplacian_points.items) |laplacian_point| {
        //     std.debug.print("lp: ", .{});
        //     for (laplacian_point.overlapping_points.slice()) |overlapping_point| {
        //         std.debug.print("({}, {} {}) ", .{ overlapping_point.global_id, overlapping_point.periodicity.data[0], overlapping_point.periodicity.data[1] });
        //     }
        //     std.debug.print("\nstencil: {any}\n", .{laplacian_point.stencil_ids.slice()});
        //     std.debug.print("rhs: {any}\n\n", .{laplacian_point.rhs.data});
        // }

        return laplacian_points;
    }

    fn appendIfUnique(overlapping_points: *std.BoundedArray(OverlappingPoint, 4), id: usize, periodicity: types.Vec2d) !void {
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

const ControlFunction = struct {
    allocator: std.mem.Allocator,
    data: []types.Vec2d,
    // algorithm: KhamaysehEtAl,
    algorithm: White,

    fn init(allocator: std.mem.Allocator, dof: usize, mesh_data: *const discrete.Mesh) !ControlFunction {
        // TODO: merge with other allocations.
        var data = try allocator.alloc(types.Vec2d, dof);
        @memset(data[0..], .{ .data = .{ 0, 0 } });
        return .{
            .allocator = allocator,
            .data = data,
            // .algorithm = try .init(allocator, mesh_data, data),
            .algorithm = .init(mesh_data, data),
        };
    }

    fn deinit(self: ControlFunction) void {
        self.allocator.free(self.data);
        self.algorithm.deinit();
    }

    const StegerSorenson = struct {
        mesh: *const discrete.Mesh,

        fn init(mesh_data: *const discrete.Mesh, control_function: []types.Vec2d) White {
            var white = White{ .mesh = mesh_data };
            white.initControlFunction(control_function);
            return white;
        }

        fn deinit(self: White) void {
            _ = self;
        }

        fn initControlFunction(self: *@This(), control_function: []types.Vec2d) void {
            var block_range_start: usize = 0;
            for (self.mesh.blocks.items[0..2]) |block| {
                const size = block.points.size;

                var local_id: usize = 0;

                {
                    // corner 0,0
                    // TODO: hard code indices to gain performance.
                    const x_0_0, const y_0_0 = block.points.data[local_id].data;
                    const x_0_1, const y_0_1 = block.points.data[local_id + 1].data;
                    const x_0_2, const y_0_2 = block.points.data[local_id + 2].data;
                    const x_1_0, const y_1_0 = block.points.data[local_id + size[1]].data;
                    const x_2_0, const y_2_0 = block.points.data[local_id + 2 * size[1]].data;

                    // forward differences
                    const x_xi = -x_0_0 + x_1_0;
                    const y_xi = -y_0_0 + y_1_0;
                    const x_xi2 = x_0_0 - 2 * x_1_0 + x_2_0;
                    const y_xi2 = y_0_0 - 2 * y_1_0 + y_2_0;

                    const x_eta = -x_0_0 + x_0_1;
                    const y_eta = -y_0_0 + y_0_1;
                    const x_eta2 = x_0_0 - 2 * x_0_1 + x_0_2;
                    const y_eta2 = y_0_0 - 2 * y_0_1 + y_0_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                    control_function[block_range_start + local_id] = .init(p, q);
                    local_id += 1;

                    for (1..block.points.size[1]) |j| {
                        const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                        control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                        local_id += 1;
                    }
                }

                // edge i_min
                {
                    // var local_id: usize = size[1];
                    for (1..block.points.size[0] - 1) |_| {
                        const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
                        const x_i_0, const y_i_0 = block.points.data[local_id].data;
                        const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
                        const x_i_2, const y_i_2 = block.points.data[local_id + 2].data;
                        const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;

                        // central differences
                        const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
                        const y_xi = 0.5 * (y_ip1_0 - y_im1_0);
                        const x_xi2 = x_ip1_0 - 2.0 * x_i_0 + x_im1_0;
                        const y_xi2 = y_ip1_0 - 2.0 * y_i_0 + y_im1_0;

                        const g11 = x_xi * x_xi + y_xi * y_xi;

                        // forward differences
                        const x_eta = -x_i_0 + x_i_1;
                        const y_eta = -y_i_0 + y_i_1;
                        const x_eta2 = x_i_0 - 2 * x_i_1 + x_i_2;
                        const y_eta2 = y_i_0 - 2 * y_i_1 + y_i_2;

                        const g22 = x_eta * x_eta + y_eta * y_eta;

                        // eq. 6.10
                        const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                        const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                        // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                        control_function[block_range_start + local_id] = .init(p, q);
                        local_id += 1;
                        // local_id += size[1];

                        for (1..block.points.size[1]) |j| {
                            const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                            control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                            local_id += 1;
                        }
                    }
                }

                {
                    // corner n,0
                    // TODO: hard code indices to gain performance.
                    // const local_id = (size[0] - 1) * size[1];
                    const x_n_0, const y_n_0 = block.points.data[local_id].data;
                    const x_n_1, const y_n_1 = block.points.data[local_id + 1].data;
                    const x_n_2, const y_n_2 = block.points.data[local_id + 2].data;
                    const x_nm1_0, const y_nm1_0 = block.points.data[local_id - size[1]].data;
                    const x_nm2_0, const y_nm2_0 = block.points.data[local_id - 2 * size[1]].data;

                    // backward differences
                    const x_xi = x_n_0 - x_nm1_0;
                    const y_xi = y_n_0 - y_nm1_0;

                    const x_xi2 = x_n_0 - 2 * x_nm1_0 + x_nm2_0;
                    const y_xi2 = y_n_0 - 2 * y_nm1_0 + y_nm2_0;

                    // forward differences
                    const x_eta = -x_n_0 + x_n_1;
                    const y_eta = -y_n_0 + y_n_1;

                    const x_eta2 = x_n_0 - 2 * x_n_1 + x_n_2;
                    const y_eta2 = y_n_0 - 2 * y_n_1 + y_n_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                    control_function[block_range_start + local_id] = .init(p, q);
                    local_id += 1;

                    for (1..block.points.size[1]) |j| {
                        const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                        control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                        local_id += 1;
                    }
                }

                const num_points = block.points.size[0] * block.points.size[1];
                block_range_start += num_points;
            }
        }

        fn update(self: White, control_function: []types.Vec2d) void {
            _ = self;
            _ = control_function;
        }
    };

    const White = struct {
        mesh: *const discrete.Mesh,

        fn init(mesh_data: *const discrete.Mesh, control_function: []types.Vec2d) White {
            var white = White{ .mesh = mesh_data };
            white.initControlFunction(control_function);
            return white;
        }

        fn deinit(self: White) void {
            _ = self;
        }

        fn initControlFunction(self: *@This(), control_function: []types.Vec2d) void {
            var block_range_start: usize = 0;
            for (self.mesh.blocks.items[0..2]) |block| {
                const size = block.points.size;

                var local_id: usize = 0;

                {
                    // corner 0,0
                    // TODO: hard code indices to gain performance.
                    const x_0_0, const y_0_0 = block.points.data[local_id].data;
                    const x_0_1, const y_0_1 = block.points.data[local_id + 1].data;
                    const x_0_2, const y_0_2 = block.points.data[local_id + 2].data;
                    const x_1_0, const y_1_0 = block.points.data[local_id + size[1]].data;
                    const x_2_0, const y_2_0 = block.points.data[local_id + 2 * size[1]].data;

                    // forward differences
                    const x_xi = -x_0_0 + x_1_0;
                    const y_xi = -y_0_0 + y_1_0;
                    const x_xi2 = x_0_0 - 2 * x_1_0 + x_2_0;
                    const y_xi2 = y_0_0 - 2 * y_1_0 + y_2_0;

                    const x_eta = -x_0_0 + x_0_1;
                    const y_eta = -y_0_0 + y_0_1;
                    const x_eta2 = x_0_0 - 2 * x_0_1 + x_0_2;
                    const y_eta2 = y_0_0 - 2 * y_0_1 + y_0_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                    control_function[block_range_start + local_id] = .init(p, q);
                    local_id += 1;

                    for (1..block.points.size[1]) |j| {
                        const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                        control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                        local_id += 1;
                    }
                }

                // edge i_min
                {
                    // var local_id: usize = size[1];
                    for (1..block.points.size[0] - 1) |_| {
                        const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
                        const x_i_0, const y_i_0 = block.points.data[local_id].data;
                        const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
                        const x_i_2, const y_i_2 = block.points.data[local_id + 2].data;
                        const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;

                        // central differences
                        const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
                        const y_xi = 0.5 * (y_ip1_0 - y_im1_0);
                        const x_xi2 = x_ip1_0 - 2.0 * x_i_0 + x_im1_0;
                        const y_xi2 = y_ip1_0 - 2.0 * y_i_0 + y_im1_0;

                        const g11 = x_xi * x_xi + y_xi * y_xi;

                        // forward differences
                        const x_eta = -x_i_0 + x_i_1;
                        const y_eta = -y_i_0 + y_i_1;
                        const x_eta2 = x_i_0 - 2 * x_i_1 + x_i_2;
                        const y_eta2 = y_i_0 - 2 * y_i_1 + y_i_2;

                        const g22 = x_eta * x_eta + y_eta * y_eta;

                        // eq. 6.10
                        const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                        const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                        // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                        control_function[block_range_start + local_id] = .init(p, q);
                        local_id += 1;
                        // local_id += size[1];

                        for (1..block.points.size[1]) |j| {
                            const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                            control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                            local_id += 1;
                        }
                    }
                }

                {
                    // corner n,0
                    // TODO: hard code indices to gain performance.
                    // const local_id = (size[0] - 1) * size[1];
                    const x_n_0, const y_n_0 = block.points.data[local_id].data;
                    const x_n_1, const y_n_1 = block.points.data[local_id + 1].data;
                    const x_n_2, const y_n_2 = block.points.data[local_id + 2].data;
                    const x_nm1_0, const y_nm1_0 = block.points.data[local_id - size[1]].data;
                    const x_nm2_0, const y_nm2_0 = block.points.data[local_id - 2 * size[1]].data;

                    // backward differences
                    const x_xi = x_n_0 - x_nm1_0;
                    const y_xi = y_n_0 - y_nm1_0;

                    const x_xi2 = x_n_0 - 2 * x_nm1_0 + x_nm2_0;
                    const y_xi2 = y_n_0 - 2 * y_nm1_0 + y_nm2_0;

                    // forward differences
                    const x_eta = -x_n_0 + x_n_1;
                    const y_eta = -y_n_0 + y_n_1;

                    const x_eta2 = x_n_0 - 2 * x_n_1 + x_n_2;
                    const y_eta2 = y_n_0 - 2 * y_n_1 + y_n_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    // std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                    control_function[block_range_start + local_id] = .init(p, q);
                    local_id += 1;

                    for (1..block.points.size[1]) |j| {
                        const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size[1])) - 1);
                        control_function[block_range_start + local_id] = .init(factor * p, factor * q);
                        local_id += 1;
                    }
                }

                const num_points = block.points.size[0] * block.points.size[1];
                block_range_start += num_points;
            }
        }

        fn update(self: White, control_function: []types.Vec2d) void {
            _ = self;
            _ = control_function;
        }
    };

    const KhamaysehEtAl = struct {
        allocator: std.mem.Allocator,
        mesh: *const discrete.Mesh,
        ghost_points: []types.Vec2d,
        g: []types.Float,
        buffer: []types.Vec2d,

        // TODO: move ghost_points into buffer to save an allocation.

        fn init(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh, control_function: []types.Vec2d) !KhamaysehEtAl {
            const n_ghost_points = blk: {
                var count: usize = 0;
                for (mesh_data.blocks.items[0..2]) |block| {
                    const size = block.points.size;
                    count += size[0] - 2;
                }
                break :blk count;
            };

            const buffer_size = types.Index2d{
                @max(mesh_data.blocks.items[0].points.size[0], mesh_data.blocks.items[1].points.size[0]),
                @max(mesh_data.blocks.items[0].points.size[1], mesh_data.blocks.items[1].points.size[1]),
            };

            var self = KhamaysehEtAl{
                .allocator = allocator,
                .mesh = mesh_data,
                .ghost_points = try allocator.alloc(types.Vec2d, n_ghost_points),
                .g = try allocator.alloc(types.Float, n_ghost_points),
                .buffer = try allocator.alloc(types.Vec2d, 2 * buffer_size[0] + 2 * buffer_size[1]),
            };

            self.computeGhostPoints(control_function);
            return self;
        }

        fn deinit(self: @This()) void {
            self.allocator.free(self.ghost_points);
            self.allocator.free(self.g);
            self.allocator.free(self.buffer);
        }

        fn computeGhostPoints(self: *@This(), control_function: []types.Vec2d) void {
            _ = control_function;
            var point_id: usize = 0;
            for (self.mesh.blocks.items[0..2], 0..2) |block, block_id| {
                const size = block.points.size;

                // @memset(self.buffer, .init(0, 0));

                // const edge_i_min = self.buffer[0..size[0]];
                // const edge_i_max = self.buffer[size[0] .. 2 * size[0]];
                // const edge_j_min = self.buffer[2 * size[0] .. 2 * size[0] + size[1]];
                // const edge_j_max = self.buffer[2 * size[0] + size[1] .. 2 * size[0] + 2 * size[1]];

                {
                    // corner 0,0
                    // TODO: hard code indices to gain performance.
                    const local_id = 0;
                    const x_0_0, const y_0_0 = block.points.data[local_id].data;
                    const x_0_1, const y_0_1 = block.points.data[local_id + 1].data;
                    const x_0_2, const y_0_2 = block.points.data[local_id + 2].data;
                    const x_1_0, const y_1_0 = block.points.data[local_id + size[1]].data;
                    const x_2_0, const y_2_0 = block.points.data[local_id + 2 * size[1]].data;

                    // forward differences
                    const x_xi = -x_0_0 + x_1_0;
                    const y_xi = -y_0_0 + y_1_0;

                    const x_xi2 = x_0_0 - 2 * x_1_0 + x_2_0;
                    const y_xi2 = y_0_0 - 2 * y_1_0 + y_2_0;

                    const x_eta = -x_0_0 + x_0_1;
                    const y_eta = -y_0_0 + y_0_1;

                    const x_eta2 = x_0_0 - 2 * x_0_1 + x_0_2;
                    const y_eta2 = y_0_0 - 2 * y_0_1 + y_0_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                    // edge_i_min[0] = .init(p, q);
                }

                // edge i_min
                {
                    var local_id: usize = size[1];
                    for (1..block.points.size[0] - 1) |_| {
                        const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
                        const x_i_0, const y_i_0 = block.points.data[local_id].data;
                        const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
                        // const x_i_2, const y_i_2 = block.points.data[local_id + 2].data;
                        const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;

                        // central differences
                        const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
                        const y_xi = 0.5 * (y_ip1_0 - y_im1_0);
                        const x_xi2 = x_ip1_0 - 2.0 * x_i_0 + x_im1_0;
                        const y_xi2 = y_ip1_0 - 2.0 * y_i_0 + y_im1_0;
                        const g11 = x_xi * x_xi + y_xi * y_xi;

                        // one sided approx normal to wall
                        const x_eta_0 = x_i_1 - x_i_0;
                        const y_eta_0 = y_i_1 - y_i_0;

                        // eq. 6.13
                        const factor = (-y_xi * x_eta_0 + x_xi * y_eta_0) / g11;
                        const x_eta_fixed = factor * y_xi;
                        const y_eta_fixed = factor * -x_xi;

                        const x_i_m1 = x_i_0 - x_eta_fixed;
                        const y_i_m1 = y_i_0 - y_eta_fixed;
                        // self.ghost_points[point_id] = .init(x_ghost, y_ghost);

                        // const g22 = x_eta_fixed * x_eta_fixed + y_eta_fixed * y_eta_fixed;
                        // self.g[point_id] = g22;

                        const x_eta = x_i_1 - x_i_0;
                        const y_eta = y_i_1 - y_i_0;
                        const x_eta2 = x_i_m1 - 2.0 * x_i_0 + x_i_1;
                        const y_eta2 = y_i_m1 - 2.0 * y_i_0 + y_i_1;

                        const g22 = x_eta * x_eta + y_eta * y_eta;

                        // // forward
                        // const x_eta = -x_i_0 + x_i_1;
                        // const y_eta = -y_i_0 + y_i_1;
                        //
                        // const x_eta2 = x_i_0 - 2 * x_i_1 + x_i_2;
                        // const y_eta2 = y_i_0 - 2 * y_i_1 + y_i_2;
                        //
                        // const g22 = x_eta * x_eta + y_eta * y_eta;

                        // eq. 6.10
                        const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                        const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                        std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });

                        point_id += 1;
                        local_id += size[1];
                    }
                }

                {
                    // corner n,0
                    // TODO: hard code indices to gain performance.
                    const local_id = (size[0] - 1) * size[1];
                    const x_n_0, const y_n_0 = block.points.data[local_id].data;
                    const x_n_1, const y_n_1 = block.points.data[local_id + 1].data;
                    const x_n_2, const y_n_2 = block.points.data[local_id + 2].data;
                    const x_nm1_0, const y_nm1_0 = block.points.data[local_id - size[1]].data;
                    const x_nm2_0, const y_nm2_0 = block.points.data[local_id - 2 * size[1]].data;

                    // backward differences
                    const x_xi = x_n_0 - x_nm1_0;
                    const y_xi = y_n_0 - y_nm1_0;

                    const x_xi2 = x_n_0 - 2 * x_nm1_0 + x_nm2_0;
                    const y_xi2 = y_n_0 - 2 * y_nm1_0 + y_nm2_0;

                    // forward differences
                    const x_eta = -x_n_0 + x_n_1;
                    const y_eta = -y_n_0 + y_n_1;

                    const x_eta2 = x_n_0 - 2 * x_n_1 + x_n_2;
                    const y_eta2 = y_n_0 - 2 * y_n_1 + y_n_2;

                    const g11 = x_xi * x_xi + y_xi * y_xi;
                    const g22 = x_eta * x_eta + y_eta * y_eta;

                    // eq. 6.10
                    const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                    const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                    std.debug.print("block {} point {} P {} Q {}\n", .{ block_id, point_id, p, q });
                }

                std.posix.exit(0);
            }
        }

        fn update(self: @This(), control_function: []types.Vec2d) !void {
            var point_id: usize = 0;
            var block_range_start: usize = 0;
            for (self.mesh.blocks.items[0..2]) |block| {
                const size = block.points.size;

                @memset(self.buffer, .init(0, 0));

                const edge_i_min = self.buffer[0..size[0]];
                const edge_i_max = self.buffer[size[0] .. 2 * size[0]];
                const edge_j_min = self.buffer[2 * size[0] .. 2 * size[0] + size[1]];
                const edge_j_max = self.buffer[2 * size[0] + size[1] .. 2 * size[0] + 2 * size[1]];

                // NOTE: hard coded for i_min right now!

                // edge i_min
                {
                    var local_id: usize = size[1];
                    for (1..block.points.size[0] - 1) |edge_id| {
                        const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
                        const x_i_0, const y_i_0 = block.points.data[local_id].data;
                        const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
                        const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;

                        // central differences
                        const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
                        const y_xi = 0.5 * (y_ip1_0 - y_im1_0);
                        const x_xi2 = x_ip1_0 - 2.0 * x_i_0 + x_im1_0;
                        const y_xi2 = y_ip1_0 - 2.0 * y_i_0 + y_im1_0;
                        const g11 = x_xi * x_xi + y_xi * y_xi;

                        const x_i_m1, const y_i_m1 = self.ghost_points[point_id].data;

                        // const x_eta_f = 0.5 * (x_i_1 - x_i_m1);
                        // const y_eta_f = 0.5 * (y_i_1 - y_i_m1);

                        const x_eta = x_i_1 - x_i_0;
                        const y_eta = y_i_1 - y_i_0;
                        const x_eta2 = x_i_m1 - 2.0 * x_i_0 + x_i_1;
                        const y_eta2 = y_i_m1 - 2.0 * y_i_0 + y_i_1;

                        // const g22 = self.g[point_id];
                        // const g22 = x_eta_f * x_eta_f + y_eta_f * y_eta_f;
                        const g22 = x_eta * x_eta + y_eta * y_eta;

                        // eq. 6.10
                        const p = -(x_xi * x_xi2 + y_xi * y_xi2) / g11 - (x_xi * x_eta2 + y_xi * y_eta2) / g22;
                        const q = -(x_eta * x_eta2 + y_eta * y_eta2) / g22 - (x_eta * x_xi2 + y_eta * y_xi2) / g11;

                        edge_i_min[edge_id] = .init(p, q);

                        local_id += size[1];
                        point_id += 1;
                    }
                }

                const num_points = block.points.size[0] * block.points.size[1];
                try tfi.linear2d(control_function[block_range_start .. block_range_start + num_points], edge_i_min, edge_i_max, edge_j_min, edge_j_max);
                block_range_start += num_points;
            }
        }
    };
};
