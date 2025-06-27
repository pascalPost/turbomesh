const std = @import("std");
const types = @import("../types.zig");
const discrete = @import("../discrete.zig");
const tfi = @import("../tfi.zig");
const smooth = @import("smooth.zig");

const Type = enum {
    laplace,
    white,
    // khamaysehEtAl,
};

pub const Algorithm = union(Type) {
    laplace: void,
    white: White,
    // khamaysehEtAl: KhamaysehEtAl,
};

pub const ControlFunction = struct {
    allocator: std.mem.Allocator,
    data: []types.Vec2d,
    algorithm: Algorithm,

    pub fn init(allocator: std.mem.Allocator, dof: usize, mesh: discrete.Mesh, algorithm: Algorithm) !ControlFunction {
        // TODO: merge with other allocations.
        var data = try allocator.alloc(types.Vec2d, dof);
        @memset(data[0..], .{ .data = .{ 0, 0 } });

        switch (algorithm) {
            .laplace => {},
            .white => White.initControlFunction(data, mesh),
        }

        return .{
            .allocator = allocator,
            .data = data,
            .algorithm = algorithm,
        };
    }

    pub fn deinit(self: ControlFunction) void {
        self.allocator.free(self.data);
    }

    pub fn update(self: ControlFunction, mesh: discrete.Mesh) void {
        switch (self.algorithm) {
            .laplace => {},
            .white => |a| a.update(self.data, mesh),
        }
    }
};

pub const White = struct {
    /// target distance to first point
    ds_target: f64,

    /// target wall angle between point on wall and first point
    theta_target: f64 = 0.5 * std.math.pi,

    pub fn init(ds_target: f64, theta_target: f64) White {
        return .{
            .ds_target = ds_target,
            .theta_target = theta_target,
        };
    }

    fn initControlFunction(control_function: []types.Vec2d, mesh: discrete.Mesh) void {
        var block_range_start: usize = 0;
        for (mesh.blocks.items[0..2]) |block| {
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
                // const x_xi = -1.5 * x_0_0 + 2 * x_1_0 - 0.5 * x_
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

            // TODO: add connections!
        }
    }

    fn computeUpdate(
        self: White,
        local_id: *usize,
        control_function: []types.Vec2d,
        x_xi: f64,
        y_xi: f64,
        x_eta: f64,
        y_eta: f64,
        size_j: usize,
        block_range_start: usize,
    ) void {
        const g11 = x_xi * x_xi + y_xi * y_xi;
        const g12 = x_xi * x_eta + y_xi * y_eta;
        const g22 = x_eta * x_eta + y_eta * y_eta;

        const ds = @sqrt(g22);
        const theta = std.math.acos(g12 / @sqrt(g11 * g22));

        const delta_ds = self.ds_target - ds;
        const delta_theta = self.theta_target - theta;

        const delta_p = -std.math.atan2(delta_theta, self.theta_target);
        const delta_q = std.math.atan2(delta_ds, self.ds_target);

        var p, var q = control_function[block_range_start + local_id.*].data;
        p += 0.1 * delta_p;
        q += 0.1 * delta_q;
        control_function[block_range_start + local_id.*] = .init(p, q);
        local_id.* += 1;

        // TODO: wouldn't it be better to do this once at the end!?
        // TODO: this is hard coded to run over j

        for (1..size_j) |j| {
            const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size_j)) - 1);
            control_function[block_range_start + local_id.*] = .init(factor * p, factor * q);
            local_id.* += 1;
        }
    }

    fn update(self: White, control_function: []types.Vec2d, mesh: discrete.Mesh) void {
        var block_range_start: usize = 0;

        // TODO: remove this hard coding!

        for (mesh.blocks.items[0..2]) |block| {
            const size = block.points.size;

            var local_id: usize = 0;

            {
                // NOTE: this is only necessary if this point is not connected!

                const x_0_0, const y_0_0 = block.points.data[local_id].data;
                const x_0_1, const y_0_1 = block.points.data[local_id + 1].data;
                const x_1_0, const y_1_0 = block.points.data[local_id + size[1]].data;

                // forward differences
                const x_xi = -x_0_0 + x_1_0;
                const y_xi = -y_0_0 + y_1_0;

                // forward differences
                const x_eta = -x_0_0 + x_0_1;
                const y_eta = -y_0_0 + y_0_1;

                self.computeUpdate(&local_id, control_function, x_xi, y_xi, x_eta, y_eta, block.points.size[1], block_range_start);
            }

            {
                for (1..block.points.size[0] - 1) |_| {
                    const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
                    const x_i_0, const y_i_0 = block.points.data[local_id].data;
                    const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
                    const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;

                    // central differences
                    const x_xi = 0.5 * (x_ip1_0 - x_im1_0);
                    const y_xi = 0.5 * (y_ip1_0 - y_im1_0);

                    // forward differences
                    const x_eta = -x_i_0 + x_i_1;
                    const y_eta = -y_i_0 + y_i_1;

                    self.computeUpdate(&local_id, control_function, x_xi, y_xi, x_eta, y_eta, block.points.size[1], block_range_start);
                }
            }

            {
                // NOTE: this is only necessary if this point is not connected!

                const x_n_0, const y_n_0 = block.points.data[local_id].data;
                const x_n_1, const y_n_1 = block.points.data[local_id + 1].data;
                const x_nm1_0, const y_nm1_0 = block.points.data[local_id - size[1]].data;

                // backward differences
                const x_xi = x_n_0 - x_nm1_0;
                const y_xi = y_n_0 - y_nm1_0;

                // forward differences
                const x_eta = -x_n_0 + x_n_1;
                const y_eta = -y_n_0 + y_n_1;

                self.computeUpdate(&local_id, control_function, x_xi, y_xi, x_eta, y_eta, block.points.size[1], block_range_start);
            }

            const num_points = block.points.size[0] * block.points.size[1];
            block_range_start += num_points;
        }

        // for (mesh.connections.items) |connection| {
        //
        //     // TODO: remove hard coding
        //     // TODO: save for which connections this needs to be done to avoid looping and chacling all
        //
        //     const range_0 = connection.ranges[0];
        //     const range_1 = connection.ranges[1];
        //
        //     if (range_0.block == 0 and range_1.block == 1) {
        //         // NOTE: not (yet) implemented... not sure if this could be needed at some point.
        //         std.debug.assert(connection.periodicity == null);
        //
        //         const point_data = [2][]types.Vec2d{
        //             self.mesh.blocks.items[connection.ranges[0].block].points.data,
        //             self.mesh.blocks.items[connection.ranges[1].block].points.data,
        //         };
        //
        //         var it = smooth.RangeFillMatrixIterator.init(connection, self.mesh);
        //
        //         const connected_points = it.next().?;
        //         const point_idx: [2]c_int = .{ @intCast(connected_points[0].value), @intCast(connected_points[1].value) };
        //
        //         // NOTE: i an j are arbitrary for connections....
        //         // TODO: make this consistent with the real directions of the block
        //
        //         const x_i_j, const y_i_j = point_data[0][@intCast(point_idx[0])];
        //         const x_i_jm1, const y_i_jm1 = point_data[0][@intCast(point_idx[0] + it.first_internal_point_shift[0])];
        //         const x_i_jp1, const y_i_jp1 = point_data[1][@intCast(point_idx[1] + it.first_internal_point_shift[1])];
        //         const x_ip1_j, const y_ip1_j = point_data[0][@intCast(point_idx[0] + it.in_connection_direction_shift[0])];
        //
        //         // const x_im1_0, const y_im1_0 = block.points.data[local_id - size[1]].data;
        //         // const x_i_0, const y_i_0 = block.points.data[local_id].data;
        //         // const x_i_1, const y_i_1 = block.points.data[local_id + 1].data;
        //         // const x_ip1_0, const y_ip1_0 = block.points.data[local_id + size[1]].data;
        //
        //         // central differences
        //         const x_xi = 0.5 * (x_i_jp1 - x_i_jm1);
        //         const y_xi = 0.5 * (y_i_jp1 - y_i_jm1);
        //
        //         // forward differences
        //         const x_eta = -x_i_j + x_ip1_j;
        //         const y_eta = -y_i_j + y_ip1_j;
        //
        //         {
        //             const g11 = x_xi * x_xi + y_xi * y_xi;
        //             const g12 = x_xi * x_eta + y_xi * y_eta;
        //             const g22 = x_eta * x_eta + y_eta * y_eta;
        //
        //             const ds = @sqrt(g22);
        //             const theta = std.math.acos(g12 / @sqrt(g11 * g22));
        //
        //             const delta_ds = self.ds_target - ds;
        //             const delta_theta = self.theta_target - theta;
        //
        //             const delta_p = -std.math.atan2(delta_theta, self.theta_target);
        //             const delta_q = std.math.atan2(delta_ds, self.ds_target);
        //
        //             var p, var q = control_function[block_range_start + local_id.*].data;
        //             p += 0.1 * delta_p;
        //             q += 0.1 * delta_q;
        //             control_function[block_range_start + local_id.*] = .init(p, q);
        //             local_id.* += 1;
        //
        //             // TODO: wouldn't it be better to do this once at the end!?
        //             // TODO: this is hard coded to run over j
        //
        //             for (1..size_j) |j| {
        //                 const factor: f64 = 1 - @as(f64, @floatFromInt(j)) / (@as(f64, @floatFromInt(size_j)) - 1);
        //                 control_function[block_range_start + local_id.*] = .init(factor * p, factor * q);
        //                 local_id.* += 1;
        //             }
        //         }
        //
        //         {
        //             // TODO: enhance this!!!
        //             var local_id = connected_points[0].value;
        //             const block = mesh.blocks.items[range_0.block];
        //             self.computeUpdate(&local_id, control_function, x_xi, y_xi, x_eta, y_eta, block.points.size[1], 0);
        //         }
        //         self.computeUpdate(&local_id, control_function, x_xi, y_xi, x_eta, y_eta, block.points.size[1], block_range_start);
        //     }
        // }
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

// fn initControlFunctionThomasAndMiddlecoff(self: *RowCompressedMatrixSystem2d) !void {
//
//     // TODO: we assume the first 2 blocks to be the O blocks
//     // TODO: for now, we just compute on the j edges
//
//     // TODO: introduce global buffer to reduce number of allocations.
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
//         // edge i_min
//         {
//             var local_id: usize = size[1];
//             for (1..block.points.size[0] - 1) |edge_id| {
//                 const im1_j = block.points.data[local_id - size[1]];
//                 const i_j = block.points.data[local_id];
//                 const ip1_j = block.points.data[local_id + size[1]];
//
//                 std.debug.assert(local_id - size[1] == block.points.index(.{ edge_id - 1, 0 }));
//                 std.debug.assert(local_id == block.points.index(.{ edge_id, 0 }));
//                 std.debug.assert(local_id + size[1] == block.points.index(.{ edge_id + 1, 0 }));
//
//                 // central differences
//                 const x_xi = types.scale(0.5, types.sub(ip1_j, im1_j));
//                 const x_xi2 = types.sub(types.add(ip1_j, im1_j), types.scale(2.0, i_j));
//
//                 // eq. 11
//                 const phi = -(x_xi.data[0] * x_xi2.data[0] + x_xi.data[1] * x_xi2.data[1]) / (x_xi.data[0] * x_xi.data[0] + x_xi.data[1] * x_xi.data[1]);
//                 const psi = 0.0;
//
//                 edge_i_min[edge_id] = .init(phi, psi);
//
//                 local_id += size[1];
//             }
//         }
//
//         // edge i_max
//         {
//             var local_id: usize = 2 * size[1] - 1;
//             for (1..block.points.size[0] - 1) |edge_id| {
//                 const im1_j = block.points.data[local_id - size[1]];
//                 const i_j = block.points.data[local_id];
//                 const ip1_j = block.points.data[local_id + size[1]];
//
//                 std.debug.assert(local_id - size[1] == block.points.index(.{ edge_id - 1, size[1] - 1 }));
//                 std.debug.assert(local_id == block.points.index(.{ edge_id, size[1] - 1 }));
//                 std.debug.assert(local_id + size[1] == block.points.index(.{ edge_id + 1, size[1] - 1 }));
//
//                 // central differences
//                 const x_xi = types.scale(0.5, types.sub(ip1_j, im1_j));
//                 const x_xi2 = types.sub(types.add(ip1_j, im1_j), types.scale(2.0, i_j));
//
//                 // eq. 11
//                 const phi = -(x_xi.data[0] * x_xi2.data[0] + x_xi.data[1] * x_xi2.data[1]) / (x_xi.data[0] * x_xi.data[0] + x_xi.data[1] * x_xi.data[1]);
//                 const psi = 0.0;
//
//                 edge_i_max[edge_id] = .init(phi, psi);
//
//                 local_id += size[1];
//             }
//         }
//
//         // edge j_min
//         {
//             var local_id: usize = 1;
//             for (1..block.points.size[1] - 1) |edge_id| {
//                 const i_jm1 = block.points.data[local_id - 1];
//                 const i_j = block.points.data[local_id];
//                 const i_jp1 = block.points.data[local_id + 1];
//
//                 std.debug.assert(local_id - 1 == block.points.index(.{ 0, edge_id - 1 }));
//                 std.debug.assert(local_id == block.points.index(.{ 0, edge_id }));
//                 std.debug.assert(local_id + 1 == block.points.index(.{ 0, edge_id + 1 }));
//
//                 // central differences
//                 const x_eta = types.scale(0.5, types.sub(i_jp1, i_jm1));
//                 const x_eta2 = types.sub(types.add(i_jp1, i_jm1), types.scale(2.0, i_j));
//
//                 // eq. 11
//                 const phi = 0.0;
//                 const psi = -(x_eta.data[0] * x_eta2.data[0] + x_eta.data[1] * x_eta2.data[1]) / (x_eta.data[0] * x_eta.data[0] + x_eta.data[1] * x_eta.data[1]);
//
//                 edge_j_min[edge_id] = .init(phi, psi);
//
//                 local_id += 1;
//             }
//         }
//
//         // edge j_max
//         {
//             var local_id = (block.points.size[0] - 1) * block.points.size[1] + 1;
//             for (1..block.points.size[1] - 1) |edge_id| {
//                 const i_jm1 = block.points.data[local_id - 1];
//                 const i_j = block.points.data[local_id];
//                 const i_jp1 = block.points.data[local_id + 1];
//
//                 std.debug.assert(local_id - 1 == block.points.index(.{ size[0] - 1, edge_id - 1 }));
//                 std.debug.assert(local_id == block.points.index(.{ size[0] - 1, edge_id }));
//                 std.debug.assert(local_id + 1 == block.points.index(.{ size[0] - 1, edge_id + 1 }));
//
//                 // central differences
//                 const x_eta = types.scale(0.5, types.sub(i_jp1, i_jm1));
//                 const x_eta2 = types.sub(types.add(i_jp1, i_jm1), types.scale(2.0, i_j));
//
//                 // eq. 11
//                 const phi = 0.0;
//                 const psi = -(x_eta.data[0] * x_eta2.data[0] + x_eta.data[1] * x_eta2.data[1]) / (x_eta.data[0] * x_eta.data[0] + x_eta.data[1] * x_eta.data[1]);
//
//                 edge_j_max[edge_id] = .init(phi, psi);
//
//                 local_id += 1;
//             }
//         }
//
//         const num_points = block.points.size[0] * block.points.size[1];
//         try tfi.linear2d(self.control_function[block_range_start .. block_range_start + num_points], edge_i_min, edge_i_max, edge_j_min, edge_j_max);
//         block_range_start += num_points;
//     }
// }
//
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
