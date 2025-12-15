// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

const blade = @import("blade.zig");
const types = @import("../types.zig");
const clustering = @import("../clustering.zig");
const discrete = @import("../discrete.zig");
const geometry = @import("../geometry.zig");
const boundary = @import("../boundary.zig");
const smooth = @import("../smoothing/smooth.zig");

const Float = types.Float;
const Vec2d = types.Vec2d;
const add = types.add;
const sub = types.sub;
const scale = types.scale;
const abs = types.abs;

/// O4H multi block template:
///
///  _______________________________________________________________________________________________
/// |               |           *                                               ** |                |
/// |               |           *               up (5)                     *****   |                |
/// |               |__________________________________________________****        |                |
/// |               |    i<-|   /               ss (0)                 \  i^ ->j   |                |
/// |               |       vj /  ____________________________________  \__________|                |
/// |               |         /  /                                    \  \         |                |
/// |  upstream (6) | IN (2) |--|* leading edge        trailing edge *|--| out (3) | downstream (7) |
/// |               |         \  \____________________________________/  / ^j      |                |
/// |               |          \                ps (1)                  /  ->i     |                |
/// |               |___________\______________________________________/___________|                |
/// | ^j            |  -> i     *                                      *           | ^j             |
/// | |             |  vj       *             down (4)                 *           | |              |
/// | ->i           |           *                                      *           | ->i            |
/// |_______________________________________________________________________________________________
pub const O4H = struct {

    // TODO: allow to determine inlet and outlet automatically from something like the default mesh size.
    inlet_axial_position: types.Float,
    outlet_axial_position: types.Float,

    ps_csv_path: []const u8,
    ss_csv_path: []const u8,
    pitch: types.Float,
    blade_clustering: clustering.Function,
    num_cells: struct {
        o_grid: types.Index,
        middle_i: types.Index,
        in_up_j: types.Index,
        in_down_j: types.Index,
        in_i: types.Index,
        out_up_j: types.Index,
        out_down_j: types.Index,
        out_i: types.Index,
        down_j: types.Index,
        bulge: types.Index,

        // TODO: merge in_i and out_i

        // TODO: remove unnecessary parameters!

        // TODO: can be made optional if computed automatically based on the average size
        upstream_i: types.Index,
        downstream_i: types.Index,
    },
    smoothing: struct {
        iterations: usize,
        backend: smooth.solver.Type = .umfpack,
        wall_control_function: smooth.control_function.Algorithm = .{ .laplace = {} },
    },

    pub fn run(self: *const O4H, allocator: std.mem.Allocator) !discrete.Mesh {

        // TODO: remove pressure side and suction side and just use up and down (to be neutral w.r.t. turbine or compressor)

        // TODO: add geometry and discrete entities to manager

        const num_cells_ss =
            self.num_cells.in_up_j + self.num_cells.middle_i + self.num_cells.bulge + self.num_cells.out_up_j + self.num_cells.out_i;
        const num_cells_ps =
            self.num_cells.in_down_j + self.num_cells.middle_i + self.num_cells.out_down_j;

        const profile = try blade.Profile.init(allocator, self.ps_csv_path, self.ss_csv_path);
        defer profile.deinit();

        var ps_edge = try discrete.Edge.init(allocator, num_cells_ps + 1, .{ .spline = profile.pressure_side }, self.blade_clustering);
        defer ps_edge.deinit();

        var ss_edge = try discrete.Edge.init(allocator, num_cells_ss + 1, .{ .spline = profile.suction_side }, self.blade_clustering);
        defer ss_edge.deinit();

        // TODO: introducing a tolerance, this should not be necessary anymore.
        // TODO: handle this by a connect points function (?)
        const leading_edge = ss_edge.points[0];
        ps_edge.points[0] = leading_edge;

        const trailing_edge = ss_edge.points[ss_edge.points.len - 1];
        ps_edge.points[ps_edge.points.len - 1] = trailing_edge;

        // o - grid for viscous computations

        // compute o-grid target by projecting the blade normal outward

        // TODO: make this runtime dependent
        // TODO: replace with a percentage value of the chord length
        const d = 0.001;

        const ps_outer_edge = discrete.Edge{ .allocator = allocator, .points = try projectNormal(allocator, ps_edge.points[0..], d), .clustering = try allocator.dupe(Float, ps_edge.clustering) };
        defer ps_outer_edge.deinit();

        const ss_outer_edge = blk: {
            var ss_outer = discrete.Edge{ .allocator = allocator, .points = try projectNormal(allocator, ss_edge.points[0..], -d), .clustering = try allocator.dupe(Float, ss_edge.clustering) };
            ss_outer.points[0] = ps_outer_edge.points[0];
            ss_outer.points[ss_outer.points.len - 1] = ps_outer_edge.points[ps_outer_edge.points.len - 1];
            break :blk ss_outer;
        };
        defer ss_outer_edge.deinit();

        var mesh = discrete.Mesh.init(allocator);
        errdefer mesh.deinit();

        // Block SS (0)
        //
        // |          /
        // |         |          ^
        // |        |           |
        // |-------< x_00      i_min
        //         LE

        const ss_i_min = ss_edge;
        const ss_i_max = ss_outer_edge;

        const ss_j_min = try discrete.Edge.init(
            allocator,
            self.num_cells.o_grid + 1,
            .{ .line = geometry.Line.init(ss_i_min.points[0], ss_i_max.points[0]) },
            .{ .single_hyperbolic_clustering = .{ .delta_s = 0.01 } },
        );
        defer ss_j_min.deinit();

        const ss_j_max = try discrete.Edge.init(
            allocator,
            self.num_cells.o_grid + 1,
            .{ .line = geometry.Line.init(ss_i_min.points[ss_edge.points.len - 1], ss_i_max.points[ss_i_max.points.len - 1]) },
            .{ .single_hyperbolic_clustering = .{ .delta_s = 0.01 } },
        );
        defer ss_j_max.deinit();

        const ss = try discrete.Block2d.init(allocator, ss_i_min, ss_i_max, ss_j_min, ss_j_max);

        try mesh.addBlock("ss", ss);
        const ss_id = mesh.blocks.items.len - 1;

        // Block PS (1)
        //
        //         LE
        // |-------< x_00      i_min
        // |        |            |
        // |         |           v
        // |          /

        const ps_i_min = ps_edge;
        const ps_i_max = ps_outer_edge;
        const ps_j_min = ss_j_min;
        const ps_j_max = ss_j_max;

        const ps = try discrete.Block2d.init(allocator, ps_i_min, ps_i_max, ps_j_min, ps_j_max);

        try mesh.addBlock("ps", ps);
        const ps_id = mesh.blocks.items.len - 1;

        // Block IN (2)
        //
        // x_01
        // |---------- x_00
        // |         |
        // |        |           |
        // |       < LE       j_min
        // |        |           |
        // |         |          V
        // |---------- x_10
        // x_11

        // TODO: remove hard coded positions for block

        const in_j_min = try discrete.Edge.combine(allocator, &.{ .{
            .edge = &ss_i_max,
            .start = self.num_cells.in_up_j,
            .end = 0,
        }, .{
            .edge = &ps_i_max,
            .start = 0,
            .end = self.num_cells.in_down_j,
        } });
        defer in_j_min.deinit();
        std.debug.assert(in_j_min.points.len == self.num_cells.in_up_j + self.num_cells.in_down_j + 1);

        const in_x_00 = in_j_min.points[0];
        const in_x_01 = in_j_min.points[in_j_min.points.len - 1];

        // TODO: remove this hardcoded distance
        const in_x_10 = sub(in_x_00, Vec2d.init(0.02, -0.001));
        const in_x_11 = sub(in_x_01, Vec2d.init(0.02, 0.02));

        const in_j_max = try discrete.Edge.init(allocator, in_j_min.points.len, .{ .line = .{ .start = in_x_10, .end = in_x_11 } }, .{ .uniform = .{} });
        defer in_j_max.deinit();
        const in_i_min = try discrete.Edge.init(allocator, self.num_cells.in_i + 1, .{ .line = .{ .start = in_x_00, .end = in_x_10 } }, .{ .uniform = .{} });
        defer in_i_min.deinit();
        const in_i_max = try discrete.Edge.init(allocator, self.num_cells.in_i + 1, .{ .line = .{ .start = in_x_01, .end = in_x_11 } }, .{ .uniform = .{} });
        defer in_i_max.deinit();

        const in = try discrete.Block2d.init(allocator, in_i_min, in_i_max, in_j_min, in_j_max);

        try mesh.addBlock("in", in);
        const in_id = mesh.blocks.items.len - 1;

        //
        // Block OUT (3)
        //

        const out_j_min = try discrete.Edge.combine(allocator, &.{ .{
            .edge = &ps_i_max,
            .start = self.num_cells.in_down_j + self.num_cells.middle_i,
            .end = ps_i_max.points.len - 1,
        }, .{
            .edge = &ss_i_max,
            .start = ss_i_max.points.len - 1,
            .end = self.num_cells.in_up_j + self.num_cells.bulge + self.num_cells.middle_i + self.num_cells.out_i,
        } });
        defer out_j_min.deinit();
        std.debug.assert(out_j_min.points.len == self.num_cells.out_down_j + self.num_cells.out_up_j + 1);

        const out_x_00 = out_j_min.points[0];
        const out_x_01 = out_j_min.points[out_j_min.points.len - 1];

        // TODO: remove hard coded coordinate
        const out_x_10 = add(out_x_00, Vec2d.init(0.01, -0.02));
        const out_x_11 = add(out_x_01, Vec2d.init(0.02, -0.01));

        const out_j_max = try discrete.Edge.init(allocator, out_j_min.points.len, .{ .line = .{ .start = out_x_10, .end = out_x_11 } }, .{ .uniform = .{} });
        defer out_j_max.deinit();

        const out_i_min = try discrete.Edge.init(allocator, self.num_cells.out_i + 1, .{ .line = .{ .start = out_x_00, .end = out_x_10 } }, .{ .uniform = .{} });
        defer out_i_min.deinit();
        const out_i_max = try discrete.Edge.init(allocator, self.num_cells.out_i + 1, .{ .line = .{ .start = out_x_01, .end = out_x_11 } }, .{ .uniform = .{} });
        defer out_i_max.deinit();

        const out = try discrete.Block2d.init(allocator, out_i_min, out_i_max, out_j_min, out_j_max);

        try mesh.addBlock("out", out);
        const out_id = mesh.blocks.items.len - 1;

        //
        // Block DOWN (4)
        //

        const down_i_min = try discrete.Edge.combine(allocator, &.{
            .{
                .edge = &in_i_max,
                .start = self.num_cells.in_i,
                .end = 0,
            },
            .{
                .edge = &ps_i_max,
                .start = self.num_cells.in_down_j,
                .end = self.num_cells.in_down_j + self.num_cells.middle_i,
            },
            .{
                .edge = &out_i_min,
                .start = 0,
                .end = self.num_cells.out_i,
            },
        });
        defer down_i_min.deinit();

        const down_x_00 = in_x_11;
        const down_x_01 = sub(leading_edge, Vec2d.init(0.0, 0.5 * self.pitch));
        const down_x_11 = sub(trailing_edge, Vec2d.init(0.0, 0.5 * self.pitch));
        const down_x_10 = out_x_10;

        const down_i_max = try discrete.Edge.init(allocator, down_i_min.points.len, .{ .line = .{ .start = down_x_01, .end = down_x_11 } }, .{ .uniform = .{} });
        defer down_i_max.deinit();

        const down_j_min = try discrete.Edge.init(allocator, self.num_cells.down_j + 1, .{ .line = .{ .start = down_x_00, .end = down_x_01 } }, .{ .uniform = .{} });
        defer down_j_min.deinit();

        const down_j_max = try discrete.Edge.init(allocator, down_j_min.points.len, .{ .line = .{ .start = down_x_10, .end = down_x_11 } }, .{ .uniform = .{} });
        defer down_j_max.deinit();

        const down = try discrete.Block2d.init(allocator, down_i_min, down_i_max, down_j_min, down_j_max);

        try mesh.addBlock("down", down);
        const down_id = mesh.blocks.items.len - 1;

        //
        // Block UP (5)
        //

        //
        // TODO: remove all the deinits w/ better design to increase efficiency!
        //

        const up_j_min = out_i_max;
        const up_i_min = try discrete.Edge.combine(allocator, &.{
            .{
                .edge = &ss_i_max,
                .start = self.num_cells.in_up_j + self.num_cells.middle_i + self.num_cells.bulge + self.num_cells.out_i,
                .end = self.num_cells.in_up_j,
            },
            .{
                .edge = &in_i_min,
                .start = 0,
                .end = self.num_cells.in_i,
            },
        });
        defer up_i_min.deinit();

        const up_x_11 = add(leading_edge, Vec2d.init(0.0, 0.5 * self.pitch));
        const up_x_i_max_middle = add(trailing_edge, Vec2d.init(0.0, 0.5 * self.pitch));
        const up_x_01 = out_x_11;
        const up_x_10 = in_x_10;

        const up_i_max_0 = try discrete.Edge.init(allocator, self.num_cells.bulge + 1, .{ .line = .{ .start = up_x_01, .end = up_x_i_max_middle } }, .{ .uniform = .{} });
        defer up_i_max_0.deinit();
        const up_i_max_1 = try discrete.Edge
            .init(allocator, up_i_min.points.len - self.num_cells.bulge, .{ .line = .{ .start = up_x_i_max_middle, .end = up_x_11 } }, .{ .uniform = .{} });
        defer up_i_max_1.deinit();

        const up_i_max = try discrete.Edge.combine(allocator, &.{
            .{
                .edge = &up_i_max_0,
                .start = 0,
                .end = self.num_cells.bulge,
            },
            .{
                .edge = &up_i_max_1,
                .start = 0,
                .end = up_i_max_1.points.len - 1,
            },
        });
        defer up_i_max.deinit();

        const up_j_max = try discrete.Edge.init(allocator, self.num_cells.out_i + 1, .{ .line = .{ .start = up_x_10, .end = up_x_11 } }, .{ .uniform = .{} });
        defer up_j_max.deinit();

        const up = try discrete.Block2d.init(allocator, up_i_min, up_i_max, up_j_min, up_j_max);

        try mesh.addBlock("up", up);
        const up_id = mesh.blocks.items.len - 1;

        //
        // Block UPSTREAM (6)
        //

        const upstream_j_max = try discrete.Edge.combine(allocator, &.{ .{
            .edge = &down_j_min,
            .start = self.num_cells.down_j,
            .end = 0,
        }, .{
            .edge = &in_j_max,
            .start = in_j_max.points.len - 1,
            .end = 0,
        }, .{
            .edge = &up_j_max,
            .start = 0,
            .end = up_j_max.points.len - 1,
        } });
        defer upstream_j_max.deinit();

        const upstream_x_10 = upstream_j_max.points[0];
        const upstream_x_11 = upstream_j_max.points[upstream_j_max.points.len - 1];

        // TODO: remove hard coding
        const upstream_x_00 = Vec2d.init(self.inlet_axial_position, leading_edge.data[1] - 0.5 * self.pitch);
        const upstream_x_01 = Vec2d.init(self.inlet_axial_position, leading_edge.data[1] + 0.5 * self.pitch);

        const upstream_j_min = try discrete.Edge.init(allocator, upstream_j_max.points.len, .{ .line = .{ .start = upstream_x_00, .end = upstream_x_01 } }, .{ .uniform = .{} });
        defer upstream_j_min.deinit();

        const upstream_i_min = try discrete.Edge.init(allocator, self.num_cells.upstream_i + 1, .{ .line = .{ .start = upstream_x_00, .end = upstream_x_10 } }, .{ .uniform = .{} });
        defer upstream_i_min.deinit();
        const upstream_i_max = try discrete.Edge.init(allocator, self.num_cells.upstream_i + 1, .{ .line = .{ .start = upstream_x_01, .end = upstream_x_11 } }, .{ .uniform = .{} });
        defer upstream_i_max.deinit();

        const upstream = try discrete.Block2d.init(allocator, upstream_i_min, upstream_i_max, upstream_j_min, upstream_j_max);

        try mesh.addBlock("upstream", upstream);
        const upstream_id = mesh.blocks.items.len - 1;

        //
        // Block DOWNSTREAM (7)
        //

        const downstream_j_min = try discrete.Edge.combine(allocator, &.{ .{
            .edge = &down_j_max,
            .start = down_j_max.points.len - 1,
            .end = 0,
        }, .{
            .edge = &out_j_max,
            .start = 0,
            .end = out_j_max.points.len - 1,
        }, .{
            .edge = &up_i_max_0,
            .start = 0,
            .end = up_i_max_0.points.len - 1,
        } });
        defer downstream_j_min.deinit();

        const downstream_x_00 = downstream_j_min.points[0];
        const downstream_x_01 = downstream_j_min.points[downstream_j_min.points.len - 1];

        // TODO: remove hard coding
        const downstream_x_10 = add(downstream_x_00, Vec2d.init(self.outlet_axial_position, 0.00));
        const downstream_x_11 = add(downstream_x_10, Vec2d.init(0.0, self.pitch));

        const downstream_j_max = try discrete.Edge.init(allocator, downstream_j_min.points.len, .{ .line = .{ .start = downstream_x_10, .end = downstream_x_11 } }, .{ .uniform = .{} });
        defer downstream_j_max.deinit();

        const downstream_i_min = try discrete.Edge.init(allocator, self.num_cells.downstream_i + 1, .{ .line = .{ .start = downstream_x_00, .end = downstream_x_10 } }, .{ .uniform = .{} });
        defer downstream_i_min.deinit();
        const downstream_i_max = try discrete.Edge.init(allocator, self.num_cells.downstream_i + 1, .{ .line = .{ .start = downstream_x_01, .end = downstream_x_11 } }, .{ .uniform = .{} });
        defer downstream_i_max.deinit();

        const downstream = try discrete.Block2d.init(allocator, downstream_i_min, downstream_i_max, downstream_j_min, downstream_j_max);

        try mesh.addBlock("downstream", downstream);
        const downstream_id = mesh.blocks.items.len - 1;

        // Connections

        try mesh.connections.appendSlice(&.{
            boundary.Connection.init(.{
                .{ .block = ss_id, .side = boundary.Side.j_min, .start = 0, .end = self.num_cells.o_grid },
                .{ .block = ps_id, .side = boundary.Side.j_min, .start = 0, .end = self.num_cells.o_grid },
            }, null),
            boundary.Connection.init(.{
                .{ .block = ss_id, .side = boundary.Side.j_max, .start = 0, .end = self.num_cells.o_grid },
                .{ .block = ps_id, .side = boundary.Side.j_max, .start = 0, .end = self.num_cells.o_grid },
            }, null),

            boundary.Connection.init(.{
                .{ .block = down_id, .side = boundary.Side.j_min, .start = self.num_cells.down_j, .end = 0 },
                .{ .block = upstream_id, .side = boundary.Side.j_max, .start = 0, .end = self.num_cells.bulge },
            }, null),
            boundary.Connection.init(.{
                .{ .block = in_id, .side = boundary.Side.j_max, .start = in_j_min.points.len - 1, .end = 0 },
                .{ .block = upstream_id, .side = boundary.Side.j_max, .start = self.num_cells.bulge, .end = self.num_cells.bulge + in_j_min.points.len - 1 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = in_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.in_i },
                .{ .block = down_id, .side = boundary.Side.i_min, .start = self.num_cells.in_i, .end = 0 },
            }, null),

            boundary.Connection.init(.{
                .{ .block = up_id, .side = boundary.Side.j_max, .start = 0, .end = self.num_cells.out_i },
                .{ .block = upstream_id, .side = boundary.Side.j_max, .start = self.num_cells.bulge + in_j_min.points.len - 1, .end = upstream_j_max.points.len - 1 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = in_id, .side = boundary.Side.i_min, .start = 0, .end = self.num_cells.in_i },
                .{ .block = up_id, .side = boundary.Side.i_min, .start = up_i_min.points.len - self.num_cells.in_i - 1, .end = up_i_min.points.len - 1 },
            }, null),

            boundary.Connection.init(.{
                .{ .block = down_id, .side = boundary.Side.j_max, .start = self.num_cells.down_j, .end = 0 },
                .{ .block = downstream_id, .side = boundary.Side.j_min, .start = 0, .end = self.num_cells.down_j },
            }, null),
            boundary.Connection.init(.{
                .{ .block = out_id, .side = boundary.Side.j_max, .start = 0, .end = out_j_max.points.len - 1 },
                .{ .block = downstream_id, .side = boundary.Side.j_min, .start = self.num_cells.down_j, .end = self.num_cells.down_j + out_j_max.points.len - 1 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = out_id, .side = boundary.Side.i_min, .start = 0, .end = self.num_cells.out_i },
                .{ .block = down_id, .side = boundary.Side.i_min, .start = down_i_min.points.len - 1 - self.num_cells.out_i, .end = down_i_min.points.len - 1 },
            }, null),

            boundary.Connection.init(.{
                .{ .block = out_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.out_i },
                .{ .block = up_id, .side = boundary.Side.j_min, .start = 0, .end = self.num_cells.out_i },
            }, null),
            boundary.Connection.init(.{
                .{ .block = up_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.bulge },
                .{ .block = downstream_id, .side = boundary.Side.j_min, .start = downstream_j_min.points.len - 1 - self.num_cells.bulge, .end = downstream_j_min.points.len - 1 },
            }, null),

            boundary.Connection.init(.{
                .{ .block = ss_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.in_up_j },
                .{ .block = in_id, .side = boundary.Side.j_min, .start = self.num_cells.in_up_j, .end = 0 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = ss_id, .side = boundary.Side.i_max, .start = self.num_cells.in_up_j, .end = self.num_cells.in_up_j + self.num_cells.middle_i + self.num_cells.bulge + self.num_cells.out_i },
                .{ .block = up_id, .side = boundary.Side.i_min, .start = up_i_min.points.len - 1 - self.num_cells.in_i, .end = 0 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = ss_id, .side = boundary.Side.i_max, .start = self.num_cells.in_up_j + self.num_cells.bulge + self.num_cells.middle_i + self.num_cells.out_i, .end = ss_i_max.points.len - 1 },
                .{ .block = out_id, .side = boundary.Side.j_min, .start = out_j_min.points.len - 1, .end = self.num_cells.out_down_j },
            }, null),

            boundary.Connection.init(.{
                .{ .block = ps_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.in_down_j },
                .{ .block = in_id, .side = boundary.Side.j_min, .start = self.num_cells.in_up_j, .end = in_j_min.points.len - 1 },
            }, null),
            boundary.Connection.init(.{
                .{ .block = ps_id, .side = boundary.Side.i_max, .start = self.num_cells.in_down_j, .end = self.num_cells.in_down_j + self.num_cells.middle_i },
                .{ .block = down_id, .side = boundary.Side.i_min, .start = self.num_cells.in_i, .end = down_i_min.points.len - 1 - self.num_cells.out_i },
            }, null),
            boundary.Connection.init(.{
                .{ .block = ps_id, .side = boundary.Side.i_max, .start = self.num_cells.in_down_j + self.num_cells.middle_i, .end = ps_i_max.points.len - 1 },
                .{ .block = out_id, .side = boundary.Side.j_min, .start = 0, .end = self.num_cells.out_down_j },
            }, null),

            boundary.Connection.init(.{
                .{ .block = upstream_id, .side = boundary.Side.i_min, .start = 0, .end = self.num_cells.upstream_i },
                .{ .block = upstream_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.upstream_i },
            }, .{ .data = .{ 0, self.pitch } }),
            boundary.Connection.init(.{
                .{ .block = down_id, .side = boundary.Side.i_max, .start = 0, .end = down_i_max.points.len - 1 },
                .{ .block = up_id, .side = boundary.Side.i_max, .start = up_i_max.points.len - 1, .end = up_i_max.points.len - down_i_max.points.len },
            }, .{ .data = .{ 0, self.pitch } }),
            boundary.Connection.init(.{
                .{ .block = downstream_id, .side = boundary.Side.i_min, .start = 0, .end = self.num_cells.downstream_i },
                .{ .block = downstream_id, .side = boundary.Side.i_max, .start = 0, .end = self.num_cells.downstream_i },
            }, .{ .data = .{ 0, self.pitch } }),
        });

        // Boundary conditions
        // TODO: add boundary conditions

        // TODO: consider removing the Side definition... Check if this eases the code.

        try smooth.mesh(allocator, &mesh, self.smoothing.iterations, self.smoothing.backend, self.smoothing.wall_control_function);

        return mesh;
    }
};

fn projectNormal(allocator: std.mem.Allocator, edge: []Vec2d, distance: Float) ![]Vec2d {
    var projected = try allocator.alloc(types.Vec2d, edge.len);

    // handle all points except the first and last
    for (1..edge.len - 1) |i| {
        const x_i = edge[i];
        const x_im1 = edge[i - 1];
        const x_ip1 = edge[i + 1];

        const x_xi = scale(0.5, sub(x_ip1, x_im1));

        const n = scale(1.0 / abs(x_xi), Vec2d.init(x_xi.data[1], -x_xi.data[0]));

        projected[i] = add(x_i, scale(distance, n));
    }

    // handle first point
    {
        const i = 0;
        const x_i = edge[i];
        const x_ip1 = edge[i + 1];

        const x_xi = sub(x_ip1, x_i);

        const n = scale(1.0 / abs(x_xi), Vec2d.init(x_xi.data[1], -x_xi.data[0]));

        projected[i] = add(x_i, scale(distance, n));
    }

    // handle last point
    {
        const i = edge.len - 1;
        const x_i = edge[i];
        const x_im1 = edge[i - 1];

        const x_xi = sub(x_i, x_im1);

        const n = scale(1.0 / abs(x_xi), Vec2d.init(x_xi.data[1], -x_xi.data[0]));

        projected[i] = add(x_i, scale(distance, n));
    }

    return projected;
}

test "O4H template" {
    const allocator = std.testing.allocator;
    const template = O4H{
        .ps_csv_path = "./examples/T106/T106_ps.dat",
        .ss_csv_path = "./examples/T106/T106_ss.dat",
        .pitch = 0.08836, // m
        .blade_clustering = .{ .roberts = .{ .alpha = 0.5, .beta = 1.03 } },
        .num_cells = .{
            .o_grid = 40,
            .in_up_j = 50,
            .in_down_j = 10,
            .in_i = 10,
            .out_up_j = 40,
            .out_down_j = 10,
            .out_i = 10,
            .middle_i = 100,
            .down_j = 40,
            .bulge = 40,
            .upstream_i = 20,
            .downstream_i = 10,
        },
        .smoothing = .{
            .iterations = 10,
            .backend = .umfpack,
            .wall_control_function = .{
                .white = .{
                    .ds_target = 1e-6,
                    .theta_target = 1.570796327,
                },
            },
        },
    };

    var mesh = try template.run(allocator);
    defer mesh.deinit();

    // try mesh.write(allocator, "o4h_linear.cgns");

    // try smooth.mesh(allocator, &mesh, 10, .umfpack, .laplace);

    try mesh.write(allocator, "o4h.cgns");
}
