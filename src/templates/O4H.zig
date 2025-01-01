const std = @import("std");

const blade = @import("blade.zig");
const types = @import("../types.zig");
const clustering = @import("../clustering.zig");
const discete = @import("../discrete.zig");
const geometry = @import("../geometry.zig");

const Float = types.Float;
const Vec2d = types.Vec2d;
const add = types.add;
const sub = types.sub;
const scale = types.scale;
const abs = types.abs;

const Turbine = struct {
    ps_csv_path: []const u8,
    ss_csv_path: []const u8,
    pitch: types.Float,
    blade_clustering: clustering.Function,
    num_cells: struct {
        o_grid: types.Index,
        in_up_j: types.Index,
        middle_i: types.Index,
        out_up_j: types.Index,
        in_down_j: types.Index,
        out_down_j: types.Index,
    },

    fn run(self: *const Turbine, allocator: std.mem.Allocator) !void {

        // TODO add geometry and discrete entities to manager

        const num_cells_ss =
            self.num_cells.in_up_j + self.num_cells.middle_i + self.num_cells.out_up_j;
        const num_cells_ps =
            self.num_cells.in_down_j + self.num_cells.middle_i + self.num_cells.out_down_j;

        const profile = try blade.Profile.init(allocator, self.ps_csv_path, self.ss_csv_path);
        defer profile.deinit();

        var ps_edge = try discete.Edge.init(allocator, num_cells_ps + 1, .{ .spline = profile.pressure_side }, self.blade_clustering);
        defer ps_edge.deinit();

        var ss_edge = try discete.Edge.init(allocator, num_cells_ss + 1, .{ .spline = profile.suction_side }, self.blade_clustering);
        defer ss_edge.deinit();

        // o - grid for viscous computations

        // compute o-grid target by projecting the blade normal outward

        // TODO make this runtime dependent
        // TODO replace with a percentage value of the chord length
        const d = 0.001;

        const ps_outer_edge = discete.Edge{ .allocator = allocator, .points = try projectNormal(allocator, ps_edge.points[0..], d), .clustering = try allocator.dupe(Float, ps_edge.clustering) };
        defer ps_outer_edge.deinit();

        const ss_outer_edge = discete.Edge{ .allocator = allocator, .points = try projectNormal(allocator, ss_edge.points[0..], d), .clustering = try allocator.dupe(Float, ss_edge.clustering) };
        defer ss_outer_edge.deinit();

        // split blade distribution

        // TODO

        var mesh = discete.Mesh.init(allocator);
        defer mesh.deinit();

        // |          /
        // |         |          ^
        // |        |           |
        // |-------< x_00      i_min
        //         LE

        const ss_j_min_edge = try discete.Edge.init(
            allocator,
            self.num_cells.o_grid + 1,
            .{ .line = geometry.Line.init(ss_edge.points[0], ss_outer_edge.points[0]) },
            .{ .uniform = .{} },
        );
        defer ss_j_min_edge.deinit();
        const ss_j_max_edge = try discete.Edge.init(
            allocator,
            self.num_cells.o_grid + 1,
            .{ .line = geometry.Line.init(ss_edge.points[ss_edge.points.len - 1], ss_outer_edge.points[ss_outer_edge.points.len - 1]) },
            .{ .uniform = .{} },
        );
        defer ss_j_max_edge.deinit();

        const ss = try discete.Block2d.init(allocator, ss_edge, ss_outer_edge, ss_j_min_edge, ss_j_max_edge);

        try mesh.addBlock("ss", ss);

        try mesh.write("o4h.cgns");
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

test "turbine template" {
    const allocator = std.testing.allocator;
    const template = Turbine{
        .ps_csv_path = "./examples/T106/T106_ps.dat",
        .ss_csv_path = "./examples/T106/T106_ss.dat",
        .pitch = 0.08836, // m
        .blade_clustering = .{ .roberts = .{ .alpha = 0.5, .beta = 1.03 } },
        .num_cells = .{
            .o_grid = 33,
            .in_up_j = 95,
            .in_down_j = 20,
            .out_up_j = 210,
            .out_down_j = 20,
            .middle_i = 120,
        },
    };

    try template.run(allocator);
}
