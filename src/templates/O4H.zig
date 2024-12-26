const std = @import("std");

const blade = @import("blade.zig");
const types = @import("../types.zig");
const clustering = @import("../clustering.zig");

const Turbine = struct {
    ps_csv_path: []const u8,
    ss_csv_path: []const u8,
    pitch: types.Float,
    blade_clustering: clustering.Function,
    num_cells_blade_half: types.Index,

    fn run(self: *const Turbine, allocator: std.mem.Allocator) !void {

        // TODO add geometry and discrete entities to manager

        const profile = try blade.Profile.init(allocator, self.ps_csv_path, self.ss_csv_path);
        defer profile.deinit();

        const blade_clustering = try clustering.create(allocator, self.blade_clustering, self.num_cells_blade_half + 1);
        defer allocator.free(blade_clustering);

        var ps_edge = try allocator.alloc(types.Vec2d, self.num_cells_blade_half + 1);
        defer allocator.free(ps_edge);
        try profile.pressure_side.interpolate(blade_clustering, types.cast(ps_edge[0..]));

        // const file = try std.fs.cwd().createFile("ps.csv", .{});
        // defer file.close();
        // const writer = file.writer();
        // for (ps_edge) |x| {
        //     try writer.print("{} {}\n", .{ x.data[0], x.data[1] });
        // }

        var ss_edge = try allocator.alloc(types.Vec2d, self.num_cells_blade_half + 1);
        defer allocator.free(ss_edge);
        try profile.suction_side.interpolate(blade_clustering, types.cast(ss_edge[0..]));

        // split blade distribution

    }
};

test "turbine template" {
    const allocator = std.testing.allocator;
    const template = Turbine{
        .ps_csv_path = "./examples/T106/T106_ps.dat",
        .ss_csv_path = "./examples/T106/T106_ss.dat",
        .pitch = 0.08836, // m
        .blade_clustering = .{ .roberts = .{ .alpha = 0.5, .beta = 1.03 } },
        .num_cells_blade_half = 30,
    };

    try template.run(allocator);
}
