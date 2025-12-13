const std = @import("std");

const csv = @import("../csv.zig");
const spline = @import("../spline.zig");
const types = @import("../types.zig");

pub const Profile = struct {
    pressure_side: spline.FittingSpline(2),
    suction_side: spline.FittingSpline(2),

    pub fn init(allocator: std.mem.Allocator, ps_csv_path: []const u8, ss_csv_path: []const u8) !Profile {
        const ps = try readSide(allocator, ps_csv_path);
        defer ps.deinit();

        const ss = try readSide(allocator, ss_csv_path);
        defer ss.deinit();

        if (!types.eql(ps.items[0], ss.items[0])) {
            std.log.err("Leading edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingLeadingEdge;
        }

        if (!types.eql(ps.items[ps.items.len - 1], ss.items[ss.items.len - 1])) {
            std.log.err("Trailing edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingTrailingEdge;
        }

        std.debug.assert(ps.items.len > 1);
        std.debug.assert(ps.items[0].data[0] < ps.items[ps.items.len - 1].data[0]);

        const pressure_side = try spline.FittingSpline(2).init(allocator, types.cast(ps.items[0..]), 3);
        const suction_side = try spline.FittingSpline(2).init(allocator, types.cast(ss.items[0..]), 3);

        return .{ .pressure_side = pressure_side, .suction_side = suction_side };
    }

    pub fn deinit(self: *const Profile) void {
        self.pressure_side.deinit();
        self.suction_side.deinit();
    }
};

fn readSide(allocator: std.mem.Allocator, csv_path: []const u8) !std.array_list.Managed(types.Vec2d) {
    const side = try csv.parseCsvIntoVec2d(allocator, csv_path);
    // TODO add check of axial direction: does the vector form origion to end point in upstream direction (?)
    if (side.items[0].data[0] > side.items[side.items.len - 1].data[0]) {
        std.mem.reverse(types.Vec2d, side.items);
    }
    return side;
}

test "read a blade profile" {
    const allocator = std.testing.allocator;
    const profile = try Profile.init(allocator, "./examples/T106/T106_ps.dat", "./examples/T106/T106_ss.dat");
    profile.deinit();
}
