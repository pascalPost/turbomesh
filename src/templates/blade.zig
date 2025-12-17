// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");

const csv = @import("../csv.zig");
const spline = @import("../spline.zig");
const types = @import("../types.zig");

pub const Profile = struct {
    down_part: spline.FittingSpline(2),
    up_part: spline.FittingSpline(2),

    pub fn init(allocator: std.mem.Allocator, down_csv_path: []const u8, up_csv_path: []const u8) !Profile {
        const down = try readSide(allocator, down_csv_path);
        defer down.deinit();

        const up = try readSide(allocator, up_csv_path);
        defer up.deinit();

        if (!types.eql(down.items[0], up.items[0])) {
            std.log.err("Leading edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingLeadingEdge;
        }

        if (!types.eql(down.items[down.items.len - 1], up.items[up.items.len - 1])) {
            std.log.err("Trailing edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingTrailingEdge;
        }

        std.debug.assert(down.items.len > 1);
        std.debug.assert(down.items[0].data[0] < down.items[down.items.len - 1].data[0]);

        const down_part = try spline.FittingSpline(2).init(allocator, types.cast(down.items[0..]), 3);
        const up_part = try spline.FittingSpline(2).init(allocator, types.cast(up.items[0..]), 3);

        return .{ .down_part = down_part, .up_part = up_part };
    }

    pub fn deinit(self: *const Profile) void {
        self.down_part.deinit();
        self.up_part.deinit();
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
