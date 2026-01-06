// Copyright (c) 2026 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const types = @import("types.zig");
const spline = @import("spline.zig");

pub const Geometry = struct {
    pitch: types.Float,
    profile: Profile,

    pub fn init(pitch: types.Float, profile: Profile) Geometry {
        return .{ .pitch = pitch, .profile = profile };
    }
};

pub const Profile = struct {
    down_part: spline.FittingSpline(2),
    up_part: spline.FittingSpline(2),

    pub fn init(allocator: std.mem.Allocator, down: []const types.Vec2d, up: []const types.Vec2d) !Profile {
        if (!types.eql(down[0], up[0])) {
            std.log.err("Leading edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingLeadingEdge;
        }

        if (!types.eql(down[down.len - 1], up[up.len - 1])) {
            std.log.err("Trailing edge of suction and pressure side must be equal.", .{});
            return error.NonMatchingTrailingEdge;
        }

        std.debug.assert(down.len > 1);
        std.debug.assert(down[0].data[0] < down[down.len - 1].data[0]);

        const down_part = try spline.FittingSpline(2).init(allocator, types.castConst(down[0..]), 3);
        const up_part = try spline.FittingSpline(2).init(allocator, types.castConst(up[0..]), 3);

        return .{ .down_part = down_part, .up_part = up_part };
    }

    pub fn deinit(self: *const Profile) void {
        self.down_part.deinit();
        self.up_part.deinit();
    }
};
