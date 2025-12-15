// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const types = @import("types.zig");
const spline = @import("spline.zig");

const Tag = enum {
    line,
    spline,
};

pub const Curve = union(Tag) {
    line: Line,
    spline: spline.FittingSpline(2),
};

pub const Line = struct {
    start: types.Vec2d,
    end: types.Vec2d,

    pub fn init(start: types.Vec2d, end: types.Vec2d) Line {
        return .{ .start = start, .end = end };
    }

    pub fn interpolate(self: *const Line, clustering: []const f64, values: []types.Vec2d) !void {
        if (clustering.len != values.len) {
            std.debug.print("Mismatch of slice length ({} != {})", .{ clustering.len, values.len });
            return error.Mismatch;
        }

        std.debug.assert(clustering[0] == 0.0);
        std.debug.assert(clustering[clustering.len - 1] == 1.0);

        const dx = types.sub(self.end, self.start);

        for (values, clustering) |*v, u| {
            v.* = types.add(self.start, types.scale(u, dx));
        }
    }
};
