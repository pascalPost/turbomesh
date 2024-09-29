const std = @import("std");
const types = @import("types.zig");
const Vec2d = types.Vec2d;

const Tag = enum(u1) { line };
const Curve = union(Tag) { line: *Line2d };

fn discretize() void {}

// TODO add Spline

pub const Line2d = struct { start: Vec2d, end: Vec2d };
