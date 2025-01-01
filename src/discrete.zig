const std = @import("std");
const types = @import("types.zig");
const cluster = @import("clustering.zig");
const geometry = @import("geometry.zig");
const tfi = @import("tfi.zig");
const cgns = @import("cgns.zig");

pub const Edge = struct {
    allocator: std.mem.Allocator,
    points: []types.Vec2d,
    clustering: []types.Float,

    pub fn init(allocator: std.mem.Allocator, n: types.Index, curve: geometry.Curve, clustering: cluster.Function) !Edge {
        const u = try cluster.create(allocator, clustering, n);

        var data = try allocator.alloc(types.Vec2d, n);
        switch (curve) {
            .line => |l| try l.interpolate(u, data[0..]),
            .spline => |s| try s.interpolate(u, types.cast(data[0..])),
        }

        return .{
            .allocator = allocator,
            .points = data,
            .clustering = u,
        };
    }

    pub fn deinit(self: Edge) void {
        self.allocator.free(self.clustering);
        self.allocator.free(self.points);
    }
};

pub const Block2d = struct {
    allocator: std.mem.Allocator,
    points: types.Mat2d,

    pub fn init(allocator: std.mem.Allocator, i_min: Edge, i_max: Edge, j_min: Edge, j_max: Edge) !Block2d {
        std.debug.assert(i_min.points.len == i_max.points.len);
        std.debug.assert(j_min.points.len == j_max.points.len);

        var points = try types.Mat2d.init(allocator, .{ i_min.points.len, j_min.points.len });
        tfi.linearBoundaryBlendedControlFunction(&points, i_min.points, i_max.points, j_min.points, j_max.points, i_min.clustering, i_max.clustering, j_min.clustering, j_max.clustering);

        // std.debug.print("{any}\n", .{points.data});

        return .{ .allocator = allocator, .points = points };
    }

    pub fn deinit(self: Block2d) void {
        self.points.deinit(self.allocator);
    }
};

// const Shell2d = struct {
//     // consists of 4 2d edges
// };

pub const Mesh = struct {
    blocks: std.ArrayList(Block2d),
    names: std.ArrayList([:0]const u8),

    // vec of boundaries

    pub fn init(allocator: std.mem.Allocator) Mesh {
        return .{
            .blocks = std.ArrayList(Block2d).init(allocator),
            .names = std.ArrayList([:0]const u8).init(allocator),
        };
    }

    pub fn addBlock(self: *Mesh, name: []const u8, block: Block2d) !void {
        try self.blocks.append(block);
        try self.names.append(try self.names.allocator.dupeZ(u8, name));
    }

    pub fn deinit(self: Mesh) void {
        for (self.blocks.items) |b| b.deinit();
        for (self.names.items) |n| self.names.allocator.free(n);
        self.blocks.deinit();
        self.names.deinit();
    }

    pub fn write(self: Mesh, allocator: std.mem.Allocator, filename: [:0]const u8) !void {

        // buffer
        var size: usize = 0;
        for (self.blocks.items) |b| {
            size = @max(size, b.points.data.len);
        }
        var buffer = try allocator.alloc(types.Float, size);
        defer allocator.free(buffer);

        // block data
        var block_points = try allocator.alloc(types.Mat2d, self.blocks.items.len);
        defer allocator.free(block_points);
        for (self.blocks.items, block_points[0..]) |b, *data| data.* = b.points;

        try cgns.write(filename, self.names.items, block_points, buffer[0..]);
    }
};
