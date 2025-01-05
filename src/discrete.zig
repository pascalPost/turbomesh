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

    pub fn combine(allocator: std.mem.Allocator, edges: []const EdgeView) !Edge {
        std.debug.assert(edges.len > 1);

        const tol = 1e-10;

        // check if edges can be merged
        {
            for (edges[0 .. edges.len - 1], edges[1..], 1..) |edge, edge_next, i| {
                if (!types.eqlApprox(edge.edge.points[edge.end], edge_next.edge.points[edge_next.start], tol)) {
                    std.debug.print("edges {} and {} cannot be combined as end points do not match: ({}, {}) and ({}, {})\n", .{
                        i,                                              i + 1,
                        edge.edge.points[edge.end].data[0],             edge.edge.points[edge.end].data[1],
                        edge_next.edge.points[edge_next.start].data[0], edge_next.edge.points[edge_next.start].data[1],
                    });
                    unreachable;
                }
            }
        }

        var n: usize = 0;
        for (edges) |e| {
            n += e.len();
        }
        n -= edges.len - 1;

        const u = try allocator.alloc(types.Float, n);
        const points = try allocator.alloc(types.Vec2d, n);

        // points
        {
            var start: usize = 0;
            for (edges) |edge| start += edge.clonePoints(points[start..]) - 1;
        }

        // cluster
        {
            // cum sum
            var start: usize = 0;
            var last_value: types.Float = 0.0;
            for (edges) |edge| {
                start += edge.cloneClustering(u[start..], last_value) - 1;
                last_value = u[start];
            }

            // rel to [0..1]
            for (u[0..]) |*v| v.* /= last_value;
        }

        return .{
            .allocator = allocator,
            .points = points,
            .clustering = u,
        };
    }
};

pub const EdgeView = struct {
    edge: *const Edge,
    start: usize,
    end: usize,

    fn len(self: EdgeView) usize {
        if (self.start > self.end) {
            return self.start - self.end + 1;
        }
        return self.end - self.start + 1;
    }

    fn clonePoints(self: EdgeView, buffer: []types.Vec2d) usize {
        if (self.start > self.end) {
            const length = self.start - self.end + 1;
            @memcpy(buffer[0..length], self.edge.points[self.end .. self.start + 1]);
            std.mem.reverse(types.Vec2d, buffer[0..length]);
            return length;
        }

        const length = self.end - self.start + 1;
        @memcpy(buffer[0..length], self.edge.points[self.start .. self.end + 1]);
        return length;
    }

    fn cloneClustering(self: EdgeView, buffer: []types.Float, initial_value: types.Float) usize {
        buffer[0] = initial_value;

        const first = if (self.start > self.end) self.end else self.start;
        const last = if (self.start > self.end) self.start else self.end;

        const last_value = self.edge.clustering[first];
        var i = first + 1;
        var i_buf: usize = 1;
        while (i <= last) : (i += 1) {
            const value = self.edge.clustering[i];
            const delta = value - last_value;
            buffer[i_buf] = initial_value + delta;
            i_buf += 1;
        }
        return i_buf; // length
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

test "combining edges" {
    const allocator = std.testing.allocator;
    const edge_1 = try Edge.init(allocator, 3, .{ .line = .{ .start = types.Vec2d.init(0.0, 0.0), .end = types.Vec2d.init(2.0, 0.0) } }, .{ .uniform = .{} });
    defer edge_1.deinit();

    const edge_2 = try Edge.init(allocator, 3, .{ .line = .{ .start = types.Vec2d.init(2.0, 0.0), .end = types.Vec2d.init(4.0, 0.0) } }, .{ .uniform = .{} });
    defer edge_2.deinit();

    {
        const edge = try Edge.combine(allocator, &.{
            .{ .edge = &edge_1, .start = 0, .end = 2 },
            .{ .edge = &edge_2, .start = 0, .end = 2 },
        });
        defer edge.deinit();

        try std.testing.expectEqualDeep(edge.points, ([_]types.Vec2d{
            types.Vec2d.init(0.0, 0.0),
            types.Vec2d.init(1.0, 0.0),
            types.Vec2d.init(2.0, 0.0),
            types.Vec2d.init(3.0, 0.0),
            types.Vec2d.init(4.0, 0.0),
        })[0..]);
        try std.testing.expectEqualDeep(edge.clustering, ([_]types.Float{ 0.0, 0.25, 0.5, 0.75, 1.0 })[0..]);
    }

    {
        const edge = try Edge.combine(allocator, &.{
            .{ .edge = &edge_1, .start = 1, .end = 2 },
            .{ .edge = &edge_2, .start = 0, .end = 1 },
        });
        defer edge.deinit();

        try std.testing.expectEqualDeep(edge.points, ([_]types.Vec2d{
            types.Vec2d.init(1.0, 0.0),
            types.Vec2d.init(2.0, 0.0),
            types.Vec2d.init(3.0, 0.0),
        })[0..]);
        try std.testing.expectEqualDeep(edge.clustering, ([_]types.Float{ 0.0, 0.5, 1.0 })[0..]);
    }

    {
        const edge = try Edge.combine(allocator, &.{
            .{ .edge = &edge_2, .start = 2, .end = 0 },
            .{ .edge = &edge_1, .start = 2, .end = 0 },
        });
        defer edge.deinit();

        try std.testing.expectEqualDeep(edge.points, ([_]types.Vec2d{
            types.Vec2d.init(4.0, 0.0),
            types.Vec2d.init(3.0, 0.0),
            types.Vec2d.init(2.0, 0.0),
            types.Vec2d.init(1.0, 0.0),
            types.Vec2d.init(0.0, 0.0),
        })[0..]);
        try std.testing.expectEqualDeep(edge.clustering, ([_]types.Float{ 0.0, 0.25, 0.5, 0.75, 1.0 })[0..]);
    }

    {
        const edge = try Edge.combine(allocator, &.{
            .{ .edge = &edge_2, .start = 1, .end = 0 },
            .{ .edge = &edge_1, .start = 2, .end = 1 },
        });
        defer edge.deinit();

        try std.testing.expectEqualDeep(edge.points, ([_]types.Vec2d{
            types.Vec2d.init(3.0, 0.0),
            types.Vec2d.init(2.0, 0.0),
            types.Vec2d.init(1.0, 0.0),
        })[0..]);
        try std.testing.expectEqualDeep(edge.clustering, ([_]types.Float{ 0.0, 0.5, 1.0 })[0..]);
    }
}
