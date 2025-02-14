const std = @import("std");
const types = @import("types.zig");
const discrete = @import("discrete.zig");

pub const Side = enum {
    i_min,
    i_max,
    j_min,
    j_max,
};

pub const Range = struct {
    block: usize,
    side: Side,
    start: usize,
    end: usize,

    pub fn len(self: Range) usize {
        if (self.start > self.end) {
            return self.start - self.end + 1;
        }
        return self.end - self.start + 1;
    }

    pub fn iterate(self: Range, mesh: *const discrete.Mesh) RangeIterator {
        const points = &mesh.blocks.items[self.block].points;
        const size = points.size;

        switch (self.side) {
            .i_min => {
                const start: isize = @intCast(points.index(.{ self.start, 0 }));
                const end: isize = @intCast(points.index(.{ self.end, 0 }));
                var step: isize = @intCast(size[1]);
                if (start > end) step = -step;
                return .{ .pos = start, .end = end + step, .step = step };
            },
            .j_max => {
                const start: isize = @intCast(points.index(.{ size[0] - 1, self.start }));
                const end: isize = @intCast(points.index(.{ size[0] - 1, self.end }));
                var step: isize = 1;
                if (start > end) step = -step;
                return .{ .pos = start, .end = end + step, .step = step };
            },
            .i_max => {
                const start: isize = @intCast(points.index(.{ self.start, size[1] - 1 }));
                const end: isize = @intCast(points.index(.{ self.end, size[1] - 1 }));
                var step: isize = @intCast(size[1]);
                if (start > end) step = -step;
                return .{ .pos = start, .end = end + step, .step = step };
            },
            .j_min => {
                const start: isize = @intCast(points.index(.{ 0, self.start }));
                const end: isize = @intCast(points.index(.{ 0, self.end }));
                var step: isize = 1;
                if (start > end) step = -step;
                return .{ .pos = start, .end = end + step, .step = step };
            },
        }
    }
};

const RangeIterator = struct {
    pos: isize,
    end: isize,
    step: isize,

    pub fn next(self: *@This()) ?usize {
        if (self.pos == self.end) return null;
        const pos = self.pos;
        self.pos += self.step;
        return @intCast(pos);
    }
};

pub const Connection = struct {
    data: [2]Range,

    pub fn len(self: Connection) usize {
        const length = self.data[0].len();
        std.debug.assert(length == self.data[1].len());
        return length;
    }

    pub fn lenInternal(self: Connection) usize {
        const length = self.len();
        std.debug.assert(length > 2);
        return length - 2;
    }

    fn internalRange(self: Connection) [2]Range {
        var new_range: [2]Range = undefined;
        @memcpy(new_range[0..], self.data);

        for (new_range) |r| {
            if (r.start < r.end) {
                r.start += 1;
                r.end -= 1;
            } else {
                r.end += 1;
                r.start -= 1;
            }
        }
        return new_range;
    }
};

pub const Periodic = struct {
    a: Range,
    b: Range,
    periodicity: types.Vec2d,
};

const ConditionTag = enum {
    wall,
    inlet,
    outlet,
};

pub const Condition = struct {
    range: Range,
    kind: ConditionTag,
};

/// PointData is a helper struct to manage the boundary points of a mesh providing a flat buffer of the block boundary points.
///
/// example block of size 6x4:
///          (0,3) (1,3) (2,3) (3,3) (4,3) (5,3)
///          13    12    11    10    09    08
/// (0,3)    x     x     x     x     x     x     (5,3)
/// (0,2) 14 x                             x  07 (5,2)
/// (0,1) 15 x                             x  06 (5,1)
/// (0,0)    x     x     x     x     x     x     (5,0)
///          00    01    02    03    04    05
///         (0,0) (1,0) (2,0) (3,0) (4,0) (5,0)
///
/// ```
pub fn PointData(comptime T: type) type {
    return struct {
        allocator: std.mem.Allocator,
        buffer: []T,
        mesh: *const discrete.Mesh,

        pub fn init(allocator: std.mem.Allocator, mesh_data: *const discrete.Mesh) !PointData(T) {
            var num_boundary_points: usize = 0;
            for (mesh_data.blocks.items) |b| {
                const points = b.points;
                const size_i = points.size[0];
                const size_j = points.size[1];

                num_boundary_points += 2 * size_i + 2 * size_j - 4;
            }

            const buffer = try allocator.alloc(T, num_boundary_points);

            return .{
                .allocator = allocator,
                .buffer = buffer,
                .mesh = mesh_data,
            };
        }

        pub fn deinit(self: @This()) void {
            self.allocator.free(self.buffer);
        }

        pub fn iterateRange(self: *@This(), range: Range) PointDataIterator(T) {
            // move to the correct block
            var start: usize = 0;

            for (self.mesh.blocks.items[0..range.block]) |b| {
                const points = b.points;
                const size_i = points.size[0];
                const size_j = points.size[1];

                const num_boundary_points = 2 * size_i + 2 * size_j - 4;

                start += num_boundary_points;
            }

            const block_size = self.mesh.blocks.items[range.block].points.size;

            var end: usize = 0;
            switch (range.side) {
                .i_min => {
                    end = start + block_size[0];
                },
                .j_max => {
                    start += block_size[0] - 1;
                    end = start + block_size[1];
                },
                .i_max => {
                    start += block_size[0] + block_size[1] - 2;
                    end = start + block_size[0];
                },
                .j_min => {
                    start += 2 * block_size[0] + block_size[1] - 3;
                    end = start + block_size[1];
                },
            }

            if (range.start > range.end) {
                return .{ .buffer = self.buffer[start..end], .end = range.start, .pos = @intCast(range.end) };
            }

            return .{ .buffer = self.buffer[start..end], .end = range.end, .pos = @intCast(range.start) };
        }
    };
}

fn PointDataIterator(comptime T: type) type {
    return struct {
        buffer: []T,
        end: usize,
        pos: isize,

        pub fn next(self: *@This()) ?T {
            if (self.pos == self.end) return null;
            const pos = self.pos;
            self.pos += 1;
            return self.buffer[pos];
        }

        pub fn nextPtr(self: *@This()) ?*T {
            if (self.pos == self.end) return null;
            const pos = self.pos;
            self.pos += 1;
            return &self.buffer[@intCast(pos)];
        }
    };
}

test "point data" {
    const allocator = std.testing.allocator;
    var mesh = discrete.Mesh.init(allocator);
    defer mesh.deinit();

    const x_00 = types.Vec2d.init(0.0, 0.0);
    const x_01 = types.Vec2d.init(0.0, 1.0);
    const x_10 = types.Vec2d.init(1.0, 0.0);
    const x_11 = types.Vec2d.init(1.0, 1.0);

    const size = .{ 6, 4 };

    const i_min = try discrete.Edge.init(allocator, size[0], .{ .line = .{ .start = x_00, .end = x_10 } }, .{ .uniform = .{} });
    defer i_min.deinit();

    const i_max = try discrete.Edge.init(allocator, size[0], .{ .line = .{ .start = x_01, .end = x_11 } }, .{ .uniform = .{} });
    defer i_max.deinit();

    const j_min = try discrete.Edge.init(allocator, size[1], .{ .line = .{ .start = x_00, .end = x_01 } }, .{ .uniform = .{} });
    defer j_min.deinit();

    const j_max = try discrete.Edge.init(allocator, size[1], .{ .line = .{ .start = x_10, .end = x_11 } }, .{ .uniform = .{} });
    defer j_max.deinit();

    const block = try discrete.Block2d.init(allocator, i_min, i_max, j_min, j_max);

    try mesh.addBlock("block_0", block);

    var boundary_point_connections = try PointData(std.BoundedArray(usize, 4)).init(allocator, &mesh);
    defer boundary_point_connections.deinit();

    {
        const range = Range{ .block = 0, .side = .i_min, .start = 0, .end = 5 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 0, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 1, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 2, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 3, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 4, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 5, 0 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const range = Range{ .block = 0, .side = .i_min, .start = 5, .end = 0 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 5, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 4, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 3, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 2, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 1, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 0 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const size_j = block.points.size[1];
        const range = Range{ .block = 0, .side = .i_max, .start = 0, .end = 5 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 0, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 1, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 2, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 3, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 4, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 5, size_j - 1 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const size_j = block.points.size[1];
        const range = Range{ .block = 0, .side = .i_max, .start = 5, .end = 0 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 5, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 4, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 3, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 2, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 1, size_j - 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, size_j - 1 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const range = Range{ .block = 0, .side = .j_min, .start = 0, .end = 3 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 0, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 2 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 3 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const range = Range{ .block = 0, .side = .j_min, .start = 3, .end = 0 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ 0, 3 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 2 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ 0, 0 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const size_i = block.points.size[0];
        const range = Range{ .block = 0, .side = .j_max, .start = 0, .end = 3 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 0 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 2 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 3 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        const size_i = block.points.size[0];
        const range = Range{ .block = 0, .side = .j_max, .start = 3, .end = 0 };
        var it = range.iterate(&mesh);
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 3 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 2 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 1 }), it.next());
        try std.testing.expectEqual(block.points.index(.{ size_i - 1, 0 }), it.next());
        try std.testing.expectEqual(null, it.next());
    }

    {
        var it = boundary_point_connections.iterateRange(.{ .block = 0, .side = .i_min, .start = 0, .end = 5 });
        var idx: usize = 0;
        while (it.nextPtr()) |*value| {
            idx += 1;
            try value.*.append(0);
        }
        try std.testing.expectEqual(idx, 5);
    }

    {
        var it = boundary_point_connections.iterateRange(.{ .block = 0, .side = .j_max, .start = 0, .end = 3 });
        var idx: usize = 0;
        while (it.nextPtr()) |*value| {
            idx += 1;
            try value.*.append(1);
        }
        try std.testing.expectEqual(idx, 3);
    }

    {
        var it = boundary_point_connections.iterateRange(.{ .block = 0, .side = .i_max, .start = 5, .end = 0 });
        var idx: usize = 0;
        while (it.nextPtr()) |*value| {
            idx += 1;
            try value.*.append(2);
        }
        try std.testing.expectEqual(idx, 5);
    }

    {
        var it = boundary_point_connections.iterateRange(.{ .block = 0, .side = .i_max, .start = 0, .end = 5 });
        var idx: usize = 0;
        while (it.nextPtr()) |*value| {
            idx += 1;
            try value.*.append(2);
        }
        try std.testing.expectEqual(idx, 5);
    }

    {
        var it = boundary_point_connections.iterateRange(.{ .block = 0, .side = .j_max, .start = 3, .end = 0 });
        var idx: usize = 0;
        while (it.nextPtr()) |*value| {
            idx += 1;
            try value.*.append(3);
        }
        try std.testing.expectEqual(idx, 3);
    }
}
