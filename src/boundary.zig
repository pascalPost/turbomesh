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
        const points = mesh.blocks.items[self.block].points;
        const size = points.size;

        var iterator: RangeIterator = undefined;

        switch (self.side) {
            .i_min => {
                iterator.idx = points.index(.{ self.start, 0 });
                iterator.incerment = @intCast(size[1]);
            },
            .j_max => {
                iterator.idx = points.index(.{ size[0] - 1, self.start });
                iterator.incerment = 1;
            },
            .i_max => {
                iterator.idx = points.index(.{ self.start, size[1] - 1 });
                iterator.incerment = @intCast(size[1]);
            },
            .j_min => {
                iterator.idx = points.index(.{ 0, self.start });
                iterator.incerment = 1;
            },
        }

        if (self.start > self.end) {
            iterator.incerment = -iterator.incerment;
            iterator.count = self.start - self.end;
        } else {
            iterator.count = self.end - self.start;
        }

        return iterator;
    }
};

const RangeIterator = struct {
    // for the last element, count will be 0 and idx will be != null
    // if no elements are left, idx will be null
    count: usize,
    idx: ?usize,
    incerment: isize,

    pub fn next(self: *@This()) ?usize {
        const idx = self.idx;
        if (self.count > 0) {
            self.count -= 1;
            self.idx = @bitCast(@as(isize, @intCast(idx.?)) + self.incerment);
        } else {
            self.idx = null;
        }
        return idx;
    }
};

pub const Connection = struct {
    ranges: [2]Range,

    pub fn init(ranges: [2]Range) Connection {
        return .{ .ranges = ranges };
    }

    pub fn len(self: Connection) usize {
        const length = self.ranges[0].len();
        std.debug.assert(length == self.ranges[1].len());
        return length;
    }

    pub fn lenInternal(self: Connection) usize {
        const length = self.len();
        std.debug.assert(length > 2);
        return length - 2;
    }

    pub fn iterate(self: Connection, mesh: *const discrete.Mesh) ConnectionIterator {
        const it_0 = self.ranges[0].iterate(mesh);
        const it_1 = self.ranges[1].iterate(mesh);
        return .{ .data = .{ it_0, it_1 } };
    }

    pub fn internalRanges(self: Connection) [2]Range {
        var new_range: [2]Range = undefined;
        @memcpy(new_range[0..], self.ranges[0..]);

        for (new_range[0..]) |*r| {
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

const ConnectionIterator = struct {
    data: [2]RangeIterator,

    pub fn next(self: *@This()) ?struct { usize, usize } {
        const n_0 = self.data[0].next() orelse {
            return null;
        };
        const n_1 = self.data[1].next() orelse {
            return null;
        };
        return .{ n_0, n_1 };
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

/// A flat array data structure for all mesh block boundary points. The block data is saved in column major ordering.
pub fn PointData(comptime T: type) type {
    return struct {
        allocator: std.mem.Allocator,
        buffer: []T,

        /// starting index of each block within buffer.
        block_range_start: []usize,

        /// returns a struct with an allocated (not yet initiallized) buffer for all boundary points
        /// of the mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const discrete.Mesh) !PointData(T) {
            var block_range_start = try allocator.alloc(usize, mesh.blocks.items.len);

            var total_boundary_points: usize = 0;
            for (mesh.blocks.items, 0..) |block, block_idx| {
                const size = block.points.size;
                block_range_start[block_idx] = total_boundary_points;
                const block_boundary_points = 2 * (size[1] + size[0] - 2);
                total_boundary_points += block_boundary_points;
            }

            return .{
                .allocator = allocator,
                .buffer = try allocator.alloc(T, total_boundary_points),
                .block_range_start = block_range_start,
            };
        }

        pub fn deinit(self: PointData(T)) void {
            self.allocator.free(self.buffer);
            self.allocator.free(self.block_range_start);
        }

        pub fn bufferIndex(self: PointData(T), boundary_point: types.MeshIndex2d, block_size: types.Index2d) !usize {

            // Here is an example for a block of size 5x7:
            //
            //                     i_max
            //             (0,6) (1,6) (2,6) (3,6) (4,6)
            //       (0,6) *[06] *[08] *[10] *[12] *[19] (4,6)
            //       (0,5) *[05]                   *[18] (4,5)
            //       (0,4) *[04]                   *[17] (4,4)
            // j_min (0,3) *[03]                   *[16] (4,3) j_max
            //       (0,2) *[02]                   *[15] (4,2)
            //       (0,1) *[01]                   *[14] (4,1)
            //       (0,0) *[00] *[07] *[09] *[11] *[13] (4,0)
            //             (0,0) (1,0) (2,0) (3,0) (4,0)
            //                     i_min

            const block_boundary_point_idx = blk: {
                if (boundary_point.point[0] == 0) {
                    // j_min
                    break :blk boundary_point.point[1];
                } else if (boundary_point.point[0] == block_size[0] - 1) {
                    // j_max
                    break :blk block_size[1] + 2 * (block_size[0] - 2) + boundary_point.point[1];
                } else if (boundary_point.point[1] == 0) {
                    // i_min
                    std.debug.assert(boundary_point.point[1] != 0);
                    break :blk block_size[1] + (boundary_point.point[0] - 1) * 2;
                } else if (boundary_point.point[1] == block_size[1] - 1) {
                    // i_max
                    std.debug.assert(boundary_point.point[1] != block_size[1] - 1);
                    break :blk block_size[1] - 1 + (boundary_point.point[0] - 1) * 2;
                } else {
                    return error.NotBoundaryIndex;
                }
            };

            return self.block_range_start[boundary_point.block] + block_boundary_point_idx;
        }
    };
}
