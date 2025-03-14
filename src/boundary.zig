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
    data: [2]Range,

    pub fn init(data: [2]Range) Connection {
        return .{ .data = data };
    }

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

    pub fn iterate(self: Connection, mesh: *const discrete.Mesh) ConnectionIterator {
        const it_0 = self.data[0].iterate(mesh);
        const it_1 = self.data[1].iterate(mesh);
        return .{ .data = .{ it_0, it_1 } };
    }

    pub fn internalRanges(self: Connection) [2]Range {
        var new_range: [2]Range = undefined;
        @memcpy(new_range[0..], self.data[0..]);

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
    data: RangeIterator[2],

    fn next(self: *@This()) ?struct { usize, usize } {
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

/// A flat array data structure for all mesh block boundary points.
/// The block data is saved in column major ordering.
pub fn PointData(comptime T: type) type {
    return struct {
        allocator: std.mem.Allocator,
        buffer: []T,

        // NOTE: a struct of arrays would be more efficient... There is something in the std to use.
        blocks: []BlockInfo,

        /// returns a struct with an allocated (not yet initiallized) buffer for all boundary points
        /// of the mesh.
        pub fn init(allocator: std.mem.Allocator, mesh: *const discrete.Mesh) PointData(T) {
            var data = PointData(T){
                .allocator = allocator,
                .buffer = undefined,
                .blocks = try allocator.alloc(BlockInfo, mesh.blocks.item.len),
            };

            var count: usize = 0;
            for (mesh.blocks.items, 0..) |block, block_idx| {
                const size = block.points.size;
                data.blocks[block_idx] = .{ .buffer_start_idx = count, .size = size };
                count += size[0] * size[1];
            }

            data.buffer = try allocator.alloc(T, count);

            return data;
        }

        pub fn deinit(self: PointData(T)) void {
            self.allocator.free(self.buffer);
            self.allocator.free(self.blocks);
        }

        pub fn get(self: PointData(T), global_point_idx: usize) T {
            return self.buffer[self.bufferIndex(global_point_idx)];
        }

        pub fn getPtr(self: *PointData(T), global_point_idx: usize) *T {
            return self.buffer[self.bufferIndex(global_point_idx)];
        }

        pub fn bufferIndex(self: PointData(T), global_point_idx: usize) usize {
            // TODO: move the transformation of a global id to a block and point id into mesh

            var block_idx: usize = self.blocks.len - 1;
            while (global_point_idx > self.blocks[block_idx].buffer_start_idx) : (block_idx -= 1) {}

            const block = self.blocks[block_idx];

            const local_point_idx = global_point_idx - block.buffer_start_idx;

            // compute 2d index
            const point_idx = blk: {
                const point_i = @divTrunc(local_point_idx, block.size[1]);
                const point_j = local_point_idx - point_i * block.size[1];
                break :blk types.Index2d{ point_i, point_j };
            };

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

            if (point_idx[1] == 0) {
                // j_min
                return point_idx[1];
            } else if (point_idx[1] == block.size[1] - 1) {
                // j_max
                return block.size[1] + 2 * (block.size[0] - 2) + point_idx[1];
            } else if (point_idx[0] == 0) {
                // i_min
                std.debug.assert(point_idx[1] != 0);
                return block.size[1] + (point_idx[0] - 1) * 2;
            } else if (point_idx[0] == block.size[0] - 1) {
                // i_max
                std.debug.assert(point_idx[1] != block.size[1] - 1);
                return block.size[1] - 1 + (point_idx[0] - 1) * 2;
            } else {
                return error.NotBoundaryIndex;
            }
        }

        /// internal type to collect the size and buffer start of each block.
        const BlockInfo = struct {
            buffer_start_idx: usize,
            size: types.Index2d,
        };
    };
}
