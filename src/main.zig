const std = @import("std");
const types = @import("types.zig");
const cgns = @import("cgns.zig");

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    var buffer: [21 * 17]Float = undefined;

    const block_names = [_][]const u8{
        "block_0",
    };

    var block_points = [_]Mat2d{
        try Mat2d.init(allocator, .{ 21, 17 }),
    };

    defer {
        for (block_points) |block| {
            block.deinit(allocator);
        }
    }

    for (block_points) |block| {
        const size = block.size;

        var idx: usize = 0;
        var i: usize = 0;
        var j: usize = 0;
        while (i < size[0]) : (i += 1) {
            while (j < size[1]) : (j += 1) {
                block.data[idx] = Vec2d{ @floatFromInt(i), @floatFromInt(j) };
                idx += 1;
            }
        }
    }

    const filename = "example.cgns";
    try cgns.write(filename, &block_names, &block_points, &buffer);
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
