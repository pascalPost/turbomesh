const types = @import("types.zig");

const Mat2d = types.Mat2d;

const Shell2d = struct {
    // consists of 4 2d edges
};

const Block2d = struct {
    points: Mat2d,
};

const Mesh = struct {
    // map of names to block ids
    // vec of blocks
    // vec of boundaries
};

test "create simple block" {}
