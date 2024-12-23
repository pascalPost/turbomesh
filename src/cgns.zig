const std = @import("std");
const types = @import("types.zig");

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

const cgns_log = std.log.scoped(.cgns);

const cgns = @cImport({
    @cInclude("cgnslib.h");
});

fn getErrorMessage() [*:0]const u8 {
    const msg = cgns.cg_get_error();
    return msg;
}

pub fn write(filename: []const u8, block_names: []const []const u8, block_points: []const Mat2d, buffer: []Float) !void {
    const n_blocks = block_names.len;
    if (n_blocks != block_points.len) {
        cgns_log.err("inconsistnet input lengths (should both equal the number of blocks) (names: {}, coordinates: {})", .{ block_names.len, block_points.len });
        return error.cgnsInconsistentInput;
    }

    var file_handle: c_int = undefined;
    var ierr: c_int = undefined;

    ierr = cgns.cg_open(filename.ptr, cgns.CG_MODE_WRITE, &file_handle);
    if (ierr != 0) {
        cgns_log.err("error opening file (filename: {s}, file_handle: {}) : {s}", .{ filename, file_handle, getErrorMessage() });
        return error.cgnsOpen;
    }

    var base_handle: c_int = undefined;
    ierr = cgns.cg_base_write(file_handle, "Base", 2, 2, &base_handle);
    if (ierr != 0) {
        cgns_log.err("error writing base (file: {}, base: {}): {s}", .{ file_handle, base_handle, getErrorMessage() });
        return error.cgnsBaseWrite;
    }

    for (block_points, 0..) |block, i_block| {
        const name = block_names[i_block];

        const num_points = block.size;
        const num_cells = Index2d{ num_points[0] - 1, num_points[1] - 1 };

        const size = [6]cgns.cgsize_t{ @intCast(num_points[0]), @intCast(num_points[1]), @intCast(num_cells[0]), @intCast(num_cells[1]), 0, 0 };

        var zone_handle: c_int = undefined;
        ierr = cgns.cg_zone_write(file_handle, base_handle, name.ptr, &size, cgns.Structured, &zone_handle);
        if (ierr != 0) {
            cgns_log.err("error writing zone (file: {}, base: {}, zone: {}): {s}", .{ file_handle, base_handle, zone_handle, getErrorMessage() });
            return error.cgnsZoneWrite;
        }

        var coord_handle: c_int = undefined;

        // transform coordinates from array of structs to struct of arrays
        if (buffer.len < num_points[0] * num_points[1]) {
            cgns_log.err("buffer error (file: {}, base: {}, zone: {}) : buffer size ({}) too small for block ({} x {} = {})", .{ file_handle, base_handle, zone_handle, buffer.len, num_points[0], num_points[1], num_points[0] * num_points[1] });
            return error.cgnsBufferInsufficient;
        }

        // for now we use a single buffer; perhaps it is faster to have a buffer twice the size and write x and y coordinates in one go
        // this might be a trade-off: mem vs speed
        {
            var idx: usize = 0;
            var j: Index = 0;
            while (j < block.size[1]) : (j += 1) {
                var i: Index = 0;
                while (i < block.size[0]) : (i += 1) {
                    buffer[idx] = block.data[block.index(.{ i, j })].data[0];
                    idx += 1;
                }
            }
        }

        ierr = cgns.cg_coord_write(file_handle, base_handle, zone_handle, cgns.RealSingle, "CoordinateX", buffer.ptr, &coord_handle);
        if (ierr != 0) {
            cgns_log.err("error writing coord (file: {}, base: {}, zone: {}, coord: {}): {s}", .{ file_handle, base_handle, zone_handle, coord_handle, getErrorMessage() });
            return error.cgnsCoordWrite;
        }

        {
            var idx: usize = 0;
            var j: Index = 0;
            while (j < block.size[1]) : (j += 1) {
                var i: Index = 0;
                while (i < block.size[0]) : (i += 1) {
                    buffer[idx] = block.data[block.index(.{ i, j })].data[1];
                    idx += 1;
                }
            }
        }

        ierr = cgns.cg_coord_write(file_handle, base_handle, zone_handle, cgns.RealSingle, "CoordinateY", buffer.ptr, &coord_handle);
        if (ierr != 0) {
            cgns_log.err("error writing coord (file: {}, base: {}, zone: {}, coord: {}): {s}", .{ file_handle, base_handle, zone_handle, coord_handle, getErrorMessage() });
            return error.cgnsCoordWrite;
        }
    }

    defer {
        ierr = cgns.cg_close(file_handle);
        if (ierr != 0) {
            cgns_log.err("error closing file: {s}", .{getErrorMessage()});
        }
    }
}

test "write a mesh to a cgns file" {
    const allocator = std.testing.allocator;

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
        while (i < size[0]) : (i += 1) {
            var j: usize = 0;
            while (j < size[1]) : (j += 1) {
                block.data[idx] = Vec2d{ .data = .{ @floatFromInt(i), @floatFromInt(j) } };
                idx += 1;
            }
        }
    }

    const filename = "example.cgns";
    try write(filename, &block_names, &block_points, &buffer);
}
