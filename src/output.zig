const std = @import("std");
const types = @import("types.zig");

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;

const cgns_log = std.log.scoped(.cgns);

const cgns = @cImport({
    @cInclude("cgnslib.h");
});

fn get_error_message() [*:0]const u8 {
    const msg = cgns.cg_get_error();
    return msg;
}

fn write(filename: [:0]const u8, block_names: []const []const u8, block_points: []const Mat2d) !void {
    const n_blocks = block_names.len;
    if (n_blocks != block_points.len) {
        cgns_log.err("inconsistnet input lengths (should both equal the number of blocks) (names: {}, coordinates: {})", .{ block_names.len, block_points.len });
        return error.cgnsInconsistentInput;
    }

    var file_handle: c_int = undefined;
    var ierr: c_int = undefined;

    ierr = cgns.cg_open(filename, cgns.CG_MODE_WRITE, &file_handle);
    if (ierr != 0) {
        cgns_log.err("error opening file (filename: {s}, file_handle: {}) : {s}", .{ filename, file_handle, get_error_message() });
        return error.cgnsOpen;
    }

    var base_handle: c_int = undefined;
    ierr = cgns.cg_base_write(file_handle, "Base", 2, 2, &base_handle);
    if (ierr != 0) {
        cgns_log.err("error writing base: {s}", .{get_error_message()});
        return error.cgnsBaseWrite;
    }

    var i_block = 0;
    while (i_block < n_blocks) : (i_block += 1) {
        const block = block_points[i_block];

        const num_points = block.size;
        const num_cells = Index2d{num_points[0]-1, num_points[1]-1};

        const size = [6]cgns.cgsize_t{num_points[0], num_points[1],num_cells[0], num_cells[1],0,0};

        var zone_handle: c_int = undefined;
        ierr = cg_zone_write(file_handle, base_handle, "block_0", cgsize_t *size,
        cgns.Structured, &zone_handle);
        if(ierr != 0){
        cgns_log.err("error writing zone: {s}", .{get_error_message()});
        return error.cgnsBaseWrite;
        }
        //
        // var coord_array_handle : c_int = undefined;
        // ierr = cg_coord_write(file_handle, base_handle, zone_handle, cgns.RealSingle,
        // "CoordinateX", void *coord_array, &coord_array_handle);
        // if(ierr != 0){
        // cgns_log.err("error writing coord: {s}", .{get_error_message()});
        // return error.cgnsBaseWrite;
        // }
    }

    defer {
        ierr = cgns.cg_close(file_handle);
        if (ierr != 0) {
            cgns_log.err("error closing file: {s}", .{get_error_message()});
        }
    }
}

test "write a mesh to a cgns file" {
const allocator = std.testing.allocator;

    // buffer

    var block_points = Mat2d.init(allocator, .{21,17});

const size = block_points.size;

    // TODO try with row-major index
    var idx = 0;
    var i = 0;
    var j = 0;
    while (i < block_points.size[0]) : (i += 1) {
        while (j < block_points.size[1]) : (j += 1) {
            block_coordinates[idx] = Vec2d{ i, j };
            idx += 1;
        }
    }

    const filename = "example.cgns";
    try write(filename);
}
