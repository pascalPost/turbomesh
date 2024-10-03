const std = @import("std");

const cgns_log = std.log.scoped(.cgns);

const cgns = @cImport({
    @cInclude("cgnslib.h");
});

fn get_error_message() [*:0]const u8 {
    const msg = cgns.cg_get_error();
    return msg;
}

fn write(filename: [:0]const u8) !void {
    var file_handle: c_int = undefined;
    var ierr: c_int = undefined;

    ierr = cgns.cg_open(filename, cgns.CG_MODE_WRITE, &file_handle);
    if (ierr != 0) {
        cgns_log.err("error opening CGNS file: {s}", .{get_error_message()});
        return error.cgnsOpen;
    }

    cgns_log.debug("file opend successfully", .{});

    // ierr = cgns.cg_base_write(file_handle, char *basename, int cell_dim,
    //     int phys_dim, int *B);

    defer {
        ierr = cgns.cg_close(file_handle);
        if (ierr != 0) {
            cgns_log.err("error closing CGNS file: {s}", .{get_error_message()});
        }
        cgns_log.debug("file closed successfully", .{});
    }
}

test "write a mesh to a cgns file" {
    const filename = "example.cgns";
    try write(filename);
}
