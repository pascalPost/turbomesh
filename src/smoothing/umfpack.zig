const std = @import("std");
const c = @cImport({
    @cInclude("umfpack.h");
});

// UMFPACK docs: https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK/Doc

pub fn solve(
    n_row: i32,
    n_col: i32,
    ap: []const i32,
    ai: []const i32,
    ax: []const f64,
    rhs_x: []const f64,
    x: []f64,
    rhs_y: []const f64,
    y: []f64,
) !void {
    var symbolic: ?*anyopaque = undefined;
    var numeric: ?*anyopaque = undefined;

    // var buffer: [1024]u8 = undefined;
    // try writeFile("Ap.txt", ap, &buffer);
    // try writeFile("Ai.txt", ai, &buffer);
    // try writeFile("Ax.txt", ax, &buffer);
    // try writeFile("Ax.txt", ax, &buffer);
    // try writeFile("rhs_x.txt", rhs_x, &buffer);
    // try writeFile("rhs_y.txt", rhs_y, &buffer);

    {
        const res = c.umfpack_di_symbolic(n_row, n_col, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, &symbolic, null, null);

        switch (res) {
            c.UMFPACK_OK => {},
            c.UMFPACK_ERROR_invalid_matrix => return error.InvalidMatrix,
            else => std.debug.panic("error return from umfpack: {}\n", .{res}),
        }
    }

    _ = c.umfpack_di_numeric(ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, symbolic, &numeric, null, null);
    c.umfpack_di_free_symbolic(&symbolic);
    _ = c.umfpack_di_solve(c.UMFPACK_Aat, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, x[0..].ptr, rhs_x[0..].ptr, numeric, null, null);
    _ = c.umfpack_di_solve(c.UMFPACK_Aat, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, y[0..].ptr, rhs_y[0..].ptr, numeric, null, null);
    c.umfpack_di_free_numeric(&numeric);
}

fn writeFile(name: []const u8, field: anytype, buffer: []u8) !void {
    var file = try std.fs.cwd().createFile(name, .{});
    defer file.close();

    var file_writer = file.writer(buffer);
    var writer = &file_writer.interface;

    for (field) |item| {
        try writer.print("{}\n", .{item});
    }

    try writer.flush();
}

test "umfpack" {
    const n: c_int = 5;
    const Ap = [_]c_int{ 0, 2, 5, 9, 10, 12 };
    const Ai = [_]c_int{ 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4 };
    const Ax = [_]f64{ 2.0, 3.0, 3.0, -1.0, 4.0, 4.0, -3.0, 1.0, 2.0, 2.0, 6.0, 1.0 };
    const b = [_]f64{ 8.0, 45.0, -3.0, 3.0, 19.0 };
    var x: [5]f64 = undefined;

    // std.debug.print("BLAS used: {s}\n", .{c.SuiteSparse_BLAS_library()});
    // std.debug.print("BLAS integer size: {} bytes\n", .{@sizeOf(c.SUITESPARSE_BLAS_INT)});

    var symbolic: ?*anyopaque = undefined;
    var numeric: ?*anyopaque = undefined;
    _ = c.umfpack_di_symbolic(n, n, Ap[0..].ptr, Ai[0..].ptr, Ax[0..].ptr, &symbolic, null, null);
    _ = c.umfpack_di_numeric(Ap[0..].ptr, Ai[0..].ptr, Ax[0..].ptr, symbolic, &numeric, null, null);
    c.umfpack_di_free_symbolic(&symbolic);
    _ = c.umfpack_di_solve(c.UMFPACK_A, Ap[0..].ptr, Ai[0..].ptr, Ax[0..].ptr, x[0..].ptr, b[0..].ptr, numeric, null, null);
    c.umfpack_di_free_numeric(&numeric);

    const tol = std.math.floatEps(f64) * 10;

    try std.testing.expectApproxEqAbs(1.0, x[0], tol);
    try std.testing.expectApproxEqAbs(2.0, x[1], tol);
    try std.testing.expectApproxEqAbs(3.0, x[2], tol);
    try std.testing.expectApproxEqAbs(4.0, x[3], tol);
    try std.testing.expectApproxEqAbs(5.0, x[4], tol);
}
