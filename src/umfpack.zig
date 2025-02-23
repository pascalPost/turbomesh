const std = @import("std");
const c = @cImport({
    @cInclude("umfpack.h");
});

// UMFPACK docs: https://github.com/PetterS/SuiteSparse/tree/master/UMFPACK/Doc

pub fn solve2(
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

    // {
    //     var file = try std.fs.cwd().createFile("Ap.txt", .{});
    //     defer file.close();

    //     const writer = file.writer();

    //     for (ap) |a| {
    //         try writer.print("{}\n", .{a});
    //     }
    // }

    // {
    //     var file = try std.fs.cwd().createFile("Ai.txt", .{});
    //     defer file.close();

    //     const writer = file.writer();

    //     for (ai) |a| {
    //         try writer.print("{}\n", .{a});
    //     }
    // }

    // {
    //     var file = try std.fs.cwd().createFile("Ax.txt", .{});
    //     defer file.close();

    //     const writer = file.writer();

    //     for (ax) |a| {
    //         try writer.print("{}\n", .{a});
    //     }
    // }

    // {
    //     var file = try std.fs.cwd().createFile("rhs_x.txt", .{});
    //     defer file.close();

    //     const writer = file.writer();

    //     for (rhs_x) |a| {
    //         try writer.print("{}\n", .{a});
    //     }
    // }

    // {
    //     var file = try std.fs.cwd().createFile("rhs_y.txt", .{});
    //     defer file.close();

    //     const writer = file.writer();

    //     for (rhs_y) |a| {
    //         try writer.print("{}\n", .{a});
    //     }
    // }

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
