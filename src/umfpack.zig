const std = @import("std");
const c = @cImport({
    @cInclude("umfpack.h");
});

test "umfpack simple" {
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
