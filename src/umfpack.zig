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

    {
        var file = try std.fs.cwd().createFile("Ap.txt", .{});
        defer file.close();

        const writer = file.writer();

        for (ap) |a| {
            try writer.print("{}\n", .{a});
        }
    }

    {
        var file = try std.fs.cwd().createFile("Ai.txt", .{});
        defer file.close();

        const writer = file.writer();

        for (ai) |a| {
            try writer.print("{}\n", .{a});
        }
    }

    {
        var file = try std.fs.cwd().createFile("Ax.txt", .{});
        defer file.close();

        const writer = file.writer();

        for (ax) |a| {
            try writer.print("{}\n", .{a});
        }
    }

    {
        const res = c.umfpack_di_symbolic(n_row, n_col, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, &symbolic, null, null);

        switch (res) {
            c.UMFPACK_ERROR_invalid_matrix => return error.InvalidMatrix,
            else => std.debug.panic("error return from umfpack: {}\n", .{res}),
        }
    }

    _ = c.umfpack_di_numeric(ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, symbolic, &numeric, null, null);
    c.umfpack_di_free_symbolic(&symbolic);
    _ = c.umfpack_di_solve(c.UMFPACK_A, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, x[0..].ptr, rhs_x[0..].ptr, numeric, null, null);
    _ = c.umfpack_di_solve(c.UMFPACK_A, ap[0..].ptr, ai[0..].ptr, ax[0..].ptr, y[0..].ptr, rhs_y[0..].ptr, numeric, null, null);
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

// test "umfpack interface" {
//     var system = try SparseSystemSquare.init(std.testing.allocator, 5);
// }

// const SparseSystemSquare = struct {
//     // NOTE can be split in sparse matrix and rest

//     dof: c_int,

//     // matrix in compressed column form
//     ap: []c_int,
//     ai: []c_int,
//     ax: []f64,

//     b: []f64, // rhs
//     x: []f64, // solution

//     fn solve(self: *SparseSystemSquare) void {
//         var symbolic: ?*anyopaque = undefined;
//         var numeric: ?*anyopaque = undefined;
//         _ = c.umfpack_di_symbolic(self.dof, self.dof, self.ap[0..].ptr, self.ai[0..].ptr, self.ax[0..].ptr, &symbolic, null, null);
//         _ = c.umfpack_di_numeric(self.ap[0..].ptr, self.ai[0..].ptr, self.ax[0..].ptr, symbolic, &numeric, null, null);
//         c.umfpack_di_free_symbolic(&symbolic);
//         _ = c.umfpack_di_solve(c.UMFPACK_A, self.ap[0..].ptr, self.ai[0..].ptr, self.ax[0..].ptr, self.x[0..].ptr, self.b[0..].ptr, numeric, null, null);
//         c.umfpack_di_free_numeric(&numeric);
//     }

//     fn init(allocator: std.mem.Allocator, dof: usize, non_zero_entries: usize) !SparseSystemSquare {
//         // allocations
//         var buffer_int = try allocator.alloc(c_int, (dof + 1) + non_zero_entries);
//         defer allocator.free(buffer_int);

//         var buffer_float = try allocator.alloc(f64, 4 * dof + non_zero_entries);
//         defer allocator.free(buffer_float);

//         // // solution vectors
//         // var x_new = buffer_float[0..dof];
//         // var y_new = buffer_float[dof .. 2 * dof];

//         // // right hand side (rhs) containing the current internal coordinates
//         // var rhs_x = buffer_float[2 * dof .. 3 * dof];
//         // var rhs_y = buffer_float[3 * dof .. 4 * dof];

//         // left hand side (lhs) in compressed column form
//         // var lhs_values = buffer_float[4 * dof ..]; // colummn-vise
//         var lhs_p = buffer_int[0 .. dof + 1]; // cum sum of non-zero entries in columns
//         lhs_p[0] = 0; // first value must be zero
//         var lhs_i = buffer_int[dof + 1 ..]; // row indicies with non-zero values
//         std.debug.assert(lhs_i.len == non_zero_entries_capacity);

//         // define non-zero matrix entries
//         {
//             var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_p[1..]));
//             var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

//             var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_i[0..]));
//             var row_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries_capacity);

//             const point_size_j: c_int = @intCast(points.size[1]);
//             const col_size = point_size_j - 2;

//             // TODO even better: assemble directly in final compressed column format

//             var idx: c_int = undefined;

//             errdefer {
//                 std.debug.print("block size: {} x {} = {}\n", .{ points.size[0], points.size[1], points.size[0] * points.size[1] });
//                 std.debug.print("all internal points => DOF: {}\n", .{dof});
//                 std.debug.print("matrix size: {} x {}\n", .{ dof, dof });
//                 std.debug.print("idx: {}\n", .{idx});
//                 std.debug.print("column: {}\n", .{row_count.items.len});
//                 std.debug.print("row count buffer size: {}\n", .{lhs_p.len});
//                 std.debug.print("row count capacity: {}\n", .{row_count.capacity});
//                 std.debug.print("non-zero entries capacity: {}, matrix buffer size: {}\n", .{ non_zero_entries_capacity, lhs_i.len });
//                 std.debug.print("{any}\n", .{row_count.items});
//             }

//             // corner A[0, 0] = (1, 1)
//             {
//                 idx = point_size_j + 1;
//                 try row_entries.append(idx); // A[i, j]
//                 try row_entries.append(idx + 1); // A[i, j+1]
//                 try row_entries.append(idx + col_size); // A[i+1, j]
//                 try row_entries.append(idx + col_size + 1); // A[i+1, j+1]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }

//             // edge A[0, j] = (1, j)
//             for (2..points.size[1] - 2) |_| {
//                 idx += 1; // points.size[1] + j
//                 try row_entries.append(idx - 1); // a[i, j-1]
//                 try row_entries.append(idx); // A[i, j]
//                 try row_entries.append(idx + 1); // A[i, j+1]
//                 try row_entries.append(idx + col_size - 1); // A[i+1, j-1]
//                 try row_entries.append(idx + col_size); // A[i+1, j]
//                 try row_entries.append(idx + col_size + 1); // A[i+1, j+1]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }

//             // corner A[0, m] = (1, J-1)
//             {
//                 idx += 1; // points.size[1] + points.size[1] - 2
//                 try row_entries.append(idx - 1); // a[i, j-1]
//                 try row_entries.append(idx); // A[i, j]
//                 try row_entries.append(idx + col_size - 1); // A[i+1, j-1]
//                 try row_entries.append(idx + col_size); // A[i+1, j]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }

//             // loop over i
//             for (2..points.size[0] - 3) |_| {

//                 // edge A[i, 0] = (i, 0)
//                 {
//                     idx += 1;
//                     try row_entries.append(idx - col_size); // A[i-1, j]
//                     try row_entries.append(idx - col_size + 1); // A[i-1, j+1]
//                     try row_entries.append(idx); // A[i, j]
//                     try row_entries.append(idx + 1); // A[i, j+1]
//                     try row_entries.append(idx + col_size); // A[i+1, j]
//                     try row_entries.append(idx + col_size + 1); // A[i+1, j+1]
//                     try row_count.append(@intCast(row_entries.items.len));
//                 }

//                 // internal points
//                 for (2..points.size[1] - 2) |_| {
//                     // loop over j
//                     idx += 1;
//                     try row_entries.append(idx - col_size - 1); // A[i-1, j-1]
//                     try row_entries.append(idx - col_size); // A[i-1, j]
//                     try row_entries.append(idx - col_size + 1); // A[i-1, j+1]
//                     try row_entries.append(idx - 1); // a[i, j-1]
//                     try row_entries.append(idx); // A[i, j]
//                     try row_entries.append(idx + 1); // A[i, j+1]
//                     try row_entries.append(idx + col_size - 1); // A[i+1, j-1]
//                     try row_entries.append(idx + col_size); // A[i+1, j]
//                     try row_entries.append(idx + col_size + 1); // A[i+1, j+1]
//                     try row_count.append(@intCast(row_entries.items.len));
//                 }

//                 // edge A[i, m] = (i, J-1)
//                 {
//                     idx += 1;
//                     try row_entries.append(idx - col_size - 1); // A[i-1, j-1]
//                     try row_entries.append(idx - col_size); // A[i-1, j]
//                     try row_entries.append(idx - 1); // a[i, j-1]
//                     try row_entries.append(idx); // A[i, j]
//                     try row_entries.append(idx + col_size - 1); // A[i+1, j-1]
//                     try row_entries.append(idx + col_size); // A[i+1, j]
//                     try row_count.append(@intCast(row_entries.items.len));
//                 }
//             }

//             // corner A[n, 0] = (I-1, 1)
//             {
//                 idx += 1;
//                 try row_entries.append(idx - col_size); // A[i-1, j]
//                 try row_entries.append(idx - col_size + 1); // A[i-1, j+1]
//                 try row_entries.append(idx); // A[i, j]
//                 try row_entries.append(idx + 1); // A[i, j+1]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }

//             // edge A[0, j] = (1, j)
//             for (2..points.size[1] - 2) |_| {
//                 idx += 1; // points.size[1] + j
//                 try row_entries.append(idx - 1); // a[i, j-1]
//                 try row_entries.append(idx); // A[i, j]
//                 try row_entries.append(idx + 1); // A[i, j+1]
//                 try row_entries.append(idx + col_size - 1); // A[i+1, j-1]
//                 try row_entries.append(idx + col_size); // A[i+1, j]
//                 try row_entries.append(idx + col_size + 1); // A[i+1, j+1]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }

//             // corner A[n, m] = (I-1, J-1)
//             {
//                 idx += 1; // points.size[1] + points.size[1] - 2
//                 try row_entries.append(idx - col_size - 1); // A[i-1, j-1]
//                 try row_entries.append(idx - col_size); // A[i-1, j]
//                 try row_entries.append(idx - 1); // a[i, j-1]
//                 try row_entries.append(idx); // A[i, j]
//                 try row_count.append(@intCast(row_entries.items.len));
//             }
//         }
//     }
// };
