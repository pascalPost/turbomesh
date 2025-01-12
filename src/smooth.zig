const std = @import("std");
const discrete = @import("discrete.zig");
const types = @import("types.zig");
const umfpack = @import("umfpack.zig");

pub fn block(allocator: std.mem.Allocator, points: *types.Mat2d, iterations: usize) !void {
    // TODO can we use a better algo for the sparse banded matrix ?

    const dof = (points.size[0] - 2) * (points.size[1] - 2);

    // allocations
    // const entries_per_dof = 9;
    // const non_zero_entries_capacity = dof * entries_per_dof; // this is an overstimation

    const non_zero_entries =
        // corners
        4 * 4 +

        // edges
        ((points.size[0] - 4) + (points.size[1] - 4)) * 2 * 6 +

        // internal
        (points.size[0] - 4) * (points.size[1] - 4) * 9;

    std.debug.print("{}\n", .{non_zero_entries});
    std.debug.assert(non_zero_entries == 24388);

    var buffer_int = try allocator.alloc(c_int, (dof + 1) + non_zero_entries);
    defer allocator.free(buffer_int);

    // left hand side (lhs) in compressed column form
    var lhs_p = buffer_int[0 .. dof + 1]; // cum sum of non-zero entries in columns
    lhs_p[0] = 0; // first value must be zero
    var lhs_i = buffer_int[dof + 1 ..]; // row indicies with non-zero values
    std.debug.assert(lhs_i.len == non_zero_entries);

    var buffer_float = try allocator.alloc(f64, 4 * dof + non_zero_entries);
    defer allocator.free(buffer_float);

    // solution vectors
    var x_new = buffer_float[0..dof];
    var y_new = buffer_float[dof .. 2 * dof];

    // right hand side (rhs) containing the current internal coordinates
    var rhs_x = buffer_float[2 * dof .. 3 * dof];
    var rhs_y = buffer_float[3 * dof .. 4 * dof];

    var lhs_values = buffer_float[4 * dof ..];

    // define non-zero matrix entries
    {
        var fba_count = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_p[1..]));
        var row_count = try std.ArrayList(c_int).initCapacity(fba_count.allocator(), dof);

        var fba_entries = std.heap.FixedBufferAllocator.init(std.mem.sliceAsBytes(lhs_i[0..]));
        var row_entries = try std.ArrayList(c_int).initCapacity(fba_entries.allocator(), non_zero_entries);

        const point_size_j: c_int = @intCast(points.size[1]);
        const col_size = point_size_j - 2;

        // TODO even better: assemble directly in final compressed column format

        var matrix_idx: c_int = 0;

        errdefer {
            std.debug.print("block size: {} x {} = {}\n", .{ points.size[0], points.size[1], points.size[0] * points.size[1] });
            std.debug.print("all internal points => DOF: {}\n", .{dof});
            std.debug.print("matrix size: {} x {}\n", .{ dof, dof });
            std.debug.print("idx: {}\n", .{matrix_idx});
            std.debug.print("column: {}\n", .{row_count.items.len});
            std.debug.print("row count buffer size: {}\n", .{lhs_p.len});
            std.debug.print("row count capacity: {}\n", .{row_count.capacity});
            std.debug.print("non-zero entries capacity: {}, matrix buffer size: {}\n", .{ non_zero_entries, lhs_i.len });
            std.debug.print("{any}\n", .{row_count.items});
        }

        // corner A[0, 0] = (1, 1)
        {
            // point_idx = point_size_j + 1;
            try row_entries.append(matrix_idx); // A[i, j]
            try row_entries.append(matrix_idx + 1); // A[i, j+1]
            try row_entries.append(matrix_idx + col_size); // A[i+1, j]
            try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
            try row_count.append(@intCast(row_entries.items.len));
        }

        // edge A[0, j] = (1, j)
        for (2..points.size[1] - 2) |_| {
            matrix_idx += 1; // points.size[1] + j
            try row_entries.append(matrix_idx - 1); // A[i, j-1]
            try row_entries.append(matrix_idx); // A[i, j]
            try row_entries.append(matrix_idx + 1); // A[i, j+1]
            try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
            try row_entries.append(matrix_idx + col_size); // A[i+1, j]
            try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
            try row_count.append(@intCast(row_entries.items.len));
        }

        // corner A[0, m] = (1, J-1)
        {
            matrix_idx += 1; // points.size[1] + points.size[1] - 2
            try row_entries.append(matrix_idx - 1); // A[i, j-1]
            try row_entries.append(matrix_idx); // A[i, j]
            try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
            try row_entries.append(matrix_idx + col_size); // A[i+1, j]
            try row_count.append(@intCast(row_entries.items.len));
        }

        // loop over i
        for (2..points.size[0] - 3) |_| {

            // edge A[i, 0] = (i, 0)
            {
                matrix_idx += 1;
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // internal points
            for (2..points.size[1] - 2) |_| {
                // loop over j
                matrix_idx += 1;
                try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + 1); // A[i, j+1]
                try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
                try row_count.append(@intCast(row_entries.items.len));
            }

            // edge A[i, m] = (i, J-1)
            {
                matrix_idx += 1;
                try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
                try row_entries.append(matrix_idx - col_size); // A[i-1, j]
                try row_entries.append(matrix_idx - 1); // A[i, j-1]
                try row_entries.append(matrix_idx); // A[i, j]
                try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
                try row_entries.append(matrix_idx + col_size); // A[i+1, j]
                try row_count.append(@intCast(row_entries.items.len));
            }
        }

        // corner A[n, 0] = (I-1, 1)
        {
            matrix_idx += 1;
            try row_entries.append(matrix_idx - col_size); // A[i-1, j]
            try row_entries.append(matrix_idx - col_size + 1); // A[i-1, j+1]
            try row_entries.append(matrix_idx); // A[i, j]
            try row_entries.append(matrix_idx + 1); // A[i, j+1]
            try row_count.append(@intCast(row_entries.items.len));
        }

        // edge A[0, j] = (1, j)
        for (2..points.size[1] - 2) |_| {
            matrix_idx += 1; // points.size[1] + j
            try row_entries.append(matrix_idx - 1); // A[i, j-1]
            try row_entries.append(matrix_idx); // A[i, j]
            try row_entries.append(matrix_idx + 1); // A[i, j+1]
            try row_entries.append(matrix_idx + col_size - 1); // A[i+1, j-1]
            try row_entries.append(matrix_idx + col_size); // A[i+1, j]
            try row_entries.append(matrix_idx + col_size + 1); // A[i+1, j+1]
            try row_count.append(@intCast(row_entries.items.len));
        }

        // corner A[n, m] = (I-1, J-1)
        {
            matrix_idx += 1; // points.size[1] + points.size[1] - 2
            try row_entries.append(matrix_idx - col_size - 1); // A[i-1, j-1]
            try row_entries.append(matrix_idx - col_size); // A[i-1, j]
            try row_entries.append(matrix_idx - 1); // A[i, j-1]
            try row_entries.append(matrix_idx); // A[i, j]
            try row_count.append(@intCast(row_entries.items.len));
        }
    }

    // iterate and fill matrix values
    for (0..iterations) |n| {
        std.debug.print("  iteration: {}\n", .{n});

        // laplace conditions (zero intialization)
        const s = 0.0;
        const t = 0.0;

        {
            const point_size_j = points.size[1];
            var points_idx = point_size_j + 1;

            // corner A[0, 0] = (1, 1)
            {
                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[0] = values.a_i_j; // A[i, j]
                lhs_values[1] = values.a_i_jp1; // A[i, j+1]
                lhs_values[2] = values.a_ip1_j; // A[i+1, j]
                lhs_values[3] = values.a_ip1_jp1; // A[i+1, j+1]

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j-1] A[i+1, j-1]

                rhs_x[0] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0] +
                    values.a_i_jm1 * i_jm1.data[0] +
                    values.a_ip1_jm1 * ip1_jm1.data[0]);

                rhs_y[0] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1] +
                    values.a_i_jm1 * i_jm1.data[1] +
                    values.a_ip1_jm1 * ip1_jm1.data[1]);
            }

            var matrix_idx: usize = 4;
            var rhs_idx: usize = 1;

            // edge A[0, j] = (1, j)
            for (2..points.size[1] - 2) |_| {
                points_idx += 1; // points.size[1] + j

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_j = points.data[points_idx + points.size[1]];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 2] = values.a_i_jp1; // A[i, j+1]
                lhs_values[matrix_idx + 3] = values.a_ip1_jm1; // A[i+1, j-1]
                lhs_values[matrix_idx + 4] = values.a_ip1_j; // A[i+1, j]
                lhs_values[matrix_idx + 5] = values.a_ip1_jp1; // A[i+1, j+1]

                matrix_idx += 6;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1]);

                rhs_idx += 1;
            }

            // corner A[0, m] = (1, J-1)
            {
                points_idx += 1; // points.size[1] + points.size[1] - 2

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 1] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 2] = values.a_ip1_jm1; // A[i+1, j-1]
                lhs_values[matrix_idx + 3] = values.a_ip1_j; // A[i+1, j]

                matrix_idx += 4;

                // RHS: A[i-1, j-1], A[i-1, j], A[i-1, j+1], A[i, j+1] A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                    values.a_im1_j * im1_j.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0] +
                    values.a_i_jp1 * i_jp1.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                    values.a_im1_j * im1_j.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1] +
                    values.a_i_jp1 * i_jp1.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1]);

                rhs_idx += 1;
            }

            // loop over i
            for (2..points.size[0] - 3) |_| {

                // edge A[i, 0] = (i, 0)
                {
                    points_idx += 1;

                    const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                    const im1_j = points.data[points_idx - points.size[1]];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                    const ip1_j = points.data[points_idx + points.size[1]];

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 1] = values.a_im1_jp1; // A[i-1, j+1]
                    lhs_values[matrix_idx + 2] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 3] = values.a_i_jp1; // A[i, j+1]
                    lhs_values[matrix_idx + 4] = values.a_ip1_j; // A[i+1, j]
                    lhs_values[matrix_idx + 5] = values.a_ip1_jp1; // A[i+1, j+1]

                    matrix_idx += 6;

                    // RHS: A[i-1, j-1], A[i, j-1], A[i+1, j-1]

                    rhs_x[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[0] +
                        values.a_i_jm1 * i_jm1.data[0] +
                        values.a_ip1_jm1 * ip1_jm1.data[0]);

                    rhs_y[rhs_idx] = -(values.a_im1_jm1 * im1_jm1.data[1] +
                        values.a_i_jm1 * i_jm1.data[1] +
                        values.a_ip1_jm1 * ip1_jm1.data[1]);

                    rhs_idx += 1;
                }

                // internal points
                for (2..points.size[1] - 2) |_| {
                    // loop over j
                    points_idx += 1;

                    const im1_j = points.data[points_idx - points.size[1]];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_j = points.data[points_idx + points.size[1]];

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.a_im1_jp1; // A[i-1, j+1]
                    lhs_values[matrix_idx + 3] = values.a_i_jm1; // A[i, j-1]
                    lhs_values[matrix_idx + 4] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 5] = values.a_i_jp1; // A[i, j+1]
                    lhs_values[matrix_idx + 6] = values.a_ip1_jm1; // A[i+1, j-1]
                    lhs_values[matrix_idx + 7] = values.a_ip1_j; // A[i+1, j]
                    lhs_values[matrix_idx + 8] = values.a_ip1_jp1; // A[i+1, j+1]

                    matrix_idx += 9;

                    rhs_x[rhs_idx] = 0;
                    rhs_y[rhs_idx] = 0;

                    rhs_idx += 1;
                }

                // edge A[i, m] = (i, J-1)
                {
                    points_idx += 1;

                    const im1_j = points.data[points_idx - points.size[1]];
                    const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                    const i_jm1 = points.data[points_idx - 1];
                    const i_jp1 = points.data[points_idx + 1];
                    const ip1_j = points.data[points_idx + points.size[1]];
                    const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                    const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                    lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                    lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                    lhs_values[matrix_idx + 2] = values.a_i_jm1; // A[i, j-1]
                    lhs_values[matrix_idx + 3] = values.a_i_j; // A[i, j]
                    lhs_values[matrix_idx + 4] = values.a_ip1_jm1; // A[i+1, j-1]
                    lhs_values[matrix_idx + 5] = values.a_ip1_j; // A[i+1, j]

                    matrix_idx += 6;

                    // RHS: A[i-1, j+1], A[i, j+1], A[i+1, j+1]

                    rhs_x[rhs_idx] = -(values.a_im1_jp1 * im1_jp1.data[0] +
                        values.a_i_jp1 * i_jp1.data[0] +
                        values.a_ip1_jp1 * ip1_jp1.data[0]);

                    rhs_y[rhs_idx] = -(values.a_im1_jp1 * im1_jp1.data[1] +
                        values.a_i_jp1 * i_jp1.data[1] +
                        values.a_ip1_jp1 * ip1_jp1.data[1]);

                    rhs_idx += 1;
                }
            }

            // corner A[n, 0] = (I-1, 1)
            {
                points_idx += 1;

                const im1_jm1 = points.data[points_idx - points.size[1] - 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 1] = values.a_im1_jp1; // A[i-1, j+1]
                lhs_values[matrix_idx + 2] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 3] = values.a_i_jp1; // A[i, j+1]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j-1], A[i-1, j-1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0] +
                    values.a_i_jm1 * i_jm1.data[0] +
                    values.a_im1_jm1 * im1_jm1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1] +
                    values.a_i_jm1 * i_jm1.data[1] +
                    values.a_im1_jm1 * im1_jm1.data[1]);

                rhs_idx += 1;
            }

            // edge A[0, j] = (1, j)
            for (2..points.size[1] - 2) |_| {
                points_idx += 1; // points.size[1] + j

                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.a_im1_jp1; // A[i-1, j+1]
                lhs_values[matrix_idx + 3] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 4] = values.a_i_j; // A[i, j]
                lhs_values[matrix_idx + 5] = values.a_i_jp1; // A[i, j+1]

                matrix_idx += 6;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1]);

                rhs_idx += 1;
            }

            // corner A[n, m] = (I-1, J-1)
            {
                points_idx += 1; // points.size[1] + points.size[1] - 2

                const im1_jp1 = points.data[points_idx - points.size[1] + 1];
                const im1_j = points.data[points_idx - points.size[1]];
                const i_jm1 = points.data[points_idx - 1];
                const i_jp1 = points.data[points_idx + 1];
                const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
                const ip1_j = points.data[points_idx + points.size[1]];
                const ip1_jp1 = points.data[points_idx + points.size[1] + 1];

                const values = computeMatrixValues(im1_j, ip1_j, i_jm1, i_jp1, s, t);

                lhs_values[matrix_idx + 0] = values.a_im1_jm1; // A[i-1, j-1]
                lhs_values[matrix_idx + 1] = values.a_im1_j; // A[i-1, j]
                lhs_values[matrix_idx + 2] = values.a_i_jm1; // A[i, j-1]
                lhs_values[matrix_idx + 3] = values.a_i_j; // A[i, j]

                matrix_idx += 4;

                // RHS: A[i+1, j-1], A[i+1, j], A[i+1, j+1], A[i, j+1], A[i-1, j+1]

                rhs_x[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[0] +
                    values.a_ip1_j * ip1_j.data[0] +
                    values.a_ip1_jp1 * ip1_jp1.data[0] +
                    values.a_i_jp1 * i_jp1.data[0] +
                    values.a_im1_jp1 * im1_jp1.data[0]);

                rhs_y[rhs_idx] = -(values.a_ip1_jm1 * ip1_jm1.data[1] +
                    values.a_ip1_j * ip1_j.data[1] +
                    values.a_ip1_jp1 * ip1_jp1.data[1] +
                    values.a_i_jp1 * i_jp1.data[1] +
                    values.a_im1_jp1 * im1_jp1.data[1]);

                rhs_idx += 1;
            }
        }

        try umfpack.solve2(
            @intCast(dof),
            @intCast(dof),
            lhs_p,
            lhs_i,
            lhs_values,
            rhs_x,
            x_new[0..],
            rhs_y,
            y_new[0..],
        );

        //     // solve
        //     let mut solver = Solver::new(config)?;
        //     solver.initialize(&lhs)?;
        //     solver.factorize()?;
        //     solver.solve(&mut x_new, &rhs_x)?;
        //     solver.solve(&mut y_new, &rhs_y)?;

        //     let mut x_norm_sqr = 0.0;
        //     let mut y_norm_sqr = 0.0;

        //     for j in 1..shape.1 - 1 {
        //         for i in 1..shape.0 - 1 {
        //             let idx = matrix_index(i, j, i_size);
        //             let Vec2d(x, y) = coords[[i, j]];
        //             let dx = x - x_new[idx];
        //             let dy = y - y_new[idx];

        //             x_norm_sqr += dx * dx;
        //             y_norm_sqr += dy * dy;
        //         }
        //     }

        //     let norm = (x_norm_sqr + y_norm_sqr).sqrt();

        //     println!("\tresidual {norm}");

        //     // copy into coordinate field
        //     // TODO consider different memory layout for higher efficiency
        //     for j in 1..shape.1 - 1 {
        //         for i in 1..shape.0 - 1 {
        //             let idx = matrix_index(i, j, i_size);
        //             coords[[i, j]] = Vec2d(x_new[idx], y_new[idx]);
        //         }
        //     }
    }
}

fn computeMatrixValues(im1_j: types.Vec2d, ip1_j: types.Vec2d, i_jm1: types.Vec2d, i_jp1: types.Vec2d, s: f64, t: f64) struct {
    a_i_j: f64,
    a_ip1_j: f64,
    a_im1_j: f64,
    a_i_jp1: f64,
    a_i_jm1: f64,
    a_ip1_jp1: f64,
    a_ip1_jm1: f64,
    a_im1_jp1: f64,
    a_im1_jm1: f64,
} {
    const x_xi = 0.5 * (ip1_j.data[0] - im1_j.data[0]);
    const x_eta = 0.5 * (i_jp1.data[0] - i_jm1.data[0]);
    const y_xi = 0.5 * (ip1_j.data[1] - im1_j.data[1]);
    const y_eta = 0.5 * (i_jp1.data[1] - i_jm1.data[1]);

    const p = x_eta * x_eta + y_eta * y_eta;
    const q = x_xi * x_eta + y_xi * y_eta;
    const r = x_xi * x_xi + y_xi * y_xi;

    return .{
        .a_i_j = -2.0 * p - 2.0 * r,
        .a_ip1_j = p + 0.5 * s,
        .a_im1_j = p - 0.5 * s,
        .a_i_jp1 = r + 0.5 * t,
        .a_i_jm1 = r - 0.5 * t,
        .a_ip1_jp1 = -0.5 * q,
        .a_ip1_jm1 = 0.5 * q,
        .a_im1_jp1 = 0.5 * q,
        .a_im1_jm1 = -0.5 * q,
    };
}

// fn computeValues(points: types.Mat2d, points_idx: usize) struct {
//     a_i_j: f64,
//     a_ip1_j: f64,
//     a_im1_j: f64,
//     a_i_jp1: f64,
//     a_i_jm1: f64,
//     a_ip1_jp1: f64,
//     a_ip1_jm1: f64,
//     a_im1_jp1: f64,
//     a_im1_jm1: f64,
// } {
//     // laplace conditions (zero intialization)
//     const s = 0.0;
//     const t = 0.0;

//     const ip1_jm1 = points.data[points_idx + points.size[1] - 1];
//     const ip1_j = points.data[points_idx + points.size[1]];
//     const ip1_jp1 = points.data[points_idx + points.size[1] + 1];
//     const im1_jm1 = points.data[points_idx - points.size[1] - 1];
//     const im1_j = points.data[points_idx - points.size[1]];
//     const im1_jp1 = points.data[points_idx - points.size[1] + 1];
//     const i_jp1 = points.data[points_idx + 1];
//     const i_jm1 = points.data[points_idx - 1];

//     const x_xi = 0.5 * (ip1_j[0] - im1_j[0]);
//     const x_eta = 0.5 * (i_jp1[0] - i_jm1[0]);
//     const y_xi = 0.5 * (ip1_j[1] - im1_j[1]);
//     const y_eta = 0.5 * (i_jp1[1] - i_jm1[1]);

//     const p = x_eta * x_eta + y_eta * y_eta;
//     const q = x_xi * x_eta + y_xi * y_eta;
//     const r = x_xi * x_xi + y_xi * y_xi;

//     const a_i_j = -2.0 * p - 2.0 * r;
//     const a_ip1_j = p + 0.5 * s;
//     const a_im1_j = p - 0.5 * s;
//     const a_i_jp1 = r + 0.5 * t;
//     const a_i_jm1 = r - 0.5 * t;
//     const a_ip1_jp1 = -0.5 * q;
//     const a_ip1_jm1 = 0.5 * q;
//     const a_im1_jp1 = 0.5 * q;
//     const a_im1_jm1 = -0.5 * q;
// }
