// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const RowCompressedMatrixSystem2d = @import("smooth.zig").RowCompressedMatrixSystem2d;
const bicgstab = @import("BiCGStab.zig");

const log = std.log.scoped(.gmres_solver);

const Preconditioner = @import("preconditioner.zig").Preconditioner;

pub const GMRESSolver = struct {
    system: RowCompressedMatrixSystem2d,
    buffer: ?[]f64 = null,
    dof: usize,
    seeded_initial_guess: bool = false,
    preconditioner: Preconditioner = .diagonal,
    ilu_values: ?[]f64 = null,
    ilu_diag_pos: ?[]c_int = null,
    ilu_marker: ?[]c_int = null,
    restart: usize = 30,
    max_iters: usize = 1000,
    rtol: f64 = 1e-6,
    atol: f64 = 1e-8,

    const Work = struct {
        restart: usize,
        diag_inv: []f64,
        r: []f64,
        w: []f64,
        z: []f64,
        v: []f64,
        h: []f64,
        cs: []f64,
        sn: []f64,
        g: []f64,
    };

    pub fn init(system: RowCompressedMatrixSystem2d, preconditioner: Preconditioner) GMRESSolver {
        return .{
            .system = system,
            .dof = system.rhs_x.len,
            .preconditioner = preconditioner,
        };
    }

    pub fn deinit(self: *GMRESSolver) void {
        if (self.buffer) |buf| {
            self.system.allocator.free(buf);
        }
        self.buffer = null;

        if (self.ilu_values) |buf| {
            self.system.allocator.free(buf);
        }
        self.ilu_values = null;

        if (self.ilu_diag_pos) |buf| {
            self.system.allocator.free(buf);
        }
        self.ilu_diag_pos = null;

        if (self.ilu_marker) |buf| {
            self.system.allocator.free(buf);
        }
        self.ilu_marker = null;
    }

    fn precondition(self: *GMRESSolver, work: Work) !void {
        switch (self.preconditioner) {
            .diagonal => self.updateDiagonalInverse(work.diag_inv),
            .ilu0 => try self.updateIlu0(),
        }
    }

    pub fn solve(self: *GMRESSolver) !void {
        if (!self.seeded_initial_guess) {
            self.seedInitialGuess();
        }

        const work = try self.ensureWorkspace();

        try self.system.fillXSpecific();
        try self.precondition(work);
        self.solveComponent(work, self.system.rhs_x, self.system.x_new, "x");

        try self.system.fillYSpecific();
        try self.precondition(work);
        self.solveComponent(work, self.system.rhs_y, self.system.y_new, "y");
    }

    fn ensureWorkspace(self: *GMRESSolver) !Work {
        const restart = @min(self.restart, self.dof);

        const v_len = (restart + 1) * self.dof;
        const h_len = (restart + 1) * restart;
        const cs_len = restart;
        const sn_len = restart;
        const g_len = restart + 1;
        const r_len = self.dof;
        const w_len = self.dof;
        const z_len = self.dof;
        const diag_len = self.dof;

        const total = v_len + h_len + cs_len + sn_len + g_len + r_len + w_len + z_len + diag_len;

        if (self.buffer == null or self.buffer.?.len != total) {
            if (self.buffer) |buf| {
                self.system.allocator.free(buf);
            }
            self.buffer = try self.system.allocator.alloc(f64, total);
        }

        const buf = self.buffer.?;
        var offset: usize = 0;

        const v = buf[offset .. offset + v_len];
        offset += v_len;

        const h = buf[offset .. offset + h_len];
        offset += h_len;

        const cs = buf[offset .. offset + cs_len];
        offset += cs_len;

        const sn = buf[offset .. offset + sn_len];
        offset += sn_len;

        const g = buf[offset .. offset + g_len];
        offset += g_len;

        const r = buf[offset .. offset + r_len];
        offset += r_len;

        const w = buf[offset .. offset + w_len];
        offset += w_len;

        const z = buf[offset .. offset + z_len];
        offset += z_len;

        const diag_inv = buf[offset .. offset + diag_len];

        return .{
            .restart = restart,
            .diag_inv = diag_inv,
            .r = r,
            .w = w,
            .z = z,
            .v = v,
            .h = h,
            .cs = cs,
            .sn = sn,
            .g = g,
        };
    }

    fn seedInitialGuess(self: *GMRESSolver) void {
        var row_block_start_id: usize = 0;
        for (self.system.mesh.blocks.items) |block| {
            var point_id: usize = 0;
            for (0..block.points.size[0]) |_| {
                for (0..block.points.size[1]) |_| {
                    const row_idx = row_block_start_id + point_id;
                    const p = block.points.data[point_id];
                    self.system.x_new[row_idx] = p.data[0];
                    self.system.y_new[row_idx] = p.data[1];
                    point_id += 1;
                }
            }
            row_block_start_id += block.points.size[0] * block.points.size[1];
        }

        self.seeded_initial_guess = true;
    }

    fn updateDiagonalInverse(self: *GMRESSolver, diag_inv: []f64) void {
        for (0..self.dof) |row| {
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);
            const row_id: c_int = @intCast(row);

            var diag: f64 = 0.0;
            for (start..end) |idx| {
                if (self.system.lhs_i[idx] == row_id) {
                    diag = self.system.lhs_values[idx];
                    break;
                }
            }

            if (diag == 0.0) {
                diag_inv[row] = 1.0;
            } else {
                diag_inv[row] = 1.0 / diag;
            }
        }
    }

    // Incomplete LU(0) factorization using the existing CSR sparsity pattern.
    fn updateIlu0(self: *GMRESSolver) !void {
        const nnz = self.system.lhs_values.len;

        if (self.ilu_values == null or self.ilu_values.?.len != nnz) {
            if (self.ilu_values) |buf| {
                self.system.allocator.free(buf);
            }
            self.ilu_values = try self.system.allocator.alloc(f64, nnz);
        }

        if (self.ilu_diag_pos == null or self.ilu_diag_pos.?.len != self.dof) {
            if (self.ilu_diag_pos) |buf| {
                self.system.allocator.free(buf);
            }
            self.ilu_diag_pos = try self.system.allocator.alloc(c_int, self.dof);
        }

        if (self.ilu_marker == null or self.ilu_marker.?.len != self.dof) {
            if (self.ilu_marker) |buf| {
                self.system.allocator.free(buf);
            }
            self.ilu_marker = try self.system.allocator.alloc(c_int, self.dof);
        }

        const lu = self.ilu_values.?;
        const diag_pos = self.ilu_diag_pos.?;
        const marker = self.ilu_marker.?;
        const invalid_pos: c_int = -1;

        @memcpy(lu, self.system.lhs_values);
        @memset(diag_pos, invalid_pos);
        @memset(marker, invalid_pos);

        for (0..self.dof) |row| {
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);
            const row_id: c_int = @intCast(row);

            for (start..end) |idx| {
                if (self.system.lhs_i[idx] == row_id) {
                    diag_pos[row] = @intCast(idx);
                    break;
                }
            }
        }

        var missing_diag_count: usize = 0;
        var zero_diag_count: usize = 0;

        for (0..self.dof) |row| {
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);

            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                marker[col] = @intCast(idx);
            }

            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                if (col >= row) continue;

                const diag_idx = diag_pos[col];
                var diag: f64 = 1.0;
                if (diag_idx >= 0) {
                    diag = lu[@intCast(diag_idx)];
                    if (diag == 0.0) {
                        diag = 1.0;
                        zero_diag_count += 1;
                    }
                } else {
                    missing_diag_count += 1;
                }

                const lij = lu[idx] / diag;
                lu[idx] = lij;

                const row_start: usize = @intCast(self.system.lhs_p[col]);
                const row_end: usize = @intCast(self.system.lhs_p[col + 1]);
                for (row_start..row_end) |row_idx| {
                    const col_k: usize = @intCast(self.system.lhs_i[row_idx]);
                    if (col_k <= col) continue;

                    const pos = marker[col_k];
                    if (pos >= 0) {
                        lu[@intCast(pos)] -= lij * lu[row_idx];
                    }
                }
            }

            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                marker[col] = invalid_pos;
            }
        }

        if (missing_diag_count > 0 or zero_diag_count > 0) {
            log.warn("ILU(0) preconditioner: {} missing diagonals, {} zero diagonals treated as 1.0", .{ missing_diag_count, zero_diag_count });
        }
    }

    fn solveSystem(self: *GMRESSolver, rhs: []const f64, x: []f64, work: Work, label: []const u8) void {
        const breakdown_eps = 1e-30;
        const restart = work.restart;
        if (restart == 0) return;

        const norm_b = norm(rhs);
        const tol = @max(self.atol, self.rtol * norm_b);

        var iter_total: usize = 0;
        var last_resid: f64 = 0.0;

        while (iter_total < self.max_iters) {
            self.matVec(x, work.w);
            for (rhs, work.w, 0..) |rhs_val, ax_val, i| {
                work.r[i] = rhs_val - ax_val;
            }

            self.applyPreconditioner(work.r, work.z, work);
            const beta = norm(work.z);
            if (beta <= tol) return;

            const v0 = self.vSlice(work, 0);
            for (0..self.dof) |i| {
                v0[i] = work.z[i] / beta;
            }

            @memset(work.h, 0);
            @memset(work.cs, 0);
            @memset(work.sn, 0);
            @memset(work.g, 0);
            work.g[0] = beta;

            var cols_used: usize = 0;
            var converged = false;
            var resid = beta;

            var j: usize = 0;
            while (j < restart and iter_total < self.max_iters) : (j += 1) {
                const vj = self.vSlice(work, j);
                self.matVec(vj, work.w);
                self.applyPreconditioner(work.w, work.z, work);

                for (0..j + 1) |i| {
                    const vi = self.vSlice(work, i);
                    const h_ij = dot(work.z, vi);
                    self.hSet(work, i, j, h_ij);
                    for (0..self.dof) |k| {
                        work.z[k] -= h_ij * vi[k];
                    }
                }

                const h_next = norm(work.z);
                self.hSet(work, j + 1, j, h_next);

                if (h_next > breakdown_eps) {
                    const vnext = self.vSlice(work, j + 1);
                    for (0..self.dof) |k| {
                        vnext[k] = work.z[k] / h_next;
                    }
                }

                for (0..j) |i| {
                    const h_i = self.hGet(work, i, j);
                    const h_ip1 = self.hGet(work, i + 1, j);
                    const temp = work.cs[i] * h_i + work.sn[i] * h_ip1;
                    self.hSet(work, i + 1, j, -work.sn[i] * h_i + work.cs[i] * h_ip1);
                    self.hSet(work, i, j, temp);
                }

                const h_j = self.hGet(work, j, j);
                const h_jp1 = self.hGet(work, j + 1, j);
                const rot = computeGivens(h_j, h_jp1);
                work.cs[j] = rot.c;
                work.sn[j] = rot.s;

                self.hSet(work, j, j, rot.r);
                self.hSet(work, j + 1, j, 0.0);

                const g_j = work.g[j];
                const g_jp1 = work.g[j + 1];
                work.g[j] = rot.c * g_j + rot.s * g_jp1;
                work.g[j + 1] = -rot.s * g_j + rot.c * g_jp1;

                resid = @abs(work.g[j + 1]);
                iter_total += 1;
                cols_used = j + 1;

                last_resid = resid;
                if (resid <= tol) {
                    converged = true;
                    break;
                }
            }

            if (cols_used == 0) break;

            var y = work.w[0..cols_used];
            var idx: usize = cols_used;
            while (idx > 0) {
                idx -= 1;
                var sum = work.g[idx];
                for (idx + 1..cols_used) |k| {
                    sum -= self.hGet(work, idx, k) * y[k];
                }

                const h_ii = self.hGet(work, idx, idx);
                if (h_ii == 0.0) break;
                y[idx] = sum / h_ii;
            }

            for (0..cols_used) |i| {
                const vi = self.vSlice(work, i);
                const yi = y[i];
                for (0..self.dof) |k| {
                    x[k] += yi * vi[k];
                }
            }

            if (converged) return;
            if (resid <= tol) return;
        }

        log.warn("gmres solve did not converge for {s}: iter={}, residual={e:.3}", .{ label, iter_total, last_resid });
    }

    fn applyPreconditioner(self: *GMRESSolver, rhs: []const f64, out: []f64, work: Work) void {
        switch (self.preconditioner) {
            .diagonal => {
                for (0..self.dof) |i| {
                    out[i] = rhs[i] * work.diag_inv[i];
                }
            },
            .ilu0 => self.applyIlu0(rhs, out),
        }
    }

    // Apply ILU(0) with forward/backward substitution on the stored LU factors.
    fn applyIlu0(self: *GMRESSolver, rhs: []const f64, out: []f64) void {
        const lu = self.ilu_values.?;
        const diag_pos = self.ilu_diag_pos.?;

        for (0..self.dof) |row| {
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);
            var sum: f64 = rhs[row];
            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                if (col < row) {
                    sum -= lu[idx] * out[col];
                }
            }
            out[row] = sum;
        }

        var row: usize = self.dof;
        while (row > 0) {
            row -= 1;
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);
            var sum: f64 = out[row];
            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                if (col > row) {
                    sum -= lu[idx] * out[col];
                }
            }

            const diag_idx = diag_pos[row];
            var diag: f64 = 1.0;
            if (diag_idx >= 0) {
                diag = lu[@intCast(diag_idx)];
                if (diag == 0.0) diag = 1.0;
            }
            out[row] = sum / diag;
        }
    }

    fn matVec(self: GMRESSolver, x: []const f64, out: []f64) void {
        for (0..self.dof) |row| {
            const start: usize = @intCast(self.system.lhs_p[row]);
            const end: usize = @intCast(self.system.lhs_p[row + 1]);
            var sum: f64 = 0.0;
            for (start..end) |idx| {
                const col: usize = @intCast(self.system.lhs_i[idx]);
                sum += self.system.lhs_values[idx] * x[col];
            }
            out[row] = sum;
        }
    }

    fn vSlice(self: *GMRESSolver, work: Work, index: usize) []f64 {
        const start = index * self.dof;
        return work.v[start .. start + self.dof];
    }

    fn hIndex(work: Work, row: usize, col: usize) usize {
        return row + (work.restart + 1) * col;
    }

    fn hGet(self: *GMRESSolver, work: Work, row: usize, col: usize) f64 {
        _ = self;
        return work.h[hIndex(work, row, col)];
    }

    fn hSet(self: *GMRESSolver, work: Work, row: usize, col: usize, value: f64) void {
        _ = self;
        work.h[hIndex(work, row, col)] = value;
    }
};

fn computeGivens(a: f64, b: f64) struct { c: f64, s: f64, r: f64 } {
    if (b == 0.0) {
        return .{ .c = 1.0, .s = 0.0, .r = a };
    }

    if (@abs(b) > @abs(a)) {
        const t = a / b;
        const s = 1.0 / std.math.sqrt(1.0 + t * t);
        return .{ .c = s * t, .s = s, .r = b / s };
    }

    const t = b / a;
    const c = 1.0 / std.math.sqrt(1.0 + t * t);
    return .{ .c = c, .s = c * t, .r = a / c };
}

fn dot(a: []const f64, b: []const f64) f64 {
    var sum: f64 = 0.0;
    for (a, b) |av, bv| {
        sum += av * bv;
    }
    return sum;
}

fn norm(a: []const f64) f64 {
    return std.math.sqrt(dot(a, a));
}
