// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const std = @import("std");
const Preconditioner = @import("preconditioner.zig").Preconditioner;
const RowCompressedMatrixSystem2d = @import("smooth.zig").RowCompressedMatrixSystem2d;

const log = std.log.scoped(.bicgstab_solver);

pub const BiCGStabSolver = struct {
    system: RowCompressedMatrixSystem2d,
    buffer: ?[]f64 = null,
    dof: usize,
    seeded_initial_guess: bool = false,
    preconditioner: Preconditioner = .diagonal,
    ilu_values: ?[]f64 = null,
    ilu_diag_pos: ?[]c_int = null,
    ilu_marker: ?[]c_int = null,
    max_iters: usize = 1000,
    rtol: f64 = 1e-6,
    atol: f64 = 1e-8,

    const Work = struct {
        diag_inv: []f64,
        r: []f64,
        r_hat: []f64,
        p: []f64,
        v: []f64,
        s: []f64,
        t: []f64,
        precond: []f64,
    };

    pub fn init(system: RowCompressedMatrixSystem2d, preconditioner: Preconditioner) BiCGStabSolver {
        return .{
            .system = system,
            .dof = system.rhs_x.len,
            .preconditioner = preconditioner,
        };
    }

    pub fn deinit(self: *BiCGStabSolver) void {
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

    fn precondition(self: *BiCGStabSolver, work: Work) !void {
        switch (self.preconditioner) {
            .diagonal => self.updateDiagonalInverse(work.diag_inv),
            .ilu0 => try self.updateIlu0(),
        }
    }

    pub fn solve(self: *BiCGStabSolver) !void {
        if (!self.seeded_initial_guess) {
            self.seedInitialGuess();
        }

        const work = try self.ensureWorkspace();

        try self.system.fillXSpecific();
        try self.precondition(work);
        self.solveSystem(self.system.rhs_x, self.system.x_new, work, "x");

        try self.system.fillYSpecific();
        try self.precondition(work);
        self.solveSystem(self.system.rhs_y, self.system.y_new, work, "y");
    }

    fn ensureWorkspace(self: *BiCGStabSolver) !Work {
        const lanes = 8;
        const total = self.dof * lanes;

        if (self.buffer == null or self.buffer.?.len != total) {
            if (self.buffer) |buf| {
                self.system.allocator.free(buf);
            }
            self.buffer = try self.system.allocator.alloc(f64, total);
        }

        const buf = self.buffer.?;
        var offset: usize = 0;

        const diag_inv = buf[offset .. offset + self.dof];
        offset += self.dof;

        const r = buf[offset .. offset + self.dof];
        offset += self.dof;

        const r_hat = buf[offset .. offset + self.dof];
        offset += self.dof;

        const p = buf[offset .. offset + self.dof];
        offset += self.dof;

        const v = buf[offset .. offset + self.dof];
        offset += self.dof;

        const s = buf[offset .. offset + self.dof];
        offset += self.dof;

        const t = buf[offset .. offset + self.dof];
        offset += self.dof;

        const precond = buf[offset .. offset + self.dof];

        return .{
            .diag_inv = diag_inv,
            .r = r,
            .r_hat = r_hat,
            .p = p,
            .v = v,
            .s = s,
            .t = t,
            .precond = precond,
        };
    }

    fn seedInitialGuess(self: *BiCGStabSolver) void {
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

    fn updateDiagonalInverse(self: *BiCGStabSolver, diag_inv: []f64) void {
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
    fn updateIlu0(self: *BiCGStabSolver) !void {
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

    fn solveSystem(self: *BiCGStabSolver, rhs: []const f64, x: []f64, work: Work, label: []const u8) void {
        const breakdown_eps = 1e-30;

        self.matVec(x, work.v);
        for (rhs, work.v, 0..) |rhs_val, ax_val, i| {
            work.r[i] = rhs_val - ax_val;
        }

        @memcpy(work.r_hat, work.r);

        const norm_b = norm(rhs);
        var norm_r = norm(work.r);
        const tol = @max(self.atol, self.rtol * norm_b);

        if (norm_r <= tol) return;

        @memset(work.p, 0);
        @memset(work.v, 0);

        var rho_old: f64 = 1.0;
        var alpha: f64 = 1.0;
        var omega: f64 = 1.0;
        var iter: usize = 0;

        while (iter < self.max_iters) : (iter += 1) {
            const rho_new = dot(work.r_hat, work.r);
            if (@abs(rho_new) < breakdown_eps) {
                break;
            }

            const beta = (rho_new / rho_old) * (alpha / omega);
            for (0..self.dof) |i| {
                work.p[i] = work.r[i] + beta * (work.p[i] - omega * work.v[i]);
            }

            self.applyPreconditioner(work.p, work.precond, work);

            self.matVec(work.precond, work.v);

            const denom = dot(work.r_hat, work.v);
            if (@abs(denom) < breakdown_eps) {
                break;
            }

            alpha = rho_new / denom;

            for (0..self.dof) |i| {
                work.s[i] = work.r[i] - alpha * work.v[i];
            }

            for (0..self.dof) |i| {
                x[i] += alpha * work.precond[i];
            }

            const norm_s = norm(work.s);
            if (norm_s <= tol) {
                return;
            }

            self.applyPreconditioner(work.s, work.precond, work);

            self.matVec(work.precond, work.t);

            const t_dot_t = dot(work.t, work.t);
            if (@abs(t_dot_t) < breakdown_eps) {
                break;
            }

            omega = dot(work.t, work.s) / t_dot_t;
            if (@abs(omega) < breakdown_eps) {
                break;
            }

            for (0..self.dof) |i| {
                x[i] += omega * work.precond[i];
            }

            for (0..self.dof) |i| {
                work.r[i] = work.s[i] - omega * work.t[i];
            }

            norm_r = norm(work.r);
            if (norm_r <= tol) {
                return;
            }

            rho_old = rho_new;
        }

        const residual = if (norm_b > 0.0) norm_r / norm_b else norm_r;
        log.warn("bicgstab solve did not converge for {s}: iter={}, residual={e:.3}", .{ label, iter, residual });
    }

    fn applyPreconditioner(self: *BiCGStabSolver, rhs: []const f64, out: []f64, work: Work) void {
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
    fn applyIlu0(self: *BiCGStabSolver, rhs: []const f64, out: []f64) void {
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

    fn matVec(self: BiCGStabSolver, x: []const f64, out: []f64) void {
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
};

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
