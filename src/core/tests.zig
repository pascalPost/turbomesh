// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

const config = @import("config");

comptime {
    _ = @import("cgns.zig");
    _ = @import("spline.zig");
    _ = @import("input.zig");
    _ = @import("types.zig");
    _ = @import("discrete.zig");
    _ = @import("smoothing/BiCGStab.zig");
    _ = @import("smoothing/GMRES.zig");
    _ = @import("boundary.zig");
    _ = @import("templates/O4H.zig");

    if (config.use_umfpack) {
        _ = @import("smoothing/umfpack.zig");
    }

    if (config.use_petsc) {
        _ = @import("smoothing/petsc.zig");
    }
}
