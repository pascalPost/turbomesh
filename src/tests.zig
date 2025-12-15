// Copyright (c) 2025 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

comptime {
    _ = @import("cgns.zig");
    _ = @import("spline.zig");
    _ = @import("templates/blade.zig");
    _ = @import("types.zig");
    _ = @import("discrete.zig");
    _ = @import("smoothing/umfpack.zig");
    _ = @import("smoothing/petsc.zig");
    _ = @import("boundary.zig");
    _ = @import("templates/O4H.zig");
}
