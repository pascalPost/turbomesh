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
