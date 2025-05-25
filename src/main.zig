const std = @import("std");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");
const platform = @import("gui/platform.zig");

const c = @cImport({
    @cInclude("GL/gl.h");
});

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    // check graphics system
    var env = try std.process.getEnvMap(allocator);
    defer env.deinit();

    // const platform = guiPlatform.init(env);
    try platform.init(env, .{ 800, 600 });
    defer platform.deinit();

    c.glClearColor(1, 1, 0.5, 1);

    while (true) {
        c.glClear(c.GL_COLOR_BUFFER_BIT);
        c.glFlush();

        if (try platform.swapBuffersAndReturnStopSignal()) break;
    }

    // // Initialize GUI
    // var gui_instance = try gui.Gui.init(allocator);
    // defer gui_instance.deinit();
    //
    // const block_points = [_]Mat2d{
    //     try Mat2d.init(allocator, .{ 21, 17 }),
    // };
    //
    // defer {
    //     for (block_points) |block| {
    //         block.deinit(allocator);
    //     }
    // }
    //
    // // Initialize mesh with some test data
    // for (block_points) |block| {
    //     const size = block.size;
    //
    //     var idx: usize = 0;
    //     var i: usize = 0;
    //     while (i < size[0]) : (i += 1) {
    //         var j: usize = 0;
    //         while (j < size[1]) : (j += 1) {
    //             block.data[idx] = Vec2d{ @floatFromInt(i), @floatFromInt(j) };
    //             idx += 1;
    //         }
    //     }
    // }
    //
    // // Set the mesh data in the App
    // app.setMesh(block_points[0]);
    //
    // // Main loop
    // while (try app.update()) {
    //     // The App update function handles all the rendering and event processing
    // }
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
