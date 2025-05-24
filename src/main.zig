const std = @import("std");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");

const c = @cImport({
    @cInclude("wayland-client.h");
    @cInclude("wayland-egl.h");

    @cDefine("WL_EGL_PLATFORM", "1");
    @cInclude("EGL/egl.h");
    @cInclude("EGL/eglext.h");
    @cUndef("WL_EGL_PLATFORM");

    // NOTE: needed as GNOME does not implement Server Side Decorations...
    // we use this lib to add Client Side Decorations that try to match
    // the used system. At least on my system, the decorations are NOT (!!)
    // matched yet.
    @cInclude("libdecor-0/libdecor.h");

    @cInclude("GL/gl.h");
});

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

const FAILURE = -1;

var compositor: ?*c.wl_compositor = null;
var running = true;
var egl_window: ?*c.wl_egl_window = null;
var egl_display: c.EGLDisplay = null;
var egl_surface: c.EGLSurface = null;

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer std.debug.assert(gpa.deinit() == .ok);
    const allocator = gpa.allocator();

    // check graphics system
    var env = try std.process.getEnvMap(allocator);
    defer env.deinit();

    const xdg_session_type = env.get("XDG_SESSION_TYPE") orelse {
        return error.UnknownSessionType;
    };

    if (std.mem.eql(u8, "x11", xdg_session_type)) {
        return error.UnsupportedX11;
    }

    if (!std.mem.eql(u8, "wayland", xdg_session_type)) {
        return error.UnsupportedSessionType;
    }

    // wayland setup

    const display = c.wl_display_connect(null);
    defer c.wl_display_disconnect(display);

    const registry = c.wl_display_get_registry(display);
    defer c.wl_registry_destroy(registry);

    var registry_listener = c.wl_registry_listener{
        .global = struct {
            fn handle(
                data: ?*anyopaque,
                registry_: ?*c.wl_registry,
                name: u32,
                interface: [*c]const u8,
                version: u32,
            ) callconv(.c) void {
                _ = data;

                if (std.mem.orderZ(u8, interface, c.wl_compositor_interface.name) == .eq) {
                    compositor = @ptrCast(c.wl_registry_bind(registry_, name, &c.wl_compositor_interface, version));
                }
            }
        }.handle,
        .global_remove = struct {
            fn handle(
                data: ?*anyopaque,
                registry_: ?*c.wl_registry,
                name: u32,
            ) callconv(.c) void {
                _ = data;
                _ = registry_;
                _ = name;
            }
        }.handle,
    };

    if (c.wl_registry_add_listener(registry, &registry_listener, null) == FAILURE) {
        return error.AddListenerFailure;
    }

    if (c.wl_display_roundtrip(display) == FAILURE) {
        return error.RoundtripFailure;
    }

    if (compositor == null) {
        return error.UndefinedCompositor;
    }

    const surface = c.wl_compositor_create_surface(compositor);
    defer c.wl_surface_destroy(surface);

    egl_display = c.eglGetPlatformDisplay(c.EGL_PLATFORM_WAYLAND_KHR, display, null);
    if (egl_display == c.EGL_NO_DISPLAY) {
        return error.EglNoDisplay;
    }

    var egl_major: c.EGLint = undefined;
    var egl_minor: c.EGLint = undefined;
    if (c.eglInitialize(egl_display, &egl_major, &egl_minor) != c.EGL_TRUE) {
        return error.FailureEglInit;
    }
    defer {
        if (c.eglTerminate(egl_display) != c.EGL_TRUE) {
            std.log.err("egl terminate error", .{});
        }
    }

    const egl_attributes = [_:c.EGL_NONE]c.EGLint{
        c.EGL_SURFACE_TYPE,    c.EGL_WINDOW_BIT,
        c.EGL_RENDERABLE_TYPE, c.EGL_OPENGL_BIT,
        c.EGL_RED_SIZE,        8,
        c.EGL_GREEN_SIZE,      8,
        c.EGL_BLUE_SIZE,       8,
        c.EGL_ALPHA_SIZE,      8,
    };
    var egl_config: c.EGLConfig = null;
    var num_config: c.EGLint = 0;
    if (c.eglChooseConfig(egl_display, &egl_attributes, &egl_config, 1, &num_config) != c.EGL_TRUE) {
        return error.FailureEglChooseConfig;
    }
    if (num_config > 1) {
        return error.UnexpectedNumEglConfigs;
    }

    if (c.eglBindAPI(c.EGL_OPENGL_API) != c.EGL_TRUE) {
        return error.FailureEglBindApi;
    }

    const context_attributes = [_:c.EGL_NONE]c.EGLint{
        c.EGL_CONTEXT_MAJOR_VERSION, 4,
        c.EGL_CONTEXT_MINOR_VERSION, 5,
    };

    const egl_context = c.eglCreateContext(egl_display, egl_config, c.EGL_NO_CONTEXT, &context_attributes);
    if (egl_context == c.EGL_NO_CONTEXT) {
        return error.FailureEglContext;
    }
    defer {
        eglMakeCurrent(egl_display, c.EGL_NO_SURFACE, c.EGL_NO_CONTEXT) catch {};

        const res = c.eglDestroyContext(egl_display, egl_context);
        switch (res) {
            c.EGL_TRUE => {},
            c.EGL_FALSE => {
                std.log.err("destruction of EGL context failed.", .{});
                eglGetError() catch {};
            },
            c.EGL_BAD_DISPLAY => std.log.err("err in eglDestroyContext: given display is not an EGL display connection.", .{}),
            c.EGL_NOT_INITIALIZED => std.log.err("err in eglDestroyContext: given display has not been intialized.", .{}),
            c.EGL_BAD_CONTEXT => std.log.err("err in eglDestroyContext: given context is not an EGL rendering context.", .{}),
            else => std.log.err("unknown return code from eglDestroyContext.", .{}),
        }
    }

    egl_window = c.wl_egl_window_create(surface, 800, 600);
    defer c.wl_egl_window_destroy(egl_window);

    egl_surface = c.eglCreatePlatformWindowSurface(egl_display, egl_config, egl_window, null);
    if (egl_surface == c.EGL_NO_SURFACE) {
        return error.FailureEglSurfaceCreate;
    }
    defer if (c.eglDestroySurface(egl_display, egl_surface) != c.EGL_TRUE) {
        std.log.err("err in eglDestroySurface", .{});
    };

    try eglMakeCurrent(egl_display, egl_surface, egl_context);

    var libdecor_iface = c.libdecor_interface{
        .@"error" = struct {
            fn handler(
                context: ?*c.libdecor,
                err: c.libdecor_error,
                message: [*c]const u8,
            ) callconv(.c) void {
                _ = context;
                std.log.err("libdecor error {} : {s}", .{ err, message });
                unreachable;
            }
        }.handler,
    };
    const libdecor_context = c.libdecor_new(display, &libdecor_iface);
    defer c.libdecor_unref(libdecor_context);

    // const Window = struct {
    //     egl_window: ?*c.wl_egl_window,
    //     egl_display: c.EGLDisplay,
    //     egl_surface: c.EGLSurface,
    // };

    var libdecor_frame_iface = c.libdecor_frame_interface{
        .configure = struct {
            fn configure(
                frame: ?*c.libdecor_frame,
                config: ?*c.libdecor_configuration,
                user_data: ?*anyopaque,
            ) callconv(.c) void {
                _ = user_data;
                // const window: *Window = @ptrCast(@alignCast(user_data));
                // const egl_window_ = window.egl_window;

                var width: c_int = undefined;
                var height: c_int = undefined;
                if (!c.libdecor_configuration_get_content_size(config, frame, &width, &height)) {
                    // NOTE: we assume this not to work as the toplevel dimensions should not be configured yet.

                    // TODO: remove
                    width = 800;
                    height = 600;
                }

                c.wl_egl_window_resize(egl_window, width, height, 0, 0);

                const state = c.libdecor_state_new(width, height);
                c.libdecor_frame_commit(frame, state, config);
                c.libdecor_state_free(state);
            }
        }.configure,
        .close = struct {
            fn close(
                frame: ?*c.libdecor_frame,
                user_data: ?*anyopaque,
            ) callconv(.c) void {
                _ = frame;
                _ = user_data;

                running = false;
            }
        }.close,
        .commit = struct {
            fn commit(
                frame: ?*c.libdecor_frame,
                user_data: ?*anyopaque,
            ) callconv(.c) void {
                _ = frame;
                _ = user_data;

                // const window: *Window = @ptrCast(@alignCast(user_data));
                // const egl_display_ = window.egl_display;
                // const egl_surface_ = window.egl_surface;

                if (c.eglSwapBuffers(egl_display, egl_surface) != c.EGL_TRUE) {
                    std.log.err("err in eglSwapBuffers", .{});
                }
            }
        }.commit,
    };

    // var window = Window{
    //     .egl_window = egl_window,
    //     .egl_display = egl_display,
    //     .egl_surface = egl_surface,
    // };

    const libdecor_frame = c.libdecor_decorate(libdecor_context, surface, &libdecor_frame_iface, null);
    defer c.libdecor_frame_unref(libdecor_frame);

    c.libdecor_frame_set_title(libdecor_frame, "TEST");
    c.libdecor_frame_set_app_id(libdecor_frame, "TEST");
    c.libdecor_frame_map(libdecor_frame);

    if (c.wl_display_roundtrip(display) == FAILURE) {
        return error.RoundtripFailure;
    }

    if (c.wl_display_roundtrip(display) == FAILURE) {
        return error.RoundtripFailure;
    }

    if (c.libdecor_dispatch(libdecor_context, 0) < 0) {
        std.log.err("err in libdecor_dispatch", .{});
    }

    c.glClearColor(1, 1, 0.5, 1);

    while (running) {
        c.glClear(c.GL_COLOR_BUFFER_BIT);
        c.glFlush();

        if (c.eglSwapBuffers(egl_display, egl_surface) != c.EGL_TRUE) {
            return error.FailedEglSwapBuffers;
        }

        if (c.wl_display_dispatch(display) == FAILURE) {
            return error.FailedDisplayDispatch;
        }
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

fn eglMakeCurrent(
    egl_display_: c.EGLDisplay,
    egl_surface_: c.EGLSurface,
    egl_context_: c.EGLContext,
) !void {
    if (c.eglMakeCurrent(egl_display_, egl_surface_, egl_surface_, egl_context_) == c.EGL_FALSE) {
        try eglGetError();
    }
}

fn eglGetError() !void {
    const err = c.eglGetError();
    if (err == c.EGL_SUCCESS) return;

    std.log.err("err in eglMakeCurrent (code: 0x{x})", .{err});
    switch (err) {
        c.EGL_BAD_ACCESS => return error.EglBadAccess,
        c.EGL_BAD_MATCH => return error.EglBadMatch,
        c.EGL_BAD_NATIVE_WINDOW => return error.EglBadNativeWindow,
        c.EGL_BAD_CONTEXT => return error.EglBadContext,
        c.EGL_BAD_ALLOC => return error.EglBadAlloc,
        c.EGL_BAD_CURRENT_SURFACE => return error.EglBadCurrentSurface,
        c.EGL_BAD_SURFACE => return error.EglBadSurface,
        else => return error.FailedEglMakeCurrent,
    }
}

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
