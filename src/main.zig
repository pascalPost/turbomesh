const std = @import("std");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");

const c = @cImport({
    @cInclude("wayland-client.h");
    @cInclude("wayland-egl.h");
    // @cInclude("xdg-shell-client-protocol.h");

    @cDefine("WL_EGL_PLATFORM", "1");
    @cInclude("EGL/egl.h");
    @cInclude("EGL/eglext.h");
    @cUndef("WL_EGL_PLATFORM");

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
// var xdg_wm_base: ?*c.xdg_wm_base = null;
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
                // else if (std.mem.orderZ(u8, interface, c.xdg_wm_base_interface.name) == .eq) {
                //     xdg_wm_base = @ptrCast(c.wl_registry_bind(registry_, name, &c.xdg_wm_base_interface, version));
                // }
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

    // if (xdg_wm_base == null) {
    //     return error.UndefinedXdgWmBase;
    // }

    // var xdg_wm_base_listener = c.xdg_wm_base_listener{
    //     .ping = struct {
    //         fn handle(
    //             data: ?*anyopaque,
    //             xdg_wm_base_: ?*c.xdg_wm_base,
    //             serial: u32,
    //         ) callconv(.c) void {
    //             _ = data;
    //             c.xdg_wm_base_pong(xdg_wm_base_, serial);
    //         }
    //     }.handle,
    // };
    //
    // if (c.xdg_wm_base_add_listener(xdg_wm_base, &xdg_wm_base_listener, null) == FAILURE) {
    //     return error.FailWmBaseAddListiner;
    // }

    const surface = c.wl_compositor_create_surface(compositor);
    defer c.wl_surface_destroy(surface);

    // const xdg_surface = c.xdg_wm_base_get_xdg_surface(xdg_wm_base, surface);
    // defer c.xdg_surface_destroy(xdg_surface);
    //
    // var xdg_surface_listener = c.xdg_surface_listener{
    //     .configure = xdgSurfaceConfigureHandler,
    // };
    //
    // if (c.xdg_surface_add_listener(xdg_surface, &xdg_surface_listener, null) == FAILURE) {
    //     return error.FailXdgSurfaceAddListener;
    // }
    //
    // const xdg_toplevel = c.xdg_surface_get_toplevel(xdg_surface);
    // defer c.xdg_toplevel_destroy(xdg_toplevel);
    //
    // c.xdg_toplevel_set_title(xdg_toplevel, "turbomesh");
    //
    // var xdg_toplevel_listener = c.xdg_toplevel_listener{
    //     .configure = toplevelConfigureHandler,
    //     .close = toplevelCloseHandler,
    //     .configure_bounds = toplevelConfigureBoundsHandler,
    //     .wm_capabilities = toplevelWmCapabilitiesHandler,
    // };
    //
    // if (c.xdg_toplevel_add_listener(xdg_toplevel, &xdg_toplevel_listener, null) == FAILURE) {
    //     return error.FailedXdgToplevelAddListener;
    // }

    // if (c.wl_display_roundtrip(display) == FAILURE) {
    //     return error.RoundtripFailure;
    // }

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
    defer if (c.eglDestroyContext(egl_display, egl_config) != c.EGL_TRUE) {
        std.log.err("error in egl destroy context", .{});
    };

    egl_window = c.wl_egl_window_create(surface, 800, 600);
    defer c.wl_egl_window_destroy(egl_window);

    egl_surface = c.eglCreatePlatformWindowSurface(egl_display, egl_config, egl_window, null);
    if (egl_surface == c.EGL_NO_SURFACE) {
        return error.FailureEglSurfaceCreate;
    }
    defer if (c.eglDestroySurface(egl_display, egl_surface) != c.EGL_TRUE) {
        std.log.err("err in eglDestroySurface", .{});
    };

    if (c.eglMakeCurrent(egl_display, egl_surface, egl_surface, egl_context) == c.EGL_FALSE) {
        const err = c.eglGetError();
        std.log.err("err in eglMakeCurrent (code: {})", .{err});
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
                    std.log.err("err in libdecor_configuration_get_content_size", .{});

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

        // if (c.libdecor_dispatch(libdecor_context, 0) < 0) {
        //     std.log.err("err in libdecor_dispatch", .{});
        // }
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

// fn xdgSurfaceConfigureHandler(
//     data: ?*anyopaque,
//     xdg_surface: ?*c.xdg_surface,
//     serial: u32,
// ) callconv(.c) void {
//     _ = data;
//
//     // the compositor confidures the surface; ack. the configure event
//     c.xdg_surface_ack_configure(xdg_surface, serial);
//
//     // if(configured) {
//     //     // if alredy configured
//     //     c.wl_surface_commit(surface);
//     // }
//
//     // configured = true;
// }

// fn toplevelConfigureHandler(
//     data: ?*anyopaque,
//     xdg_toplevel: ?*c.xdg_toplevel,
//     width: i32,
//     height: i32,
//     states: [*c]c.wl_array,
// ) callconv(.c) void {
//     _ = data;
//     _ = xdg_toplevel;
//     _ = width;
//     _ = height;
//     _ = states;
// }
//
// fn toplevelCloseHandler(
//     data: ?*anyopaque,
//     xdg_toplevel: ?*c.xdg_toplevel,
// ) callconv(.c) void {
//     _ = data;
//     _ = xdg_toplevel;
//     running = false;
// }
//
// fn toplevelConfigureBoundsHandler(
//     data: ?*anyopaque,
//     xdg_toplevel: ?*c.xdg_toplevel,
//     width: i32,
//     height: i32,
// ) callconv(.c) void {
//     _ = data;
//     _ = xdg_toplevel;
//     _ = width;
//     _ = height;
// }
//
// fn toplevelWmCapabilitiesHandler(
//     data: ?*anyopaque,
//     xdg_toplevel: ?*c.xdg_toplevel,
//     capabilities: [*c]c.wl_array,
// ) callconv(.c) void {
//     _ = data;
//     _ = xdg_toplevel;
//     _ = capabilities;
// }

test "simple test" {
    var list = std.ArrayList(i32).init(std.testing.allocator);
    defer list.deinit(); // try commenting this out and see if zig detects the memory leak!
    try list.append(42);
    try std.testing.expectEqual(@as(i32, 42), list.pop());
}
