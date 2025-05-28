const std = @import("std");
pub const gl = @import("gl");
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

const FAILURE = -1;

var width: usize = 800;
var height: usize = 600;

var display: ?*c.wl_display = null;
var surface: ?*c.wl_surface = null;
var compositor: ?*c.wl_compositor = null;
pub var signal_stop = false;
var egl_window: ?*c.wl_egl_window = null;
var egl_display: c.EGLDisplay = null;
var egl_surface: c.EGLSurface = null;
var egl_context: c.EGLContext = null;

var libdecor_context: ?*c.libdecor = null;
var libdecor_frame: ?*c.libdecor_frame = null;

var libdecor_frame_iface = c.libdecor_frame_interface{
    .configure = struct {
        fn configure(
            frame: ?*c.libdecor_frame,
            config: ?*c.libdecor_configuration,
            user_data: ?*anyopaque,
        ) callconv(.c) void {
            _ = user_data;

            {
                var width_c: c_int = undefined;
                var height_c: c_int = undefined;
                if (c.libdecor_configuration_get_content_size(config, frame, &width_c, &height_c)) {
                    // NOTE: we assume this not to work as the toplevel dimensions should not be configured yet.
                    width = @intCast(width_c);
                    height = @intCast(height_c);
                }
            }

            c.wl_egl_window_resize(egl_window, @intCast(width), @intCast(height), 0, 0);

            const state = c.libdecor_state_new(@intCast(width), @intCast(height));
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

            signal_stop = true;
        }
    }.close,
    .commit = struct {
        fn commit(
            frame: ?*c.libdecor_frame,
            user_data: ?*anyopaque,
        ) callconv(.c) void {
            _ = frame;
            _ = user_data;

            if (c.eglSwapBuffers(egl_display, egl_surface) != c.EGL_TRUE) {
                std.log.err("err in eglSwapBuffers", .{});
            }
        }
    }.commit,
};

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

pub fn swapBuffers() !void {
    if (c.eglSwapBuffers(egl_display, egl_surface) != c.EGL_TRUE) {
        return error.FailedEglSwapBuffers;
    }

    if (c.wl_display_dispatch(display) == FAILURE) {
        return error.FailedDisplayDispatch;
    }
}

pub fn deinit() void {
    c.libdecor_frame_unref(libdecor_frame);
    c.libdecor_unref(libdecor_context);

    eglMakeCurrent(egl_display, c.EGL_NO_SURFACE, c.EGL_NO_CONTEXT) catch {
        std.log.err("err in wayland deinit: eglMakeCurrent to no_surface and no_context failed.", .{});
    };

    c.wl_egl_window_destroy(egl_window);

    if (c.eglDestroySurface(egl_display, egl_surface) != c.EGL_TRUE) {
        std.log.err("err in eglDestroySurface", .{});
    }

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

    if (c.eglTerminate(egl_display) != c.EGL_TRUE) {
        std.log.err("egl terminate error", .{});
    }

    c.wl_surface_destroy(surface);
    c.wl_display_disconnect(display);
}

pub fn init(env: std.process.EnvMap, size: ?struct { usize, usize }, gl_proc_table: *gl.ProcTable) !void {
    if (size) |s| {
        width = s[0];
        height = s[1];
    }

    const xdg_session_type = env.get("XDG_SESSION_TYPE") orelse {
        return error.UnknownSessionType;
    };

    if (std.mem.eql(u8, "x11", xdg_session_type)) {
        return error.UnsupportedX11;
    }

    if (!std.mem.eql(u8, "wayland", xdg_session_type)) {
        return error.UnsupportedSessionType;
    }

    display = c.wl_display_connect(null);

    {
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
    }

    if (compositor == null) {
        return error.UndefinedCompositor;
    }

    surface = c.wl_compositor_create_surface(compositor);

    egl_display = c.eglGetPlatformDisplay(c.EGL_PLATFORM_WAYLAND_KHR, display, null);
    if (egl_display == c.EGL_NO_DISPLAY) {
        return error.EglNoDisplay;
    }

    var egl_major: c.EGLint = undefined;
    var egl_minor: c.EGLint = undefined;
    if (c.eglInitialize(egl_display, &egl_major, &egl_minor) != c.EGL_TRUE) {
        return error.FailureEglInit;
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

    egl_context = c.eglCreateContext(egl_display, egl_config, c.EGL_NO_CONTEXT, &context_attributes);
    if (egl_context == c.EGL_NO_CONTEXT) {
        return error.FailureEglContext;
    }

    egl_window = c.wl_egl_window_create(surface, @intCast(width), @intCast(height));

    egl_surface = c.eglCreatePlatformWindowSurface(egl_display, egl_config, egl_window, null);
    if (egl_surface == c.EGL_NO_SURFACE) {
        return error.FailureEglSurfaceCreate;
    }

    try eglMakeCurrent(egl_display, egl_surface, egl_context);

    libdecor_context = c.libdecor_new(display, &libdecor_iface);

    libdecor_frame = c.libdecor_decorate(libdecor_context, surface, &libdecor_frame_iface, null);

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

    // gl init
    if (!gl_proc_table.init(c.eglGetProcAddress)) {
        return error.GlInitFailed;
    }
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
