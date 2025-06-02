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

    @cInclude("linux/input-event-codes.h");
});

const log = std.log.scoped(.wayland);

const FAILURE = -1;

var width: usize = 800;
var height: usize = 600;

var display: ?*c.wl_display = null;
var surface: ?*c.wl_surface = null;
var compositor: ?*c.wl_compositor = null;
var seat: ?*c.wl_seat = null;
var pointer: ?*c.wl_pointer = null;
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

var seat_listener = c.wl_seat_listener{
    .capabilities = struct {
        fn capabilities(data: ?*anyopaque, seat_: ?*c.wl_seat, capabilities_: u32) callconv(.c) void {
            _ = data;
            if (capabilities_ & c.WL_SEAT_CAPABILITY_POINTER != 0) {
                pointer = c.wl_seat_get_pointer(seat_);
            } else {
                std.log.err("no pointer for main sear.", .{});
                unreachable;
            }
        }
    }.capabilities,
    .name = struct {
        fn name(user_data: ?*anyopaque, seat_: ?*c.wl_seat, name_: [*c]const u8) callconv(.c) void {
            _ = user_data;
            _ = seat_;
            _ = name_;
        }
    }.name,
};

var pointer_state: struct {
    dragging: bool,
    last: struct { f64, f64 },
    offset: struct { f64, f64 },
} = .{ .dragging = false, .last = .{ 0, 0 }, .offset = .{ 0, 0 } };

const PointerEventMask = struct {
    const enter: usize = 0;
    const leave: usize = 1 << 1;
    const motion: usize = 1 << 2;
    const button: usize = 1 << 3;
    const axis: usize = 1 << 4;
    const axis_source: usize = 1 << 5;
    const axis_stop: usize = 1 << 6;
    const axis_discrete: usize = 1 << 7;

    value: u8 = 0,
};

const PointerEvent = struct {
    mask: PointerEventMask = .{},
    surface_x: c.wl_fixed_t = 0,
    surface_y: c.wl_fixed_t = 0,
    button: u32 = 0,
    state: u32 = 0,
    time: u32 = 0,
    serial: u32 = 0,
    axis: [2]struct {
        valid: bool = false,
        value: ?c.wl_fixed_t = null,
        discrete: i32 = 0,
    } = .{ .{}, .{} },
    axis_source: u32 = 0,
};

var pointer_event = PointerEvent{};

// TODO: change debug print to log with reasonable level
var pointer_listener = c.wl_pointer_listener{
    .enter = struct {
        fn enter(
            user_data: ?*anyopaque,
            pointer_: ?*c.wl_pointer,
            serial: u32,
            surface_: ?*c.wl_surface,
            surface_x: c.wl_fixed_t,
            surface_y: c.wl_fixed_t,
        ) callconv(.c) void {
            _ = user_data;
            _ = pointer_;
            _ = surface_;

            log.debug("enter at {} {}\n", .{ surface_x, surface_y });

            pointer_event.mask.value |= PointerEventMask.enter;
            pointer_event.serial = serial;
            pointer_event.surface_x = surface_x;
            pointer_event.surface_y = surface_y;
        }
    }.enter,

    .leave = struct {
        fn leave(user_data: ?*anyopaque, pointer_: ?*c.wl_pointer, serial: u32, surface_: ?*c.wl_surface) callconv(.c) void {
            _ = user_data;
            _ = pointer_;
            _ = surface_;

            log.debug("leave\n", .{});

            pointer_event.mask.value |= PointerEventMask.leave;
            pointer_event.serial = serial;
        }
    }.leave,

    .motion = struct {
        fn motion(user_data: ?*anyopaque, pointer_: ?*c.wl_pointer, time: u32, surface_x: c.wl_fixed_t, surface_y: c.wl_fixed_t) callconv(.c) void {
            _ = user_data;
            _ = pointer_;

            log.debug("motion {} {}\n", .{ surface_x, surface_y });

            pointer_event.mask.value |= PointerEventMask.motion;
            pointer_event.time = time;
            pointer_event.surface_x = surface_x;
            pointer_event.surface_y = surface_y;
        }
    }.motion,

    .button = struct {
        fn button(user_data: ?*anyopaque, pointer_: ?*c.wl_pointer, serial: u32, time: u32, button_: u32, state: u32) callconv(.c) void {
            _ = user_data;
            _ = pointer_;

            log.debug("button {} , time {} , state {}\n", .{ button_, time, state });

            pointer_event.mask.value |= PointerEventMask.button;
            pointer_event.serial = serial;
            pointer_event.time = time;
            pointer_event.button = button_;
            pointer_event.state = state;
        }
    }.button,

    .axis = struct {
        fn axis(user_data: ?*anyopaque, pointer_: ?*c.wl_pointer, time: u32, axis_: u32, value: c.wl_fixed_t) callconv(.c) void {
            _ = user_data;
            _ = pointer_;

            pointer_event.mask.value |= PointerEventMask.axis;
            pointer_event.time = time;
            pointer_event.axis[axis_].valid = true;
            pointer_event.axis[axis_].value = value;
        }
    }.axis,

    .frame = struct {
        fn frame(user_data: ?*anyopaque, pointer_: ?*c.wl_pointer) callconv(.c) void {
            _ = user_data;
            _ = pointer_;

            log.debug("frame\n", .{});

            // const time = pointer_event.time;
            // std.log.debug("pointer frame @ {}\n", .{time});

            const x = c.wl_fixed_to_double(pointer_event.surface_x);
            const y = c.wl_fixed_to_double(pointer_event.surface_y);

            if (pointer_event.button == c.BTN_LEFT and pointer_event.state == c.WL_POINTER_BUTTON_STATE_PRESSED) {
                pointer_state.dragging = true;

                pointer_state.offset[0] += x - pointer_state.last[0];
                pointer_state.offset[1] += y - pointer_state.last[1];

                std.debug.print("dragging offset {} {}\n", .{ pointer_state.offset[0], pointer_state.offset[1] });
            } else {
                pointer_state.dragging = false;
            }

            pointer_state.last = .{ x, y };
        }
    }.frame,
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
                    } else if (std.mem.orderZ(u8, interface, c.wl_seat_interface.name) == .eq) {
                        seat = @ptrCast(c.wl_registry_bind(registry_, name, &c.wl_seat_interface, version));
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

    if (seat == null) {
        return error.UndefinedSeat;
    }

    if (c.wl_seat_add_listener(seat, &seat_listener, null) == FAILURE) {
        return error.SeatAddListenerFailed;
    }

    pointer = c.wl_seat_get_pointer(seat);

    if (c.wl_pointer_add_listener(pointer, &pointer_listener, null) == FAILURE) {
        return error.PointerAddListenerFailed;
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
