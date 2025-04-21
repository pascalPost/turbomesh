// x11 + opengl doc: https://github.com/gamedevtech/X11OpenGLWindow

const std = @import("std");
const c = @cImport({
    @cInclude("X11/Xlib.h");
});

pub fn main() !void {
    // connect fo X server on default display
    const display = c.XOpenDisplay(null);
    if (display == null) {
        return error.XServerConnectionNoSuccess;
    }

    const screen_id = c.XDefaultScreen(display);
    const root_window = c.XDefaultRootWindow(display);
    const window = c.XCreateSimpleWindow(display, root_window, 0, 0, 320, 200, 1, c.BlackPixel(display, screen_id), c.WhitePixel(display, screen_id));

    // handle delete window
    var atom_wm_delete_window = c.XInternAtom(display, "WM_DELETE_WINDOW", c.False);
    if (c.XSetWMProtocols(display, window, &atom_wm_delete_window, 1) == 0) {
        return error.XHandleWindowDeleteNoSuccess;
    }

    _ = c.XSelectInput(display, window, c.StructureNotifyMask | c.KeyPressMask | c.KeyReleaseMask | c.KeymapStateMask);

    _ = c.XClearWindow(display, window);
    _ = c.XMapRaised(display, window);

    // message loop
    var event: c.XEvent = undefined;
    while (true) {
        _ = c.XNextEvent(display, &event);
        switch (event.type) {
            c.ClientMessage => if (event.xclient.data.l[0] == atom_wm_delete_window) break,
            c.DestroyNotify => break,
            c.KeymapNotify => _ = c.XRefreshKeyboardMapping(&event.xmapping),
            c.KeyPress => {},
            c.KeyRelease => {},
            // c.ButtonPress => {},
            // c.ButtonRelease => {},
            // c.MotionNotify => {},
            // c.EnterNotify => {
            //     // mouse enters window
            // },
            // c.LeaveNotify => {
            //     // mouse leaves window
            // },
            else => {},
        }
    }

    // cleanup
    _ = c.XDestroyWindow(display, window);
    _ = c.XCloseDisplay(display);
}
