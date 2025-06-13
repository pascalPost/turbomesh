const std = @import("std");
const builtin = @import("builtin");

const State = @import("state.zig").State;

var abs_path_buffer: [2048]u8 = undefined;
var abs_path: []u8 = undefined;

const lib_path = switch (builtin.os.tag) {
    .linux => "zig-out/lib/libturbomesh.so",
    .macos => "zig-out/lib/libturbomesh.dylib",
    else => @compileError("unknown OS"),
};

const InitFn = *const fn (*State) callconv(.c) void;
// const DeinitFn = *const fn (*anyopaque) callconv(.c) void;
const UpdateFn = *const fn (*State) callconv(.c) void;

pub const Lib = struct {
    lib_ptr: std.DynLib,
    mod_time: i128,

    init_fn: InitFn,
    // deinit_fn: DeinitFn,
    update_fn: UpdateFn,

    pub fn open(state: *State) !Lib {
        const cwd = std.fs.cwd();
        const stat = try cwd.statFile(lib_path);
        abs_path = try cwd.realpath(lib_path, &abs_path_buffer);
        const res = try load();

        res.init_fn(state);

        return .{
            .lib_ptr = res.lib_ptr,
            .mod_time = stat.mtime,
            .init_fn = res.init_fn,
            // .deinit_fn = res.deinit_fn,
            .update_fn = res.update_fn,
        };
    }

    pub fn close(self: *Lib) void {
        self.lib_ptr.close();
    }

    pub fn reload(self: *Lib, state: *State) !void {
        const cwd = std.fs.cwd();
        const stat = try cwd.statFile(lib_path);
        if (stat.mtime <= self.mod_time) return; // return if no loading is needed.

        self.lib_ptr.close();

        const res = try load();
        self.lib_ptr = res.lib_ptr;
        self.init_fn = res.init_fn;
        self.update_fn = res.update_fn;
        self.mod_time = stat.mtime;

        self.init_fn(state);
    }

    // pub fn init(self: Lib, state: *State) void {
    //     return self.init_fn(state);
    // }

    // pub fn deinit(self: Lib, state: *anyopaque) void {
    //     return self.deinit_fn(state);
    // }

    pub fn update(self: Lib, state: *State) void {
        return self.update_fn(state);
    }
};

fn load() !struct {
    lib_ptr: std.DynLib,
    init_fn: InitFn,
    // deinit_fn: DeinitFn,
    update_fn: UpdateFn,
} {
    var lib = try std.DynLib.open(abs_path);

    const init_fn = lib.lookup(InitFn, "init") orelse return error.LookupInitFnFailed;
    // const deinit_fn: DeinitFn = lib.lookup(DeinitFn, "deinit") orelse return error.LookupDeinitFnFailed;
    const update_fn: UpdateFn = lib.lookup(UpdateFn, "update") orelse return error.LookupUpdateFnFailed;

    std.log.debug("lib loaded successfully.", .{});

    return .{
        .lib_ptr = lib,
        .init_fn = init_fn,
        // .deinit_fn = deinit_fn,
        .update_fn = update_fn,
    };
}
