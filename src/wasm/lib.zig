const std = @import("std");
const core = @import("core");

extern "env" fn console_log(ptr: [*]const u8, len: usize) void;

fn wasmLog(
    comptime level: std.log.Level,
    comptime scope: @TypeOf(.enum_literal),
    comptime format: []const u8,
    args: anytype,
) void {
    var buf: [2048]u8 = undefined;

    const level_txt = comptime level.asText();
    const prefix2 = if (scope == .default) ": " else "(" ++ @tagName(scope) ++ "): ";
    const full_fmt = level_txt ++ prefix2 ++ format ++ "\n";

    const msg = std.fmt.bufPrint(&buf, full_fmt, args) catch {
        const short = "[log message too long]\n";
        console_log(short.ptr, short.len);
        return;
    };

    console_log(msg.ptr, msg.len);
}

pub const std_options: std.Options = .{
    // .log_level = .debug,
    .logFn = wasmLog,
};

// state
var mesh_global: ?core.discrete.Mesh = null;

fn runFromInput(allocator: std.mem.Allocator, input: core.input.Input) !void {
    freeMesh();

    const profile = try core.input.create_profile(allocator, input.geometry.profile, input.geometry.scale);
    defer profile.deinit();

    const geometry = core.machine.Geometry.init(input.geometry.scale * input.geometry.pitch, profile);

    var mesh = try input.template.run(allocator, geometry);
    errdefer mesh.deinit();

    try core.smoothing.smooth.mesh(
        allocator,
        &mesh,
        input.smoothing.iterations,
        input.smoothing.solver,
        input.smoothing.wall_control_function,
    );

    mesh_global = mesh;
}

export fn freeMesh() void {
    if (mesh_global) |m| {
        m.deinit();
        mesh_global = null;
    }
}

export fn alloc(size: usize) usize {
    const allocator = std.heap.wasm_allocator;
    const buffer = allocator.alloc(u8, size) catch return 0;
    return @intFromPtr(buffer.ptr);
}

export fn dealloc(ptr: usize, size: usize) void {
    if (ptr == 0 or size == 0) return;
    const allocator = std.heap.wasm_allocator;
    const buffer = @as([*]u8, @ptrFromInt(ptr))[0..size];
    allocator.free(buffer);
}

export fn run(input_ptr: [*]const u8, input_len: usize) void {
    const allocator = std.heap.wasm_allocator;

    if (input_len == 0) {
        std.log.err("encountered 0 input.", .{});
        return;
    }

    const input_json = input_ptr[0..input_len];
    var parsed = std.json.parseFromSlice(core.input.Input, allocator, input_json, .{}) catch |err| {
        std.log.err("failed to parse input: {s}", .{@errorName(err)});
        return;
    };
    defer parsed.deinit();
    runFromInput(allocator, parsed.value) catch |err| {
        std.log.err("failed to run input: {s}", .{@errorName(err)});
    };
    return;
}

export fn blocksCount() usize {
    return mesh_global.?.blocks.items.len;
}

export fn blockSizeI(block_idx: usize) u32 {
    const sz = mesh_global.?.blocks.items[block_idx].points.size;
    return @intCast(sz[0]); // i
}

export fn blockSizeJ(block_idx: usize) u32 {
    const sz = mesh_global.?.blocks.items[block_idx].points.size;
    return @intCast(sz[1]); // j
}

export fn blockSize(block_idx: usize, out: [*]u32) void {
    const sz = mesh_global.?.blocks.items[block_idx].points.size;
    out[0] = @intCast(sz[0]); // i
    out[1] = @intCast(sz[1]); // j
}

export fn blockPointsPtr(block_idx: usize) [*]const f64 {
    const pts = mesh_global.?.blocks.items[block_idx].points.data;
    return @ptrCast(pts.ptr); // points packed as x0,y0,x1,y1,...
}

export fn blockPointsLen(block_idx: usize) usize {
    return mesh_global.?.blocks.items[block_idx].points.data.len * 2; // f64 entries
}
