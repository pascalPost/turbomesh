const std = @import("std");
const types = @import("types.zig");
const templates = @import("templates/templates.zig");
const smoothing = @import("smoothing/smoothing.zig");
const machine = @import("machine.zig");

const native_arch = @import("builtin").target.cpu.arch;

pub const ProfileInputTag = enum {
    data, // inline data
    csv, // csv file
};

pub const ProfileInput = union(ProfileInputTag) {
    data: struct {
        down: []const [2]f64,
        up: []const [2]f64,
    },
    csv: struct {
        down_csv_path: []const u8,
        up_csv_path: []const u8,
    },
};

const Input = struct {
    template: templates.Template,
    smoothing: struct {
        iterations: usize = 0,
        solver: smoothing.solver.Option,
        wall_control_function: smoothing.wall_control_function.Algorithm = .{ .laplace = {} },
    },
    geometry: struct {
        pitch: types.Float,
        profile: ProfileInput,
    },
};

pub fn create_profile(allocator: std.mem.Allocator, input: ProfileInput, scale: types.Float) !machine.Profile {
    switch (input) {
        .data => |d| {
            const up = try allocVec2dFromTuples(allocator, d.up);
            defer allocator.free(up);
            errdefer allocator.free(up);
            const down = try allocVec2dFromTuples(allocator, d.down);
            defer allocator.free(down);
            errdefer allocator.free(down);

            if (scale != 1.0) {
                for (down) |*point| {
                    point.data[0] *= scale;
                    point.data[1] *= scale;
                }
                for (up) |*point| {
                    point.data[0] *= scale;
                    point.data[1] *= scale;
                }
            }

            return try machine.Profile.init(allocator, down, up);
        },
        .csv => |csv| {
            if (native_arch == .wasm32) {
                return error.ExternalFilesNotSupported;
            } else {
                const down = try readSide(allocator, csv.down_csv_path);
                defer down.deinit();
                const up = try readSide(allocator, csv.up_csv_path);
                defer up.deinit();

                if (scale != 1.0) {
                    for (down.items) |*point| {
                        point.data[0] *= scale;
                        point.data[1] *= scale;
                    }
                    for (up.items) |*point| {
                        point.data[0] *= scale;
                        point.data[1] *= scale;
                    }
                }

                return try machine.Profile.init(allocator, down.items, up.items);
            }
        },
    }
}

pub fn allocVec2dFromTuples(allocator: std.mem.Allocator, tuples: []const [2]f64) ![]types.Vec2d {
    const out = try allocator.alloc(types.Vec2d, tuples.len);
    for (tuples, out) |source, *dest| {
        dest.* = types.Vec2d.init(source[0], source[1]);
    }
    return out;
}

fn readSide(allocator: std.mem.Allocator, csv_path: []const u8) !std.array_list.Managed(types.Vec2d) { //     const side = try csv.parseCsvIntoVec2d(allocator, csv_path);
    const csv = @import("csv.zig");
    const side = try csv.parseCsvIntoVec2d(allocator, csv_path);
    // TODO add check of axial direction: does the vector form origion to end point in upstream direction (?)
    if (side.items[0].data[0] > side.items[side.items.len - 1].data[0]) {
        std.mem.reverse(types.Vec2d, side.items);
    }
    return side;
}

test "read a blade profile from csv files" {
    const allocator = std.testing.allocator;
    var profile = try create_profile(allocator, .{ .csv = .{ .down_csv_path = "./examples/T106/T106_ps.dat", .up_csv_path = "./examples/T106/T106_ss.dat" } });
    profile.deinit();
}
