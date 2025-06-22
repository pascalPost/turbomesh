const std = @import("std");
const types = @import("types.zig");

const Float = types.Float;
const Vec2d = types.Vec2d;

pub fn parseCsvIntoVec2d(allocator: std.mem.Allocator, file_path: []const u8) !std.ArrayList(Vec2d) {
    std.log.info("reading csv: {s}", .{file_path});
    var file = try std.fs.cwd().openFile(file_path, .{});
    defer file.close();

    var buf_reader = std.io.bufferedReader(file.reader());
    var in_stream = buf_reader.reader();

    // we assume that the csv is not bigger than 1000 lines (if it is bigger, reallocations will happen)
    const assumed_csv_lines: usize = 1000;
    var list = try std.ArrayList(Vec2d).initCapacity(allocator, assumed_csv_lines);

    // iterate over all lines
    var i_line: usize = 0;
    var buf: [1024]u8 = undefined;
    while (try in_stream.readUntilDelimiterOrEof(&buf, '\n')) |line| {

        // ignore lines starting with a comment character #
        if (line[0] == '#') continue;

        const expected_n_floats = 2;
        var line_data: Vec2d = undefined;

        // try to split the line into floats
        var i_float: u2 = 0;
        var it = std.mem.tokenizeScalar(u8, line, ' '); // TODO: extend to different seperators
        while (it.next()) |entry| {
            if (i_float == expected_n_floats) {
                std.log.err("csv parsing error: too many entries in line detected\n  file: {s}\n  line {}: {s}", .{ file_path, i_line, line });
                return error.line;
            }

            const float_value = try std.fmt.parseFloat(Float, entry);
            line_data.data[i_float] = float_value;
            i_float += 1;
        }

        if (i_float != expected_n_floats or it.next() != null) {
            std.log.err("csv read error", .{});
        }

        try list.append(line_data);
        i_line += 1;
    }

    list.shrinkAndFree(list.items.len);
    return list;
}

test "parse csv" {
    const allocator = std.testing.allocator;
    const file_path = "examples/T106/T106_ps.dat";
    const data = try parseCsvIntoVec2d(allocator, file_path);
    defer data.deinit();

    try std.testing.expect(types.eql(data.items[0], .{ .data = [_]Float{ 1.127030384, -0.047185256 } }));
    try std.testing.expect(types.eql(data.items[data.items.len - 1], .{ .data = [_]Float{ 1.047805900, 0.000076595 } }));
}
