const std = @import("std");

const help_message =
    \\turbomesh - Turbomachinery mesh processing tool
    \\
    \\USAGE:
    \\    turbomesh <config.json>
    \\
    \\ARGUMENTS:
    \\    <config.json>    Path to the JSON configuration file
    \\
    \\OPTIONS:
    \\    -h, --help       Print this help message and exit
    \\    -v, --version    Show version information
    \\
    \\DESCRIPTION:
    \\    turbomesh is a mesh processing tool for turbomachinery configurations
    \\    based on a user-defined JSON configuration.
    \\
    \\    The configuration file should specify all necessary input/output paths,
    \\    transformation rules, and processing options.
    \\
    \\EXAMPLE:
    \\    turbomesh config.json
;

pub fn parseArgs() !std.fs.File {
    // TODO: add tests
    // TODO: add test for relative file path
    // TODO: add test for absolute file path
    var args = std.process.args();

    if (args.inner.count <= 1) {
        _ = try std.io.getStdErr().write(
            \\Error: missing config file
            \\
            \\Usage: turbomesh <config.json>
            \\Try 'turbomesh --help' for more information.
        );
        std.process.exit(64);
    }

    var config_file_path: [:0]const u8 = undefined;
    while (args.next()) |arg| {
        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            _ = try std.io.getStdOut().write(help_message);
            std.process.exit(0);
        } else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-v")) {
            // TODO: retrieve the version from the build.zig
            _ = try std.io.getStdOut().writer().print("v{s}\n", .{"0.1.0"});
        } else {
            config_file_path = arg;
        }
    }

    const config_file = std.fs.cwd().openFile(config_file_path, .{}) catch {
        _ = try std.io.getStdErr().writer().print(
            \\Error: could not open file "{s}": file not found
            \\
            \\Usage: turbomesh <config.json>
            \\Try 'turbomesh --help' for more information.
        , .{config_file_path});
        std.process.exit(66);
    };

    std.log.info("Config file: {s}", .{config_file_path});

    return config_file;
}
