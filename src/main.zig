const std = @import("std");
const types = @import("types.zig");
const cgns = @import("cgns.zig");
const spline = @import("spline.zig");
const glfw = @import("zglfw");
pub const gl = @import("gl");

var gl_proc_table: gl.ProcTable = undefined;

const Mat2d = types.Mat2d;
const Index2d = types.Index2d;
const Vec2d = types.Vec2d;
const Float = types.Float;
const Index = types.Index;

pub fn main() !void {
    // var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    // defer std.debug.assert(gpa.deinit() == .ok);
    // const allocator = gpa.allocator();
    //
    // // check graphics system
    // var env = try std.process.getEnvMap(allocator);
    // defer env.deinit();

    try glfw.init();
    defer glfw.terminate();

    const gl_major = 4;
    const gl_minor = 5;
    glfw.windowHint(.context_version_major, gl_major);
    glfw.windowHint(.context_version_minor, gl_minor);
    glfw.windowHint(.opengl_profile, .opengl_core_profile);
    glfw.windowHint(.opengl_forward_compat, true);
    glfw.windowHint(.client_api, .opengl_api);
    glfw.windowHint(.doublebuffer, true);

    const window = try glfw.createWindow(800, 600, "zig-gamedev: minimal_glfw_gl", null);
    defer glfw.destroyWindow(window);

    glfw.makeContextCurrent(window);

    glfw.swapInterval(1);

    gl.makeProcTableCurrent(&gl_proc_table);
    if (!gl_proc_table.init(glfw.getProcAddress)) {
        return error.LoadGlAddressesFailed;
    }
    defer gl.makeProcTableCurrent(null);

    // TODO: add dark mode detection
    // TODO: add mode selection
    gl.ClearColor(1, 1, 1, 1); // white bg

    {
        var width: c_int = undefined;
        var height: c_int = undefined;
        glfw.getWindowSize(window, &width, &height);

        var xscale: f32 = undefined;
        var yscale: f32 = undefined;
        glfw.getWindowContentScale(window, &xscale, &yscale);

        const w: c_int = @intFromFloat(@as(f32, @floatFromInt(width)) * xscale);
        const h: c_int = @intFromFloat(@as(f32, @floatFromInt(height)) * yscale);
        gl.Viewport(0, 0, w, h);
    }

    const vertex_shader_source =
        \\#version 330 core
        \\layout (location = 0) in vec2 aPos;
        \\uniform vec2 uOffset;
        \\
        \\void main()
        \\{
        \\  gl_Position = vec4(aPos + uOffset, 0.0, 1.0);
        \\}
    ;

    const fragment_shader_source =
        \\#version 330 core
        \\out vec4 FragColor;
        \\
        \\void main()
        \\{
        \\  FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);
        \\}
    ;

    const vertex_shader = gl.CreateShader(gl.VERTEX_SHADER);
    gl.ShaderSource(vertex_shader, 1, &.{vertex_shader_source}, null);
    gl.CompileShader(vertex_shader);

    {
        var info_log: [512:0]u8 = undefined;
        var success: c_int = undefined;
        gl.GetShaderiv(vertex_shader, gl.COMPILE_STATUS, &success);
        if (success == gl.FALSE) {
            gl.GetShaderInfoLog(vertex_shader, info_log.len, null, &info_log);
            std.debug.print("vertex shader compilation error:\n{str}\n", .{info_log});
        }
    }

    const fragment_shader = gl.CreateShader(gl.FRAGMENT_SHADER);
    gl.ShaderSource(fragment_shader, 1, &.{fragment_shader_source}, null);
    gl.CompileShader(fragment_shader);

    {
        var info_log: [512:0]u8 = undefined;
        var success: c_int = undefined;
        gl.GetShaderiv(fragment_shader, gl.COMPILE_STATUS, &success);
        if (success == gl.FALSE) {
            gl.GetShaderInfoLog(fragment_shader, info_log.len, null, &info_log);
            std.debug.print("fragment shader compilation error:\n{str}\n", .{info_log});
        }
    }

    const program = gl.CreateProgram();
    if (program == 0) {
        return error.CreateProgramFailed;
    }

    gl.AttachShader(program, vertex_shader);
    gl.AttachShader(program, fragment_shader);
    gl.LinkProgram(program);

    {
        var success: c_int = undefined;
        gl.GetProgramiv(program, gl.LINK_STATUS, &success);
        if (success == gl.FALSE) {
            var info_log: [512:0]u8 = undefined;
            gl.GetProgramInfoLog(program, info_log.len, null, &info_log);
            std.debug.print("program linkage error:\n{str}\n", .{info_log});
        }
    }

    const offset_location = gl.GetUniformLocation(program, "uOffset");
    if (offset_location == -1) {
        return error.GetUniformLocationFailed;
    }

    gl.UseProgram(program);
    gl.DeleteShader(vertex_shader);
    gl.DeleteShader(fragment_shader);

    const vertices = [_]f32{
        -0.5, -0.5,
        0.5,  -0.5,
        0.0,  0.5,
    };

    // TODO: allow to zoom with a scaling option

    // TODO: fill this with an offset
    gl.Uniform2f(offset_location, 0.0, 0.0);

    var vao: gl.uint = undefined;
    gl.GenVertexArrays(1, (&vao)[0..1]);
    gl.BindVertexArray(vao);

    var vbo: gl.uint = undefined;
    gl.GenBuffers(1, (&vbo)[0..1]);

    gl.BindBuffer(gl.ARRAY_BUFFER, vbo);
    gl.BufferData(gl.ARRAY_BUFFER, @sizeOf(@TypeOf(vertices)), &vertices, gl.STATIC_DRAW);

    gl.VertexAttribPointer(0, 2, gl.FLOAT, gl.FALSE, 2 * @sizeOf(gl.float), 0);
    gl.EnableVertexAttribArray(0);

    gl.PointSize(10);

    while (!window.shouldClose()) {
        glfw.pollEvents();

        gl.Clear(gl.COLOR_BUFFER_BIT);

        gl.UseProgram(program);
        gl.BindVertexArray(vao);
        gl.DrawArrays(gl.POINTS, 0, 3);

        window.swapBuffers();
    }
}
