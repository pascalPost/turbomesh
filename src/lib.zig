const gl = @import("gl");
const State = @import("state.zig").State;

export fn update(state: *State) callconv(.c) void {
    gl.makeProcTableCurrent(state.gl_proc_table_ptr);

    gl.Clear(gl.COLOR_BUFFER_BIT);

    gl.UseProgram(state.program);
    gl.Uniform1f(state.scale_location, state.scale);
    gl.Uniform2f(state.center_location, state.center[0], state.center[1]);
    gl.Uniform2f(state.offset_location, @floatCast(state.offset[0]), @floatCast(state.offset[1]));
    gl.BindVertexArray(state.vao);
    gl.Uniform4f(state.color_location, 1.0, 1.0, 1.0, 1.0);
    // gl.DrawArrays(gl.POINTS, 0, @intCast(point_buffer.len / 2));
    gl.DrawElements(gl.LINES, @intCast(state.element_buffer.len), gl.UNSIGNED_INT, 0);
}
