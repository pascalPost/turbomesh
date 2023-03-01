// Copyright (c) 2023 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

use glow::HasContext;
use std::mem::size_of;

pub fn create_vao_vbo(
    gl: &glow::Context,
    vertices: &[[f32; 2]],
) -> (glow::VertexArray, glow::Buffer) {
    unsafe {
        // data represents a flat array of f32s that are to be interpreted as vec2s.
        // let triangle_vertices = [0.5f32, 1.0f32, 0.0f32, 0.0f32, 1.0f32, 0.0f32];
        let data_u8 = std::slice::from_raw_parts(
            vertices.as_ptr() as *const u8,
            vertices.len() * 2 * size_of::<f32>(),
        );

        // construct the vertex buffer object (vbo) and upload the data
        let vbo = gl
            .create_buffer()
            .expect("Cannot create vertex buffer object.");
        gl.bind_buffer(glow::ARRAY_BUFFER, Some(vbo));
        gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, data_u8, glow::STATIC_DRAW);

        // construct a vertex array object (vao) to describe the format of the input vbo
        let vao = gl
            .create_vertex_array()
            .expect("Cannot create vertex array");
        gl.bind_vertex_array(Some(vao));
        gl.enable_vertex_attrib_array(0);
        gl.vertex_attrib_pointer_f32(0, 2, glow::FLOAT, false, 2 * size_of::<f32>() as i32, 0);
        assert_eq!(8, 2 * size_of::<f32>());

        (vao, vbo)
    }
}
