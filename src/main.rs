// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use glow::HasContext;

use eframe::egui;
use egui::mutex::Mutex;
use egui::*;
use std::mem::size_of;
use std::sync::Arc;

fn main() {
    let options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1024.0, 1480.0)),
        multisampling: 0,
        renderer: eframe::Renderer::Glow,
        ..Default::default()
    };
    eframe::run_native(
        "turbomesh",
        options,
        Box::new(|cc| Box::new(TurbomeshApp::new(cc))),
    );
}

struct Point {
    name: String,
    pos: [f32; 2],
}

struct TurbomeshApp {
    /// Behind an `Arc<Mutex<…>>` so we can pass it to [`egui::PaintCallback`] and paint later.
    point_rendering: Arc<Mutex<PointRendering>>,
    offset: Vec2,
    // TODO add option to define variables
    points: Vec<Point>,
}

impl TurbomeshApp {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let gl = cc
            .gl
            .as_ref()
            .expect("You need to run eframe with the glow backend");
        Self {
            point_rendering: Arc::new(Mutex::new(PointRendering::new(gl))),
            offset: Vec2::new(0.0, 0.0),
            points: vec![],
        }
    }
}

impl eframe::App for TurbomeshApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::SidePanel::left("my_side_panel")
            .resizable(false)
            .show(ctx, |ui| {
                CollapsingHeader::new("Geometry")
                    .default_open(false)
                    .show(ui, |ui| {
                        CollapsingHeader::new("Points")
                            .default_open(false)
                            .show(ui, |ui| {
                                self.points
                                    .iter()
                                    .for_each(|p| if ui.button(p.name.clone()).clicked() {});
                                if ui.button("+").clicked() {
                                    self.points.push(Point {
                                        name: String::from("test"),
                                        pos: [0.0, 0.0],
                                    })
                                }
                            });
                        CollapsingHeader::new("Curves")
                            .default_open(false)
                            .show(ui, |ui| if ui.button("+").clicked() {});
                    });
                CollapsingHeader::new("Mesh")
                    .default_open(false)
                    .show(ui, |_| {});
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::Frame::canvas(ui.style()).show(ui, |ui| {
                if ui.button("Reset View").clicked() {
                    self.offset = Vec2::new(0.0, 0.0);
                }
                self.custom_painting(ui);
            });
        });
    }

    fn on_exit(&mut self, gl: Option<&glow::Context>) {
        if let Some(gl) = gl {
            self.point_rendering.lock().destroy(gl);
        }
    }
}

impl TurbomeshApp {
    fn custom_painting(&mut self, ui: &mut egui::Ui) {
        let (rect, response) =
            // ui.allocate_exact_size(egui::Vec2::splat(300.0), egui::Sense::drag());
            // ui.all
            ui.allocate_exact_size(ui.available_size(), egui::Sense::drag());

        const MOUSE_SENSITIVITY: f32 = 0.5 * 1e-2;
        const INVERT_MOUSE: f32 = -1.0; // to not invert set to 1.0

        self.offset += response.drag_delta() * Vec2::new(1.0, INVERT_MOUSE) * MOUSE_SENSITIVITY;

        // Clone locals so we can move them into the paint callback:
        let offset = self.offset;
        let point_rendering = self.point_rendering.clone();

        let callback = egui::PaintCallback {
            rect,
            callback: std::sync::Arc::new(egui_glow::CallbackFn::new(move |_info, painter| {
                point_rendering.lock().render(painter.gl(), offset);
            })),
        };
        ui.painter().add(callback);
    }
}

struct PointRendering {
    program: glow::Program,
    vao_vbo_vec: Vec<(glow::VertexArray, glow::Buffer)>,
}

fn create_vao_vbo(gl: &glow::Context, data: &[[f32; 2]]) -> (glow::VertexArray, glow::Buffer) {
    unsafe {
        // data represents a flat array of f32s that are to be interpreted as vec2s.
        // let triangle_vertices = [0.5f32, 1.0f32, 0.0f32, 0.0f32, 1.0f32, 0.0f32];
        let data_u8 = std::slice::from_raw_parts(
            data.as_ptr() as *const u8,
            data.len() * 2 * size_of::<f32>(),
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

impl PointRendering {
    fn new(gl: &glow::Context) -> Self {
        use glow::HasContext as _;

        let shader_version = if cfg!(target_arch = "wasm32") {
            "#version 300 es"
        } else {
            "#version 330"
        };

        unsafe {
            gl.enable(glow::PROGRAM_POINT_SIZE);
        }

        let mut vao_vbo_vec: Vec<(glow::VertexArray, glow::Buffer)> = vec![];
        vao_vbo_vec.push(create_vao_vbo(
            gl,
            [[0.5f32, 1.0f32], [0.0f32, 0.0f32], [1.0f32, 0.0f32]].as_slice(),
        ));
        vao_vbo_vec.push(create_vao_vbo(
            gl,
            [[0.0f32, 0.5f32], [-0.5f32, -0.5f32], [0.5f32, -0.5f32]].as_slice(),
        ));

        unsafe {
            let program = gl.create_program().expect("Cannot create program");

            let (vertex_shader_source, fragment_shader_source) = (
                r#"
                    in vec2 in_position;
                    out vec4 v_color;
                    uniform vec2 u_offset_xy;
                    void main() {
                        gl_PointSize = 10;
                        gl_Position = vec4(in_position, 0.0, 1.0);
                        gl_Position += vec4(u_offset_xy, 0.0, 0.0);
                    }
                "#,
                r#"
                    precision mediump float;
                    in vec2 position;
                    out vec4 out_color;
                    uniform float blue;
                    void main() {
                        out_color = vec4(1.0, 1.0, 1.0, 1.0);
                    }
                "#,
            );

            let shader_sources = [
                (glow::VERTEX_SHADER, vertex_shader_source),
                (glow::FRAGMENT_SHADER, fragment_shader_source),
            ];

            let shaders: Vec<_> = shader_sources
                .iter()
                .map(|(shader_type, shader_source)| {
                    let shader = gl
                        .create_shader(*shader_type)
                        .expect("Cannot create shader");
                    gl.shader_source(shader, &format!("{}\n{}", shader_version, shader_source));
                    gl.compile_shader(shader);
                    if !gl.get_shader_compile_status(shader) {
                        panic!("{}", gl.get_shader_info_log(shader));
                    }
                    gl.attach_shader(program, shader);
                    shader
                })
                .collect();

            gl.link_program(program);
            if !gl.get_program_link_status(program) {
                panic!("{}", gl.get_program_info_log(program));
            }

            for shader in shaders {
                gl.detach_shader(program, shader);
                gl.delete_shader(shader);
            }

            Self {
                program,
                vao_vbo_vec,
            }
        }
    }

    fn destroy(&self, gl: &glow::Context) {
        use glow::HasContext as _;
        unsafe {
            gl.delete_program(self.program);
            for (vao, vbo) in self.vao_vbo_vec.iter() {
                gl.delete_vertex_array(*vao);
                gl.delete_buffer(*vbo);
            }
        }
    }

    fn render(&self, gl: &glow::Context, offset: Vec2) {
        use glow::HasContext as _;
        unsafe {
            gl.use_program(Some(self.program));
            gl.uniform_2_f32(
                gl.get_uniform_location(self.program, "u_offset_xy")
                    .as_ref(),
                offset.x,
                offset.y,
            );

            for (vao, _) in self.vao_vbo_vec.iter() {
                gl.bind_vertex_array(Some(*vao));
                gl.draw_arrays(glow::POINTS, 0, 3);
            }
        }
    }
}
