// Copyright (c) 2022 Pascal Post
// This code is licensed under AGPL license (see LICENSE.txt for details)

#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")] // hide console window on Windows in release

use eframe::egui;
use egui::mutex::Mutex;
use egui::*;
use glow::{HasContext, NativeShader};
use std::sync::Arc;
use turbomesh::gl::create_vao_vbo;
use turbomesh::mesh::Mesh;

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

struct Points {
    name: String,
    coords: Vec<[f32; 2]>,
}

struct TurbineTemplateData {
    ps_csv_path: Option<String>,
    ss_csv_path: Option<String>,
}

impl TurbineTemplateData {
    fn new() -> Self {
        Self {
            ps_csv_path: None,
            ss_csv_path: None,
        }
    }
}

struct Shader {
    gl: Arc<glow::Context>,
    shader: NativeShader,
}

impl Shader {
    fn from_source(gl: &Arc<glow::Context>, shader_type: u32, shader_source: &str) -> Self {
        let shader_version = if cfg!(target_arch = "wasm32") {
            "#version 300 es"
        } else {
            "#version 330"
        };

        unsafe {
            let shader = gl.create_shader(shader_type).expect("Cannot create shader");

            gl.shader_source(shader, &format!("{}\n{}", shader_version, shader_source));
            gl.compile_shader(shader);
            if !gl.get_shader_compile_status(shader) {
                panic!("{}", gl.get_shader_info_log(shader));
            }

            Self {
                gl: gl.clone(),
                shader,
            }
        }
    }

    pub fn vertex_shader_from_source(gl: &Arc<glow::Context>, source: &str) -> Self {
        Self::from_source(gl, glow::VERTEX_SHADER, source)
    }

    pub fn fragment_shader_from_source(gl: &Arc<glow::Context>, source: &str) -> Self {
        Self::from_source(gl, glow::FRAGMENT_SHADER, source)
    }
}

impl Drop for Shader {
    fn drop(&mut self) {
        unsafe {
            self.gl.delete_shader(self.shader);
        }
    }
}

struct Program {
    gl: Arc<glow::Context>,
    program: glow::NativeProgram,
}

impl Program {
    fn new(
        gl: &Arc<glow::Context>,
        vertex_shader_source: &str,
        fragment_shader_source: &str,
    ) -> Self {
        unsafe {
            let program = gl.create_program().expect("Cannot create program");

            let vertex_shader = Shader::vertex_shader_from_source(gl, vertex_shader_source);
            let fragment_shader = Shader::fragment_shader_from_source(gl, fragment_shader_source);

            gl.attach_shader(program, vertex_shader.shader);
            gl.attach_shader(program, fragment_shader.shader);

            gl.link_program(program);
            if !gl.get_program_link_status(program) {
                panic!("{}", gl.get_program_info_log(program));
            }

            gl.detach_shader(program, vertex_shader.shader);
            gl.detach_shader(program, fragment_shader.shader);

            Self {
                gl: gl.clone(),
                program,
            }
        }
    }
}

impl Drop for Program {
    fn drop(&mut self) {
        unsafe {
            self.gl.delete_program(self.program);
        }
    }
}

struct PointsViewRenderer {
    gl_program: Arc<Program>,

    /// VAO and VBO to render points
    gl_data: (glow::VertexArray, glow::Buffer),

    num_points: usize,
}

impl PointsViewRenderer {
    fn new(gl_program: &Arc<Program>, points: &Points) -> Self {
        let gl_data = create_vao_vbo(gl_program.gl.as_ref(), &points.coords.as_slice());

        Self {
            gl_program: gl_program.clone(),
            gl_data,
            num_points: points.coords.len(),
        }
    }

    fn render(&self, offset: Vec2, scale: Vec2) {
        use glow::HasContext as _;
        unsafe {
            self.gl_program
                .gl
                .use_program(Some(self.gl_program.program));
            self.gl_program.gl.uniform_2_f32(
                self.gl_program
                    .gl
                    .get_uniform_location(self.gl_program.program, "u_offset_xy")
                    .as_ref(),
                offset.x,
                offset.y,
            );
            self.gl_program.gl.uniform_2_f32(
                self.gl_program
                    .gl
                    .get_uniform_location(self.gl_program.program, "u_scale_xy")
                    .as_ref(),
                scale.x,
                scale.y,
            );

            self.gl_program.gl.bind_vertex_array(Some(self.gl_data.0));
            self.gl_program
                .gl
                .draw_arrays(glow::POINTS, 0, self.num_points as i32);
        }
    }
}

struct MeshViewRenderer {
    gl_program: Arc<Program>,

    /// VAOs, VBOs and number of points to render the blocks of the mesh
    gl_data: Vec<(glow::VertexArray, glow::Buffer, usize)>,
    // /// vector for the visible of the blocks
    // block_visibility: Vec<bool>,
}

impl MeshViewRenderer {
    fn new(gl_program: &Arc<Program>, mesh: &Mesh) -> Self {
        let mut gl_data = Vec::<(glow::VertexArray, glow::Buffer, usize)>::new();
        gl_data.reserve(mesh.blocks.len());

        for block in mesh.blocks.iter() {
            let mut coords = Vec::<[f32; 2]>::new();
            coords.reserve(block.coords.size());
            block
                .coords
                .as_slice()
                .iter()
                .for_each(|v| coords.push([v.0 as f32, v.1 as f32]));

            let (vao, vbo) = create_vao_vbo(gl_program.gl.as_ref(), coords.as_slice());

            gl_data.push((vao, vbo, coords.len()));
        }

        // let block_visibility = vec![true; mesh.blocks.len()];

        Self {
            gl_program: gl_program.clone(),
            gl_data,
            // block_visibility,
        }
    }

    fn render(&self, offset: Vec2, scale: Vec2) {
        use glow::HasContext as _;
        unsafe {
            self.gl_program
                .gl
                .use_program(Some(self.gl_program.program));
            self.gl_program.gl.uniform_2_f32(
                self.gl_program
                    .gl
                    .get_uniform_location(self.gl_program.program, "u_offset_xy")
                    .as_ref(),
                offset.x,
                offset.y,
            );
            self.gl_program.gl.uniform_2_f32(
                self.gl_program
                    .gl
                    .get_uniform_location(self.gl_program.program, "u_scale_xy")
                    .as_ref(),
                scale.x,
                scale.y,
            );

            for (vao, _, count) in self.gl_data.iter() {
                self.gl_program.gl.bind_vertex_array(Some(*vao));
                self.gl_program
                    .gl
                    .draw_arrays(glow::POINTS, 0, *count as i32);
            }
        }
    }
}

const VERTEX_SHADER_SOURCE: &str = r#"
in vec2 in_position;
out vec4 v_color;
uniform vec2 u_offset_xy;
uniform vec2 u_scale_xy;
void main() {
    gl_PointSize = 10;
    gl_Position = vec4(in_position, 0.0, 1.0);
    gl_Position += vec4(u_offset_xy, 0.0, 0.0);
    gl_Position *= vec4(u_scale_xy, 1.0, 1.0);
}
"#;

const FRAGMENT_SHADER_SOURCE: &str = r#"
precision mediump float;
in vec2 position;
out vec4 out_color;
uniform float blue;
void main() {
    out_color = vec4(1.0, 1.0, 1.0, 1.0);
}
"#;

struct MeshView {
    mesh: Mesh,
    renderer: Arc<Mutex<MeshViewRenderer>>,
}

impl MeshView {
    fn new(gl_program: &Arc<Program>, mesh: Mesh) -> Self {
        let renderer = MeshViewRenderer::new(gl_program, &mesh);
        Self {
            mesh: mesh,
            renderer: Arc::new(Mutex::new(renderer)),
        }
    }
}

struct PointsView {
    points: Points,
    renderer: Arc<Mutex<PointsViewRenderer>>,
}

impl PointsView {
    fn new(gl_program: &Arc<Program>, points: Points) -> Self {
        let renderer = PointsViewRenderer::new(gl_program, &points);
        Self {
            points: points,
            renderer: Arc::new(Mutex::new(renderer)),
        }
    }
}

struct TurbomeshApp {
    // gl: Arc<glow::Context>,
    program: Arc<Program>,
    /// Behind an `Arc<Mutex<…>>` so we can pass it to [`egui::PaintCallback`] and paint later.
    // point_rendering: Arc<Mutex<PointRendering>>,
    offset: Vec2,

    scale: Vec2,

    // TODO add option to define variables
    points_views: Vec<PointsView>,

    mesh_view: MeshView,

    show_turbine_template: bool,
    turbine_template_data: Option<TurbineTemplateData>,
}

impl TurbomeshApp {
    fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let gl = cc
            .gl
            .as_ref()
            .expect("You need to run eframe with the glow backend");

        let program = Arc::new(Program::new(
            gl,
            VERTEX_SHADER_SOURCE,
            FRAGMENT_SHADER_SOURCE,
        ));

        Self {
            // gl: gl.clone(),
            program: program.clone(),
            // point_rendering: Arc::new(Mutex::new(PointRendering::new(gl.as_ref()))),
            offset: Vec2::new(0.0, 0.0),
            scale: Vec2::new(1.0, 1.0),
            points_views: vec![],
            mesh_view: MeshView::new(&program, Mesh::new()),
            show_turbine_template: false,
            turbine_template_data: None,
        }
    }
}

impl eframe::App for TurbomeshApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::SidePanel::left("side_panel")
            .resizable(true)
            .show(ctx, |ui| {
                CollapsingHeader::new("Templates")
                    .default_open(false)
                    .show(ui, |ui| {
                        if ui.button("Turbine").clicked() {
                            self.show_turbine_template = true;
                        }
                    });

                CollapsingHeader::new("Geometry")
                    .default_open(false)
                    .show(ui, |ui| {
                        CollapsingHeader::new("Points")
                            .default_open(false)
                            .show(ui, |ui| {
                                self.points_views
                                    .iter()
                                    .for_each(|p| if ui.button(p.points.name.clone()).clicked() {});
                                if ui.button("+").clicked() {
                                    self.points_views.push(PointsView::new(
                                        &self.program,
                                        Points {
                                            name: String::from("test"),
                                            coords: vec![[0.0, 0.0]; 1],
                                        },
                                    ))
                                }
                            });
                        CollapsingHeader::new("Curves")
                            .default_open(false)
                            .show(ui, |ui| if ui.button("+").clicked() {});
                    });

                CollapsingHeader::new("Mesh")
                    .default_open(false)
                    .show(ui, |ui| {
                        for block in self.mesh_view.mesh.blocks.iter() {
                            ui.horizontal(|ui| {
                                if ui.button(block.name.as_str()).clicked() {}

                                let mut visible = true;
                                ui.checkbox(&mut visible, "");

                                if visible {}
                            });
                        }
                    });

                ui.separator();

                if self.show_turbine_template {
                    if self.turbine_template_data.is_none() {
                        self.turbine_template_data = Some(TurbineTemplateData::new());
                    }

                    let turbine_template_data = self
                        .turbine_template_data
                        .as_mut()
                        .expect("Object should exist.");

                    ui.heading("Turbine Template");

                    ui.horizontal_wrapped(|ui| {
                        ui.label("Pressure Side");

                        if ui
                            .button(
                                turbine_template_data
                                    .ps_csv_path
                                    .as_ref()
                                    .map_or("Open csv file...", |path| path.as_str()),
                            )
                            .clicked()
                        {
                            if let Some(path) = rfd::FileDialog::new().pick_file() {
                                turbine_template_data.ps_csv_path =
                                    Some(path.display().to_string());
                            }
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        ui.label("Suction Side");

                        if ui
                            .button(
                                turbine_template_data
                                    .ss_csv_path
                                    .as_ref()
                                    .map_or("Open csv file...", |path| path.as_str()),
                            )
                            .clicked()
                        {
                            if let Some(path) = rfd::FileDialog::new().pick_file() {
                                turbine_template_data.ss_csv_path =
                                    Some(path.display().to_string());
                            }
                        }
                    });

                    ui.horizontal_wrapped(|ui| {
                        if ui.button("Reset").clicked() {
                            self.show_turbine_template = true;
                            self.turbine_template_data = None;
                        }
                        if ui.button("Apply").clicked() {
                            // TODO add red highlight of missing data

                            if self
                                .turbine_template_data
                                .as_ref()
                                .unwrap()
                                .ps_csv_path
                                .is_some()
                                && self
                                    .turbine_template_data
                                    .as_ref()
                                    .unwrap()
                                    .ss_csv_path
                                    .is_some()
                            {
                                let (_, mesh) = turbomesh::turbine::run_turbine_template(
                                    self.turbine_template_data
                                        .as_ref()
                                        .unwrap()
                                        .ps_csv_path
                                        .as_ref()
                                        .unwrap()
                                        .as_str(),
                                    self.turbine_template_data
                                        .as_ref()
                                        .unwrap()
                                        .ss_csv_path
                                        .as_ref()
                                        .unwrap()
                                        .as_str(),
                                );

                                self.mesh_view = MeshView::new(&self.program, mesh);

                                self.show_turbine_template = false;
                            }
                        }
                        if ui.button("Close").clicked() {
                            self.show_turbine_template = false;
                        }
                    });
                }
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::Frame::canvas(ui.style()).show(ui, |ui| {
                if ui.button("Reset View").clicked() {
                    self.offset = Vec2::new(0.0, 0.0);
                    self.scale = Vec2::new(1.0, 1.0);
                }
                if ui.button("+").clicked() {
                    self.scale += Vec2::new(1.0, 1.0);
                }
                if ui.button("-").clicked() {
                    self.scale -= Vec2::new(1.0, 1.0);
                }
                self.custom_painting(ui);
            });
        });
    }
}

impl TurbomeshApp {
    fn custom_painting(&mut self, ui: &mut egui::Ui) {
        let (rect, response) =
            ui.allocate_exact_size(ui.available_size_before_wrap(), egui::Sense::drag());

        const MOUSE_SENSITIVITY: f32 = 0.5 * 1e-2;
        const INVERT_MOUSE: f32 = -1.0; // to not invert set to 1.0

        self.offset += response.drag_delta() * Vec2::new(1.0, INVERT_MOUSE) * MOUSE_SENSITIVITY;

        // Clone locals so we can move them into the paint callback:
        let offset = self.offset;
        let scale = self.scale;

        let points_views_renderes: Vec<Arc<Mutex<PointsViewRenderer>>> = self
            .points_views
            .iter()
            .map(|p| p.renderer.clone())
            .collect();

        let mesh_view_renderer = self.mesh_view.renderer.clone();

        let callback = egui::PaintCallback {
            rect,
            callback: Arc::new(egui_glow::CallbackFn::new(move |_info, painter| {
                unsafe {
                    painter.gl().enable(glow::PROGRAM_POINT_SIZE);
                }

                mesh_view_renderer.lock().render(offset, scale);
                points_views_renderes
                    .iter()
                    .for_each(|pvr| pvr.lock().render(offset, scale));
            })),
        };

        ui.painter().add(callback);
    }
}
