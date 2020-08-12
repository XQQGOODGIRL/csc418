#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#define IMGUI_DEFINE_MATH_OPERATORS

#include <igl/unproject.h>
#include <pick_nearest_vertices.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

//stl
#include <vector>
#include <array>
#include <deque>

//Eigen
#include <Eigen/Dense>

#include <Util.h>

namespace Visualize {
    
    //custom phase space plot
    bool plot_phase_space(const char *label, ImVec2 q_bounds, ImVec2 q_dot_bounds, const Eigen::VectorXd &q, const Eigen::VectorXd &q_dot);
    
    
    
    void setup(const Util::MarkerParticle &marker_particles, bool ps_plot = false);
    
    
    void update_particle_positions(Eigen::MatrixXd pos);
    
    void set_picking_tolerance(double);
    
    //UI methods
    bool mouse_down(igl::opengl::glfw::Viewer &viewer, int x, int y);
    
    bool mouse_up(igl::opengl::glfw::Viewer &viewer, int x, int y);
    
    bool mouse_move(igl::opengl::glfw::Viewer &viewer, int x, int y);
    
    igl::opengl::glfw::Viewer & viewer();
    
    igl::opengl::glfw::imgui::ImGuiMenu & viewer_menu();
    
    const Eigen::Vector3d & mouse_world();
    
    const Eigen::Vector3d & mouse_drag_world();
    
    
    const std::vector<unsigned int> & picked_vertices();
    
    bool is_mouse_dragging();
    
}


#endif
