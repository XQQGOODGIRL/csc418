#include <iostream>
#include <thread>
#include "assignment_setup.h"

#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <visualization.h>

//Simulation State
Util::MacGrid mac_grid;
Util::MarkerParticle marker_particles;

bool simulating = true;

bool simulation_callback() {
    
    while(simulating) {
        simulate(mac_grid, marker_particles);
    }
    
    return false;
}

bool draw_callback(igl::opengl::glfw::Viewer &viewer) {
    
    draw(marker_particles);
    
    return false;
}

int main(int argc, char **argv) {

    std::cout<<"Start Fluid Simulating\n";
    
    // GRID SETUP
    init_state(mac_grid, marker_particles);
    
    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulation_callback);
    simulation_thread.detach();
    
    // VISULIZATION SETUP
    Visualize::setup(marker_particles, true);
    Visualize::viewer().callback_post_draw = &draw_callback;
    Visualize::viewer().launch();
}
