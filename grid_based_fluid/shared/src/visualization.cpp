#include <visualization.h>

//libigl viewer
namespace Visualize {
    
    igl::opengl::glfw::Viewer g_viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    
    Util::MarkerParticle const *g_marker_particles;
    
    //cache for phase space data
    std::deque<std::pair<float, float> > g_state;
    //picking variables
    std::vector<unsigned int > g_picked_vertices;
    
    //mouse UI state variables
    bool g_mouse_dragging = false;
    double g_picking_tol = 0.001;
    Eigen::Vector3d g_mouse_win; //mouse window coordinates
    Eigen::Vector3d g_mouse_drag; //last mouse drag vector
    Eigen::Vector3d g_mouse_world;
    Eigen::Vector3d g_mouse_drag_world; //mouse drag delta in the world space
    
}

igl::opengl::glfw::imgui::ImGuiMenu & Visualize::viewer_menu() { return menu; }

igl::opengl::glfw::Viewer & Visualize::viewer() { return g_viewer; }

bool Visualize::plot_phase_space(const char *label, ImVec2 q_bounds, ImVec2 q_dot_bounds, const Eigen::VectorXd &q, const Eigen::VectorXd &q_dot) {
    
    using namespace ImGui;
    
    unsigned int cache_size = 10000;
    unsigned int num_lines = 5; //always odd number because I want lines to cross at 0,0
    
    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    
    if(g_state.size() > cache_size) {
        g_state.pop_front();
    }
    
    //update plotting cache
    g_state.push_back(std::make_pair(q(0), q_dot(0)));
    
    //ImGUI stuff that I don't understand (taken from code example here: https://github.com/ocornut/imgui/issues/786)
    const ImGuiStyle& Style = GetStyle();
    const ImGuiIO& IO = GetIO();
    ImDrawList* DrawList = GetWindowDrawList();
    ImGuiWindow* Window = GetCurrentWindow();
    
    if (Window->SkipItems)
        return false;
    
    // header and spacing
    int hovered = IsItemActive() || IsItemHovered(); // IsItemDragged() ?
    Dummy(ImVec2(0,3));
    
    // prepare canvas
    ImVec2 avail = GetContentRegionAvail();
    ImVec2 Canvas(ImMin(avail.x, avail.y), ImMin(avail.x, avail.y));
    
    Canvas = CalcItemSize(Canvas, style.FramePadding.x * 2.0f, style.FramePadding.y * 2.0f);
    ImRect bb(Window->DC.CursorPos, Window->DC.CursorPos + Canvas);
    
    const ImGuiID id = Window->GetID(label);
    
    RenderFrame(bb.Min, bb.Max, GetColorU32(ImGuiCol_FrameBg, 1), true, Style.FrameRounding);
    
    //local grid coordinates are -1,1 in both directions
    auto pix_to_normalized =  [&bb, &Canvas](ImVec2 pix) {  return ImVec2((pix.x-bb.Min.x)/Canvas.x,(pix.y-bb.Min.y)/Canvas.y); };
    auto normalized_to_pix =  [&bb, &Canvas] (ImVec2 norm) {  return ImVec2(norm.x*Canvas.x + bb.Min.x, norm.y*Canvas.y + bb.Min.y); };
    auto data_to_normalized =  [&q_bounds, &q_dot_bounds] (ImVec2 state) {  return ImVec2((state.x - q_bounds.x)/(q_bounds.y - q_bounds.x), (state.y - q_dot_bounds.x)/(q_dot_bounds.y - q_dot_bounds.x)); };
    
    //background grid centered on origin
    for (float i = 0.f; i <= 1.f; i+= 1.f/static_cast<float>(num_lines-1)) {
        DrawList->AddLine(
                          normalized_to_pix(ImVec2(i, 0.f)),
                          normalized_to_pix(ImVec2(i, 1.f)),
                          GetColorU32(ImGuiCol_TextDisabled), 1.2);
    }
    
    for (float i = 0.f; i <= 1.f; i+= 1.f/static_cast<float>(num_lines-1)) {
        DrawList->AddLine(
                          normalized_to_pix(ImVec2(0.f, i)),
                          normalized_to_pix(ImVec2(1.f, i)),
                          GetColorU32(ImGuiCol_TextDisabled), 1.2);
    }
    
    //plot phase space trajectory
    bool clip_p1;
    bool clip_p2;
    for(unsigned int i=0; i<g_state.size()-1; ++i) {
        
        clip_p1 = false;
        clip_p2 = false;
        
        ImVec2 p1 = data_to_normalized(ImVec2(g_state[i].first, g_state[i].second));
        ImVec2 p2 = data_to_normalized(ImVec2(g_state[i+1].first, g_state[i+1].second));
        
        if(p1.x < 0.f || p1.x > 1.f || p1.y < 0.f || p1.y > 1.f) {
            clip_p1 = true;
        }
        
        if(p2.x < 0.f || p2.x > 1.f || p2.y < 0.f || p2.y > 1.f) {
            clip_p2 = true;
        }
        
        p1.x = ImMin(ImMax(p1.x,0.f),1.f);
        p1.y = ImMin(ImMax(p1.y,0.f),1.f);
        p2.x = ImMin(ImMax(p2.x,0.f),1.f);
        p2.y = ImMin(ImMax(p2.y,0.f),1.f);
        
        if(!clip_p1 || !clip_p2) {
            DrawList->AddLine(
                              normalized_to_pix(p1),
                              normalized_to_pix(p2),
                              GetColorU32(ImGuiCol_ColumnActive), 2);
        }
        //std::cout<<data_to_normalized(ImVec2(g_state[i].first, g_state[i].second)).x<<"\n";
    }
    
    //label axes
    
    return true;
    
}





void Visualize::setup(const Util::MarkerParticle &marker_particles, bool ps_plot) {
    
    g_marker_particles = &marker_particles;
    
    //add new menu for phase space plotting
    Visualize::g_viewer.plugins.push_back(&menu);
    
    menu.callback_draw_viewer_menu = [&](){};
    
    Visualize::g_viewer.callback_mouse_down = mouse_down;
    Visualize::g_viewer.callback_mouse_up = mouse_up;
    Visualize::g_viewer.callback_mouse_move = mouse_move;
    
    Visualize::g_viewer.core().background_color.setConstant(1.0);
    Visualize::g_viewer.core().is_animating = true;
}




void Visualize::update_particle_positions(Eigen::MatrixXd pos) {
    //g_viewer.data().clear_points();
    g_viewer.data().set_points(pos,Eigen::RowVector3d(0,0,1));
    g_viewer.data().point_size = 3.;
    
}

const std::vector<unsigned int> & Visualize::picked_vertices() {
    return g_picked_vertices;
}

void Visualize::set_picking_tolerance(double r) {
    g_picking_tol = r;
}

bool Visualize::mouse_down(igl::opengl::glfw::Viewer &viewer, int x, int y) {
    
    g_mouse_win = Eigen::Vector3d(g_viewer.current_mouse_x,viewer.core().viewport(3) - g_viewer.current_mouse_y,0.);
    igl::unproject(
                   g_mouse_win,
                   g_viewer.core().view,
                   g_viewer.core().proj,
                   g_viewer.core().viewport,
                   g_mouse_world);
    
    //if you click on the mesh select the vertex, otherwise do nothing
    
    
    return false;
}

bool Visualize::mouse_up(igl::opengl::glfw::Viewer &viewer, int x, int y) {
    
    g_mouse_dragging = false;
    g_picked_vertices.clear();
    g_mouse_drag_world.setZero();
    return false;
}

const Eigen::Vector3d & Visualize::mouse_world() {
    return g_mouse_world;
}

const Eigen::Vector3d & Visualize::mouse_drag_world() {
    return g_mouse_drag_world;
}

bool Visualize::mouse_move(igl::opengl::glfw::Viewer &viewer, int x, int y) {
    
    g_mouse_drag = Eigen::Vector3d(g_viewer.current_mouse_x,viewer.core().viewport(3) - g_viewer.current_mouse_y,0.) - g_mouse_win;
    g_mouse_win = Eigen::Vector3d(g_viewer.current_mouse_x,viewer.core().viewport(3) - g_viewer.current_mouse_y,0.);
    
    igl::unproject(
                   g_mouse_win,
                   g_viewer.core().view,
                   g_viewer.core().proj,
                   g_viewer.core().viewport,
                   g_mouse_drag_world);
    
    
    g_mouse_drag_world -= g_mouse_world;
    
    //std::cout<<"Test: "<<g_mouse_drag_world.transpose()<<"\n";
    igl::unproject(
                   g_mouse_win,
                   g_viewer.core().view,
                   g_viewer.core().proj,
                   g_viewer.core().viewport,
                   g_mouse_world);
    
    
    if(g_mouse_dragging && g_picked_vertices.size() > 0 ) {
        return true;
    }
    
    return false;
}

bool Visualize::is_mouse_dragging() {
    return g_mouse_dragging;
}


