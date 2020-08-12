#include <pick_nearest_vertices.h>
#include <iostream>
#include <igl/Hit.h>
#include <igl/ray_mesh_intersect.h>
#include <igl/unproject.h>
#include <igl/unproject_ray.h>

bool pick_nearest_vertices(std::vector<unsigned int> &verts, Eigen::Ref<const Eigen::Vector3d> win,
                           Eigen::Ref<const Eigen::Matrix44f> view, Eigen::Ref<const Eigen::Matrix44f> proj, Eigen::Vector4f viewport,
                           Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, double radius) {
    
    verts.clear();
    
    // Source, destination and direction in world
    Eigen::Vector3f start,dir;
    igl::Hit hit;
    
    //compute start and direction in the world to check for picked vertex
    //YOUR CODE HERE
    Eigen::Vector3f des;
    Eigen::Vector3f win_s = Eigen::Vector3f(win[0], win[1], -1.);
    igl::unproject(win.cast<float>(),view,proj,viewport,start);
    igl::unproject(win_s.cast<float>(),view,proj,viewport,des);
    dir = start-des;
    
    const auto & shoot_ray = [&V,&F](const Eigen::Vector3f& s, const Eigen::Vector3f& dir, igl::Hit & hit)->bool
    {
        std::vector<igl::Hit> hits;
        
        if(!igl::ray_mesh_intersect(s,dir,V,F,hits))
        {
            return false;
        }
        hit = hits[0];
        return true;
    };
    
    if(!shoot_ray(start,dir,hit))
    {
        return false;
    }
    
    //check if any of the hit vertices are within the selection radius
    //YOUR CODE HERE
    
    Eigen::Vector3d q;
    if (hit.u >= hit.v){
        if (hit.u >= (1-hit.u-hit.v)){
            q = V.row(F(hit.id, 1));
        }else{
            q = V.row(F(hit.id, 0));
        }
    }else{
        if (hit.v >= (1-hit.u-hit.v)){
            q = V.row(F(hit.id, 2));
        }else{
            q = V.row(F(hit.id, 0));
        }
    }

    for (int i = 0; i < V.rows(); i++){
        Eigen::Vector3d q0 = V.row(i);
        double d = (q.cast<double>() - q0).norm();
        if (d <= radius){
	    std::cout<<i<<"\n\n\n\n\n";
            verts.push_back(i);
        }
    }
    return (verts.size() == 0 ? false : true);
}
