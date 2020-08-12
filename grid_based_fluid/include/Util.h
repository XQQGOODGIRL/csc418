#ifndef Util_h
#define Util_h
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <vector>

namespace Util{
    class MacGrid
    {
    public:
        MacGrid();
        //MacGrid(double x, double y, double z);
        ~MacGrid();
        int n_cellx;
        int n_celly;
        int n_cellz;
        double g_cell_size;
        double g_x;
        double g_z;
        double g_y;
        Eigen::Vector3d min_vert;
        Eigen::Vector3d max_vert;
        Eigen::VectorXd vx;
        Eigen::VectorXd vy;
        Eigen::VectorXd vz;
        Eigen::VectorXd p;
        Eigen::VectorXi label;
        void init(double l, double w, double h, double grid_size);
    private:
        //void resize(double x, double y, double z);
        
        
    };
    
    class MarkerParticle
    {
    public:
        MarkerParticle();
        //MarkerParticle(Eigen::Vector3d pos, Eigen::Vector3d vel);
        ~MarkerParticle();
        Eigen::MatrixXd pos;
        std::vector<Eigen::Vector3d> vel;
        std::vector<int> gridIndex;
        void init(MacGrid &mg, std::vector<Eigen::Vector3i> indices);
        
    };

};

#endif /* Util_h */
