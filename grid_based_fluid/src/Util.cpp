#include "Util.h"
#include <iostream>
#include <stdlib.h>

namespace Util{
    // mac grid class implementation
    MacGrid::MacGrid(){
        
    }
    
    MacGrid::~MacGrid(){
        
    }
    
    void MacGrid::init(double l, double w, double h, double grid_size){
        g_cell_size = grid_size;
        min_vert = Eigen::Vector3d(-l, -h, -w);
        max_vert = Eigen::Vector3d(l, h, w);
        
        n_cellx = 2*l / grid_size;  // # of cells
        n_celly = 4*h / grid_size;
        n_cellz = 2*w / grid_size;

        p.resize(n_cellx * n_celly * n_cellz);
        vx.resize((n_cellx + 1) * n_celly * n_cellz);
        vy.resize(n_cellx * (n_celly + 1) * n_cellz);
        vz.resize(n_cellx * n_celly * (n_cellz + 1));
        p.setZero();
        vx.setZero();
        vy.setZero();
        vz.setZero();
        // label: 0 air; 1 solid; 2 fluid
        label.resize(n_cellx * n_celly * n_cellz);
        label.setZero();
        for (int i = 0; i < n_cellx; i++){
            for (int j = 0; j < n_celly; j++){
                for (int k = 0; k < n_cellz; k++){
                    if (i == 0 || i == n_cellx - 1 || j == 0 || k == 0 || k == n_cellz - 1){
                        // set boundary as solid
                        label[(n_celly * n_cellz) * i + n_cellz * j + k] = 1;
                    }
                }
            }
        }
    }

    // marker particle class implementation
    MarkerParticle::MarkerParticle(){
        
    }
    
    MarkerParticle::~MarkerParticle(){
        
    }
    
    void MarkerParticle::init(MacGrid &mg, std::vector<Eigen::Vector3i> cells){
        pos.resize(8 * cells.size(), 3);
        pos.setZero();
        vel.clear();
        gridIndex.clear();
        double dx = mg.g_cell_size;
        int n_x = mg.n_cellx;
        int n_y = mg.n_celly;
        int n_z = mg.n_cellz;
        Eigen::Vector3d min_vert = mg.min_vert;
        
        for (int n = 0; n < cells.size(); n++){
            int a = cells[n][0];
            int b = cells[n][1];
            int c = cells[n][2];
            double pxl = (a - 1) * dx + dx / 4.;
            double pxr = a * dx - dx / 4.;
            double pyb = (b - 1) * dx + dx / 4.;
            double pyf = b * dx - dx / 4.;
            double pzb = (c - 1) * dx + dx / 4.;
            double pzt = c * dx - dx / 4.;
            
            pos.row(8*n) = min_vert + Eigen::Vector3d(pxl, pyb, pzt);
            pos.row(8*n+1) = min_vert + Eigen::Vector3d(pxr, pyb, pzt);
            pos.row(8*n+2) = min_vert + Eigen::Vector3d(pxl, pyf, pzt);
            pos.row(8*n+3) = min_vert + Eigen::Vector3d(pxr, pyf, pzt);
            pos.row(8*n+4) = min_vert + Eigen::Vector3d(pxl, pyb, pzb);
            pos.row(8*n+5) = min_vert + Eigen::Vector3d(pxr, pyb, pzb);
            pos.row(8*n+6) = min_vert + Eigen::Vector3d(pxl, pyf, pzb);
            pos.row(8*n+7) = min_vert + Eigen::Vector3d(pxr, pyf, pzb);
            for (int i = 0; i < 8; i++){
                // jitter=0.02  dx/4=0.025
                double jitterX = ((rand() % 41) - 20) / 100. * dx;
                double jitterY = ((rand() % 41) - 20) / 100. * dx;
                double jitterZ = ((rand() % 41) - 20) / 100. * dx;
                pos.row(8*n+i) += Eigen::Vector3d(jitterX, jitterY, jitterZ);
                
                gridIndex.push_back((n_z*n_y)*a+n_z*b+c);
                
                vel.push_back(Eigen::Vector3d(0., 0., 0.));
            }
        }
    }
    
}
