## Fluid Simulation

## Implemented
* 3D Mac grid
* Semi-Lagrangian advection algorithmd
* Marker particle visualization
* Using pressure solver from Bridson and Müller-Fischer 

An example of an ellipsoid dropping into water
<img src="https://github.com/XQQquxixi/graphics/blob/master/grid_based_fluid/data/0.png" width="250" height="200"> <img src="https://github.com/XQQquxixi/graphics/blob/master/grid_based_fluid/data/1.png" width="250" height="200"> <img src="https://github.com/XQQquxixi/graphics/blob/master/grid_based_fluid/data/final.png" width="250" height="180"> 


## TODO
* Make it faaster and prettier and more interactive!!
* Other advection algorithm
* Level set and Signed distance 
* A fluid bunny
* Curved solid boudary (other weird collision)
* Marching Cube 

## To Run
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
./fluid
```
Type S key to start.

