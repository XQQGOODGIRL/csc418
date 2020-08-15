## Fluid Simulation

## Implemented
3D Mac grid
Semi-Lagrangian advection algorithmd
Marker particle visualization
Using pressure solver from Bridson and Müller-Fischer
![Ellipsoid fluid example](data/ellipsoid_ex.jpg)

## TODO
Other advection algorithm
Level set and Signed distance 
Marching Cube 

## To Run
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
./fluid
```

