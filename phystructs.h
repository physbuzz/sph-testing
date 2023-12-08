#ifndef PHYSTRUCTS_H
#define PHYSTRUCTS_H
/** PointMass
dumb structure that stores the data of a gravitating particle.
*/

template<typename Float> 
struct Particle {
    Float x,y,vx,vy,rho,ax,ay;
    uint idx,idy;
};

#endif
