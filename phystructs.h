#ifndef PHYSTRUCTS_H
#define PHYSTRUCTS_H

#include "VectorND.h"

template<typename Float, int DIM> 
struct Particle {
    VectorND<Float,DIM> pos;
    VectorND<Float,DIM> vel;
    VectorND<Float,DIM> acc;
    Float rho, h, m;
};

#endif
