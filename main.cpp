#include <tgmath.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "easytime.h"
#include "ImageUtil.h"

#define EPS 0.000001
template<typename Float> 
struct Particle {
    Float x,y,vx,vy,rho,ax,ay;
    uint idx,idy;
};

enum BoundaryType{
    CLAMP,
    PERIODIC,
    SPRING
};
enum EOSType{
    LINEAR,
    IDEALGAS,
    FLUID
};

template <typename Float>
class SPHSim {

public:
    std::vector<Particle<Float> > p;
    std::vector<std::vector<uint> > idarr;

    Float domainW;
    Float domainH;
    uint numCellsX;
    uint numCellsY;

    //each particle is radius pr
    //In SPH terms, I take h=2*pr. 
    Float pr; 
    Float pm; //particle mass

    Float pressureConst; //constant in p(rho).
    Float totalEnergy;



    inline Float eosLinear(Float arg);
    inline Float eosIdealGas(Float arg);
    inline Float intEnergyIdealGas(Float arg);
    inline Float eosFluid(Float arg);
    inline bool valididQ(int x,int y);

    void recalculateIDs(); 
    void stepPositionsAndVelocities(Float dt);
    void calculatePressures();
    Float pressureAtPoint(Float x,Float y);
    void calculateForces();

    SPHSim &operator=(const SPHSim&) = delete;
    SPHSim(const SPHSim&) = delete;
    SPHSim();
    Float timestep(Float dt);
};

int main() {
    SPHSim<float> sph;

    sph.domainW=1.0;
    sph.domainH=1.0;
    sph.pm=1.0;
    sph.pr=0.0005;

    sph.numCellsX=std::ceil(sph.domainW/(2*sph.pr));
    sph.numCellsY=std::ceil(sph.domainH/(2*sph.pr));
    
    int nparticles=1000000;
    srand(0);
    for(int n=0;n<nparticles;n++){
        float x=(rand()*1.0f)/RAND_MAX,y=(2.5f/3.0f*rand())/RAND_MAX;
        sph.p.push_back(Particle<float>{x,y,0.f,0.f,0.f,0.f,0.f,0,0});
        //std::cout<<x<<std::endl;
    }
    
    sph.recalculateIDs();

    for(int a=0;a<3;a++){
        for(int b=0;b<3;b++){
            std::cout<<sph.idarr[a*3+b].size()<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<sph.pressureAtPoint(0,0)<<std::endl;
    DoubleImage dimg(400,400);
    for(int a=0;a<400;a++){
        for(int b=0;b<400;b++){
            dimg.put(a,b,sph.pressureAtPoint(a*1.0/400.0*0.01,b*1.0/400.0*0.01));
        }
    }
    std::cout<<"Done, saving."<<std::endl;
    dimg.unitStretch();
    Image outimg(dimg.getData(),400,400);
    outimg.save("pressure.bmp");
    std::cout<<"Done."<<std::endl;
    std::cout<<sizeof(int)<<std::endl;

    return 0;
}


template <typename Float>
SPHSim<Float>::SPHSim() {

}

template <typename Float>
inline Float SPHSim<Float>::eosLinear(Float arg) {
    return 0;
}
template <typename Float>
inline Float SPHSim<Float>::eosIdealGas(Float arg) {
    return Float(2.0/3.0)*pressureConst*std::pow(arg,Float(5.0/3.0));
}
template <typename Float>
inline Float SPHSim<Float>::intEnergyIdealGas(Float arg) {
    return pressureConst*std::pow(arg,Float(2.0/3.0));
}
template <typename Float>
inline Float eosFluid(Float arg) {
    return 0;
}

template <typename Float>
inline bool SPHSim<Float>::valididQ(int x,int y){
    return x>=0 && x<numCellsX && y>=0 && y<numCellsY;
}

template <typename Float>
void SPHSim<Float>::recalculateIDs() {
    idarr=std::vector<std::vector<uint> >(numCellsX*numCellsY);
    Float boxWidth=domainW/numCellsX;
    Float boxHeight=domainH/numCellsY;
    for(uint n=0;n<(uint)p.size();n++) {
        if(p[n].x<0)
            p[n].x=0;
        if(p[n].x>=domainW)
            p[n].x=domainW-Float(EPS);
        if(p[n].y<0)
            p[n].y=0;
        if(p[n].y>=domainH)
            p[n].y=domainH-Float(EPS);
        p[n].idx=std::floor(p[n].x/boxWidth);
        p[n].idy=std::floor(p[n].y/boxHeight);
        idarr[p[n].idy*numCellsX+p[n].idx].push_back(n);
    }
}

template <typename Float>
Float SPHSim<Float>::pressureAtPoint(Float x, Float y) {
    const Float kernelConst=5.0f/(M_PI*(2*pr)*(2*pr));
    Float boxWidth=domainW/numCellsX;
    Float boxHeight=domainH/numCellsY;
    int cx=std::floor(x/boxWidth);
    int cy=std::floor(y/boxHeight);
    Float ret=0;
    for(int a=-1;a<=1;a++){
        for(int b=-1;b<=1;b++){
            int midx=cx+b;
            int midy=cy+a;
            if(valididQ(midx,midy)){
                for(size_t k=0;k<idarr[midy*numCellsX+midx].size();k++){
                    size_t m=idarr[midy*numCellsX+midx][k];

        Float dx=p[m].x-x;
        Float dy=p[m].y-y;
        Float d=std::sqrt(dx*dx+dy*dy);
        Float R=d/(2*pr);
        
        if(R<Float(1.0)){
            Float wab=std::sqrt(1-R*R);//(1+3*R)*(1-R)*(1-R)*(1-R);
            ret+=pm*wab;
        }
                    
                }
            }
        }
    }
    return ret;
}

template <typename Float>
void SPHSim<Float>::calculatePressures() {
    totalEnergy=0;
    const Float kernelConst=5.0f/(M_PI*(2*pr)*(2*pr));
    for(size_t n=0;n<p.size();n++){
        p[n].rho=0;
        for(int a=-1;a<=1;a++){
            for(int b=-1;b<=1;b++){
                int midx=p[n].idx+b;
                int midy=p[n].idy+a;
                if(valididQ(midx,midy)){
                    for(size_t k=0;k<idarr[midy*numCellsX+midx].size();k++){
                        size_t m=idarr[midy*numCellsX+midx][k];

            Float dx=p[m].x-p[n].x;
            Float dy=p[m].y-p[n].y;
            Float d=std::sqrt(dx*dx+dy*dy);
            Float R=d/(2*pr);
            
            if(d<2*pr){
                Float wab=(1+3*R)*(1-R)*(1-R)*(1-R);
                p[n].rho+=pm*wab;
            }
                        
                    }
                }
            }
        }
    }
}
template <typename Float>
void SPHSim<Float>::stepPositionsAndVelocities(Float dt) {

}

template <typename Float>
void SPHSim<Float>::calculateForces() {
}

