#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "VectorND.h"
#include "phystructs.h"
#include "ParticleList.h"
#include "PGrid.h"
#include "ImageUtil.h"
#include "easytime.h"

#define EPS 0.000001

//Some string manipulation functions for saving files. pad_int(1234,5) returns "01234".
std::string pad_int(int arg, int padcount) {
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(padcount) << arg;
    return ss.str();
}

//Returns a file name in the form of "prefix00###suffix". For example "image0032.bmp"
std::string getFilename(std::string prefix, int num, int padcount, std::string suffix) {
    return prefix + pad_int(num, padcount) + suffix;
}

struct ImageParams {
    int imgw;
    int imgh;
    float realsize;
    float cx;
    float cy;
};


//Ideal gas internal energy U(rho)
template <typename Float>
inline Float intEnergyIdealGas(Float density) {
    return std::pow(density,Float(2.0/3.0));
}

//Ideal gas EOS P(rho)=(rho**2)*U'(rho)
template <typename Float>
inline Float eosIdealGas(Float density) {
    return Float(2.0/3.0)*std::pow(density,Float(5.0/3.0));
}

// Smoothing kernel Wab
template <typename Float, int DIM>
inline Float smoothingKernelQuartic(Float h, Float r) {
    if(DIM==2)
        return Float(5)/(M_PI*h*h)*(1+3*r)*(1-r)*(1-r)*(1-r);
    else
        return 0;
}

template <typename Float, int DIM>
inline VectorND<Float,DIM> smoothingKernelQuarticGradient(Float h, Float r);

// Smoothing kernel Wab
//template <typename Float>
//inline Float smoothingKernelQuarticGradient<Float,2>(Float h, Float r) {
//    return Float(5)/(M_PI*h*h)*(1+3*r)*(1-r)*(1-r)*(1-r);
//}



template<typename Float, int DIM>
class SPHSim {
public:
    //Grid structure for pressure calculations
    PGrid<Float,DIM> pgrid;
    //Maximum radius of any of the particles
    Float maxH;

    SPHSim(ParticleList<Float,DIM> &pl, VectorND<Float,DIM> domainSize, Float maxH) :
        pgrid(&pl.plist,domainSize,maxH),
        maxH(maxH)
    { }

    Float densityAtPoint(VectorND<Float,DIM> pos) {
        Float rho=0;
        for(auto p2 : pgrid.nearbyLoop(pos,maxH)){

            VectorND<Float,DIM> deltapos=p2->pos - pos;
            Float d=deltapos.length()/p2->h;
            if(d<1.0){
                Float wab=Float(5)/(M_PI*p2->h*p2->h)*(1+3*d)*(1-d)*(1-d)*(1-d);
                rho+=p2->m*wab;
            }
        }
        return rho;
    }

    /*
    void updateOnce(Float radius,Float dt) {
        //Update the density
        for(Particle<Float,DIM> &p : *s.plist){
            p.acc=VectorND<Float,DIM>();
            p.rho=densityAtPoint(p.pos);
        }

        //Next, 
        Float radiusDIM=radius*radius;

        //TODO: make this generic
        for(Particle<Float,DIM> *p1 : s.updateLoop()){
            if(p1->collision<0) {
                for(Particle<Float,DIM> *pDIM : s.nearbyLoop(p1->pos,maxH)){
                    assert(pDIM!=nullptr);
                    if(p1==pDIM || pDIM->collision==0)
                        continue;
                    auto x1=p1->pos;
                    auto xDIM=pDIM->pos;
                    auto v1=p1->vel;
                    auto vDIM=pDIM->vel;
                    auto dx=xDIM-x1;
                    auto dv=vDIM-v1;
                    Float inner=dx.dot(dv);
                    if(inner>=0)
                        continue;
                    Float dvDIM=dv.lengthDIM();
                    Float d=inner*inner-dvDIM*(dx.lengthDIM()-4.0f*radiusDIM);
                    if(d<=0)
                        continue;
                    Float t1=(-inner-sqrt(d))/dvDIM;
                    if(t1<EPSDIM || t1>dt)
                        continue;
                    //If we get to this point, there's a collision within time dt.
                    //If the other particle has already undergone a collision, also ignore
                    //but set a flag for it.
                    //
                    //Basically, we only ever handle collisions for two pairs of particles
                    //with their collision flags set to -1.
                    if(pDIM->collision>0){
                        //This is not a perfect count of the number of double collisions. 
                        //Need to double check that it has some bearing to ground truth!
                        stats.nTwoOrMore+=1;
                        stats.nOne-=1;
                        continue;
                    }

                    Float tDIM=dt-t1;
                    p1->collision=1;
                    pDIM->collision=1;
                    stats.nOne+=DIM;

                    p1->posnew=p1->pos+p1->vel*t1;
                    pDIM->posnew=pDIM->pos+pDIM->vel*t1;
                    //Collision of equal masses; we reverse each relative velocity along the direction
                    //of their collision vector.
                    //
                    //Start with velocity v1. Put it in CM frame:
                    //v1_CM=v1-(v1+vDIM)/DIM
                    //reverse it along the rhat direction:
                    //v1_CM -> v1_CM-DIM*rhat(rhat.dot(v1_CM))
                    //add (v1+vDIM)/DIM to put it back in non-CM and simplify:
                    //v1_new = v1+DIM*rhat*rhat.dot((vDIM-v1)/DIM)
                    //       = v1+rhat*rhat.dot(vDIM-v1)
                    //cout<<(pDIM->posnew-p1->posnew).length()<<endl;
                    auto dxDIM=(pDIM->posnew-p1->posnew).normalized();
                    p1->velnew=p1->vel+dxDIM*dxDIM.dot(dv);
                    pDIM->velnew=pDIM->vel-dxDIM*dxDIM.dot(dv);

                    //time evolve the rest of the way.
                    p1->posnew+=p1->velnew*tDIM;
                    pDIM->posnew+=pDIM->velnew*tDIM;
                    break;
                }
            }

            if(p1->collision<0){
                p1->posnew=p1->pos+dt*p1->vel;
                p1->velnew=p1->vel;
                p1->collision=0;
                stats.nZero+=1;
            }
            p1->pos=p1->posnew;
            p1->vel=p1->velnew;
            if(p1->posnew.x[0]<0){
                p1->vel.x[0]=-p1->vel.x[0];
                p1->pos.x[0]=-p1->pos.x[0];
            }
            if(p1->posnew.x[0]>s.domainSize[0]){
                p1->vel.x[0]=-p1->vel.x[0];
                p1->pos.x[0]=DIM*s.domainSize[0]-p1->pos.x[0];
            }
            if(p1->posnew.x[1]<0){
                p1->vel.x[1]=-p1->vel.x[1];
                p1->pos.x[1]=-p1->pos.x[1];
            }
            if(p1->posnew.x[1]>s.domainSize[1]){
                p1->vel.x[1]=-p1->vel.x[1];
                p1->pos.x[1]=DIM*s.domainSize[1]-p1->pos.x[1];
            }
        } 
    }
    //Return physical values after sampling particles within rmax of pos.
    //c is the calculated value from radialwc(rmax,p)
    PhysicsQueryStruct querySimulation(VectorND<Float,DIM> pos,Float rmax, Float p,Float c=-1.0f){
        if(c<=0.0f)
            c=radialwc(rmax,p);

        PhysicsQueryStruct ret{0.0f,0.0f,0.0f,0.0f,0.0f,0.0f};

        for(Particle<Float,DIM> *pDIM : s.nearbyLoop(pos,rmax)){
            Float r=(pDIM->pos-pos).length();
            Float w=radialw(r,p,rmax);
            ret.n+=w;
            ret.px+=pDIM->vel.x[0]*w;
            ret.py+=pDIM->vel.x[1]*w;
            ret.e+=0.5f*pDIM->vel.lengthDIM()*w;
        }
        //If we didn't pick up any particles, just return zero.
        if(ret.n<=EPSDIM)
            return ret;
        //calculate the expected momentum and expected energy
        ret.px/=ret.n;
        ret.py/=ret.n;
        ret.e/=ret.n;
        Float h=0.0f;
        int nparticles=0;
        for(Particle<Float,DIM> *pDIM : s.nearbyLoop(pos,rmax)){
            Float r=(pDIM->pos-pos).length();
            Float w=radialw(r,p,rmax);
            if(w>0){
                nparticles++;
            }
            Float p1x=pDIM->vel.x[0]-ret.px;
            Float p1y=pDIM->vel.x[1]-ret.py;
            h+=0.5f*(p1x*p1x+p1y*p1y)*w;
        }
        ret.beta=ret.n/h;
        ret.n/=c;
        Float z11=100000.0f;
        ret.s=ret.n*(DIM.0f-log(ret.n*ret.beta/z11));

        if(nparticles<=1){
            //Can't estimate beta and s if there's only one particle
            ret.beta=0.0f;
            ret.s=0.0f;
        }
        return ret;
    }*/
    void saveImage( ImageParams ip, std::string prefix, int fnamei, int padcount ){
        int imw=ip.imgw;
        int imh=ip.imgh;

        Image outimg(imw,imh);
        Float realsize=ip.realsize;
        Float cx=ip.cx;
        Float cy=ip.cy;
        Float aspect=Float(imh)/imw;

        for(int a=0;a<imw;a++){
            for(int b=0;b<imh;b++){
                Float x=cx+(Float(a)/imw-0.5f)*realsize;
                Float y=cy+(Float(b)/imh-0.5f)*realsize*aspect;
                VectorND<Float,2> pos({x,y});

                //Draw gridlines
                //Float x2=cx+(Float(a+1)/imw-0.5f)*realsize;
                //Float y2=cy+(Float(b+1)/imh-0.5f)*realsize*aspect;
                //VectorND<Float,2> pos2({x2,y2});
                //if(!(pgrid.positionToIntvec(pos)==pgrid.positionToIntvec(pos2))) {
                //    outimg.put(a,b,intToRGB(255,255,255));
                //    continue;
                //}
                Float rho=densityAtPoint(pos);
                //Float rho=0;
                //auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                outimg.put(a,b,intToRGB(rho*255,rho*255,rho*255));
            }
        }
        //cout<<"Checked "<<(Float(nParticlesChecked)/(imw*imh))<<" particles per pixel"<<endl;
        outimg.save(getFilename(prefix,fnamei,padcount,".bmp"));
    }
};

using namespace std;


int main() {
    ImageParams ip;
    ip.imgw=640;
    ip.imgh=480;
    ip.realsize=2;
    ip.cx=1.0;
    ip.cy=1.0;

    float L=2.0f;
    //Expected velocities are sqrt(2T/m)
    //time to cross a boundary ~= dx/sqrt(2T/m)
    int nparticles=10000;

    //The average number of overlapping circles at any given point = expectedN
    // expectedN == nparticles*(M_PI*r*r)/(L*L);
    float expectedN=2.5;
    // expectedRho = nparticles * mass/L*L;
    float expectedRho=0.5;

    //Invert the two above expressions to find the particle sizes
    float particleR=std::sqrt(expectedN*L*L/(nparticles*M_PI));
    float particleM=expectedRho*L*L/nparticles;


    cout<<particleR<<" "<<particleM<<endl;


    VectorND<float,2> domainSize({L,L});

    float maxH=2*particleR;
    //float dt=maxH/(6.0f*std::sqrt(2.0f*temperature));
    //float timeelapsed=0.0f;

    ParticleList<float,2> particleList;

    for ( int i = 0; i < nparticles; i++ ) 
    {
        //float vmag=std::sqrt(2*temperature);
        float vmag=0.0;
        float theta=(rand()*2.0f*M_PI)/RAND_MAX;
        VectorND<float,2> pnew({(rand()*L)/RAND_MAX,(rand()*L)/RAND_MAX});
        VectorND<float,2> vnew({vmag*std::cos(theta),vmag*std::sin(theta)});
        particleList.plist.push_back(Particle<float,2>{pnew,vnew,VectorND<float,2>(),0, particleR,particleM});
    }
    SPHSim<float,2> sph(particleList,domainSize,maxH);

    sph.pgrid.rebuildGrid();
    
    sph.saveImage(ip,"test",0,1);
//    SPHSim<float> sph;
//
//    sph.domainW=1.0;
//    sph.domainH=1.0;
//    sph.pm=1.0;
//    sph.pr=0.0005;
//
//    sph.numCellsX=std::ceil(sph.domainW/(2*sph.pr));
//    sph.numCellsY=std::ceil(sph.domainH/(2*sph.pr));
//    
//    int nparticles=1000000;
//    srand(0);
//    for(int n=0;n<nparticles;n++) {
//        float x=(rand()*1.0f)/RAND_MAX,y=(2.5f/3.0f*rand())/RAND_MAX;
//        sph.p.push_back(Particle<float>{x,y,0.f,0.f,0.f,0.f,0.f,0,0});
//        //std::cout<<x<<std::endl;
//    }
//    
//    sph.recalculateIDs();
//
//    for(int a=0;a<3;a++){
//        for(int b=0;b<3;b++){
//            std::cout<<sph.idarr[a*3+b].size()<<" ";
//        }
//        std::cout<<std::endl;
//    }
//
//    std::cout<<sph.pressureAtPoint(0,0)<<std::endl;
//    DoubleImage dimg(400,400);
//    for(int a=0;a<400;a++){
//        for(int b=0;b<400;b++){
//            dimg.put(a,b,sph.pressureAtPoint(a*1.0/400.0*0.01,b*1.0/400.0*0.01));
//        }
//    }
//    std::cout<<"Done, saving."<<std::endl;
//    dimg.unitStretch();
//    Image outimg(dimg.getData(),400,400);
//    outimg.save("pressure.bmp");
//    std::cout<<"Done."<<std::endl;
//    std::cout<<sizeof(int)<<std::endl;
//
    return 0;
}



