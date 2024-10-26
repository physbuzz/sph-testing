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

template<typename Float, int DIM>
struct ImageParamsND {
    int imgw, imgh;
    Float realsize;
    VectorND<Float,DIM> center;
    VectorND<Float,DIM> xhat;
    VectorND<Float,DIM> yhat;
    ImageParamsND() : imgw(640), imgh(480), realsize(1), center(), xhat(), yhat() {
        xhat[0]=1;
        yhat[1]=1;
    }
};

//Ideal gas internal energy U(rho)
template <typename Float>
static inline Float intEnergyIdealGas(Float density) {
    return std::pow(density,Float(2.0/3.0));
}

//Ideal gas EOS P(rho)=(rho**2)*U'(rho)
template <typename Float>
static inline Float eosIdealGas(Float density) {
    return Float(2.0/3.0)*std::pow(density,Float(5.0/3.0));
}

template<typename Float, int DIM>
static inline constexpr Float smoothingKernelQuarticConstant() {
    static_assert(0<DIM && DIM<=5, "Generic DIM not implemented");
    if constexpr (DIM == 1){
        return Float(5)/4;
    } else if constexpr (DIM == 2){
        return Float(5)/(M_PI);
    } else if constexpr (DIM == 3){
        return Float(105)/(16*M_PI);
    } else if constexpr (DIM == 4){
        return Float(28)/(M_PI*M_PI);
    } else if constexpr (DIM == 5){
        return Float(315)/(8*M_PI*M_PI);
    }
    //The full normalization constant (for h=1) is 
    //2  Pi^(n/2)/Gamma[n/2]  Integrate[(1 + 3 r) (1 - r)^3  r^(n - 1), {r, 0, 1}]
    //= (48 Pi^(n/2))/((24 n + 26 n^2 + 9 n^3 + n^4) Gamma[n/2])
    return 0;
}
template <typename Float, int exponent>
static inline Float templatePower(Float arg) {
    static_assert(exponent>0,"Exponent must be positive");
    if constexpr (exponent==1) {
        return arg;
    } else {
        return arg*templatePower<Float,exponent-1>(arg);
    }
}

template <typename Float, int DIM>
static inline Float smoothingKernelQuartic(Float h, Float r) { 
    r=r/h;
    return smoothingKernelQuarticConstant<Float,DIM>()/templatePower<Float,DIM>(h)*(1+3*r)*(1-r)*(1-r)*(1-r);
}

template <typename Float, int DIM>
static inline VectorND<Float,DIM> smoothingKernelQuartic(Float h, VectorND<Float, DIM> rvec) { 
    float r=rvec.length();
    r=r/h;
    return smoothingKernelQuarticConstant<Float,DIM>()*(-12)/templatePower<Float,DIM>(h)*(1-r)*(1-r)*(rvec/r);
}

template <typename Float, int DIM>
inline VectorND<Float,DIM> smoothingKernelQuarticGradient(Float h, Float r) { 
}

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

    Float gravity;

    SPHSim(ParticleList<Float,DIM> &pl, VectorND<Float,DIM> domainSize, Float maxH) :
        pgrid(&pl.plist,domainSize,maxH),
        maxH(maxH), gravity(0)
    { 
        recalculateDensities();
    }


    void recalculateDensities(){ 
        pgrid.rebuildGrid(); 
        for(auto &p : *pgrid.plist){
            p.rho=densityAtPoint(p.pos);
        }
    }

    Float densityAtPoint(VectorND<Float,DIM> pos) {
        Float rho=0;
        for(auto p2 : pgrid.nearbyLoop(pos,maxH)){

            VectorND<Float,DIM> deltapos=p2->pos - pos;
            Float d=deltapos.length();
            if(d<p2->h){
                Float wab=smoothingKernelQuartic<Float,DIM>(p2->h,d);
                rho+=p2->m*wab;
            }
        }
        return rho;
    }

    // Precondition: we need to have already calculated p1->rho
    void annealPositionsBadly(Float temperature, Float stepsize){
        for(Particle<Float,DIM> *p1 : pgrid.updateLoop()){
            //The energy of the whole system is the sum of kinetic energies + mass[i]*intEnergyIdealGas(rho[i]).
            //(note that the internal energy function is really the energy per unit mass)
            auto proposed_position=p1->pos+VectorND<Float,DIM>::randomGaussian()*stepsize/std::sqrt(DIM);
            //Float current_energy=p1->m*(gravity*p1->pos[0] + intEnergyIdealGas(p1->rho));
            Float rho_new=0;
            for(auto p2 : pgrid.nearbyLoop(proposed_position,maxH)){
                if(p1==p2)
                    continue;
                VectorND<Float,DIM> deltapos=p2->pos - proposed_position;
                Float d=deltapos.length();
                if(d<p2->h){
                    Float wab=smoothingKernelQuartic<Float,DIM>(p2->h,d);
                    rho_new+=p2->m*wab;
                }
            }
            //Account for the self-contribution to density
            rho_new+=p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,0);
            //Float new_energy=p1->m*(gravity*proposed_position[0]+intEnergyIdealGas(rho_new));
            //std::cout<<"Rho_new = "<<rho_new<<", rho_old = "<<p1->rho<<std::endl;

            Float pAccept=std::exp(-(rho_new-p1->rho)/temperature);
            if( pAccept>1 || Float(rand())/RAND_MAX < pAccept){
                // Accept the move of the particle
                // The problem is that moving the particle changes rho for all nearby particles
                // So we have to first subtract off the density contribution from the old position,
                // then add back the new contribution.
                //
                // Each particle 2 has its density updated by p2->rho = p2 -> rho +
                //      - p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,(p2->pos-p1->pos).length())
                //      + p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,(p2->pos-proposed_position).length())
                for(auto p2 : pgrid.nearbyLoop(p1->pos,maxH)){
                    if(p1==p2)
                        continue;

                    VectorND<Float,DIM> deltapos=p2->pos - p1->pos;
                    Float d=deltapos.length();
                    if(d<p1->h)
                        p2->rho-=p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,d);
                }
                for(auto p2 : pgrid.nearbyLoop(proposed_position,maxH)){
                    if(p1==p2)
                        continue;

                    VectorND<Float,DIM> deltapos=p2->pos - proposed_position;
                    Float d=deltapos.length();
                    if(d<p1->h)
                        p2->rho+=p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,d);
                }
                p1->pos=proposed_position;
                p1->rho=rho_new;
            }
        }
    }
    void annealPositions(Float temperature, Float stepsize){
        for(Particle<Float,DIM> *p1 : pgrid.updateLoop()){
            //The energy of the whole system is the sum of kinetic energies + mass[i]*intEnergyIdealGas(rho[i]).
            //(note that the internal energy function is really the energy per unit mass)
            auto proposed_position=p1->pos+VectorND<Float,DIM>::randomGaussian()*stepsize/std::sqrt(DIM);
            Float current_energy=p1->m*(gravity*p1->pos[0] + intEnergyIdealGas(p1->rho));

            Float rho_new=0;
            for(auto p2 : pgrid.nearbyLoop(proposed_position,maxH)){
                if(p1==p2)
                    continue;
                VectorND<Float,DIM> deltapos=p2->pos - proposed_position;
                Float d=deltapos.length();
                if(d<p2->h){
                    Float wab=smoothingKernelQuartic<Float,DIM>(p2->h,d);
                    rho_new+=p2->m*wab;
                }
            }
            //Account for the self-contribution to density
            rho_new+=p1->m*smoothingKernelQuartic<Float,DIM>(p1->h,0);
            //Float new_energy=p1->m*(gravity*proposed_position[0]+intEnergyIdealGas(rho_new));
            //Float pAccept=std::exp(-(rho_new-p1->rho)/temperature);
            //Float pAccept=std::exp(-(new_energy-current_energy)/temperature);
            //if( pAccept>1 || Float(rand())/RAND_MAX < pAccept){
            //std::cout<<"Rho_new = "<<rho_new<<", rho_old = "<<p1->rho<<std::endl;

            if(rho_new<p1->rho) {
                p1->pos=proposed_position;
                p1->rho=rho_new;
            }
        }
    }

    //Specialized case for 1D
    Float densityAtPoint1D(VectorND<Float,2> pos) {
        static_assert(DIM==1, "Can only use densityAtPoint1D in one dimension");
        Float rho=0;
        for(auto p2 : pgrid.nearbyLoop(VectorND<Float,1>{pos[0]},maxH)){

            //VectorND<Float,DIM> deltapos=p2->pos - pos;
            //Float d=deltapos.length();
            Float d=std::sqrt((p2->pos[0]-pos[0])*(p2->pos[0]-pos[0])+pos[1]*pos[1]);
            if(d<p2->h){
                Float wab=smoothingKernelQuartic<Float,DIM>(p2->h,d);
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
};
template<typename Float,int DIM>
void saveSPHImage( SPHSim<Float,DIM> &sph, ImageParamsND<Float,DIM> ip, std::string prefix, int fnamei, int padcount,
        bool gridlines=false){
    Image outimg(ip.imgw,ip.imgh);
    Float aspect=Float(ip.imgh)/ip.imgw;

    if constexpr(DIM==1) {
        for(int a=0;a<ip.imgw;a++){
            for(int b=0;b<ip.imgh;b++){
                Float x=(Float(a)/ip.imgw-0.5f)*ip.realsize;
                Float y=(Float(b)/ip.imgh-0.5f)*ip.realsize*aspect;
                auto pos=ip.center+ip.xhat*x;

                /*if(gridlines){
                    Float x2=(Float(a+1)/ip.imgw-0.5f)*ip.realsize;
                    Float y2=(Float(b+1)/ip.imgh-0.5f)*ip.realsize*aspect;
                    auto pos2=ip.center+ip.xhat*x;
                    if(!(sph.pgrid.positionToIntvec(pos)==sph.pgrid.positionToIntvec(pos2))) {
                        outimg.put(a,b,intToRGB(255,255,255));
                        continue;
                    }
                }*/
                //SPHSim<Float,1> *sph2=static_cast<SPHSim<Float,1>*>(&sph);
                Float rho=sph.densityAtPoint1D({pos[0],y});
                //auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                outimg.put(a,b,intToRGB(rho*255,rho*255,rho*255));
            }
        }
    } else {
        for(int a=0;a<ip.imgw;a++){
            for(int b=0;b<ip.imgh;b++){
                Float x=(Float(a)/ip.imgw-0.5f)*ip.realsize;
                Float y=(Float(b)/ip.imgh-0.5f)*ip.realsize*aspect;
                auto pos=ip.center+ip.xhat*x+ip.yhat*y;

                if(gridlines){
                    Float x2=(Float(a+1)/ip.imgw-0.5f)*ip.realsize;
                    Float y2=(Float(b+1)/ip.imgh-0.5f)*ip.realsize*aspect;
                    auto pos2=ip.center+ip.xhat*x+ip.yhat*y;
                    if(!(sph.pgrid.positionToIntvec(pos)==sph.pgrid.positionToIntvec(pos2))) {
                        outimg.put(a,b,intToRGB(255,255,255));
                        continue;
                    }
                }
                Float rho=sph.densityAtPoint(pos);
                //auto rgb=hsl2rgb(0.75*c*c+0.25*s,m*0.5f+0.25f,m);
                outimg.put(a,b,intToRGB(rho*255,rho*255,rho*255));
            }
        }
    }
    //cout<<"Checked "<<(Float(nParticlesChecked)/(imw*imh))<<" particles per pixel"<<endl;
    outimg.save(getFilename(prefix,fnamei,padcount,".bmp"));
}

template<typename Float,int DIM>
static inline constexpr Float hypersphereVolume(){
    static_assert(0<DIM && DIM<=5, "Dimension higher than 5 not implemented");
    if constexpr(DIM==1){
        return 2;
    } else if constexpr(DIM==2) {
        return M_PI;
    } else if constexpr(DIM==3) {
        return Float(4)*M_PI/3;
    } else if constexpr(DIM==4) {
        return Float(1)*M_PI*M_PI/2;
    } else if constexpr(DIM==5) {
        return Float(8)*M_PI*M_PI/15;
    }
}


// Simple test that image plotting works in all dimensions.
template <typename Float, int DIM>
void initializationTest(std::string fname, bool verbose=false, int nparticles=10000){
    srand(1234567);
    ImageParamsND<Float,DIM> ip;
    ip.imgw=640;
    ip.imgh=480;
    ip.realsize=2;
    ip.center=VectorND<Float,DIM>(1);

    Float L=2.0f;
    //int nparticles=10000;

    //The average number of overlapping circles at any given point = expectedN
    // in two dimensions:
    // expectedN == nparticles*(M_PI*r*r)/(L*L);
    // in n-dimensions:
    // expectedN == nparticles*hyersphereVolume<Float,DIM>()*std::pow(r,DIM)/std::pow(L,DIM);
    Float expectedN=2.5;
    // expectedRho = nparticles * mass/std::pow(L,DIM);
    Float expectedRho=0.5;

    //Invert the two above expressions to find the particle sizes
    Float particleR=L*std::pow(Float(expectedN)/(nparticles*hypersphereVolume<Float,DIM>()),Float(1)/DIM);
    Float particleM=expectedRho*std::pow(L,DIM)/nparticles;

    if(verbose){
        std::cout<<"Running configuration "<<fname<<" for dimension "<<DIM<<"."<<std::endl;
        std::cout<<"Particle number:\t"<<nparticles<<std::endl;
        std::cout<<"Particle mass:\t"<<particleM<<std::endl;
        std::cout<<"Particle radius:\t"<<particleR<<std::endl;
    }

    VectorND<Float,DIM> domainSize(L);

    Float maxH=particleR;
    //Float dt=maxH/(6.0f*std::sqrt(2.0f*temperature));
    //Float timeelapsed=0.0f;

    ParticleList<Float,DIM> particleList;

    for ( int i = 0; i < nparticles; i++ ) 
    {
        VectorND<Float,DIM> vnew;
        VectorND<Float,DIM> pnew;
        for(int j=0;j<DIM;j++){
            pnew[j]=(rand()*L)/RAND_MAX;
        }
        particleList.plist.push_back(Particle<Float,DIM>{pnew,vnew,VectorND<Float,DIM>(),0,particleR,particleM});
    }
    SPHSim<Float,DIM> sph(particleList,domainSize,maxH);

    saveSPHImage(sph,ip,fname,0,1);
}

void configuration_drawingTest(){ 
    initializationTest<float,1>("1dtestf",true,100);
    initializationTest<float,2>("2dtestf",true);
    initializationTest<float,3>("3dtestf",true);
    initializationTest<float,4>("4dtestf",true);
    initializationTest<float,5>("5dtestf",true);
    initializationTest<double,1>("1dtestd",true,100);
    initializationTest<double,2>("2dtestd",true);
    initializationTest<double,3>("3dtestd",true);
    initializationTest<double,4>("4dtestd",true);
    initializationTest<double,5>("5dtestd",true);
    initializationTest<long double,1>("1dtestl",true,100);
    initializationTest<long double,2>("2dtestl",true);
    initializationTest<long double,3>("3dtestl",true);
    initializationTest<long double,4>("4dtestl",true);
    initializationTest<long double,5>("5dtestl",true);
}

template <typename Float, int DIM>
void annealingTest(std::string fname, bool verbose=false, int nparticles=10000, int nannealingsteps=100){
    srand(1234567);
    ImageParamsND<Float,DIM> ip;
    ip.imgw=640;
    ip.imgh=480;
    ip.realsize=2;
    ip.center=VectorND<Float,DIM>(1);

    Float L=2.0f;
    //int nparticles=10000;

    //The average number of overlapping circles at any given point = expectedN
    // in two dimensions:
    // expectedN == nparticles*(M_PI*r*r)/(L*L);
    // in n-dimensions:
    // expectedN == nparticles*hyersphereVolume<Float,DIM>()*std::pow(r,DIM)/std::pow(L,DIM);
    Float expectedN=10;
    // expectedRho = nparticles * mass/std::pow(L,DIM);
    Float expectedRho=0.5;

    //Invert the two above expressions to find the particle sizes
    Float particleR=L*std::pow(Float(expectedN)/(nparticles*hypersphereVolume<Float,DIM>()),Float(1)/DIM);
    Float particleM=expectedRho*std::pow(L,DIM)/nparticles;

    if(verbose){
        std::cout<<"Running configuration "<<fname<<" for dimension "<<DIM<<"."<<std::endl;
        std::cout<<"Particle number:\t"<<nparticles<<std::endl;
        std::cout<<"Particle mass:\t"<<particleM<<std::endl;
        std::cout<<"Particle radius:\t"<<particleR<<std::endl;
    }

    VectorND<Float,DIM> domainSize(L);

    Float maxH=particleR;
    //Float dt=maxH/(6.0f*std::sqrt(2.0f*temperature));
    //Float timeelapsed=0.0f;

    ParticleList<Float,DIM> particleList;

    for ( int i = 0; i < nparticles; i++ ) 
    {
        VectorND<Float,DIM> vnew;
        VectorND<Float,DIM> pnew;
        for(int j=0;j<DIM;j++){
            pnew[j]=(rand()*L)/RAND_MAX;
        }
        particleList.plist.push_back(Particle<Float,DIM>{pnew,vnew,VectorND<Float,DIM>(),0,particleR,particleM});
    }
    SPHSim<Float,DIM> sph(particleList,domainSize,maxH);

    Float temperature=0.00001;
    for(int i=0;i<nannealingsteps;i++){
        sph.annealPositions(temperature,particleR/80);
        temperature *= (1-Float(1)/nannealingsteps);
        if(i%10==0){
            saveSPHImage(sph,ip,fname,i/10,3);
        }
    }

}

using namespace std;

int main() {
    annealingTest<float,2>("out/annealing", true,100000,3000);

    /*
    ImageParamsND<float,2> ip;
    ip.imgw=640;
    ip.imgh=480;
    ip.realsize=2;
    //ip.center=VectorND<float,2>({1.0,1.0});
    ip.center={0.0,1.0};

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

    saveSPHImage(sph,ip,"test",0,1);
    */
    
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



