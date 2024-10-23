#ifndef PARTICLEMAN_H
#define PARTICLEMAN_H
#include <string>
#include <vector>
#include "phystructs.h"





class ParticleManager {
    double x,y,vx,vy,pR,pM;
public:
    std::vector<PointMass> masses;

    ParticleManager();

    void save(std::string filename,const std::vector<PointMass> *data=nullptr);
    void load(std::string filename,std::vector<PointMass> *data=nullptr);
    
    void setCursorPosition(double x, double y);
    void setCursorVelocity(double vx, double vy);
    void setCursorParticleRadius(double pR);
    void setCursorParticleMass(double pR);
    void placeBall(int N,double rmult=1, double angvel=0);
    void placeGaussian(int N, double sdev, double rot1=1000, double rot2=0);
};
#endif //PARTICLEMAN_H
