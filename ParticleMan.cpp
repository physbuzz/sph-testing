#include "ParticleMan.h"
#include "glm/glm.hpp"
#include "glm/gtc/random.hpp"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <sstream>
#include <fstream>
#include <iostream>

//std::vector<PointMass> masses;

ParticleManager::ParticleManager() : x(0),y(0),vx(0),vy(0),pR(4),pM(10) { }

void ParticleManager::save(std::string filename,const std::vector<PointMass> *data){
    if(data==nullptr)
        data=&masses;
    std::ofstream pfile;
    pfile.open (filename.c_str());
    if(pfile.is_open()){
        for(int i=0;i<data->size();i++){
            const PointMass &p=data->at(i);
            pfile<<p.position.x<<" "<<p.position.y<<" "<<p.velocity.x<<" "<<p.velocity.y<<" "<<p.mass<<" "<<p.radius<<"\n";
        }
        pfile.close();
    } else {std::cerr<<"Could not open file "<<filename<<" for particle data writing";}
}

void  ParticleManager::load(std::string filename,std::vector<PointMass> *data){
    if(data==nullptr){
        data=&masses;
    }
    std::ifstream pfile;
    pfile.open (filename.c_str());
    if(pfile.is_open()){
        *data=std::vector<PointMass>();

        for(std::string line; std::getline(pfile, line); )   //read stream line by line
        {
            std::istringstream in(line);      //make a stream for the line itself
            PointMass p;
            in>>p.position.x>>p.position.y>>p.velocity.x>>p.velocity.y>>p.mass>>p.radius;
            p.gravforce.x=p.gravforce.y=p.springforce.x=p.springforce.y=0;
            (*data).push_back(p);
        }

        pfile.close();
    } else {std::cerr<<"Could not open file "<<filename<<" for particle data writing";}
}

    
void ParticleManager::setCursorPosition(double x, double y) { this->x=x; this->y=y;}
void ParticleManager::setCursorVelocity(double vx, double vy) { this->vx=vx; this->vy=vy;}
void ParticleManager::setCursorParticleRadius(double pR) { this->pR=pR;}
void ParticleManager::setCursorParticleMass(double pM) { this->pM=pM;}
void ParticleManager::placeBall(int N, double mult, double angvel){
    double golden=2.39996322972865332;

    for(int i = 0; i<N; i++){
        double theta=i*golden;
        double r=mult*pR*sqrt(theta);
        glm::dvec2 position = glm::dvec2(x+r*cos(theta),y+r*sin(theta));

        PointMass p;
        p.position=position+glm::dvec2(x,y);
        p.velocity=glm::dvec2(vx+r*cos(theta+1.57079633)*angvel, vy+r*sin(theta+1.57079633)*angvel);
        p.gravforce=glm::dvec2();
        p.springforce=glm::dvec2();
        p.mass=pM;
        p.radius=pR;
        masses.push_back(p);
    }

}

void ParticleManager::placeGaussian(int N, double sdev, double rot1, double rot2){

    typedef boost::mt19937                     ENG;    // Mersenne Twister
    typedef boost::normal_distribution<double> DIST;   // Normal Distribution
    typedef boost::variate_generator<ENG, DIST> GEN;    // Variate generator
    ENG eng;
    DIST dist(0, sdev);
    GEN gen(eng, dist);


    for(int i = 0; i<N; i++){
        glm::dvec2 position = glm::dvec2(gen(), gen());
        double theta = atan2(position.y, position.x) + 1.57079633;
        double dist = sqrt(position.x*position.x + position.y*position.y);
        double r = dist;
        double mag = exp(-rot1-rot2*r*r)*r;
        PointMass p;
        p.position=position+glm::dvec2(x,y);
        p.velocity=glm::dvec2(vx+cos(theta)*mag, vy+sin(theta)*mag);
        p.gravforce=glm::dvec2();
        p.springforce=glm::dvec2();
        p.mass=pM;
        p.radius=pR;
        masses.push_back(p);
    }
}

