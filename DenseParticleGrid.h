#ifndef DENSEPARTICLEGRID_H 
#define DENSEPARTICLEGRID_H

#include <vector>
#include <iostream>
#include <cmath>
#include "VectorND.h"





/* PType must have a VectorND<Float,DIM> pos for position, and a SizeType id
 * The idea with SizeType is that we don't want to waste 4 bytes with an 8 bit unsigned integer
 * if our array is small and we have millions of ids. But I can already see I'm being inconsistent
 * about it vs size_t vs int, so I need to revisit this.
 * */
template<typename Float, int DIM, typename PType, typename SizeType>
class DenseParticleGrid {
    static_assert(1<=DIM);
    std::vector<std::vector<SizeType> > idarr;
    VectorND<int,DIM> numCells;
    VectorND<Float,DIM> boxWidth;
    VectorND<Float,DIM> domainSize;

    VectorND<Float,DIM> domainMax;//iota smaller than domainSize.
    Float maxRadius;
    bool needsRebuild;
    SizeType productOfSizes;

public:
    //This should really be a private, but if we're pushing the boundaries of the RAM 
    //I really don't want to make a deep copy!
    std::vector<PType> *p;
    DenseParticleGrid() : idarr(), numCells(), boxWidth(), domainSize(),domainMax(),maxRadius(0),needsRebuild(true),productOfSizes(1), p(nullptr) {
        for(size_t i=0;i<DIM;i++){
            numCells[i]=1;
            boxWidth[i]=Float(1);
            domainSize[i]=Float(1);
            domainMax[i]=std::nextafter(domainSize[i],Float(0));
        }
    }
    void setParameters(VectorND<Float,DIM> domainSize,VectorND<int,DIM> numCells, Float maxRadius) {
        //I'm imagining all sorts of floating point issues;
        //what matters is that boxWidth[i]*numCells[i] is the exclusive rightmost part of the domain. 
        //So I want to set domainMax[i] to this number minus iota, so setting a particle's position to 
        //domainMax[i] will let particle.position[i]<boxWidth[i]*numCells[i].
        //Well, my reasoning is sound but I should really check that being this careful does anything

        this->maxRadius=maxRadius;
        this->numCells=numCells;
        for(size_t i=0;i<DIM;i++){
            if(numCells[i]<=0){
                std::cerr<<"Sanity Check Failed: DenseParticleGrid::setParameters called with nonsense number of cells "<<numCells[i]<<"!"<<std::endl;
            }
            if(domainSize[i]<=0){
                std::cerr<<"Sanity Check Failed: DenseParticleGrid::setParameters called with nonsense domain size "<<domainSize[i]<<"!"<<std::endl;
            }
            this->boxWidth[i]=domainSize[i]/numCells[i];
            this->domainSize[i]=this->boxWidth[i]*numCells[i];
            this->domainMax[i]=std::nextafter(domainSize[i],Float(0));
        }
        needsRebuild=true;
    }
    
    void rebuildAndClamp(){

        idarr=std::vector<std::vector<uint> >(numCellsX*numCellsY);
        for(size_t n=0;n<p.size();n++) {
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
        needsRebuild=false;
    }
    


};

/** GridHandler
 * takes in a bunch of Point Masses and calculates short-distance forces 
 * */

class GridHandler
{
	gridmaptype *gridmap;
	double multiplier;
	int expectedParticles;
public:
	GridHandler(double positionMultiplier, int expectedParticles);
	~GridHandler();
	void clearPoints();
	void addPointMass(PointMass *arg);
	void calculateForces(PointMass *arg, const double& K, const double& damping) const;
};
//	std::unordered_map<std::pair<int,int>,std::vector<PointMass*>> *gridmap;
//	double multiplier;
GridHandler::GridHandler(double positionMultiplier,int expectedParticles) :
	gridmap(NULL),
	multiplier(positionMultiplier),
	expectedParticles(expectedParticles) {


}
GridHandler::~GridHandler(){
	clearPoints();
}
void GridHandler::clearPoints(){
	if(gridmap){
		delete gridmap;
		gridmap=NULL;
	}
}

void GridHandler::addPointMass(PointMass *arg) {
	if(!gridmap)
		gridmap=new gridmaptype(expectedParticles);

	//if the map is instantiated, insert an element with a key of type std::pair<int,int>, representing integer location coordinates.
	//(the particle's (x,y) position is mapped to (floor(x*c),floor(y*c)).
	//If no list has been created to store particles at the grid location, then
	std::pair<int,int> key=std::make_pair((int)floor(arg->position.x*multiplier),(int)floor(arg->position.y*multiplier));
	auto f=gridmap->find(key);
	if(f==gridmap->end()){
		gridmap->insert(std::make_pair(key,std::vector<PointMass*>(1,arg)));
	} else {
		f->second.push_back(arg);
	}
}
void GridHandler::calculateForces(PointMass *arg, const double& K, const double& damping) const{
	int x0=(int)floor(arg->position.x*multiplier);
	int y0=(int)floor(arg->position.y*multiplier);

	//iterate over the 9x9 grid around the argument's position
	for(int x=-1;x<=1;x++){
		for(int y=-1;y<=1;y++){
			auto f=gridmap->find(std::make_pair(x0+x,y0+y));
			if(f!=gridmap->end()){
				//if there's a vector, iterate over it and calculate each resulting force.
				for(size_t i=0; i<f->second.size(); i++) {
					PointMass *val=f->second.at(i);

					glm::dvec2 diff=val->position-arg->position;
					double d=diff.x*diff.x+diff.y*diff.y;
					if(d<.00001) continue;

					d=sqrt(d);
					diff=diff*1.0/d; //normalize difference

					double overlap=arg->radius+val->radius-d;
					if(overlap>0) //collision
					{
						double force=-K*overlap; //spring force

						//the dot product gives the relative velocity in the direction from arg to val.
						glm::dvec2 relvel=val->velocity-arg->velocity;
						force+=damping*glm::dot(diff,relvel);

						arg->springforce+=diff*force;
					}
				}
			}
		}
	}
}
#endif
