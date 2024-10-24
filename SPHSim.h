
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
