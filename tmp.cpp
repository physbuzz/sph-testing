
var npart=1000;
var allRadius=14;
var allMass=326.6;
var kforce=1.9; //spring repulsion force
var mu=0.008; //velocity damping in collisions
var gforce=0.01;



//Particle array
var p = {x:[],y:[],vx:[],vy:[],rho:[],ax:[],ay:[],r:[],m:[],col:[],idx:[],idy:[]};
randomSeed(0);
//Create all the particles
for(var n=0;n<npart;n++) {
    p[n].x=random(0,400);
    //weird initial condition to concentrate particles at the bottom
    //var ytmp=random(0,sqrt(280));
    //p[n].y=400-allRadius-ytmp*ytmp;
    p[n].y=400-random(allRadius,400);
    p[n].vx=random(-1,1);
    p[n].vy=random(-1,1);
    p[n].rho=1.0;
    p.ax[n]=0;
    p.ay[n]=0;
    p.r[n]=allRadius;
    p.m[n]=allMass;
    p.col[n]=color(random(0,255),random(0,255),random(0,255));
    p.col[n]=color(0,0,0,30);
}

if(npart===400){
    for(var i=0;i<20;i++){
        for(var j=0;j<20;j++){
            var mult=0.9;
            p.x[j*20+i]=0+mult*allRadius*i;
            p.y[j*20+i]=371-mult*allRadius*j;
            p.vx[j*20+i]=random(-0.1,0.1);
            p.vy[j*20+i]=random(-0.1,0.1);
        }
    }
}
var boxwidth=2*allRadius+1;
var boxN=Math.floor(400/boxwidth);
var idarr=[];

var valididQ=function(x,y){
    return x>=0 && x<=boxN && y>=0 && y<=boxN;
};

var initializeids=function(){
    idarr=[];
    for(var i=0;i<=boxN;i++){
        idarr[i]=[];
        for(var j=0;j<=boxN;j++){
            idarr[i][j]=[];
        }
    }
    for(var n=0;n<npart;n++) {
        p[n].idx=Math.floor(p[n].x/boxwidth);
        p[n].idy=Math.floor(p[n].y/boxwidth);
        if(valididQ(p[n].idx,p[n].idy)){
            idarr[p[n].idy][p[n].idx].push(n);
            //println(idarr[p[n].idy][p[n].idx].length);
        } else {
            println("Error in initializeids(), particle out of bounds!");
        }
    }
};
initializeids();


var debugDrawBoxes=function(){
    for(var i=0;i<=boxN;i++){
        for(var j=0;j<=boxN;j++){
            noStroke();
            fill(255, 0, 0,50*idarr[j][i].length);
            rect(i*boxwidth,j*boxwidth,boxwidth,boxwidth);
            strokeWeight(1);
            stroke(0, 0, 0);
            line(i*boxwidth,j*boxwidth,(i+1)*boxwidth,j*boxwidth);
            line(i*boxwidth,j*boxwidth,i*boxwidth,(j+1)*boxwidth);
        }
    }
};

//Loop through all particles & call above function
var drawAllParticles=function(){
    noStroke();
    for(var n=0;n<npart;n++){
        fill(p.col[n]);
        ellipse(p[n].x,p[n].y,2*p.r[n],2*p.r[n]);
    }
};

//Lock the particles to the screen.
var handleParticleBounds=function(){
    for(var n=0;n<npart;n++){
    }
    
};
//Saving for posterity! Not used.
var slowForceCalc=function(dt){
    
    for(var n=0;n<npart;n++){
        for(var m=0;m<npart;m++){
            if(m===n){
                continue;
            }
            var dx=p[m].x-p[n].x;
            var dy=p[m].y-p[n].y;
            var d=dx*dx+dy*dy;
            var rn=p.r[n];
            var rm=p.r[m];
            if(d<(rn+rm)*(rn+rm)){
                d=sqrt(d);
                if(d>0.01){
                    p.ax[n]+= -dx/d * (kforce*(rn+rm-d));
                    p.ay[n]+= -dy/d * (kforce*(rn+rm-d));
                }
                var dvx=p[m].vx-p[n].vx;
                var dvy=p[m].vy-p[n].vy;
                var dv=sqrt(dvx*dvx+dvy*dvy);
                if(dv>0.01){
                    p.ax[n]+= mu*Math.abs(dvx)*dvx/dv*(rn+rm-d)/(rn+rm);
                    p.ay[n]+= mu*Math.abs(dvy)*dvy/dv*(rn+rm-d)/(rn+rm);
                }
            }
        }
        p.ay[n]+=gforce;
        if(p[n].y>400-p.r[n]){
            p.ay[n]-=kforce*(p[n].y+p.r[n]-400);
        
        }
    }
};
var pconst=1000.3;
var intEnergy=function(rho){
    //return pconst*pow(rho,2.0/3.0);
    //This should be the integral of pressureEOS d rho.
    return 0;
};
var pressureEOS=function(rho){
    //return 2.0/3.0*pconst*pow(rho,5.0/3.0);
return pconst*(rho/2.0 - 1.0)/max(rho*rho, 0.001);

};
var totalEnergy=0.0;
var fastForceCalc=function(dt){
    totalEnergy=0.0;
    for(var n=0;n<npart;n++){
        p[n].rho=0;
        for(var a=-1;a<=1;a++){
            for(var b=-1;b<=1;b++){
                var midx=p[n].idx+a;
                var midy=p[n].idy+b;
                if(valididQ(midx,midy)){
                    for(var k=0;k<idarr[midy][midx].length;k++){
                        var m=idarr[midy][midx][k];
            var dx=p[m].x-p[n].x;
            var dy=p[m].y-p[n].y;
            var d=Math.sqrt(dx*dx+dy*dy);
            var R=d/(2*allRadius);
            if(d<2*allRadius){
                var wab=5.0/(Math.PI*(2*allRadius)*(2*allRadius))*(1+3*R)*(1-R)*(1-R)*(1-R);
                p[n].rho+=p.m[m]*wab;
            }
                        
                    }
                }
            }
        }
    }
    for(var n=0;n<npart;n++){
        if(n===0){
            //println(p[n].rho);
        
        }
        totalEnergy+=0.5*p.m[n]*(p[n].vx*p[n].vx+p[n].vy*p[n].vy);
        totalEnergy+=p.m[n]*intEnergy(p[n].rho);
        for(var a=-1;a<=1;a++){
            for(var b=-1;b<=1;b++){
                var midx=p[n].idx+a;
                var midy=p[n].idy+b;
                if(valididQ(midx,midy)){
                    for(var k=0;k<idarr[midy][midx].length;k++){
                        var m=idarr[midy][midx][k];
                        if(m===n){
                            continue;
                        }
            var dx=p[m].x-p[n].x;
            var dy=p[m].y-p[n].y;
            var d=dx*dx+dy*dy;
            var rn=p.r[n];
            var rm=p.r[m];
            if(d<(rn+rm)*(rn+rm)){
                d=sqrt(d);
                var dvx=p[m].vx-p[n].vx;
                var dvy=p[m].vy-p[n].vy;
                var dv=sqrt(dvx*dvx+dvy*dvy);
                
                if(d>0.001){
                    var sphConst=5*12/(Math.PI*(2*rn)*(2*rn)*(2*rn));
                    var R=d/(rn+rm);
                    var gradwabmag=sphConst*R*(1-R)*(1-R);
                    var gradwabx=-dx/d*gradwabmag;
                    var gradwaby=-dy/d*gradwabmag;
                    var fmag=pressureEOS(p.rho[m])/(p.rho[m]*p.rho[m])+pressureEOS(p[n].rho)/(p[n].rho*p[n].rho);
                    p.ax[n]+=p.m[m]*gradwabx*fmag;
                    p.ay[n]+=p.m[m]*gradwaby*fmag;
                    
                    //p.drho[n]+= p.m[m]/p.rho[m]*(dvx*gradwabx+dvy*gradwaby);
                }
                
                if(dv>0.01){
                    p.ax[n]+= mu*Math.abs(dvx)*dvx/dv*(rn+rm-d)/(rn+rm);
                    p.ay[n]+= mu*Math.abs(dvy)*dvy/dv*(rn+rm-d)/(rn+rm);
                }
            }
                            
                            
                            
                        
                    }
                }
            }
        }
        
        p.ay[n]+=gforce;
        totalEnergy+= -gforce*p.m[n]*(p[n].y-400);
        
        if(p[n].y>400-p.r[n]){
            p.ay[n]-=kforce*(p[n].y+p.r[n]-400);
            totalEnergy+= 0.5*kforce*p.m[n]*(p[n].y+p.r[n]-400)*(p[n].y+p.r[n]-400);
        }
        
        if(p[n].x<=350 && p[n].y<=350 && p[n].x>= 250 && p[n].y>=250){
            p.ay[n]-=0.05;
        
        }
        
        if(p[n].x<=150 && p[n].y<=350 && p[n].x>= 050 && p[n].y>=250){
            p.ay[n]-=0.05;
        
        }
    }
};
var handleParticlePhysics = function(dt){
    
    for(var n=0;n<npart;n++){
        p[n].vx+=p.ax[n]*dt;
        p[n].vy+=p.ay[n]*dt;
            var s=(Math.sin(3*p[n].vx)+1.0)/2.0;
            var c=(Math.sin(3*p[n].vy)+1.0)/2.0;
            p.col[n]=color(0, s*100, c*100,80);
        p.ax[n]=0;
        p.ay[n]=0;
        p[n].x+=p[n].vx*dt;
        p[n].y+=p[n].vy*dt;
        if(p[n].x<0){
            p[n].x=0;
            p[n].vx=-p[n].vx;
        }
        if(p[n].y<0){
            p[n].y=0;
            p[n].vy=-p[n].vy;
        }
        if(p[n].x>400){
            p[n].x=400;
            p[n].vx=-p[n].vx;
        }
        if(p[n].y>400){
            p[n].y=400;
            p[n].vy=-p[n].vy;
        }
        var newidx=Math.floor(p[n].x/boxwidth);
        var newidy=Math.floor(p[n].y/boxwidth);
        if(newidx!==p[n].idx || newidy!==p[n].idy){
            if(!(newidx>=0 && newidx<=boxN && newidy>=0 && newidy<=boxN)){
                println("Particle out of bounds after bounds checking! This should never happen!");
                continue;
            }
            var nlocation=idarr[p[n].idy][p[n].idx].indexOf(n);
            if (nlocation > -1) {
                idarr[p[n].idy][p[n].idx].splice(nlocation, 1);
            } else {
                println("Particle not found in idarr! idarr got out of sync?");
                continue;
            }
            idarr[newidy][newidx].push(n);
            p[n].idx=newidx;
            p[n].idy=newidy;
        }
    }
    fastForceCalc();
};
    
handleParticlePhysics(0.0);
var dt=0.6;
//println("Simulating with dt="+dt);
//println("Energy at t=0: "+totalEnergy);
var ttot=0;
var lastt=0;
//called once per frame
var draw= function() {
    background(255, 255, 255,2);
    drawAllParticles();
    for(var k=0;k<1;k++){
        ttot+=dt;
        handleParticlePhysics(dt);
    }
    //println(totalEnergy);

    lastt=ttot;
    //debugDrawBoxes();
};
