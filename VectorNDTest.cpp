#include "VectorND.h"
#include <iostream>
#include <cassert>
#include <sstream>


/*
class VectorND {
public:
    std::array<Float,DIM> x;
    VectorND();
    VectorND(std::array<Float,DIM> y);

    //vector addition and subtraction
    VectorND operator+(const VectorND& b) const;
    VectorND& operator+= (const VectorND& b);
    VectorND operator- (const VectorND& b) const;
    VectorND& operator-= (const VectorND& b);
    VectorND operator* (Float b) const;
    VectorND& operator*= (Float b);
    VectorND operator/ (Float b) const;
    VectorND& operator/= (Float b);
    //return the length squared of the vector
    Float length2() const;
    //return the length of the vector
    Float length() const;
    //calculate length and mutate to (x/length,y/length).
    void normalize();
    //return normalized version of the vector.
    VectorND normalized() const;
};

template<typename Float,int DIM>
inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c);
template<typename Float,int DIM>
inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c);
//Linearly interpolate between vec1 and vec2.
template<typename Float,int DIM>
inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2);
//cubic interpolation between v1 and v2 with control points v0 v3.
template<typename Float,int DIM>
inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3);
*/
int main() {
    //empty constructor VectorND();
    VectorND<double,3> a;
    assert(a.x[0]==double(0)&&a.x[1]==double(0)&&a.x[2]==double(0));
    a.x[0]=1;
    //copy constructor VectorND(std::array<Float,DIM> y);
    VectorND<double,3> b(a);
    assert(b.x[0]==double(1)&&b.x[1]==double(0)&&b.x[2]==double(0));

    //VectorND operator+(const VectorND& b) const;
    assert((a+a+a).x[0]==double(3)&&(a+a+a).x[1]==double(0)&&(a+a+a).x[2]==double(0));

    //VectorND& operator+= (const VectorND& b);
    b+=a;
    assert(b.x[0]==double(2)&&b.x[1]==double(0)&&b.x[2]==double(0));

    //VectorND operator- (const VectorND& b) const;
    assert((a+a-a).x[0]==double(1)&&(a+a-a).x[1]==double(0)&&(a+a-a).x[2]==double(0));
    
    //VectorND& operator-= (const VectorND& b);
    b-=a;
    assert(b.x[0]==double(1)&&b.x[1]==double(0)&&b.x[2]==double(0));
   
    //VectorND operator* (Float b) const;
    VectorND<double,2> c({1,2});
    double scalar=5;
    assert((c*scalar).x[0]==double(5)&&(c*scalar).x[1]==double(10));
  
    //inline VectorND<Float,DIM> operator* (Float b,const VectorND<Float,DIM>& c);
    assert((scalar*c).x[0]==double(5)&&(scalar*c).x[1]==double(10));

    //VectorND& operator*= (Float b);
    c*=5;
    assert(c.x[0]==double(5)&&c.x[1]==double(10));
 
    //VectorND operator/ (Float b) const;
    assert((c/5).x[0]==double(1)&&(c/5).x[1]==double(2));

    //VectorND& operator/= (Float b);
    c/=5;
    assert(c.x[0]==double(1)&&c.x[1]==double(2));

    //Float length2() const;
    VectorND<int,2> d({3,4});
    assert(d.length2()==25);

    //Float length() const;
    assert(d.length()==5);
    
    //void normalize();
    VectorND<double,3> e({1,1,1});
    e.normalize();
    assert(e.x[2]==1.0/std::sqrt(3.0));
    
    //VectorND normalized() const;
    e=VectorND<double,3>({1,1,1});
    assert(e.normalized().x[2]==1.0/std::sqrt(3.0));

    //inline std::ostream& operator<<(std::ostream& out,const VectorND<Float,DIM>& c);
    //std::cout<<VectorND<double,3>({1,2,3})<<std::endl;
    std::stringstream o;
    o<<VectorND<double,3>({1,2,3});
    assert(o.str().compare("(1,2,3)")==0);

    //inline VectorND<Float,DIM> lerp(Float t,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2);
    VectorND<double,2> v0({0,0});
    VectorND<double,2> v1({1,1});
    assert(lerp(0.5,v0,v1).x[0]==0.5&&lerp(0.5,v0,v1).x[1]==0.5);

    //inline VectorND<Float,DIM> cerp(Float t,const VectorND<Float,DIM>& vec0,const VectorND<Float,DIM>& vec1,const VectorND<Float,DIM>& vec2,const VectorND<Float,DIM>& vec3);
    VectorND<double,2> w0({1,0});
    VectorND<double,2> w1({0,1});
    VectorND<double,2> w2({2,0});
    VectorND<double,2> w3({4,1});
    assert(cerp(0.5,w0,w1,w2,w3).x[0]==0.625 && cerp(0.5,w0,w1,w2,w3).x[1]==0.5);
    //okay, this ^ is pushing double floating point equality a bit far, but it works. 
    std::cout<<"All asserts passed, woohoo."<<std::endl;

    /*
    b.x[0]+=1;
    VectorND<double,3> c(b);
    c.x[1]-=1;
    VectorND<double,3> d(c);
    d.x[2]-=1;
    VectorND<double,3> e(d);
    VectorND<double,4> f({1,2,3,4});

    std::cout<<"a: "<<a<<std::endl;
    std::cout<<"b: "<<b<<std::endl;
    std::cout<<"c: "<<c<<std::endl;
    std::cout<<"c+b: "<<(c+b)<<std::endl;
    d-=VectorND<double,3>({1,2,3});
    std::cout<<"d-(1,2,3): "<<d<<std::endl;
    std::cout<<"Vector: "<<f<<std::endl;
    std::cout<<"mynormalized: "<<(VectorND<double,6>({1,1,1,1,1,1})).normalized()<<std::endl;
    std::cout<<"myint: "<<(VectorND<int,6>({3,3,-3})).length()<<std::endl;*/
    return 0;
}

