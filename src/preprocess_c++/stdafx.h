// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

//#pragma once


// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//


// 1.Containing the all the header files(#include <> or #include "" ) we need to be used in application.
// syntax format: #include <> or #include ""

// #include "targetver.h"

#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
using namespace std;



// TODO: reference additional headers your program requires here




// 2.Containing the all the constant variables we need to be used in application.
// syntax format: constant 

double const  pi=3.1415926535897932384626433832795028841971693993751;


//====================================================================================================================
// 3.Containing the all the basic classes we need to be used in application.
// syntax format: constant


//3.1 3D vector class -----------------------------------------------------------252

#define V3D Vector3D
#define tnT template <typename nT>

template <class nTc>
class Vector3D 
{
public:
//constructor
    V3D() : x(0),y(0),z(0) {}
	tnT
    V3D( nT t) {x=t;y=t;z=t;}
	tnT
    V3D(nT vx, nT vy, nT vz) : x(vx),y(vy),z(vz) {}
	tnT
    V3D(const V3D<nT> &vt) {x=vt.x;y=vt.y;z=vt.z;}                  // need to add const ????
//special function for Vector    
    double VNorm() {return sqrt(x*x+y*y+z*z);}
    V3D VDirection() {return V3D(x/VNorm(),y/VNorm(),z/VNorm());}
	V3D VReciprocal() {return V3D(1/x,1/y,1/z);}
    
// Vector3D computer + - * /.
    // operator with the other Vector3D object.

    V3D operator + (V3D &vt) {return V3D(x+vt.x,y+vt.y,z+vt.z);}

    V3D operator - (V3D &vt) {return V3D(x-vt.x,y-vt.y,z-vt.z);}

    V3D operator / (V3D &vt) {return V3D(x/vt.x,y/vt.y,z/vt.z);}
    
    //operator '*' reload to be dot product function.
    double operator * (V3D &vt) {return (x*vt.x+y*vt.y+z*vt.z);}
    //operator '/' reload to be cross product function.
    V3D operator % (V3D &vt) {return V3D(y*vt.z-z*vt.y,-(x*vt.z-z*vt.x),x*vt.y-y*vt.x);}
	//operator += -=
	V3D operator += (V3D &vt) {x+=vt.x;y+=vt.y;z+=vt.z;return V3D(x,y,z);}
	V3D operator -= (V3D &vt) {x-=vt.x;y-=vt.y;z-=vt.z;return V3D(x,y,z);}
    // operator with the number (int float double and so on.).
	// "+"
	V3D operator + (int t) {return V3D(x+t,y+t,z+t);}
	V3D operator + (float t) {return V3D(x+t,y+t,z+t);}
	V3D operator + (double t) {return V3D(x+t,y+t,z+t);}
	V3D operator + (long double t) {return V3D(x+t,y+t,z+t);}
	// "-"
	V3D operator - (int t) {return V3D(x-t,y-t,z-t);}
	V3D operator - (float t) {return V3D(x-t,y-t,z-t);}
	V3D operator - (double t) {return V3D(x-t,y-t,z-t);}
	V3D operator - (long double t) {return V3D(x-t,y-t,z-t);}
    // "*"
	V3D operator * (int t) {return V3D(x*t,y*t,z*t);}
	V3D operator * (float t) {return V3D(x*t,y*t,z*t);}
	V3D operator * (double t) {return V3D(x*t,y*t,z*t);}
	V3D operator * (long double t) {return V3D(x*t,y*t,z*t);}
	// "/"	
	V3D operator / (int t) {return V3D(x/t,y/t,z/t);}
	V3D operator / (float t) {return V3D(x/t,y/t,z/t);}
	V3D operator / (double t) {return V3D(x/t,y/t,z/t);}
	V3D operator / (long double t) {return V3D(x/t,y/t,z/t);}

//  input and output  >> <<
    friend istream& operator >> (istream &input,V3D &vt)
    {
        input >> vt.x >> vt.y>>vt.z;
        return input;
    }
    friend ostream& operator << (ostream &output,V3D &vt)
    {
        output << '\t'<< vt.x<<'\t'<< vt.y<<'\t'<<vt.z;
        return output;
    }
//logical opertor ==
    V3D operator == (V3D &vt) {return (x==vt.x && y==vt.y && z==vt.z);}
    tnT
    V3D operator == (nT t) {return (x==t && y==t && z==t);}
 //variable  
    nTc x,y,z;
 };


//3.2 mass point(single particle) class ------------------------------------------------
#define MP MassPoint
#define Pv Point_velocity
#define Pr Point_coordinate
#define Pa Point_accelation
#define Pm Point_mass
#define PE_k Point_kinetic_Energy

class MassPoint
{
public:
//constructor
    MP(): Point_mass(1.0){}
    MP(const MP &mpt) : Pv(mpt.Pv),Pr(mpt.Pr),Pa(mpt.Pa),Pm(mpt.Pm){}
    MP(const V3D<double> &vt,const V3D<double> &rt) : Pv(vt),Pr(rt){}
    MP( V3D<double> &ft,double m) {Pm=m;Pa=ft/Pm;}

//  input and output  >> <<
    friend istream& operator >> (istream &input,MP &mpt)
    {
        input >> mpt.Pr>> mpt.Pv>>mpt.Pa>>mpt.Pm;
        return input;
    }
	friend ostream& operator << (ostream &output,MP &mpt)
    {
        output << mpt.Pr<<'\t'<<mpt.Pv<<'\t'<<mpt.Pa<<'\t'<< mpt.Pm ;
        return output;
    }
// couple MassPoint interaction physical quantity
    V3D<double> DistanceTwoPionts(MP &mpt) {return mpt.Pr-Pr;}
    V3D<double> RelativeVelocity(MP &mpt) {return mpt.Pv-Pv;}
// computing function
    int StepCompute(V3D<double> &ft,double dt) 
    {
        V3D<double> Pv0(Pv);
        V3D<double> Pr0(Pr);
        Pa=ft/Pm;
        Pv=Pa*dt;
		Pr=Pr+Pv0*dt+Pa*(dt*dt/2.0);
		PE_k = Pv*Pv;
        return 0;
    }

//variable
    double Point_mass;
    V3D<double> Point_coordinate;
    V3D<double> Point_velocity;
    V3D<double> Point_accelation;
protected:
	V3D<double> Point_kinetic_Energy;
};



//3.3 DPDparticle class --------------------------------------------------------------------
#define DPDp DPDparticle
#define P_ID Particle_Identification
#define P_nCB Particle_nCrossBoundary
#define P_nC Particle_nCell
#define P_rc Particle_Radius_Cutoff
#define P_af Particle_Force_amplitude
#define P_E_k Particle_kinetic_Energy
#define P_E_p Particle_potential_Energy
#define P_E_t Particle_total_Energy

class DPDparticle : public MassPoint
{	
public:
//constructor
	DPDp(): P_ID(-1),P_nCB(0),P_nC(0),P_E_k(0.0),P_E_p(0.0),P_E_t(0.0){};

//special function for DPD particles
	V3D<int> PnCellCompute(V3D<double> &region,V3D<double> &lengthCell)		//compute the cell number to which the particle belongs
	{
		V3D<double> temp_V3D;
		V3D<int> nCell;

		temp_V3D=Pr/lengthCell;
		nCell.x=ceil(temp_V3D.x);
		nCell.y=ceil(temp_V3D.y);
		nCell.z=ceil(temp_V3D.z);
		return nCell;
	}
//variable
	int Particle_Identification;
	V3D<int> Particle_nCrossBoundary;
	V3D<int> Particle_nCell;
	double Particle_Radius_Cutoff[3],Particle_Force_amplitude[3];
	double Particle_kinetic_Energy, Particle_potential_Energy,Particle_total_Energy;
};


//3.4 MDPDparticle class -------------------------------------------------------------------
#define MDPDp MDPDparticle

class MDPDparticle : public DPDparticle
{
public:
//constructor
	MDPDp(): P_ID(1){};
//variable
	int Particle_Identification;
};


//3.5 Single Cell class -----------------------------------------------------------
#define SG SingleGrid
#define G_l Grid_length
#define G_Sx Grid_Stress_x
#define G_Sy Grid_Stress_y
#define G_Sz Grid_Stress_z
#define G_P Grid_Press
#define G_rho Grid_Density
#define G_E_k Grid_kinetic_Energy
#define G_E_p Grid_potential_Energy
#define G_E_t Grid_total_Energy

class SingleGrid : public MassPoint
{
public:
	SG(): G_P(0.0),G_rho(0.0),G_E_k(0.0),G_E_p(0.0),G_E_t(0.0){};
//variable
	V3D<double> Grid_length;
	V3D<double> Grid_Stress_x,Grid_Stress_y,Grid_Stress_z;
	double Grid_Press,Grid_Density;
	double Grid_kinetic_Energy;
	double Grid_potential_Energy;
	double Grid_total_Energy;
};