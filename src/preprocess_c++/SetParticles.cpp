/*-----------------------------------------------------------------
SetParticles.cpp
This file contain the function to set all kinds of particles 
coordination.
----------------------------------------------------------------*/

#include "stdafx.h"
#include <numeric>   

#define P_rho Particles_denstiy
#define Cc Centre_coordinate
#define W_th Wall_thickness
#define W_ly Wall_layer
extern vector<double> P_rho;
extern vector<V3D<double> > Cc,gap;
extern int nChain,nDp;
extern double ChainLen,RdsDp;
extern V3D<double> region,initCell,W_th,W_ly;

#define P_r Particle_coordinate
#define P_n Particles_number
#define P_tn Particle_type_number
list<V3D<double> > P_r;
vector<int> P_n;
int P_tn;

//local function 
int setWallParticles();
int setFluidParticles();
int setSphereParticles();
int setChainParticles();

int setParticles()
{
	setWallParticles();
	setFluidParticles();
	if (nDp > 0) setSphereParticles();
	if (nChain > 0) setChainParticles();
	P_tn=P_n.size();

	return 0;
}

int setWallParticles()
{
	#define TV temp_V3D	
	#define S_f Switch_fact
	V3D<double> temp_V3D;
	V3D<int> Switch_fact(1);
	int i,j,k;

	TV=((W_ly+1)/gap[0].VReciprocal())*(-0.5);
	
	for (i=0;i<W_ly.x;i++){
		TV.x+=(gap[0].x*S_f.x);
		for (j=0;j<W_ly.y;j++){
			TV.y+=(gap[0].y*S_f.y);
			for (k=0;k<W_ly.z;k++){
				TV.z+=(gap[0].z*S_f.z);
				if((fabs(TV.x)>(region.x-gap[1].x)/2)||(fabs(TV.y)> (region.y-gap[1].y)/2)||(fabs(TV.z)>(region.z-gap[1].z)/2))
					P_r.push_back(TV);
			}TV.z+=(gap[0].z*S_f.z);S_f.z*=-1;			
		}TV.y+=(gap[0].y*S_f.y);S_f.y*=-1;
	}

	P_n.push_back(P_r.size());
	cout <<"Wall Particles Number:"<< P_n.back()<<endl;
	return 0;
}

int setFluidParticles()
{
	#define TV temp_V3D
	#define S_f Switch_fact
	V3D<double> temp_V3D,temp;
	V3D<int> Switch_fact(1);
	int i,j,k,m,n;

	TV=(region+gap[1])*(-0.5);
	
	for (i=0;i<initCell.x;i++){
		TV.x+=(gap[1].x*S_f.x);
		for (j=0;j<initCell.y;j++){
			TV.y+=(gap[1].y*S_f.y);
			for (k=0;k<initCell.z;k++){
				TV.z+=(gap[1].z*S_f.z);
				if(nDp == 0){ 
					P_r.push_back(TV);}
				else{
					m=0;
					for (n=0;n<nDp;n++) {temp = TV-Cc[n]; if(temp_V3D.VNorm()>RdsDp) m++;}
					if (m == nDp) P_r.push_back(TV);}
			}TV.z+=(gap[1].z*S_f.z);S_f.z*=-1;
		}TV.y+=(gap[1].y*S_f.y);S_f.y*=-1;
	}

	P_n.push_back(P_r.size()-accumulate(P_n.begin(), P_n.end(),0));
	cout <<"Fluid Particles Number:"<< P_n.back()<<endl;

	return 0;
}

int setSphereParticles()
{
	#define TV temp_V3D
	#define S_f Switch_fact
	#define S_ly Sphere_layer
	V3D<double> temp_V3D(0);
	V3D<int> Switch_fact(1),temp(1);
	double Sphere_layer;
	int i,j,k,m,n;
	for (n=0;n<nDp;n++)
	{
		TV+=Cc[n];
		S_ly=ceil(RdsDp/gap[2].x);
		P_r.push_back(TV);
		for (m=0;m<8;m++)
		{
			TV=gap[2]*(-0.5);

			switch (m)
			{
			case 1: S_f.x*=-1;break;
			case 2: S_f.y*=-1;break;
			case 3: S_f.z*=-1;break;
			case 4: S_f.x*=-1;break;
			case 5: S_f.y*=-1;break;
			case 6: S_f.x*=-1;break;
			case 7: S_f=S_f*(-1);break;
			default: break;
			}
			TV.x*=S_f.x;TV.y*=S_f.y;TV.z*=S_f.z;
			temp=S_f;

			for (i=0;i<S_ly;i++){
				TV.x+=(gap[2].x*S_f.x);
				for (j=0;j<S_ly;j++){
					TV.y+=(gap[2].y*S_f.y);
					for (k=0;k<S_ly;k++){
						TV.z+=(gap[2].z*S_f.z);
						if(TV.VNorm()<RdsDp) P_r.push_back(TV);
					}TV.z+=(gap[2].z*S_f.z);S_f.z*=-1;
				}TV.y+=(gap[2].y*S_f.y);S_f.y*=-1;
			}
			S_f=temp;
		}
		

		P_n.push_back(P_r.size()-accumulate(P_n.begin(), P_n.end(),0));
		cout <<"No. "<< n+1<<" Sphere Part Particles Number:"<< P_n.back()<<endl;
		TV-=TV;
	}
	return 0;
}

int setChainParticles()
{

	return 0;
}