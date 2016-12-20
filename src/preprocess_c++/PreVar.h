
#include "stdafx.h"

#define P_r Particle_coordinate
#define P_n Particles_number

extern list<V3D<double> > Particle_coordinate;


extern vector<int> Particles_number;


//IOdata.cpp
#define P_rho Particles_denstiy
#define Cc Centre_coordinate

extern vector<double> Particles_denstiy;
extern vector<V3D<double> > Centre_coordinate;
extern int xuc,zuc;
extern vector<int> ncontrlp;
extern int nChain,nDp;
extern double ChainLen,RdsDp;
extern vector<V3D<double> > gap;
extern V3D<double> region,regionH,initCell;