/*-------------------------------------------------------------------
IOdata.cpp
This part uses to contain all the function about I/O data.
-------------------------------------------------------------------*/


#include <fstream>
#include <string> 
#include <iomanip>
#include "stdafx.h"

#define P_rho Particles_denstiy
#define Cc Centre_coordinate
#define W_th Wall_thickness
#define W_ly Wall_layer
vector<double> P_rho;
vector<V3D<double> > Cc,gap;
int nChain,nDp;
double ChainLen,RdsDp;
V3D<double> region,initCell,W_th,W_ly;

int inptconf(string filelocation)	//input the inital dat from file or terminal
{	
	#define TV temp_V3D
	#define Tdb temp_double
	double Tdb;
	V3D<double> TV;
	int i;	
	ifstream inspd;
	inspd.open(filelocation,ifstream::in);
	if(!inspd) 
	{
		cerr << "error: unable to open input file: " <<filelocation <<endl;
		system("pause");
//input the dat from the terminal
		cout << "Let's input the inital data from the terminal"<<endl;
		cout << "The wall and fluid particle density:"<<endl;
		cin >> Tdb;
		P_rho.push_back(Tdb);
		cin >> Tdb;
		P_rho.push_back(Tdb);
		cout << "The computing region scale and wall thickness:"<<endl;
		cin >>region>>W_th;
		cout << "The chain length and number:"<<endl;
		cin >> ChainLen>>nChain;
		cout << "The sphere part number:"<<endl;
		cin >> nDp;
		if(nDp > 0) 
		{
			cout << "The sphere part radius and density:"<<endl;
			cin >> RdsDp >> Tdb;
			P_rho.push_back(Tdb);
			cout <<nDp<< "  sphere parts centre coordination(= x y z): "<<endl;
			for(i = 0;i<nDp;i++) {cin >>TV ;Cc.push_back(TV);}
		}
		else{nDp = 0;RdsDp = 0.;}
	}
	else
	{
		cout <<"read in data from given Data file ..." <<endl;
		inspd >>Tdb;
		P_rho.push_back(Tdb);
		inspd >>Tdb;
		P_rho.push_back(Tdb);
		inspd >> region >> W_th;
		inspd >> ChainLen>>nChain;
		inspd >> nDp;
		if(nDp > 0) 
		{
			inspd >> RdsDp>> Tdb;
			P_rho.push_back(Tdb);
			for(i = 0;i<nDp;i++) {inspd >>TV ;Cc.push_back(TV);}

		}
		else{nDp = 0;RdsDp = 0.;}
	}
	
	for(i=0;i<3;i++) gap.push_back(pow(1.0/P_rho[i],(1./3.)));
	initCell.x=ceil(region.x/gap[1].x);
	initCell.y=ceil(region.y/gap[1].y);
	initCell.z=ceil(region.z/gap[1].z);
	region.x=gap[1].x*initCell.x;
	region.y=gap[1].y*initCell.y;
	region.z=gap[1].z*initCell.z;
	
	W_ly.x=ceil((region.x+W_th.x)/gap[0].x);
	W_ly.y=ceil((region.y+W_th.y)/gap[0].y);
	W_ly.z=ceil((region.x+W_th.z)/gap[0].z);
	W_th=W_ly/(gap[0].VReciprocal())-region;

//rewirte to the input file -->filelocation
	ofstream otspd;
	otspd.open(filelocation,ofstream::out|ofstream::trunc);
	if(!otspd){cerr << "error: unable to open output file: " <<filelocation <<endl;system("pause");}
	else{
		otspd <<P_rho[0]<<'\t'<<P_rho[1]<<endl;
		otspd <<region <<endl;
		otspd <<W_th <<endl;
		otspd << ChainLen<<'\t'<<nChain<<endl;
		otspd << nDp<<endl;
		if(nDp > 0) 
		{
			otspd << RdsDp<<'\t'<<P_rho[2]<<endl;
			for(i = 0;i<nDp;i++) {otspd <<Cc[i]<<endl;}

		}
		otspd.close();
	}

	return 0;

}


#include <numeric>
#define P_r Particle_coordinate
#define P_n Particles_number
#define P_tn Particle_type_number
extern list<V3D<double> > P_r;
extern vector<int> P_n;
extern int P_tn;

int optconf(string filelocation)
{
	#define TV temp_V3D
	int k=0,j;
	V3D<double> TV;
	list<V3D<double> > ::iterator i;	
	fstream otspd;

	otspd.open("./data/initconf.plt",ofstream::out);
	if(!otspd) {cerr << "error: unable to open file: ./data/initconf.plt"<<endl;system("pause");}
	otspd <<"ZONE"<<endl;
	otspd <<"F = POINT"<<endl;
	for (i=P_r.begin(); i!=P_r.end(); ++i) 
	{
		otspd <<*i <<endl;
		k++;
		for (j=1;j <P_tn;j++){
			if (k==accumulate(P_n.begin(),P_n.begin()+j,0)) {
				otspd <<"ZONE"<<endl;
				otspd <<"F = POINT"<<endl;
			}
		}
	}
	otspd <<"TEXT x=1 y=90 T=\" " ;
	for (j=0;j < P_tn;j++) otspd <<P_n[j]<<'\t';
	otspd <<"\""<<endl;	
	otspd.close();

// Write to binary file
	otspd.open(filelocation,fstream::out|fstream::binary|fstream::trunc);
	if(!otspd) {cerr << "error: unable to open file: "<< filelocation<<endl; system("pause");}
	else{
		otspd.write((char *)&region,sizeof(V3D<double>));
		otspd.write((char *)&W_th,sizeof(V3D<double>));
		otspd.write((char *)&nChain,sizeof(int));
		otspd.write((char *)&ChainLen,sizeof(double));
		otspd.write((char *)&nDp,sizeof(int));
		otspd.write((char *)&RdsDp,sizeof(double));
		otspd.write((char *)&P_tn,sizeof(int));
		for (j=0;j < P_tn;j++) otspd.write((char *)&P_n[j],sizeof(int));
		for (i=P_r.begin(); i!=P_r.end(); ++i) otspd.write((char *)&*i,sizeof(V3D<double>));
	otspd.close();}

	return 0;
}