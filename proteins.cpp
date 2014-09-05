//implementation file for proteins class

#include "headers/proteins.h"

using namespace std;

//constructors

//Default: constructs only non-crowder parameter proteins
protein::protein()
{
	for(int i=0;i<dim;i++){coordinates.push_back(0);}
	mass=M;
	radius=r;
	crowder=0;
}

//specific constructor. 
protein::protein(double x,double y, double z, double ms, double rad, bool crow)
{	
	crowder=crow;
	for(int i=0;i<dim;i++){coordinates.push_back(0); newcoords.push_back(0);vel.push_back(0);}
	coordinates[0]=x; newcoords[0]=x;
	coordinates[1]=y; newcoords[1]=y;
	coordinates[2]=z; newcoords[2]=z;
	mass=ms;
	radius=rad;
}

// member functions



void protein::setNNs(vector<int> nns, vector<int> nnns)
{
	NNs=nns;
	notNNs=nnns;
}


vector<double> protein::getv(const char* a)
{
	if(a="coords"){return coordinates;}
	else if(a="new coords"){return newcoords;}
	else if(a="velocity"){return vel;}
	else{cout<<a<<" is not a thing"<<endl; return vector<double>(3);}
	
	
}

double protein::getp(const char* a)
{
	if(a="mass"){return mass;}
	else if(a="radius"){return radius;}
	else{cout<<a<<" is not a thing"<<endl; return 0;}
}




vector<int> protein::getNNs(bool nn)
{
	if(nn){return NNs;}
	else{return notNNs;}
}


void protein::move(mt19937& gen, normal_distribution<> distro)
{
	double xx,yy,zz;
	xx=distro(gen);
	yy=distro(gen);
	zz=distro(gen);
	newcoords[0]+=xx; vel[0]=xx/h;
	newcoords[1]+=yy; vel[1]=yy/h;
	newcoords[2]+=zz; vel[2]=zz/h;
	cout<<"move vel"<<vel[0]<<endl;	
}

void protein::resetpos()
{
	for(int i=0;i<dim;i++){newcoords[i]=coordinates[i];}
}

void protein::newpos(vector<double> pos)
{
	newcoords[0]+=pos[0];
	newcoords[1]+=pos[1];
	newcoords[2]+=pos[2];
}

void protein::newpos(vector<double> nvel, double dt)
{
	newcoords[0]+=nvel[0]*dt;
	newcoords[1]+=nvel[1]*dt;
	newcoords[2]+=nvel[2]*dt;
}

void protein::update()
{
	for(int i=0;i<dim;i++)
	{
		coordinates[i]=newcoords[i];
		if(coordinates[i]>L){coordinates[i]-=2*L;}
		else if(coordinates[i]<-1.0*L){coordinates[i]+=2*L;}
	}

}

int protein::chkreac()
{
	int reac=2;
	double xx,yy,zz,mag2,mag;
	xx=newcoords[0];
	yy=newcoords[1];
	zz=newcoords[2];
	mag2=xx*xx+yy*yy+zz*zz;
	mag=pow(mag2,0.5);
	if(mag<r+radius){reac=1;}
	else if(mag>q){reac=0;}
	return reac;
}

bool protein::chkcntr()
{
	bool col=false;
	double xx,yy,zz,mag2,mag;
	xx=newcoords[0];
	yy=newcoords[1];
	zz=newcoords[2];
	mag2=xx*xx+yy*yy+zz*zz;
	mag=pow(mag2,0.5);
	if(mag<r+radius){col=true;}
	return col;
}

	























	
