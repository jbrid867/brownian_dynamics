//implementation file for proteins class

#include "headers/dem_funcs.h"


using namespace std;

//constructors

//Default: constructs only non-crowder parameter proteins
protein::protein()
{
	for(int i=0;i<dim;i++){coordinates.push_back(0);}
	mass=M;
	radius=r;
	nearcenter=false;
	//crowder=0;
}

mainP::mainP()
	: protein()
{
	escape=false;
}

//specific constructors. 
protein::protein(double x,double y, double z, double ms, double rad)
{	
	double mag2=0, mag;
	//crowder=crow;
	for(int i=0;i<dim;i++)
		{
			coordinates.push_back(0);
		 	newcoords.push_back(0);
		 	vel.push_back(0);
		 	colCoords.push_back(0);
		 }
	coordinates[0]=x; newcoords[0]=x; colCoords[0]=x;
	coordinates[1]=y; newcoords[1]=y; colCoords[1]=y;
	coordinates[2]=z; newcoords[2]=z; colCoords[2]=z;
	mass=ms;
	radius=rad;
	for(int i=0;i<3;i++)
	{
		mag2+=x*x+y*y+z*z;
	}
	cdist=pow(mag2,0.5);
	if(cdist<4*radius)
	{
		nearcenter=true;
	}
	else
	{
		nearcenter=false;
	}
}

mainP::mainP(double x,double y, double z, double ms, double rad, double centerrad)
	: protein(x,y,z,ms,rad)
{ 
	crad=centerrad;
	if(cdist>q-4*radius)
	{
		escape=true;
	}
	else
	{
		escape=false;
	}
}

///////////////////////////////////////////////////////////////////////////////////
// Protein member functions
///////////////////////////////////////////////////////////////////////////////////

void protein::setNNs(vector<int> nns, vector<int> nnns)
{
	NNs=nns;
	notNNs=nnns;
}


vector<double> protein::getv(string a)
{
	if(a=="coords"){return coordinates;}
	else if(a=="new coords"){return newcoords;}
	else if(a=="velocity"){return vel;}
	else{cout<<a<<" is not a thing"<<endl; return vector<double>(3);}
	
	
}

double protein::getp(string a)
{
	if(a=="mass"){return mass;}
	else if(a=="radius"){return radius;}
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
	newcoords[0]=coordinates[0]+xx; vel[0]=xx/h;
	newcoords[1]=coordinates[1]+yy; vel[1]=yy/h;
	newcoords[2]=coordinates[2]+zz; vel[2]=zz/h;
		
}

void protein::setpos(vector<double> pos)
{
	for(int i=0;i<3;i++){coordinates[i]=pos[i];}
}




void protein::newvel(vector<double> v)
{
	for(int i=0;i<3;i++){vel[i]=v[i];}
}

vector<double> protein::PBCswitch(int crowds, int index)
{

	if(index>=crowds && index<2*crowds){colCoords[0]=coordinates[0]-2*L;}
	else if(index>=2*crowds && index<3*crowds){colCoords[0]=coordinates[0]+2*L;}
	else if(index>=3*crowds && index<4*crowds){colCoords[1]=coordinates[1]-2*L;}
	else if(index>=4*crowds && index<5*crowds){colCoords[1]=coordinates[1]+2*L;}
	else if(index>=5*crowds && index<6*crowds){colCoords[2]=coordinates[2]-2*L;}
	else if(index>=6*crowds){colCoords[2]=coordinates[2]+2*L;}
		
	return colCoords;
}

void protein::update()
{
	for(int i=0; i<3; i++)
	{
		coordinates[i]=newcoords[i];
	}
}

vector<double> protein::mvVel(double t)
{
	for(int i=0;i<3;i++)
	{
		newcoords[i]=coordinates[i]+vel[i]*t;
	}
	return newcoords;
}
void protein::nudge(double t)
{
	for(int i=0;i<3;i++)
	{
		coordinates[i]+=vel[i]*t;
	}
}

bool protein::nearcntr()
{
	return nearcenter;
}




void protein::energy()
{
	double vx,vy,vz;
	vx=vel[0]; vy=vel[1]; vz=vel[2];
	double E;
	E=0.5*mass*(vx*vx+vy*vy+vz*vz);
	cout<<E<<endl;
}

// mainP functions

bool mainP::nearEsc()
{
	return escape;
}