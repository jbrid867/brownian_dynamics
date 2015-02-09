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
	for(int i=0;i<6;i++){Boundaries.push_back(false);}
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
		mag2+=coordinates[i]*coordinates[i];
		if(abs(coordinates[i]-L)<3*pow(10,-10)){Boundaries[i]=true;}
		else if(abs(coordinates[i]+L)<3*pow(10,-10)){Boundaries[i+1]=true;}
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
bool protein::IsNN(protein& other, int index, bool remover, int rindex)
{
	bool isNN=false;
	vector<double> pos2=other.getv("coords");
	double mag2=0, mag;
	double xx,yy,zz;
	for(int i=0;i<3;i++)
	{
		mag2+=(coordinates[i]-pos2[i])*(coordinates[i]-pos2[i]);
	}
	if(mag2<cut){if(!remover){NNs.push_back(index);} isNN=true;}
	if(Boundaries[0]&&other.NearBound(1)) // near  +x
	{
		xx=pos2[0]+2*L;
		mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+N);} 
			isNN=true;
		}
		if(Boundaries[2]&&other.NearBound(3)) // near +x and +y
		{
			yy=pos2[1]+2*L;
			mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+2*N);} isNN=true;
			}		
			if(Boundaries[4]&&other.NearBound(5)) // near +x, +y, +z
			{
				zz=pos2[2]+2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+3*N);} isNN=true;
				}
			}
			else if(Boundaries[5]&&other.NearBound(4)) // +x, +y, -z
			{
				zz=pos2[2]-2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+4*N);} isNN=true;
				}
			}
		}
		else if(Boundaries[3]&&other.NearBound(2)) // +x, -y
		{
			yy=pos2[1]+2*L;
			mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+5*N);} isNN=true;
			}		
			if(Boundaries[4]&&other.NearBound(5)) // near +x, -y, +z
			{
				zz=pos2[2]+2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+6*N);} isNN=true;
				}
			}
			else if(Boundaries[5]&&other.NearBound(4)) // +x, -y, -z
			{
				zz=pos2[2]-2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+7*N);} isNN=true;
				}
			}
		}
		else if(Boundaries[4]&&other.NearBound(5)) //+x +z
		{
			zz=pos2[2]+2*L;
			mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+8*N);} isNN=true;
			}
		}
		else if(Boundaries[5]&&other.NearBound(4)) //+x, -z
		{
			zz=pos2[2]-2*L;
			mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+9*N);} isNN=true;
			}
		}
	} // end + x boundary, never need +x again
	else if(Boundaries[0]&&other.NearBound(1)) // near  -x
	{
		xx=pos2[0]-2*L;
		mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+10*N);} isNN=true;
		}
		if(Boundaries[2]&&other.NearBound(3)) // near -x and +y
		{
			yy=pos2[1]+2*L;
			mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+11*N);} isNN=true;
			}		
			if(Boundaries[4]&&other.NearBound(5)) // near -x, +y, +z
			{
				zz=pos2[2]+2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+12*N);} isNN=true;
				}
			}
			else if(Boundaries[5]&&other.NearBound(4)) // -x, +y, -z
			{
				zz=pos2[2]-2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+13*N);} isNN=true;
				}
			}
		}
		else if(Boundaries[3]&&other.NearBound(2)) // -x, -y
		{
			yy=pos2[1]+2*L;
			mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+14*N);} isNN=true;
			}		
			if(Boundaries[4]&&other.NearBound(5)) // near -x, -y, +z
			{
				zz=pos2[2]+2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+15*N);} isNN=true;
				}
			}
			else if(Boundaries[5]&&other.NearBound(4)) // -x, -y, -z
			{
				zz=pos2[2]-2*L;
				mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
				mag=pow(mag2,0.5);
				if(mag<cut)
				{
					if(!remover){NNs.push_back(index+16*N);} isNN=true;
				}
			}
		}
		else if(Boundaries[4]&&other.NearBound(5)) //-x +z
		{
			zz=pos2[2]+2*L;
			mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+17*N);} isNN=true;
			}
		}
		else if(Boundaries[5]&&other.NearBound(4)) //-x, -z
		{
			zz=pos2[2]-2*L;
			mag2=mag2=(coordinates[0]-xx)*(coordinates[0]-xx)+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+18*N);} isNN=true;
			}
		}
	} // end - x boundary, never need -x again
	if(Boundaries[2]&&other.NearBound(3)) // +y
	{
		yy=pos2[1]+2*L;
		mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+19*N);} isNN=true;
		}
		if(Boundaries[4]&&other.NearBound(5)) // +y, +z
		{
			zz=pos2[2]+2*L;
			mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+20*N);} isNN=true;
			}
		}
		else if(Boundaries[5]&&other.NearBound(4)) // +y, -z
		{
			zz=pos2[2]-2*L;
			mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+21*N);} isNN=true;
			}
		}
	} // end +y
	else if(Boundaries[3]&&other.NearBound(2)) // -y
	{
		yy=pos2[1]-2*L;
		mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-pos2[2])*(coordinates[2]-pos2[2]);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+22*N);} isNN=true;
		}
		if(Boundaries[4]&&other.NearBound(5)) // -y, +z
		{
			zz=pos2[2]+2*L;
			mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+23*N);} isNN=true;
			}
		}
		else if(Boundaries[5]&&other.NearBound(4)) // -y, -z
		{
			zz=pos2[2]-2*L;
			mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-yy)*(coordinates[1]-yy)+(coordinates[2]-zz)*(coordinates[2]-zz);
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				if(!remover){NNs.push_back(index+24*N);} isNN=true;
			}
		}
	}
	if(Boundaries[4]&&other.NearBound(5)) //+z
	{
		zz=pos2[2]+2*L;
		mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+25*N);} isNN=true;
		}
	}
	else if(Boundaries[5]&&other.NearBound(4)) //-z
	{
		zz=pos2[2]-2*L;
		mag2=(coordinates[0]-pos2[0])*(coordinates[0]-pos2[0])+(coordinates[1]-pos2[1])*(coordinates[1]-pos2[1])+(coordinates[2]-zz)*(coordinates[2]-zz);
		mag=pow(mag2,0.5);
		if(mag<cut)
		{
			if(!remover){NNs.push_back(index+26*N);} isNN=true;
		}
	}
	if(!isNN && remover)
	{
		NNs.erase(NNs.begin()+rindex);
	}
	return isNN;
}

bool protein::NearBound(int i)
{
	return Boundaries[i];
}

void protein::setNNs(vector<int> nnns)
{
	notNNs=nnns;
}
void protein::setNNs(vector<int> nns, vector<int> nnns)
{
	NNs=nns;
	notNNs=nnns;
}


vector<double> protein::getv(string a)
{
	vector<double> blank(3);
	if(a=="coords"){return coordinates;}
	else if(a=="new coords"){return newcoords;}
	else if(a=="velocity"){return vel;}
	else{cout<<a<<" is not a thing"<<endl; return blank;}
	
	
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
	double mag2=xx*xx+yy*yy+zz*zz;
	double mag=pow(mag2,0.5);
	//cout<<"step magnitude = "<<mag<<endl;
}



void protein::setpos(vector<double> pos)
{
	for(int i=0;i<3;i++)
	{
		if(pos[i]<-1.0*L){pos[i]+=2*L;}
		else if(pos[i]>L){pos[i]-=2*L;}
		coordinates[i]=pos[i];
	}
}




void protein::newvel(vector<double> v)
{
	for(int i=0;i<3;i++){vel[i]=v[i];}
}

vector<double> protein::PBCswitch(int crowds, int index)
{

	if(index>=crowds && index<2*crowds) // through +x
	{
	colCoords[0]=coordinates[0]+2*L;
} 
	else if(index>=2*crowds && index<3*crowds) // +x +y
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]+2*L;
	} 
	else if(index>=3*crowds && index<4*crowds) // +x, +y +z
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=4*crowds && index<5*crowds) //x, y, -z
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=5*crowds && index<6*crowds) //x -y
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]-2*L;
	}
	else if(index>=6*crowds && index<7*crowds) // x -y z
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=7*crowds && index<8*crowds) // x -y -z
	{
		colCoords[0]=coordinates[0]+2*L;colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=8*crowds && index<9*crowds) // x z
	{
		colCoords[0]=coordinates[0]+2*L; colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=9*crowds && index<10*crowds) // x -z
	{
		colCoords[0]=coordinates[0]+2*L; colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=10*crowds && index<11*crowds) // -x
	{
		colCoords[0]=coordinates[0]-2*L;
	}
	else if(index>=11*crowds && index<12*crowds) // -x +y
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]+2*L;
	}
	else if(index>=12*crowds && index<13*crowds) // -x +y +z
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=13*crowds && index<14*crowds) // -x +y -z
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=14*crowds && index<15*crowds) //-x -y
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]-2*L;
	}
	else if(index>=15*crowds && index<16*crowds) //-x-y+z
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=16*crowds && index<17*crowds) //-x-y-z
	{
		colCoords[0]=coordinates[0]-2*L; colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=17*crowds && index<18*crowds) //-x+z
	{
		colCoords[0]=coordinates[0]-2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=18*crowds && index<19*crowds) //-x-z
	{
		colCoords[0]=coordinates[0]-2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=19*crowds && index<20*crowds) //+y
	{
		colCoords[1]=coordinates[1]+2*L;
	}
	else if(index>=20*crowds && index<21*crowds) //+y+z
	{
		colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=21*crowds && index<22*crowds) //+y-z
	{
		colCoords[1]=coordinates[1]+2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=22*crowds && index<23*crowds) //-y
	{
		colCoords[1]=coordinates[1]-2*L;
	}
	else if(index>=23*crowds && index<24*crowds) //-y+z
	{
		colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=24*crowds && index<25*crowds) //-y-z
	{
		colCoords[1]=coordinates[1]-2*L;colCoords[2]=coordinates[2]-2*L;
	}
	else if(index>=25*crowds && index<26*crowds) //+z
	{
		colCoords[2]=coordinates[2]+2*L;
	}
	else if(index>=26*crowds) //-z
	{
		colCoords[2]=coordinates[2]-2*L;
	}
		
	return colCoords;
}

void protein::update()
{
	double mag2=0,mag;
	for(int i=0; i<3; i++)
	{
		coordinates[i]=newcoords[i];
		if(coordinates[i]>L)
		{
			coordinates[i]-=2*L;
		}
		else if(coordinates[i]<-1.0*L)
		{
			coordinates[i]+=2*L;
		}
		// find way to change NN indices
		colCoords[i]=coordinates[i];
		mag2+=coordinates[i]*coordinates[i];
	}
	if(pow(mag2,0.5)<4*radius)
	{
		nearcenter=true;
	}
	else
	{
		nearcenter=false;
	}
}

vector<double> protein::mvVel(double t)
{
	for(int i=0;i<3;i++)
	{
		newcoords[i]=coordinates[i]+vel[i]*t;
		colCoords[i]=newcoords[i];
	}
	return newcoords;
}
void protein::nudge(double t)
{
	for(int i=0;i<3;i++)
	{
		coordinates[i]+=vel[i]*t;
		colCoords[i]=coordinates[i];
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

void protein::shift(vector<double> cpos, vector<protein>& crowders)
{
	vector<bool> pbcs(6);
	vector<int> nns;
	pbcs[0]=false;pbcs[1]=false;pbcs[2]=false;
	pbcs[3]=false;pbcs[4]=false;pbcs[5]=false;
	for(int i=0;i<3;i++)
	{
		coordinates[i]-=cpos[i];
		if(coordinates[i]>L)
		{
			coordinates[i]-=2*L;
			pbcs[2*i]=true;
		}
		else if(coordinates[i]<-1.0*L)
		{
			coordinates[i]+=2*L;
			pbcs[2*i+1]=true;
		}
	}
	/*for(int i=0;i<6;i++)
	{
		if(pbcs[i])
		{
			nns=
		}
	}*/
}

// mainP functions

bool mainP::nearEsc()
{
	return escape;
}

void mainP::upmain()
{
	double mag2=0,mag;
	for(int i=0; i<3; i++)
	{
		coordinates[i]=newcoords[i];
		mag2+=coordinates[i]*coordinates[i];
	}
	mag=pow(mag2,0.5);
	if(mag<5*radius)
	{
		nearcenter=true;
	}
	else
	{
		nearcenter=false;
	}
	if(mag>q-5*radius)
	{
		escape=true;
	}
	else
	{
		escape=false;
	}
}