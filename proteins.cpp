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
	collision=false;
	//crowder=0;
}

mainP::mainP()
	: protein()
{
	escape=false; reaction = false;
	collision=false;
}

//specific constructors. 
protein::protein(double x,double y, double z, double ms, double rad)
{	
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
	collision=false;
}

mainP::mainP(double x,double y, double z, double ms, double rad, double centerrad)
	: protein(x,y,z,ms,rad)
{
	
	escape=false; reaction = false; collision=false;
	crad=centerrad;
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


void protein::newpos(vector<double> pos)
{
	
	newcoords[0]=pos[0];
	newcoords[1]=pos[1];
	newcoords[2]=pos[2];
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

bool protein::colCheck(int& index, double& t, vector<protein> crowds) // returns [0] = index, [1] = time, [2]=bool
{
	// vectors
	vector<double> pos2;
	//vars
	double R2, R, dr2_o=0, dr2=0, v2=0, vdr=0, rad2, rad, tm, tp;
	int Ncr=crowds.size();
	bool col=false;
	for(int k=0;k<3;k++)
	{
		v2+=vel[k]*vel[k];
	}

	for(int i=0; i<NNs.size();i++)
	{
		dr2_o=0, dr2=0, vdr=0;
		pos2=crowds[NNs[i]%Ncr].PBCswitch(Ncr, NNs[i]);
		R=crowds[NNs[i]%Ncr].getp("radius")+radius;
		for(int k=0; k<3; k++)
		{
			dr2_o+=(pos2[k]-coordinates[k])*(pos2[k]-coordinates[k]);
			dr2+=(pos2[k]-newcoords[k])*(pos2[k]-newcoords[k]);
			vdr+=vel[k]*(pos2[k]-newcoords[k]);
		}
		if(pow(dr2,0.5)<R)
		{
			rad2=vdr*vdr - v2*(dr2_o-R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tp=(vdr+rad)/v2;
				tm=(vdr-rad)/v2;
				if(tm>0 && tm<t)
				{
					t=tm;
					col=true;
					index=NNs[i];
				}
				else if(tp>0 && tp<t)
				{
					t=tp;
					col=true;
					index=NNs[i];
				}
				else if(tp<0 && tm<0)
				{
					cout<<"err: negative times"<<endl;
				}
			}
			else
			{
				cout<<"err: negative rad"<<endl;
			}
		}
	}
	return col;
}


/*void proteins::resolve(protein& two, double t_el, int index)
{
	vector<double> pos2=two.PBCswitch(N, index), v1f(3), v2f(3);
	double 
}*/

void protein::energy()
{
	double vx,vy,vz;
	vx=vel[0]; vy=vel[1]; vz=vel[2];
	double E;
	E=0.5*mass*(vx*vx+vy*vy+vz*vz);
	cout<<E<<endl;
}


////////////////////////////////////////////////////////////////////////////
// MainP member functions
////////////////////////////////////////////////////////////////////////////

void mainP::EscReac(double& t_el, vector<bool>& events)
{
	double v2=0, vdr=0, mag_r=0,rf=0,ro=0, ro2=0, rf2=0, rad=0, rad2=0;
	double tp,tm,esc_t=0;
	for(int i=0;i<3;i++)
	{
		ro2+=coordinates[i]*coordinates[i];
		rf2+=newcoords[i]*newcoords[i];
		v2+=vel[i]*vel[i];
		vdr+=vel[i]*coordinates[i];
	}

	rf=pow(rf2,0.5);
	
	if(rf>q)
	{
		cout<<"potential escape"<<endl;
		/*cout<<"velocity "<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<endl;
		cout<<"position "<<newcoords[0]<<" "<<newcoords[1]<<" "<<newcoords[2]<<endl;*/
		rad2=vdr*vdr-v2*(ro2-q*q);
		if(rad2>0)
		{
			rad=pow(rad2,0.5);
			tp=(-1.0*vdr + rad)/v2;
			tm=(-1.0*vdr - rad)/v2;
			if(tm>0 && tm<h)
			{
				t_el=tm;
				
				events[0]=true;
			}
			else if(tp>0 && tp<h)
			{
				t_el=tp;
				
				events[0]=true;
			}
			else
			{
				cout<<"err: escape time is less than 0 or something"<<endl;
			}
		}
		else
		{
			cout<<"err: negative in the radical"<<endl;
		}
	}
	else if(rf<radius+crad)
	{
		cout<<"reaction"<<endl<<rf<<endl;
		rad2=vdr*vdr-v2*(ro2-(radius+crad)*(radius+crad));
		if(rad2>0)
		{
			rad=pow(rad2,0.5);
			tp=(-1.0*vdr + rad)/v2;
			tm=(-1.0*vdr - rad)/v2;
			if(tm>0 && tm<h)
			{
				t_el=tm;
				events[1]=true;
			}
			else if(tp>0 && tp<h)
			{
				t_el=tp;
				events[1]=true;
			}
			else
			{
				cout<<"err: reaction time is less than 0"<<endl;
			}
		}
		else
		{
			cout<<"err: negative in radical"<<endl;
		}
	}
	else{cout<<"nothing happened"<<endl;} //temporary
}


void mainP::collisions(double& t_el, vector<protein>& crowds, vector<bool>& events)
{
	vector<double> pos1(3),pos2(3), npos(3), dr(3), dr_o(3), vi(3);
	vector<int> moving_particles, nns;
	double t_rem=h;
	double tp, tm, t, mag2, mag, mag2_o, mag_o, r2, vdr, R, rad2, rad, v2;
	int Ncr=crowds.size();
	int index1, index2;
	bool ismain=true, mainCheck=false, resolved=false, iscenter=false;

	//Collision resolution stuff
	vector<double> vp(3), vt(3), v2f(3), v1i(3), v2i(3), v1f(3);
	double vp2=0, vpi, m1, m2, v1p, v2p, vperp;

	while(!resolved)
	{
		if(ismain) // set position of particle 1 and get NNs
		{
			
			if(mainCheck)
			{
				EscReac(t_el, events);
			}
		}
		else if(iscenter) // work on this later
		{
			pos1[0]=0;pos1[1]=0;pos1[2]=0;
		}
		else
		{
			npos=crowds[index1].getv("new coords");
			pos1=crowds[index1].getv("coords");
			nns=crowds[index1].getNNs(true);
		}

		for(int i=0;i<nns.size();i++)
		{
			r2=0; v2=0; vdr=0;
			pos2=crowds[nns[i]%Ncr].PBCswitch(Ncr, nns[i]);
			for(int k=0; k<3; k++)
			{
				
			}
		}



	}


}

















	
