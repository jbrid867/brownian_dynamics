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
	for(int i=0;i<dim;i++){coordinates.push_back(0); newcoords.push_back(0);vel.push_back(0);}
	coordinates[0]=x; newcoords[0]=x;
	coordinates[1]=y; newcoords[1]=y;
	coordinates[2]=z; newcoords[2]=z;
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
	newcoords[0]+=xx; vel[0]=xx/h;
	newcoords[1]+=yy; vel[1]=yy/h;
	newcoords[2]+=zz; vel[2]=zz/h;
		
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
			if(tm>0)
			{
				t_el+=tm;
				
				events[0]=true;
			}
			else if(tp>0)
			{
				t_el+=tp;
				
				events[0]=true;
			}
			else
			{
				cout<<"err: escape time is less than 0"<<endl;
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
			if(tm>0)
			{
				t_el+=tm;
				events[1]=true;
			}
			else if(tp>0)
			{
				t_el+=tp;
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
	vector<double> cpos(3), dr(3), PBCvec(3);
	double tp, tm, t, mag2, mag, r2, vdr, R, rad2, rad, v2;
	int Ncr=crowds.size();

	for(int i=0;i<NNs.size();i++)
	{
		vdr=0;r2=0;mag2=0;v2=0;
		cpos=crowds[NNs[i]%Ncr].getv("velocity");
		//PBCvec=PBCswitch(Ncr, NNs[i]);
		for(int k=0;k<3;k++)
		{
			//cpos[k]+=L*PBCvec[k];
			dr[k]=cpos[k]-coordinates[k];
			mag2+=dr[k]*dr[k];
			vdr+=vel[k]*dr[k];
			v2+=vel[k]*vel[k];
			r2+=dr[k]*dr[k];	
		}
		mag=pow(mag2,0.5);
		R=radius+crowds[NNs[i]%Ncr].getp("radius");
		if(mag<R)
		{
			rad2=vdr*vdr - v2*(r2-R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tp=(vdr + rad)/v2;
				tm=(vdr - rad)/v2;
				if(tm<0 && tp<t_el)
				{
					t_el=tp;
					events[2]=true;
					events[1]=false;
					events[0]=false;
				}
				else if(tm>0 && tm<t_el)
				{
					t_el=tm;
					events[2]=true;
					events[1]=false;
					events[0]=false;
				}
			}
		}
	}// end crowder loop
}

















	
