#include <vector>
#include "proteins.h"

using namespace std;

vector<double> PBCswitch(int crowds, int index)
{
	vector<double> b(3);
	for(int i=0;i<3;i++){b[i]=0;}

	if(index>=crowds && index<2*crowds){b[0]=-2;}
	else if(index>=2*crowds && index<3*crowds){b[0]=2;}
	else if(index>=3*crowds && index<4*crowds){b[1]=-2;}
	else if(index>=4*crowds && index<5*crowds){b[1]=2;}
	else if(index>=5*crowds && index<6*crowds){b[2]=-2;}
	else if(index>=6*crowds){b[2]=2;}
		
	return b;
}

vector<double> PDF(vector<protein> crowds, int Ncr) // For all crowders
{
	int bins=1000;
	int count=0;
	double cut=L/2.0;
	double space=cut/100.0;
	double mag, mag2;
	int spot;

	vector<double> pos1(3), pos2(3);
	vector<double> distr(100);

	for(int i=0;i<Ncr;i++)
	{
		pos1=crowds[i].getv("coords");
		for (int j=i+1;j<Ncr;j++)
		{
			mag=0; mag2=0;
			pos2=crowds[j].getv("coords");

			for(int k=0;k<3;k++)
			{
				mag2+=(pos1[k]-pos2[k])*(pos1[k]-pos2[k]);
			}	
			mag=pow(mag2,0.5);
			if(mag<cut)
			{
				spot = floor(mag/space);
				distr[spot]+=1;
				count+=1;
			}
		}
	}

	for(int i=0;i<bins;i++){distr[i]=distr[i]/(count*phi); distr[i]/(4*pi*((i+1)*space)*((i+1)*space)*space);}

	return distr;
}


vector<double> PDF1(vector<protein> crowds, int Ncr) //g(r) for central particle
{
	int bins=100;
	int count=0;
	//double cut=L-3*pow(10,-10); gone for PBC stuff
	//double space=cut/bins;
	double mag, mag2, mag2p;

	int spot;


	// PBC box stuff
	vector<double> pbcmag2(3);
	vector<bool> pbcCheck(3);
	double R = pow(6.0/pi,1/3.0)*L;
	double space=R/bins;
	int x,y,z;


	vector<double> posc(3);
	vector<double> distr(bins);
	//cout<<"size of crowder vector is "<<crowds.size()<<endl;
	for(int i=0;i<Ncr;i++)
	{
		for(int j=0;j<3;j++){pbcmag2[j]=0;pbcCheck[j]=false;}
		mag2=0;
		cout<<"accessing crowder "<<i<<endl;

		posc=crowds[i].getv("coords");
		//cout<<"check1.9"<<endl;
		for(int k=0;k<3;k++)
		{
			mag2+=(posc[k])*(posc[k]);
			if(posc[k]>L-R)
			{
				pbcCheck[k]=true;
				//cout<<k<<" pbs pos = "<<posc[k]<<endl;
				pbcmag2[k]=(posc[k]-2*L)*(posc[k]-2*L)-posc[k]*posc[k];
			}
			else if(posc[k]<R-L)
			{
				pbcCheck[k]=true;
				//cout<<k<<" pbcpos = "<<posc[k]<<endl;
				pbcmag2[k]=(posc[k]+2*L)*(posc[k]+2*L)-posc[k]*posc[k];
			}
		}
		mag=pow(mag2,0.5);
		if(mag<R)
		{
			spot=floor(mag/space);
			cout<<"bin number = "<<spot<<endl;
			distr[spot]+=1;
			count+=1;
			
			for(int el=0;el<3;el++)
			{
				if(pbcCheck[el])
				{
					mag2p=mag2+pbcmag2[el];
					mag=pow(mag2p,0.5);
					if(mag<R){
						spot=floor(mag/space);
						cout<<"bin number = "<<spot<<endl;
						distr[spot]+=1;
						count++;
					}
				}
			}
		}
	}

	for(int i=0;i<bins;i++)
	{
		distr[i]=distr[i]/(count);
		distr[i]=8*L*L*L*distr[i]/(4*pi*((i+1)*space)*((i+1)*space)*space);
	}
	cout<<"space = "<<space<<endl;
	cout<<"count = "<<count<<endl;
	cout<<"L = "<<L<<endl;
	cout<<"R = "<<R<<endl;
	return distr;
}

vector<double> EscCheck(protein main, double t)
{
	// should make use of proteins internal storage for velocities and whatnot
}