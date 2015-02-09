#include "headers/brownsys.h"


using namespace std;

diffusion_sys::diffusion_sys(int s, int track)
	: brownsys()
{
	vector<double> Pos(3);
	double n=ceil(pow(N,1.0/3.0));
	int Np=n*n*n; // number of lattice points
	int dN=Np-N; // number of lattice points which should eventually not contain crowders
	double l=2*L/n; // cubic lattice spacing
	vector< int > rem(dN);

	for(int i=0;i<s;i++)
	{
		stepvec.push_back(0);
	}

	track_index=track;
	steps=s;
	R2=0;

	//initiate random int generator
	int ran;
	random_device rd;
	srand (time(NULL));
	default_random_engine gen(rd());
	uniform_int_distribution<int> dist(1,Np); //random ints between 0 and Np

	int count=0;
	while(count<dN)
	{
		ran=dist(gen);	
		if (find(rem.begin(), rem.end(), ran) == rem.end()){rem[count]=ran;count++;}	
	}
	int numb=0;
	count=0;
	while(numb<N)
	{
		for(double k=0;k<n;k++)
		{// k for 
			for(double j=0;j<n;j++)
			{// j for
				for(double i=0;i<n;i++)
				{// i for
					if (find(rem.begin(), rem.end(), count) == rem.end())
					{//check	
						Pos[0]=l/2.0 + i*l-L;
						Pos[1]=l/2.0 + j*l-L;
						Pos[2]=l/2.0 + k*l-L;
						if(numb<N){crowders[numb].setup(Pos);}
						numb++;
					}//end check
					count++;
				}
			}
		}
	}
}

void diffusion_sys::updateNNs()
{	
	vector<int> NNs;
	vector<int> nNNs;
	int NNnum;
	int nNNnum;
	for(int j=0; j<N; j++)
	{ //STARTS REMOVAL LOOP

		NNs=crowders[j].getNNs(true);
		nNNs=crowders[j].getNNs(false);
		NNnum=NNs.size();
		nNNnum=nNNs.size();

		for(int i=0; i<NNnum; i++) // crowder NN removal loop
		{	
			if(!crowders[j].IsNN(crowders[NNs[i]%Ncr],NNs[i]%Ncr,true,i))
			{
				nNNs.push_back(NNs[i]%Ncr);
				NNs.erase(NNs.begin()+i);
				NNnum-=1;i-=1;
			}
		} // ends the i removal loops	

		for(int i=0; i<nNNnum; i++) //CAN MAKE THIS BETTER BY LINKING PARTICLES FOR WHICH IVE ALREADY FOUND NNS
		{
			if(crowders[j].IsNN(crowders[nNNs[i]],nNNs[i],false,0))
			{
				nNNs.erase(nNNs.begin()+i);
				nNNnum-=1;i-=1;
			}
		} // ends i addition loop
		crowders[j].setNNs(nNNs);
	} // ends j loop for crowder nearest neighbor update	
	cout<<"updated"<<endl;
} // ends NNupdate





void diffusion_sys::dmove(mt19937& gen, normal_distribution<> distro)
{
	vector<double> disp(3), vel(3), pos1(3), pos2(3);
	for(int i=0;i<steps;i++)
	{
		for(int j=0;j<crowders.size();j++)
		{
			crowders[j].move(gen,distro);

			if(j==track_index)
			{
				//track_check();
			}
			else
			{

			}
		}

	}
}

bool track_check()
{
	
}






void diffusion_sys::startNNs()
{
	vector<int> nns;
	vector<int> nnns;
	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2;
	
	bool NN;	
	
	
	
	for(int u=0; u<N; u++) //can improve by linking. UNFINISHED
	{
		vector<int> nns;
		vector<int> nnns;
		dumby=crowders[u].getv("coords");
		xx=dumby[0];yy=dumby[1];zz=dumby[2];
		for(int v=0; v<N; v++){NN=false;
			if(u!=v)
			{
				if(!crowders[u].IsNN(crowders[v],v,false,0))
				{
					nnns.push_back(v);
				}
			}// end if u!=v			
		}// end for v
	
	crowders[u].setNNs(nnns);
	
	}// end for u
	cout<<"NN list created"<<endl;
}// ends function