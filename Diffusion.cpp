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







void diffusion_sys::dmove(mt19937& gen, normal_distribution<> distro)
{
	/*vector<double> disp(3), vel(3), pos1(3), pos2(3);
	double t;
	for(int i=0;i<steps;i++)
	{
		for(int j=0;j<crowders.size();j++)
		{
			crowders[j].move(gen,distro);

			if(j==track_index && track_check(t))
			{
				//track_check();
			}
			else
			{

			}
		}

	}*/
}

/*bool diffusion_sys::track_check(double& t)
{
	vector<double> pos, pos1, pos2, vel;
	vector<int> nns=crowders[track_index].getNNs(true);
	double del_t, mag2=0, mag, r1,r2;
	int num=nns.size();

	pos=crowders[track_index].getv("new coords");
	pos1=crowders[track_index].getv("coords");
	vel=crowders[track_index].getv("velocity");
	r1=crowders[track_index].getp("radius");
	for(int j=0;j<num;j++)
	{
		pos2=crowders[index].getv("coords");
		r2=crowders[index].getp("radius");

		for(int i=0;i<dim;i++)
		{
			mag2+=(pos2[i]-pos[i])*(pos2[i]-pos[i]);
		}
		mag=pow(mag2,0.5);
		if(mag<r1+r2)
		{
			
		}
	}
}*/






void diffusion_sys::startNNs()
{
	vector<double> pos1,pos2;
	vector<int> nns;
	vector<int> nnns;
	vector<int> empty(N);

	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2=0, delta2, d2, dp, dp2, d, costh;
	vector<vector< bool > >  conditions(2, vector<bool>(6,false));
	int switcher=0;	
	double root2=pow(2,0.5);
	double root3=pow(3,0.5);

	vector<bool> cross(3);
	bool isNN;
	

	for(int i=0;i<N;i++)
	{
		pos1=crowders[i].getv("coords");
		nearest.push_back(empty);
		for(int k=0;k<6;k++)
		{
			conditions[0][k]=crowders[i].NearBound(k);
		}
		for(int j=0;j<N;j++)
		{
			isNN=false;
			pos2=crowders[j].getv("coords");
			mag2=0, d2=0;
			for(int k=0;k<6;k++)
			{
				conditions[1][k]=crowders[j].NearBound(k);
			}
			cross[0]=((conditions[0][0]==conditions[1][1])||(conditions[0][1]==conditions[1][0]));
			cross[1]=((conditions[0][2]==conditions[1][3])||(conditions[0][3]==conditions[1][2]));
			cross[2]=((conditions[0][4]==conditions[1][5])||(conditions[0][5]==conditions[1][4]));

			for(int k=0;k<3;k++)
			{
				d2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
			}
			d=pow(d2,0.5);
			if(!cross[0]&&!cross[1]&&!cross[2]) // no boundary
			{
				if(pow(d2,0.5)<cut)
				{
					nearest[i][j]=1;
					isNN=true;
				}
			}
			else if(cross[0]&&!cross[1]&&!cross[2]) // across x
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1])+(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=2;
					isNN=true;
				}

			}
			else if(!cross[0]&&cross[1]&&!cross[2]) // y
			{
				delta2=(pos2[0]-pos1[0])*(pos2[0]-pos1[0])+(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=3;
					isNN=true;
				}
			}
			else if(!cross[0]&&!cross[1]&&cross[2]) //z
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1])+(pos2[0]-pos1[0])*(pos2[0]-pos1[0]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=4;
					isNN=true;
				}
			}
			else if(cross[0]&&cross[1]&&!cross[2]) //xy
			{
				delta2=(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=5;
					isNN=true;
				}
			}
			else if(cross[0]&&!cross[1]&&cross[2]) //xz
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=6;
					isNN=true;
				}
			}
			else if(!cross[0]&&cross[1]&&cross[2]) //yz
			{
				delta2=(pos2[0]-pos1[0])*(pos2[0]-pos1[0]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=7;
					isNN=true;
				}
			}
			else 
			{
				dp=root3*2*L-d;
				if(dp<cut)
				{
					nearest[i][j]=8;
					isNN=true;
				}
			}
		}
		if(isNN)
		{
			crowders[i].setNNs(nnns);
			cout<<"there were nearest neighbors"<<endl;
		}
	}
	
cout<<"NN matrix created"<<endl;
for(int i=0;i<10;i++)
{
	for(int j=0;j<10;j++)
	{
		cout<<nearest[i][j]<<" ";
	}
	cout<<endl;
}
}// ends function

void diffusion_sys::updateNNs()
{	
	vector<double> pos1,pos2;
	vector<int> nnns;

	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2=0, delta2, d2, dp, dp2, d, costh;
	vector<vector< bool > >  conditions(2, vector<bool>(6,false));
	double root2=pow(2,0.5);
	double root3=pow(3,0.5);

	vector<bool> cross(3);
	bool isNN, anyNNs;
	

	for(int i=0;i<N;i++)
	{
		anyNNs=false;
		pos1=crowders[i].getv("coords");
		for(int k=0;k<6;k++)
		{
			conditions[0][k]=crowders[i].NearBound(k);
		}
		for(int j=0;j<N;j++)
		{
			isNN=false;
			pos2=crowders[j].getv("coords");
			mag2=0, d2=0;
			for(int k=0;k<6;k++)
			{
				conditions[1][k]=crowders[j].NearBound(k);
			}
			cross[0]=((conditions[0][0]==conditions[1][1])||(conditions[0][1]==conditions[1][0]));
			cross[1]=((conditions[0][2]==conditions[1][3])||(conditions[0][3]==conditions[1][2]));
			cross[2]=((conditions[0][4]==conditions[1][5])||(conditions[0][5]==conditions[1][4]));

			for(int k=0;k<3;k++)
			{
				d2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
			}
			d=pow(d2,0.5);
			if(!cross[0]&&!cross[1]&&!cross[2]) // no boundary
			{
				if(pow(d2,0.5)<cut)
				{
					nearest[i][j]=1;
					isNN=true;
				}
			}
			else if(cross[0]&&!cross[1]&&!cross[2]) // across x
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1])+(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=2;
					isNN=true;
				}

			}
			else if(!cross[0]&&cross[1]&&!cross[2]) // y
			{
				delta2=(pos2[0]-pos1[0])*(pos2[0]-pos1[0])+(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=3;
					isNN=true;
				}
			}
			else if(!cross[0]&&!cross[1]&&cross[2]) //z
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1])+(pos2[0]-pos1[0])*(pos2[0]-pos1[0]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+4*L*L-4*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=4;
					isNN=true;
				}
			}
			else if(cross[0]&&cross[1]&&!cross[2]) //xy
			{
				delta2=(pos2[2]-pos1[2])*(pos2[2]-pos1[2]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=5;
					isNN=true;
				}
			}
			else if(cross[0]&&!cross[1]&&cross[2]) //xz
			{
				delta2=(pos2[1]-pos1[1])*(pos2[1]-pos1[1]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=6;
					isNN=true;
				}
			}
			else if(!cross[0]&&cross[1]&&cross[2]) //yz
			{
				delta2=(pos2[0]-pos1[0])*(pos2[0]-pos1[0]);
				costh=pow(1-(delta2/d2),0.5);
				dp2=d2+8*L*L-4*root2*L*d*costh;
				if(pow(dp2,0.5)<cut)
				{
					nearest[i][j]=7;
					isNN=true;
				}
			}
			else 
			{
				dp=root3*2*L-d;
				if(dp<cut)
				{
					nearest[i][j]=8;
					isNN=true;
				}
			}
			if(!isNN)
			{
				nearest[i][j]=0;
			}
			else
			{
				anyNNs=true;
			}
		}
		if(anyNNs)
		{
			crowders[i].setNNs(nnns);
			cout<<"there were nearest neighbors"<<endl;
		}

	}
	
cout<<"NN matrix updated"<<endl;
cout<<nearest.size()<<endl;
} // ends NNupdate