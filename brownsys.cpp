// implemention files for brownsys class

#include "headers/brownsys.h"


using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
///// CONSTRUCTORS  CONSTRUCTORS  CONSTRUCTORS CONSTRUCTORS........
////////////////////////////////////////////////////////////////////////////////////////

//default, one protein system constructor
brownsys::brownsys()
{

	Ncr=0;
	double xx,yy,zz,mag,mag2; // dummies
	srand (time(NULL));
	//Generate random point on surface of sphere radius b, centered at origin
	xx=2*(( (double) rand() / RAND_MAX + 1)-1.5); 
	yy=2*(( (double) rand() / RAND_MAX + 1)-1.5);
	zz=2*(( (double) rand() / RAND_MAX + 1)-1.5);
	mag2=xx*xx+yy*yy+zz*zz;
	mag=pow(mag2,0.5);
	xx=b*xx/mag;
	yy=b*yy/mag;
	zz=b*zz/mag;
	main = mainP(xx,yy,zz,M,r,rc);
	for(int i=0;i<dim;i++){cvel.push_back(0);cpos.push_back(0);}
	crad=r;
	cmass=M;


}

//N protein system constructor
brownsys::brownsys(int num)
{
	
	Ncr=num;
	double xx,yy,zz,mag,mag2; // dummies
	int u,v,w,ran;
	double n=ceil(pow(num,1.0/3.0));
	int Np=n*n*n; // number of lattice points
	int dN=Np-num; // number of lattice points which should eventually not contain crowders
	double l=2*L/n; // cubic lattice spacing
	vector< int > rem(dN-1);
	events.push_back(false); events.push_back(false); centered=true;

	for(int i=0;i<dim;i++){cvel.push_back(0);cpos.push_back(0);ncpos.push_back(0);}
	crad=r;
	cmass=M;

	
	//initiate random int generator
	random_device rd;
	srand (time(NULL));
	default_random_engine gen(rd());
	uniform_int_distribution<int> dist(1,Np); //random ints between 0 and Np
	//done


	//Generate random point on surface of sphere radius b, centered at origin
	xx=2*(( (double) rand() / RAND_MAX + 1)-1.5); 
	yy=2*(( (double) rand() / RAND_MAX + 1)-1.5);
	zz=2*(( (double) rand() / RAND_MAX + 1)-1.5);
	mag2=xx*xx+yy*yy+zz*zz;
	mag=pow(mag2,0.5);
	xx=b*xx/mag;
	yy=b*yy/mag;
	zz=b*zz/mag;
	main = mainP(xx,yy,zz,1,r,rc); // set mass to one for scale effect
	

	//Find indices of lattice point nearest to protein	
	u=round((xx+L)/l - 0.5);	
	v=round((yy+L)/l - 0.5);
	w=round((zz+L)/l - 0.5);	
	//Indices found

	//Generate indices of lattice points to be removed
	int count=0;
	while(count<dN-1){
		ran=dist(gen);
		
		if (find(rem.begin(), rem.end(), ran) == rem.end()){rem[count]=ran;count++;}
		
	}
	
	//Generated

	int numb=1;
	count=1;
	// START COMPLICATED IF PROCESS TO AVOID PROTEIN
		while(numb<N){
			for(double k=0;k<n;k++){// k for
				for(double j=0;j<n;j++){// j for
	
					for(double i=0;i<n;i++){// i for
						if (find(rem.begin(), rem.end(), count) == rem.end()){//check
						if(k!=w){// k if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,1,rc));
						numb++;}// end k if
						else if(j!=v){// j if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,1,rc));
						numb++;}// end j if
						else if(i!=u){// i if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,1,rc));
						numb++;}// end i if
					
						}//end check
						count++;
					
				
}}}}
		


}

brownsys::~brownsys()
{
	cout<<"system removed"<<endl;
}



///////////////////////////////////////////////////////////////////////////////////////////
////// MEMBERS MEMBERS MEMBERS MEMBERS MEMBERS////////////
////////////////////////////////////////////////////////////////////////////////////////

// need getr function for protein class
void brownsys::startNNs()
{
	
	vector<int> nns;
	vector<int> nnns;
	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2;
	dumby=main.getv("coords");
	xx=dumby[0];yy=dumby[1];zz=dumby[2];
	bool NN;	
	//NNs of main
	for(int u=0; u<Ncr; u++)
	{	
		dumby=crowders[u].getv("coords");
		x=dumby[0];y=dumby[1];z=dumby[2]; 
		mag2=x*x+y*y+z*z;
		mag=pow(mag2,0.5);
		if(mag<cut){nearcntr.push_back(u);} //NNs of central protein		
		

		
		NN=false;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		
		
		// CAN IMPROVE THIS BY NESTING ELSE [DONE, ELSES NESTED]		

		if(mag<cut){nns.push_back(u);NN=true;} // main will never start near a boundary
	}
	
	//sort(nearcntr.begin(),nearcntr.end()); // this is already sorted
	main.setNNs(nns, nnns);
	
	
	for(int u=0; u<Ncr; u++) //can improve by linking. UNFINISHED
	{
		vector<int> nns;
		vector<int> nnns;
		dumby=crowders[u].getv("coords");
		xx=dumby[0];yy=dumby[1];zz=dumby[2];
		for(int v=0; v<Ncr; v++){NN=false;
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


void brownsys::updateNNs()
{
	vector<double> pos1, pos2;
	double dx,dy,dz,mag2,mag,xx,yy,zz;
	vector<int> NNs;
	vector<int> nNNs;
	int index;
	int index2;
	//update main
	pos1=main.getv("coords");
	NNs=main.getNNs(true);
	nNNs=main.getNNs(false);
	int NNnum=NNs.size();
	int cnum=nearcntr.size();
	int nNNnum=nNNs.size();
	
	if(cnum!=0)
	{int i=0;
	while(i<cnum) //Updates near center list
	{	
		pos2=crowders[nearcntr[i]].getv("coords");
		xx=pos2[0];yy=pos2[1];zz=pos2[2];
		mag2=xx*xx+yy*yy+zz*zz;
		mag=pow(mag2,0.5);
		if(mag>cut){nearcntr.erase(nearcntr.begin()+i);cnum-=1;}
		else{i++;}
	}
	sort (nearcntr.begin(),nearcntr.end());
	}
	
	cnum=nearcntr.size();
	
	int count=0;
	for(int i=0;i<Ncr;i++) 
	{
		if(cnum==0)
		{
			
			pos2=crowders[i].getv("coords");
			mag2=pos2[0]*pos2[0]+pos2[1]*pos2[1]+pos2[2]*pos2[2];
			mag=pow(mag2,0.5);
			if(mag<cut){nearcntr.push_back(i);}
		}
		else if(i != nearcntr[count]) //assumes nearcntr is sorted
		{
			pos2=crowders[i].getv("coords");
			mag2=pos2[0]*pos2[0]+pos2[1]*pos2[1]+pos2[2]*pos2[2];
			mag=pow(mag2,0.5);
			if(mag<cut){nearcntr.push_back(i);}
		}
		else
		{
			count++;
			if(count==cnum){cnum=0;}
		}	
	}
	sort(nearcntr.begin(),nearcntr.end());

	
	for(int i=0; i<NNnum; i++) // Main removal loop
	{
		if(!main.IsNN(crowders[NNs[i]%Ncr], NNs[i]%Ncr, true, i))
		{
			nNNs.push_back(NNs[i]%Ncr);
			NNs.erase(NNs.begin()+i);
			NNnum-=1; i-=1;
		}
	} // ends the NN removal for loop		

	for(int i=0; i<nNNnum; i++) //CAN MAKE THIS BETTER BY LINKING PARTICLES FOR WHICH IVE ALREADY FOUND NNS
	{
		if(main.IsNN(crowders[nNNs[i]],nNNs[i],false,0))
		{
			nNNs.erase(nNNs.begin()+i);
			i-=1; nNNnum-=1;
		}

	} // ends NN addition loop for MAIN

	main.setNNs(nNNs);

	//////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// NN update for every crowder
	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////


	for(int j=0; j<Ncr; j++){ //STARTS REMOVAL LOOP

	pos1=crowders[j].getv("coords");
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
} // ends NNupdate

void brownsys::moveall(mt19937& gen, normal_distribution<> distro, int& betacount)
{
	// parameters
	double t_el=h; // elapsed time
	double t_rem=h; // remaining time

	// OUTS
	vector<double> vi_out, vf_out, posi_out, posf_out;

	// Vectors
	vector<int> empty;
	
	int colIndex, counter=0; 
	
	while(!events[0] && !events[1])
	{
		if((counter)%10==0)
		{
			updateNNs();
			cout<<"ten steps"<<endl;
		}
		// Move main and check escape/reaction/collisions
		main.move(gen, distro);
		if(mainColCheck(colIndex, t_el))
		{
			cout<<"collision"<<endl;
			mainRes(t_el, colIndex, 0, empty);	
		}
		if(!events[0] && !events[1])
		{
			for(int j=0; j<Ncr; j++) // move all crowders
			{
				moveCrowd(j, gen, distro);
			}
			for(int j=0;j<3;j++) // move center
			{
				ncpos[j]=distro(gen);
				cvel[j]=ncpos[j]/h;	
			}
			centered=false;
			cColCheck();
		}
		counter++;
	}
	if(events[0])
	{
		cout<<"reaction"<<endl;
		cout<<counter<<endl;
		betacount++;
	}
	if(events[1])
	{
		cout<<"escape"<<endl;
		cout<<counter<<endl;
	}	
	
}// end moveall

/*void brownsys::mainCheck(int index, double& t_el)
{
	vector<double> pos1=crowders[index].getv
}*/

void brownsys::cColCheck()
{
	vector<double> pos2;
	double t_el=h, mag2, mag, R, vdr, r2, rad2, rad, v2,tm,tp;
	bool col=false, esc=false, reac=false;
	vector<int> moving_particles(1);
	for(int i=0; i<nearcntr.size(); i++)
	{
		vdr=0, r2=0, v2=0;
		pos2=crowders[nearcntr[i]].getv("coords");
		R=crad+crowders[nearcntr[i]].getp("radius");
		for(int k=0;k<3;k++)
		{
			mag2+=(pos2[k]-ncpos[k])*(pos2[k]-ncpos[k]);
		}
		if(pow(mag2,0.5)<R)
		{
			for(int k=0;k<3;k++)
			{
				r2+=(pos2[k]-cpos[k])*(pos2[k]-cpos[k]);
				vdr+=(pos2[k]-cpos[k])*cvel[k];
				v2+=cvel[k]*cvel[k];
			}
			rad2=vdr*vdr-v2*(r2-R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tm=(vdr-rad)/v2;
				tp=(vdr+rad)/v2;
				if(tm>0 && tm<t_el)
				{
					col=true;
					esc=false;
					reac=false;
					t_el=tm;
					moving_particles[0]=nearcntr[i];
				}
				else if(tp>0 && tp<t_el)
				{
					col=true;
					esc=false;
					reac=false;
					t_el=tp;
					moving_particles[0]=nearcntr[i];
				}
			}

		}
	}
	if(col)
	{
		mainRes(t_el, moving_particles[0], 4, moving_particles);
	}
}


void brownsys::moveCrowd(int j ,mt19937& gen, normal_distribution<> distro)
{
	crowders[j].move(gen, distro);
	vector<double> pos1=crowders[j].getv("new coords"), pos2;
	vector<int> nns=crowders[j].getNNs(true);
	double R1=crowders[j].getp("radius"), R2, mag2=0, mag;
	bool clear=true;
	bool collision=false;

	for(int i=0;i<nns.size();i++)
	{
		pos2=crowders[nns[i]%Ncr].PBCswitch(Ncr, nns[i]);
		R2=crowders[nns[i]%Ncr].getp("radius");
		for(int k=0;k<3;k++)
		{
			mag2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
		}
		mag=pow(mag2,0.5);
		if(mag<R1+R2)
		{
			clear=false;
		}
	}
	if(clear) //check main & center
	{
		crColCheck(j);
	}

}

void brownsys::crColCheck(int j)
{
	bool mcol=false, ccol=false, col=false;
	vector<double> pos1=crowders[j].getv("new coords");
	vector<double> pos=crowders[j].getv("coords");
	vector<double> mpos=main.getv("coords");
	double R1=crowders[j].getp("radius"), Rm=main.getp("radius");
	double rm2=0, rc2=0, rm, Rc, t_el=h, t_rem=h, vdr=0, Rsq=0, v2=0, dr2=0;
	double tmain, tcen, rad2, rad, R, tm, tp;
	vector<double> v1;
	vector<int>moving_particles;
	moving_particles.push_back(j);

	for(int k=0;k<3;k++)
	{
		rm2+=(mpos[k]-pos1[k])*(mpos[k]-pos1[k]);
		rc2+=(cpos[k]-pos1[k])*(cpos[k]-pos1[k]);
	}
	if(pow(rm2,0.5)<R1+Rm)
	{	
		v1=crowders[j].getv("velocity");
		mcol=true;
		R=Rm+R1;
		for(int k=0; k<3;k++)
		{
			vdr+=v1[k]*(mpos[k]-pos[k]);
			v2+=v1[k]*v1[k];
			dr2+=(mpos[k]-pos[k])*(mpos[k]-pos[k]);
		}
		rad2=vdr*vdr - v2*(dr2-R*R);
		if(rad2>0)
		{
			rad=pow(rad2,0.5);
			tp=(vdr+rad)/v2;
			tm=(vdr-rad)/v2;
			if(tm>0 && tm<t_el)
			{
				t_el=tm;
				col=true;
			}
			else if(tp>0 && tp<t_el)
			{
				t_el=tp;
				col=true;
			}
		}
	}
	if(pow(rc2,0.5)<R1+crad)
	{
		v1=crowders[j].getv("velocity");
		ccol=true;
		R=crad+R1;
		for(int k=0; k<3;k++)
		{
			vdr+=v1[k]*(mpos[k]-pos[k]);
			v2+=v1[k]*v1[k];
			dr2+=(mpos[k]-pos[k])*(mpos[k]-pos[k]);
		}
		rad2=vdr*vdr - v2*(dr2-R*R);
		if(rad2>0)
		{
			rad=pow(rad2,0.5);
			tp=(vdr+rad)/v2;
			tm=(vdr-rad)/v2;
			if(tm>0 && tm<t_el)
			{
				t_el=tm;
				mcol=false;
				col=true;
			}
			else if(tp>0 && tp<t_el)
			{
				t_el=tp;
				mcol=false;
				col=true;
			}
		}
	}
	if(mcol)
	{
		mainRes(t_el, -1, 2, moving_particles);
	}
	else if(ccol)
	{
		mainRes(t_el, -1, 3, moving_particles);
	}
}

void brownsys::mainRes(double& t_el, int index, int colSwitch, vector<int> moving_particles)
{
	// first move initial particle to its pre-collision location
	

	// some initial params
	vector<double> pos1(3), pos2(3), vi(3), v1f(3), v2f(3), vt(3), rhat(3);
	double m1, m2, R1, R2, R;
	double t_rem=h, vni, v1nf, v2nf;
	int index2;

	vector<bool> condition(4);
	if(colSwitch==0)
	{
		condition[0]=false;condition[1]=true;
		condition[2]=true;condition[3]=false;
		main.nudge(t_el);
	}
	else if(colSwitch==2)
	{
		condition[0]=false;condition[1]=true;
		condition[2]=true;condition[3]=false;
		crowders[moving_particles[0]].nudge(t_el);
		index2=moving_particles[0];
	}
	else if(colSwitch==3)
	{
		condition[0]=false;condition[1]=false;
		condition[2]=true;condition[3]=true;
		crowders[moving_particles[0]].nudge(t_el);
		index2=moving_particles[0];
	}
	else if(colSwitch==4)
	{
		condition[0]=false;condition[1]=false;
		condition[2]=true;condition[3]=true;
		for(int i=0;i<3;i++){cpos[i]+=cvel[i]*t_el;} // nudge center
	}
	int maincount=0;

	///////////////////////////////////////////////////////////
	// Tracking conditions
	// [0] - Resolved, initial = false
	// [1] - Tracking main, initial = true
	// [2] - Tracking crowders, initial = true
	// [3] - Tracking center, initial = false
	///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////
	// colSwitch values
	// 0=main collision -> INITIAL
	// 1=crowder collisions
	// 2=crowder collision with main
	// 3=crowder collision with center
	// 4=center collision with crowder
	///////////////////////////////////////////////////////////

	while(!condition[0])
	{
		vni=0;
		if(!events[0] && !events[1])
		{
			switch (colSwitch)
			{
				case 0:  // main collision	
					//cout<<"m-cr"<<endl;
					maincount++;
					pos1=main.getv("coords");
					pos2=crowders[index%Ncr].PBCswitch(Ncr, index);
					vi=main.getv("velocity");
					m1=main.getp("mass"), m2=crowders[index%Ncr].getp("mass");
					R1=main.getp("radius"), R2=crowders[index%Ncr].getp("radius");
					R=R1+R2;

					/////Velocity Calculation
					for(int k=0; k<3; k++)
					{
						rhat[k]=(pos2[k]-pos1[k])/R;
						vni+=rhat[k]*vi[k];
					}
					v1nf = vni*(m1-m2)/(m1+m2);
					v2nf = vni*2*m1/(m1+m2);
					for(int k=0;k<3;k++)
					{
						vt[k]=vi[k]-vni*rhat[k];
						v1f[k]=vt[k]+v1nf*rhat[k];
						v2f[k]=v2nf*rhat[k];
					}

					////// Set new velocities
					main.newvel(v1f);
					crowders[index%Ncr].newvel(v2f);
					t_rem-=t_el;
					t_el=t_rem;
					moving_particles.push_back(index%Ncr);	

					break;
				
				case 1: // crowd-crowd collision
					//cout<<"cr-cr"<<endl;
					pos1=crowders[index2].getv("coords");
					pos2=crowders[index%Ncr].PBCswitch(Ncr, index);
					vi=crowders[index2].getv("velocity");
					m1=crowders[index2].getp("mass"), m2=crowders[index%Ncr].getp("mass");
					R1=crowders[index2].getp("radius"), R2=crowders[index%Ncr].getp("radius");
					R=R1+R2;

					/////Velocity Calculation
					for(int k=0; k<3; k++)
					{
						rhat[k]=(pos2[k]-pos1[k])/R;
						vni+=rhat[k]*vi[k];
					}
					v1nf = vni*(m1-m2)/(m1+m2);
					v2nf = vni*2*m1/(m1+m2);
					for(int k=0;k<3;k++)
					{
						vt[k]=vi[k]-vni*rhat[k];
						v1f[k]=vt[k]+v1nf*rhat[k];
						v2f[k]=v2nf*rhat[k];
					}

					////// Set new velocities
					crowders[index2].newvel(v1f);
					crowders[index%Ncr].newvel(v2f);
					t_rem-=t_el;
					t_el=t_rem;
					moving_particles.push_back(index%Ncr);	

					break;
					
				case 2: // crowd-main collision
					//cout<<"cr-m"<<endl;
					pos1=crowders[index2].getv("coords");
					pos2=main.getv("coords");
					vi=crowders[index2].getv("velocity");
					m1=crowders[index2].getp("mass"), m2=main.getp("mass");
					R1=crowders[index2].getp("radius"), R2=main.getp("radius");
					R=R1+R2;

					/////Velocity Calculation
					for(int k=0; k<3; k++)
					{
						rhat[k]=(pos2[k]-pos1[k])/R;
						vni+=rhat[k]*vi[k];
					}
					v1nf = vni*(m1-m2)/(m1+m2);
					v2nf = vni*2*m1/(m1+m2);
					for(int k=0;k<3;k++)
					{
						vt[k]=vi[k]-vni*rhat[k];
						v1f[k]=vt[k]+v1nf*rhat[k];
						v2f[k]=v2nf*rhat[k];
					}

					////// Set new velocities
					crowders[index2].newvel(v1f);
					main.newvel(v2f);
					t_rem-=t_el;
					t_el=t_rem;	

					break;
				
				case 3: // crowd-center collision
					//cout<<"cr-cn"<<endl;
					pos1=crowders[index2].getv("coords");
					pos2=cpos;
					vi=crowders[index2].getv("velocity");
					m1=crowders[index2].getp("mass"), m2=cmass;
					R1=crowders[index2].getp("radius"), R2=crad;
					R=R1+R2;

					/////Velocity Calculation
					for(int k=0; k<3; k++)
					{
						rhat[k]=(pos2[k]-pos1[k])/R;
						vni+=rhat[k]*vi[k];
					}
					v1nf = vni*(m1-m2)/(m1+m2);
					v2nf = vni*2*m1/(m1+m2);
					for(int k=0;k<3;k++)
					{
						vt[k]=vi[k]-vni*rhat[k];
						v1f[k]=vt[k]+v1nf*rhat[k];
						v2f[k]=v2nf*rhat[k];
					}

					////// Set new velocities
					crowders[index2].newvel(v1f);
					cvel=v2f;
					t_rem-=t_el;
					t_el=t_rem;

					break;

				case 4: // center-crowd collision
					//cout<<"cn-cr"<<endl;
					pos1=cpos;
					pos2=crowders[index].getv("coords");
					vi=cvel;
					m1=cmass, m2=crowders[index].getp("mass");
					R1=crad, R2=crowders[index].getp("radius");
					R=R1+R2;

					/////Velocity Calculation
					for(int k=0; k<3; k++)
					{
						rhat[k]=(pos2[k]-pos1[k])/R;
						vni+=rhat[k]*vi[k];
					}
					v1nf = vni*(m1-m2)/(m1+m2);
					v2nf = vni*2*m1/(m1+m2);
					for(int k=0;k<3;k++)
					{
						vt[k]=vi[k]-vni*rhat[k];
						v1f[k]=vt[k]+v1nf*rhat[k];
						v2f[k]=v2nf*rhat[k];
					}

					////// Set new velocities
					cvel=v1f;
					crowders[index].newvel(v2f);
					t_rem-=t_el;
					t_el=t_rem;	

					break;
			}// end switch
			
			// advance newcoords of things
			if(condition[1])
			{
				main.mvVel(t_rem);
			}
			if(condition[2])
			{
				for(int k=0;k<moving_particles.size();k++)
				{
					crowders[moving_particles[k]].mvVel(t_rem);
				}
			}
			if(condition[3])
			{
				for(int op=0;op<3;op++)
				{
					ncpos[op]=cpos[op]+cvel[op]*t_rem;
				}
				centered=false;
			}

			if(ColCheckAll(colSwitch, index, index2, moving_particles, t_el, condition))
			{
				// nudge all tracked particles to t_el
				if(condition[1]) // advance main
				{
					main.nudge(t_el);
				} 
				if(condition[2]) // advance crowders
				{
					for(int k=0;k<moving_particles.size();k++)
					{
						crowders[moving_particles[k]].nudge(t_el);
					}
				}
				if(condition[3]) // advance center
				{
					for(int i=0;i<3;i++)
					{
						cpos[i]+=cvel[i]*t_el;
					}
				}
			}
			else if(!centered)
			{
				shftcntr();
			}
		}
		else
		{
			condition[0]=true;
		}
	}
}

void brownsys::upall(vector<int> moving_particles)
{
	main.upmain();
			
	for(int k=0;k<moving_particles.size();k++)
	{
		crowders[moving_particles[k]].update();
		// make update fix PBC coordinates
	}

	/*
	if(!centered)
	{
		center that bitch
	}
	*/
	
}

void brownsys::shftcntr()
{
	for(int i=0;i<crowders.size();i++)
	{
		crowders[i].shift(cpos, crowders);
	}
	main.shift(cpos, crowders);
	for(int j=0;j<3;j++){cpos[j]=0;ncpos[j]=0;}
	centered=true;
}

bool brownsys::ColCheckAll(int& switcher, int& index1, int& index2, vector<int> moving_particles, double& t_el, vector<bool>& condition)
{
	vector<int> nns;
	vector<double> pos2, pos1_o, pos1;
	vector<double> vel;

	double v2, r2, r2_o, vdr, R1, R2, R, rad, rad2, tm, tp;
	// dumby bools for collisions
	bool maincol=false, cCol=false, crCol=false, crCol2=false, col=false;

	int num=moving_particles.size();
	

	// check for reaction/escape
	if(condition[1])
	{
		vel=main.getv("velocity");
		pos1_o=main.getv("coords");
		pos1=main.getv("new coords");
		nns=main.getNNs(true);
		R1=main.getp("radius");
		if(main.nearcntr())
		{
			v2=0, r2=0, r2_o=0, vdr=0;
			R2=crad;
			R=R1+R2;
			for(int k=0;k<3;k++)
			{
				r2_o+=(cpos[k]-pos1_o[k])*(cpos[k]-pos1_o[k]);
				r2+=(cpos[k]-pos1[k])*(cpos[k]-pos1[k]);
				v2+=vel[k]*vel[k];
				vdr+=(cpos[k]-pos1_o[k])*vel[k];
			}
			if(pow(r2,0.5)<R)
			{
				rad2=vdr*vdr - v2*(r2_o - R*R);
				if(rad2>0)
				{
					rad=pow(rad2,0.5);
					tm=(vdr-rad)/v2;
					tp=(vdr+rad)/v2;
					if(tm>0 && tm<t_el)
					{
						t_el=tm;
						events[0]=true;
					}
					else if(tp<t_el)
					{
						t_el=tp;
						events[0]=true;
					}
				}
			}
		} //reac check
		else if(main.nearEsc())
		{
			v2=0, r2=0, r2_o=0, vdr=0;
			R=q;
			for(int k=0;k<3;k++)
			{
				r2_o+=(pos1_o[k])*(pos1_o[k]);
				r2+=(pos1[k])*(pos1[k]);
				v2+=vel[k]*vel[k];
				vdr+=(pos1_o[k])*vel[k];
			}
			if(pow(r2,0.5)>R)
			{
				rad2=vdr*vdr - v2*(r2_o - R*R);
				if(rad2>0)
				{
					rad=pow(rad2,0.5);
					tm=(vdr-rad)/v2;
					tp=(vdr+rad)/v2;
					if(tm>0 && tm<t_el)
					{
						t_el=tm;
						events[1]=true;
					}
					else if(tp<t_el)
					{
						t_el=tp;
						events[1]=true;
					}
				}
			}
		} // esc check
		
		for(int i=0; i<nns.size();i++) // main-cr collisions
		{	
			v2=0, r2=0, r2_o=0, vdr=0;
			R2=crowders[nns[i]%Ncr].getp("radius");
			R=R1+R2;
			pos2=crowders[nns[i]%Ncr].PBCswitch(Ncr, nns[i]);
			for(int k=0;k<3;k++)
			{
				r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
				r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
				v2+=vel[k]*vel[k];
				vdr+=(pos2[k]-pos1_o[k])*vel[k];
			}

			if(r2<R*R)
			{
				rad2=vdr*vdr - v2*(r2_o - R*R);
				if(rad2>0)
				{
					rad=pow(rad2,0.5);
					tm=(vdr-rad)/v2;
					tp=(vdr+rad)/v2;
					if(tm>0 && tm<t_el)
					{
						t_el=tm;
						col=true;
						maincol=true;
						index1=nns[i];
						events[0]=false;
						events[1]=false;
						switcher=0;
					}
					else if(tp<t_el)
					{
						t_el=tp;
						col=true;
						maincol=true;
						index1=nns[i];
						events[0]=false;
						events[1]=false;
						switcher=0;
					}
				}
			}
		}
		if(!maincol)
		{
			condition[1]=false;
			main.upmain();
		}
	} // end tracking main check
	if(!events[0] && !events[1]) //only if there is no escape or reaction
	{
		if(condition[2]) // tracking crowders
		{
			for(int i=0; i<num; i++)
			{
				crCol=false;
				vel=crowders[moving_particles[i]].getv("velocity");
				pos1_o=crowders[moving_particles[i]].getv("coords");
				pos1=crowders[moving_particles[i]].getv("new coords");
				nns=crowders[moving_particles[i]].getNNs(true);
				R1=crowders[moving_particles[i]].getp("radius");
				// check for collisions with other crowders
				for(int j=0; j<nns.size();j++)
				{	
					v2=0, r2=0, r2_o=0, vdr=0;
					R2=crowders[nns[j]%Ncr].getp("radius");
					R=R1+R2;
					pos2=crowders[nns[j]%Ncr].PBCswitch(Ncr, nns[j]);
					for(int k=0;k<3;k++)
					{
						r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
						r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
						v2+=vel[k]*vel[k];
						vdr+=(pos2[k]-pos1_o[k])*vel[k];
					}

					if(r2<R*R)
					{
						rad2=vdr*vdr - v2*(r2_o - R*R);
						if(rad2>0)
						{
							rad=pow(rad2,0.5);
							tm=(vdr-rad)/v2;
							tp=(vdr+rad)/v2;
							if(tm>0 && tm<t_el)
							{
								t_el=tm;
								col=true;
								crCol=true;
								index1=nns[j];
								index2=moving_particles[i];
								switcher=1; //cr-cr collision
							}
							else if(tp<t_el)
							{
								t_el=tp;
								col=true;
								crCol=true;
								index1=nns[j];
								index2=moving_particles[i];
								switcher=1;
							}
						}
					} // ends collision if
				} // ends nns loop

				//check for collision with main (this may be totally pointless)
				//also, this would mean i need a complete final vel algorithm
				//could throw in a bool for near main but whatever
				v2=0, r2=0, r2_o=0, vdr=0;
				R2=main.getp("radius");
				R=R1+R2;
				pos2=main.getv("coords");
				for(int k=0;k<3;k++)
				{
					r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
					r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
					v2+=vel[k]*vel[k];
					vdr+=(pos2[k]-pos1_o[k])*vel[k];
				}

				if(r2<R*R)
				{
					rad2=vdr*vdr - v2*(r2_o - R*R);
					if(rad2>0)
					{
						rad=pow(rad2,0.5);
						tm=(vdr-rad)/v2;
						tp=(vdr+rad)/v2;
						if(tm>0 && tm<t_el)
						{
							t_el=tm;
							col=true;
							crCol=true;
							index2=moving_particles[i];
							condition[1]=true; // definitely tracking main again... not necessary?
							switcher=2; //cr-main collision
						}
						else if(tp<t_el)
						{
							t_el=tp;
							col=true;
							crCol=true;
							index2=moving_particles[i];
							condition[1]=true;
							switcher=2;
						}
					}
				} // ends main collision if

				// check for collisions with center
				if(crowders[moving_particles[i]].nearcntr())
				{
					v2=0, r2=0, r2_o=0, vdr=0;
					R2=crad;
					R=R1+R2;
					pos2=cpos;
					for(int k=0;k<3;k++)
					{
						r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
						r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
						v2+=vel[k]*vel[k];
						vdr+=(pos2[k]-pos1_o[k])*vel[k];
					}

					if(r2<R*R)
					{
						rad2=vdr*vdr - v2*(r2_o - R*R);
						if(rad2>0)
						{
							rad=pow(rad2,0.5);
							tm=(vdr-rad)/v2;
							tp=(vdr+rad)/v2;
							if(tm>0 && tm<t_el)
							{
								t_el=tm;
								col=true;
								crCol=true;
								index2=moving_particles[i];
								condition[3]=true; // tracking center
								switcher=3; // cr-center collision
							}
							else if(tp<t_el)
							{
								t_el=tp;
								col=true;
								crCol=true;
								index2=moving_particles[i];
								condition[3]=true;
								switcher=3;
							}
						}
					} // ends collision if
				} // ends center collision if
				if(!crCol)
				{
					crowders[moving_particles[i]].update(); // set coords=newcoords;
					moving_particles.erase(moving_particles.begin()+i); // remove from list
					i-=1; // back up i 
					num-=1; // size has decreased
				}
				else
				{
					crCol2=true;
					// if main hits crowder or any crowder hits crowder
				}
			}// ends moving_particles loop
		}// ends tracking crowders loop

		if(condition[3]) // tracking center
		{
			vel=cvel;
			pos1_o=cpos;
			pos1=ncpos;
			R1=crad;
			for(int i=0;i<nearcntr.size();i++)
			{
				v2=0, r2=0, r2_o=0, vdr=0;
				R2=crowders[nearcntr[i]].getp("radius");
				R=R1+R2;
				pos2=crowders[nearcntr[i]].getv("coords");
				for(int k=0;k<3;k++)
				{
					r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
					r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
					v2+=vel[k]*vel[k];
					vdr+=(pos2[k]-pos1_o[k])*vel[k];
				}

				if(r2<R*R)
				{
					rad2=vdr*vdr - v2*(r2_o - R*R);
					if(rad2>0)
					{
						rad=pow(rad2,0.5);
						tm=(vdr-rad)/v2;
						tp=(vdr+rad)/v2;
						if(tm>0 && tm<t_el)
						{
							t_el=tm;
							col=true;
							cCol=true;
							index1=nearcntr[i];
							condition[3]=true; // tracking center
							condition[2]=true; // tracking new crowder
							switcher=4; // cr-center collision
						}
						else if(tp<t_el)
						{
							t_el=tp;
							col=true;
							cCol=true;
							index1=nearcntr[i];
							condition[3]=true;
							condition[2]=true;
							switcher=4;
						}
					}
				} // ends collision if
			} // ends center nns loop
		} // ends tracking center
		if(!cCol)
		{
			cpos=ncpos;
			// update every particles position and check for escape
			condition[3]=false; // stop tracking center
		}
	} // ends if no events

	if(!col) // no collisions at all
	{
		condition[0]=true; // resolved!
	}

	/*condition[0]=resolved;
	condition[1]=cond1;
	condition[2]=cond2;
	condition[3]=cond3;*/
	return col;
}

bool brownsys::mainColCheck(int& index, double& t_el)
{

	vector<int> nns=main.getNNs(true);
	vector<double> pos2(3), pos1_o=main.getv("coords"), pos1=main.getv("new coords");
	vector<double> vel=main.getv("velocity");

	double v2, r2, r2_o, vdr, R1=main.getp("radius"), R2, R, rad, rad2, tm, tp;
	bool col=false;

	if(main.nearcntr())
	{
		
		v2=0, r2=0, r2_o=0, vdr=0;
		R2=crad;
		R=R1+R2;
		for(int k=0;k<3;k++)
		{
			r2_o+=(cpos[k]-pos1_o[k])*(cpos[k]-pos1_o[k]);
			r2+=(cpos[k]-pos1[k])*(cpos[k]-pos1[k]);
			v2+=vel[k]*vel[k];
			vdr+=(cpos[k]-pos1_o[k])*vel[k];
		}
		if(pow(r2,0.5)<R)
		{
			rad2=vdr*vdr - v2*(r2_o - R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tm=(vdr-rad)/v2;
				tp=(vdr+rad)/v2;
				if(tm>0 && tm<t_el)
				{
					t_el=tm;
					events[0]=true;
				}
				else if(tp<t_el)
				{
					t_el=tp;
					events[0]=true;
				}
			}
		}
	}
	else if(main.nearEsc())
	{

		v2=0, r2=0, r2_o=0, vdr=0;
		R=q*q;
		for(int k=0;k<3;k++)
		{
			r2_o+=(pos1_o[k])*(pos1_o[k]);
			r2+=(pos1[k])*(pos1[k]);
			v2+=vel[k]*vel[k];
			vdr+=(pos1_o[k])*vel[k];
		}
		if(pow(r2,0.5)<R)
		{
			rad2=vdr*vdr - v2*(r2_o - R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tm=(vdr-rad)/v2;
				tp=(vdr+rad)/v2;
				if(tm>0 && tm<t_el)
				{
					t_el=tm;
					events[1]=true;
				}
				else if(tp<t_el)
				{
					t_el=tp;
					events[1]=true;
				}
			}
		}
	}

	for(int i=0; i<nns.size();i++)
	{	
		v2=0, r2=0, r2_o=0, vdr=0;
		R2=crowders[nns[i]%Ncr].getp("radius");
		R=R1+R2;
		pos2=crowders[nns[i]%Ncr].PBCswitch(Ncr, nns[i]);
		for(int k=0;k<3;k++)
		{
			r2_o+=(pos2[k]-pos1_o[k])*(pos2[k]-pos1_o[k]);
			r2+=(pos2[k]-pos1[k])*(pos2[k]-pos1[k]);
			v2+=vel[k]*vel[k];
			vdr+=(pos2[k]-pos1_o[k])*vel[k];
		}

		if(pow(r2,0.5)<R)
		{
			rad2=vdr*vdr - v2*(r2_o - R*R);
			if(rad2>0)
			{
				rad=pow(rad2,0.5);
				tm=(vdr-rad)/v2;
				tp=(vdr+rad)/v2;
				if(tm>0 && tm<t_el)
				{
					t_el=tm;
					col=true;
					index=nns[i];
					events[0]=false;
					events[1]=false;
				}
				else if(tp<t_el)
				{
					t_el=tp;
					col=true;
					index=nns[i];
					events[0]=false;
					events[1]=false;
				}
			}
		}
	}
	if(!col)
	{
		main.upmain();
	}
	return col;
}


void brownsys::equilibrate(mt19937& gen, normal_distribution<> distro, int eqsteps)
{
	vector<double> pos1(3), pos2(3), posm(3), disp(3), PBCvec(3);
	vector<int> nns;
	int numnns, count, clashcount=0, counter=0;
	bool clash;
	double mag, mag2, r1, r2;

	posm=main.getv("coords");


	for(int i=0;i<eqsteps;i++) // should be while !equilibrated
	{
	for(int k=0;k<Ncr;k++){counter++;
		
		nns=crowders[k].getNNs(true);
		numnns=nns.size();
		pos1=crowders[k].getv("coords");
				

		for(int j=0;j<3;j++){disp[j]=distro(gen); pos1[j]+=disp[j];} //get displacement
		
		count=0; mag=0; mag2=0;
		clash=false;
		
		//clash with main
		for(int dumb=0;dumb<3;dumb++){mag2+=(pos1[i]-posm[i])*(pos1[i]-posm[i]);}
		mag=pow(mag2,0.5);
		if(mag<crowders[k].getp("radius")+main.getp("radius")){clash=true;}	

		//clash with center
		mag=0;mag2=0;
		for(int dumb=0;dumb<3;dumb++){mag2+=pos1[i]*pos1[i];}
		mag=pow(mag2,0.5);
		if(mag<crowders[k].getp("radius")+crad){clash=true;}	

		//clash with crowder
		while(count<numnns && clash==false)
		{	
						
			mag2=0; mag=0;
			pos2=crowders[nns[count]%Ncr].PBCswitch(Ncr, nns[count]);
			//PBCvec = pbc(Ncr, nns[count]); // removed for new PBCswitch function

			for(int dumb=0;dumb<3;dumb++){mag2+=(pos1[dumb]-pos2[dumb])*(pos1[dumb]-pos2[dumb]);} //
			mag=pow(mag2,0.5);
			if(mag<crowders[k].getp("radius")+crowders[nns[count]%Ncr].getp("radius"))
			{clash=true; clashcount++;}
			else{count++;}
		}
		if(!clash){crowders[k].setpos(pos1);}


	}}//ends eqsteps and Ncr for loops
cout<<"clash count = "<<clashcount<<endl;
cout<<counter<<endl;
}// ends equilibrate function



void brownsys::printcoords(protein test)
{
	vector<double> pos(3);
	pos=test.getv("coords");
	cout<<"x= "<<pos[0]<<endl;
	cout<<"y= "<<pos[1]<<endl;
	cout<<"z= "<<pos[2]<<endl;	
}









	
