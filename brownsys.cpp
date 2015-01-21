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

	for(int i=0;i<dim;i++){cvel.push_back(0);cpos.push_back(0);}
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
	main = mainP(xx,yy,zz,M,r,rc);
	

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
						crowders.push_back(protein(xx,yy,zz,Mc,rc));
						numb++;}// end k if
						else if(j!=v){// j if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,Mc,rc));
						numb++;}// end j if
						else if(i!=u){// i if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,Mc,rc));
						numb++;}// end i if
					
						}//end check
						count++;
					
				
}}}}
		


}



///////////////////////////////////////////////////////////////////////////////////////////
////// MEMBERS MEMBERS MEMBERS MEMBERS MEMBERS////////////
////////////////////////////////////////////////////////////////////////////////////////

// need getr function for protein class
void brownsys::startNNs(double cut)
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

		if(mag<cut){nns.push_back(u);NN=true;}
		else{
			x-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+Ncr);NN=true;}
			else{
			x+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+2*Ncr);NN=true;}
			else{
			x-=2*L;
			y-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+3*Ncr);NN=true;}
			else{
			y+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+4*Ncr);NN=true;}
			else{
			y-=2*L;
			z-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+5*Ncr);NN=true;}
			else{
			z+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(u+6*Ncr);NN=true;}
			if(!NN){nnns.push_back(u);}
			}}}}}}//ends all 6 elses
	}
	
	sort(nearcntr.begin(),nearcntr.end());
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
			dumby=crowders[v].getv("coords");
			x=dumby[0];y=dumby[1];z=dumby[2];
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);

			if(mag<cut){nns.push_back(v);NN=true;}
			else{
			x-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+Ncr);NN=true;}
			else{
			x+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+2*Ncr);NN=true;}
			else{
			x-=2*L;
			y-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+3*Ncr);NN=true;}
			else{
			y+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+4*Ncr);NN=true;}
			else{
			y-=2*L;
			z-=2*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+5*Ncr);NN=true;}
			else{
			z+=4*L;
			mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
			mag=pow(mag2,0.5);
			if(mag<cut){nns.push_back(v+6*Ncr);NN=true;}
			}}}}}} //end all 6 elses
			if(!NN){nnns.push_back(u);}
			
		}// end if u!=v			
		}// end for v
	
	crowders[u].setNNs(nns,nnns);
	
	}// end for u
	
}// ends function


void brownsys::updateNNs(double cut)
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
		
		index=NNs[i];
		
		if(index<Ncr)
		{
			pos2=crowders[index].getv("coords");
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=Ncr && index<2*Ncr) //-x border
		{
			index2=index-Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=2*Ncr && index<3*Ncr) //+x border
		{
			index2=index-2*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=3*Ncr && index<4*Ncr) //-y border
		{
			index2=index-3*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=4*Ncr && index<5*Ncr) //+y border
		{
			index2=index-4*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=5*Ncr && index<6*Ncr) //-z border
		{
			index2=index-5*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=6*Ncr && index<7*Ncr) //+z border
		{
			index2=index-6*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
	
	
	} // ends the NN removal for loop		

	for(int i=0; i<nNNnum; i++) //CAN MAKE THIS BETTER BY LINKING PARTICLES FOR WHICH IVE ALREADY FOUND NNS
	{
		index=nNNs[i];
		pos2=crowders[index].getv("coords");
		dx=pos1[0]-pos2[0];
		dy=pos1[1]-pos2[1];
		dz=pos1[2]-pos2[2];
		if(dx<cut)
		{if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//all ifs ended
		else if(dx<cut-2*L) // -x boundary
		{dx+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-x ifs ended
		else if(dx>2*L-cut) // +x boundary
		{dx-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+2*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+x ifs ended
		else if(dy<cut-2*L) // -y boundary
		{dy+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+3*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-y ifs ended
		else if(dy>2*L-cut) // +y boundary
		{dy-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+4*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+y ifs ended
		else if(dz<cut-2*L) // -z boundary
		{dz+=2*L; if(pow(dz*dz+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+5*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-z ifs ended
		else if(dz>2*L-cut) // +z boundary
		{dz-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+6*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+z ifs ended

	} // ends NN addition loop for MAIN

main.setNNs(NNs, nNNs);

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

	for(int i=0; i<NNnum; i++) // Main removal loop
	{
		
		index=NNs[i];
		
		if(index<Ncr)
		{
			pos2=crowders[index].getv("coords");
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=Ncr && index<2*Ncr) //-x border
		{
			index2=index-Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=2*Ncr && index<3*Ncr) //+x border
		{
			index2=index-2*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=3*Ncr && index<4*Ncr) //-y border
		{
			index2=index-3*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=4*Ncr && index<5*Ncr) //+y border
		{
			index2=index-4*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=5*Ncr && index<6*Ncr) //-z border
		{
			index2=index-5*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		else if(index>=6*Ncr && index<7*Ncr) //+z border
		{
			index2=index-6*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;i-=1;}
		}
		
		
	} // ends the i removal loops	




for(int i=0; i<nNNnum; i++) //CAN MAKE THIS BETTER BY LINKING PARTICLES FOR WHICH IVE ALREADY FOUND NNS
	{
		index=nNNs[i];
		pos2=crowders[index].getv("coords");
		dx=pos1[0]-pos2[0];
		dy=pos1[1]-pos2[1];
		dz=pos1[2]-pos2[2];
		if(dx<cut)
		{if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//all ifs ended
		else if(dx<cut-2*L) // -x boundary
		{dx+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-x ifs ended
		else if(dx>2*L-cut) // +x boundary
		{dx-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+2*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+x ifs ended
		else if(dy<cut-2*L) // -y boundary
		{dy+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+3*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-y ifs ended
		else if(dy>2*L-cut) // +y boundary
		{dy-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+4*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+y ifs ended
		else if(dz<cut-2*L) // -z boundary
		{dz+=2*L; if(pow(dz*dz+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+5*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//-z ifs ended
		else if(dz>2*L-cut) // +z boundary
		{dz-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+6*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;i-=1;}
		}}//+z ifs ended

	} // ends i addition loop

} // ends j loop for crowder nearest neighbor update







	
}

void brownsys::moveall(mt19937& gen, normal_distribution<> distro, int steps)
{

	// parameters
	double t_el=h; // elapsed time
	double t_rem=h; // remaining time

	// OUTS
	vector<double> vi_out, vf_out, posi_out, posf_out;

	// Vectors
	vector<double> pos1=main.getv("coords"), pos2;
	int colIndex; 
	for(int i=0; i<steps; i++)
	{
		// Move main and check escape/reaction/collisions
		main.move(gen, distro);

		if(mainColCheck(colIndex, t_el))
		{
			cout<<"COLLISION"<<endl;
			cout<<"E_o = ";
			main.energy();
			mainRes(t_el, colIndex);
			cout<<"E_f1 = ";
			main.energy();
			cout<<"E_f2 = ";
			crowders[colIndex%Ncr].energy();
		}
		else
		{
			main.update();
		}
	}
	pos2=main.getv("coords");
	cout<<pos1[0]<<" "<<pos1[1]<<" "<<pos1[2]<<" "<<endl;
	cout<<pos2[0]<<" "<<pos2[1]<<" "<<pos2[2]<<" "<<endl;

		
	/* 
		* index = collision check -> also moves
		* collision resolve, including secondary
		

	// move center and check same

	// move crowders and check for collisions with main and center
	*/
	
	

	
	
}// end moveall


void brownsys::mainRes(double& t_el, int index)
{
	// OUTS
	vector<double> vi_out, vf_out,v2_out, pos1_out, pos2_out;
	int indexout=index;
	vi_out=main.getv("velocity");
	//v2_out=crowders[index%Ncr].PBCswitch(Ncr, index);
	pos1_out=main.getv("coords");
	// first move main to its pre-collision location
	main.nudge(t_el);
	pos1_out=main.getv("coords");
	pos2_out=crowders[index%Ncr].PBCswitch(Ncr, index);
	cout<<"r1 = ("<<pos1_out[0]<<","<<pos1_out[1]<<","<<pos1_out[2]<<")"<<endl;
	cout<<"r2 = ("<<pos2_out[0]<<","<<pos2_out[1]<<","<<pos2_out[2]<<")"<<endl;
	
	


	// some params
	vector<double> pos1=main.getv("coords");
	vector<double> pos2=crowders[index%Ncr].PBCswitch(Ncr, index);
	vector<double> vi=main.getv("velocity");
	cout<<"v1i = ("<<vi[0]<<","<<vi[1]<<","<<vi[2]<<")"<<endl;
	vector<double> v1f(3), v2f(3), vt(3), rhat(3);
	vector<int> moving_particles;
	double m1=main.getp("mass"), m2=crowders[index%Ncr].getp("mass");
	double R1=main.getp("radius"), R2=crowders[index%Ncr].getp("radius");
	double R=R1+R2;
	double t_rem=h;
	double vni, v1nf, v2nf; // the normal component of velocity
	int index2;

	//0resolved=false, 1 ismain=true, 2secmain=false, 3secondary=false;
	vector<bool> condition(4);
	condition[0]=false;condition[1]=true;condition[2]=false;condition[3]=false;

	while(!condition[0])
	{
		vni=0;
		for(int k=0; k<3; k++)
		{
			rhat[k]=(pos2[k]-pos1[k])/R;
			vt[k]=vi[k]-(vi[k]*rhat[k])*rhat[k];
			vni+=rhat[k]*vi[k];
		}
		v1nf = vni*(m1-m2)/(m1+m2);
		v2nf = vni*2*m1/(m1+m2);
		for(int k=0;k<3;k++)
		{
			v1f[k]=vt[k]+v1nf*rhat[k];
			v2f[k]=v2nf*rhat[k];
		}
		if(condition[1])
		{
			main.newvel(v1f);
			crowders[index%Ncr].newvel(v2f);
			cout<<"v1f = ("<<v1f[0]<<","<<v1f[1]<<","<<v1f[2]<<")"<<endl;
			cout<<"v2f = ("<<v2f[0]<<","<<v2f[1]<<","<<v2f[2]<<")"<<endl;
			t_rem-=t_el;
			t_el=t_rem;
			moving_particles.push_back(index%Ncr);

			// colcheck(main=true)
		}
		else if(condition[2])
		{
			crowders[index2].newvel(v1f);
			main.newvel(v2f);
			t_rem-=t_el;
			t_el=t_rem;
		}
		else if(condition[3])
		{
			crowders[index2].newvel(v1f);
			crowders[index%Ncr].newvel(v2f);
			t_rem-=t_el;
			t_el=t_rem;
			moving_particles.push_back(index%Ncr);
			// colcheck(main=false);
		}
		/*else // not sure about this
		{
			condition[0]=true;
		}*/
		
			// advance newcoords of things
		main.mvVel(t_rem);
		for(int k=0;k<moving_particles.size();k++)
		{
			crowders[moving_particles[k]].mvVel(t_rem);
		}

		if(ColCheckAll(index, index2, moving_particles, t_el, condition))
		{
			// nudge everything to t_el
			main.nudge(t_el);
			for(int k=0;k<moving_particles.size();k++)
			{
				crowders[moving_particles[k]].nudge(t_el);
			}
			if(condition[1])
			{
				pos1=main.getv("coords");
				pos2=crowders[index%Ncr].PBCswitch(Ncr, index);
				vi=main.getv("velocity");
				m1=main.getp("mass"), m2=crowders[index%Ncr].getp("mass");
				R1=main.getp("radius"), R2=crowders[index%Ncr].getp("radius");
				R=R1+R2;
			}
			else if(condition[2])
			{
				// i think this is unnecessary
			}
			else if(condition[3])
			{
				cout<<"SECONDARY"<<endl;
				pos1=crowders[index2].getv("coords");
				pos2=crowders[index%Ncr].PBCswitch(Ncr, index);
				vi=crowders[index2].getv("velocity");
				m1=crowders[index2].getp("mass"), m2=crowders[index%Ncr].getp("mass");
				R1=crowders[index2].getp("radius"), R2=crowders[index%Ncr].getp("radius");
				R=R1+R2;
			}
		}
		
		else
		{
			cout<<"RESOLVED"<<endl;
			main.update();
			//vf_out=main.getv("velocity");
			//v2_out=crowders[indexout%Ncr].getv("velocity");
			//cout<<"v1i = ("<<vi_out[0]<<","<<vi_out[1]<<","<<vi_out[2]<<")"<<endl;
			//cout<<"v1f = ("<<vf_out[0]<<","<<vf_out[1]<<","<<vf_out[2]<<")"<<endl;
			//cout<<"v2f = ("<<v2_out[0]<<","<<v2_out[1]<<","<<v2_out[2]<<")"<<endl;
			for(int k=0;k<moving_particles.size();k++)
			{
				crowders[moving_particles[k]].update();
			}
			condition[0]=true;
		}
	}

}

bool brownsys::ColCheckAll(int& index1, int& index2, vector<int> moving_particles, double& t_el, vector<bool>& condition)
{
	vector<int> nns;
	vector<double> pos2, pos1_o, pos1;
	vector<double> vel;

	double v2, r2, r2_o, vdr, R1, R2, R, rad, rad2, tm, tp;
	bool col=false;

	/*if(main.nearCenter())
	{
		// check stuff and set events
	}
	else if(main.nearReac())
	{
		// check other stuff and set events
	}*/
	if(condition[1])
	{
		vel=main.getv("velocity");
		pos1_o=main.getv("coords");
		pos1=main.getv("new coords");
		nns=main.getNNs(true);
		R1=main.getp("radius");
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
						index1=nns[i];
						events[0]=false;
						events[1]=false;
						condition[2]=false;
						condition[3]=false;
					}
					else if(tp<t_el)
					{
						t_el=tp;
						col=true;
						index1=nns[i];
						events[0]=false;
						events[1]=false;
						condition[2]=false;
						condition[3]=false;
					}
				}
			}
		}
		if(!col)
		{
			condition[1]=false;
		}
	}

	for(int i=0; i<moving_particles.size(); i++)
	{
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
						index1=nns[j];
						index2=moving_particles[i];
						events[0]=false;
						events[1]=false;
						condition[1]=false;
						condition[2]=false;
						condition[3]=true;
					}
					else if(tp<t_el)
					{
						t_el=tp;
						col=true;
						index1=nns[j];
						index2=moving_particles[i];
						events[0]=false;
						events[1]=false;
						condition[1]=false;
						condition[2]=false;
						condition[3]=true;
					}
				}
			}
		}
		//later, check for collision with center and main
	}
	return col;
}

bool brownsys::mainColCheck(int& index, double& t_el)
{

	vector<int> nns=main.getNNs(true);
	vector<double> pos2(3), pos1_o=main.getv("coords"), pos1=main.getv("new coords");
	vector<double> vel=main.getv("velocity");

	double v2, r2, r2_o, vdr, R1=main.getp("radius"), R2, R, rad, rad2, tm, tp;
	bool col=false;

	/*if(main.nearCenter())
	{
		// check stuff and set events
	}
	else if(main.nearReac())
	{
		// check other stuff and set events
	}*/

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

	/*vector<double> PDFunc=PDF1(crowders, Ncr);
	ofstream pdf;
	pdf.open("pdfbeg.txt");
	for(int i=0;i<PDFunc.size();i++)
	{
		pdf << PDFunc[i] << endl;
	}
	pdf.close();*/



	for(int i=0;i<eqsteps;i++) // should be while !equilibrated
	{cout<<"equilibration step # "<<i<<endl;
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
/*
PDFunc=PDF1(crowders, Ncr);
ofstream pdfe;
pdfe.open("./outs/pdfend10.txt");
for(int i=0;i<PDFunc.size();i++)
{
	pdfe << PDFunc[i] << endl;
}
pdfe.close();
*/



}// ends equilibrate function

//////////////// Useful functions





//////////////////////////////////////////////////
	








double brownsys::coltime(protein two)
{
	double xy,yz,xz,xy2,xz2,yx2,yz2,zx2,zy2,rdcl2,rdcl,tm,tp,mag,t;
	double vdr, v2, mag2, v1x, v2x, v0x = 0;
	double r=two.getp("radius");
	double m=two.getp("mass");
	double R2=(r+crad)*(r+crad);
	vector<double> pos=two.getv("coords");
	vector<double> rhat(3),dr0(3),newvelc(3),newvel(3);
	for(int i=0;i<dim;i++){vdr+=cvel[i]*pos[i]; v2+=cvel[i]*cvel[i];}
	

	// doesnt work for dim~=3
	xy=2*pos[0]*pos[1]*cvel[0]*cvel[1];
	yz=2*pos[2]*pos[1]*cvel[2]*cvel[1];
	xz=2*pos[0]*pos[2]*cvel[0]*cvel[2];
	xy2=cvel[0]*cvel[0]*pos[1]*pos[1];
	xz2=cvel[0]*cvel[0]*pos[2]*pos[2];
	yz2=cvel[1]*cvel[1]*pos[2]*pos[2];
	yx2=cvel[1]*cvel[1]*pos[0]*pos[0];
	zx2=cvel[2]*cvel[2]*pos[0]*pos[0];
	zy2=cvel[2]*cvel[2]*pos[1]*pos[1];
	/*for(int i=0;i<3;i++)
	{
		cout<<i<<"coord="<<pos[i]<<", vel="<<cvel[i]<<endl;
	}*/	
	
	rdcl2=xy+xz+yz-xy2-xz2-yx2-yz2-zx2-zy2+R2*v2;
	rdcl=pow(rdcl2,0.5);
	if(rdcl2>0)
	{
		tm=(vdr-rdcl)/v2;
		tp=(vdr+rdcl)/v2;
		if(tm>0){t=tm;}
		else if(tp>0){t=tp;}
		else{cout<<"negative time"<<endl;}
		if(t>del_t){cout<<"too long"<<endl;}	
		
	}
	else{cout<<"complex time"<<endl;}
	
	return t;
}




/*void brownsys::upall()
{
	main.update();
	for(int i=0;i<Ncr;i++){crowders[i].update();}
}*/


void brownsys::something()
{
	vector<int> test;
	test=main.getNNs(0);
	cout<<test.size()<<endl;
	vector<double> test1;
	 
	test1=crowders[1].getv("coords");
	cout<<test1[0]<<" "<<test1[1]<<" "<<test1[2]<<endl;
	

}	


void brownsys::printcoords(protein test)
{
	vector<double> pos(3);
	pos=test.getv("coords");
	cout<<"x= "<<pos[0]<<endl;
	cout<<"y= "<<pos[1]<<endl;
	cout<<"z= "<<pos[2]<<endl;	
}

void brownsys::NCout()
{
	int num=nearcntr.size();
	for(int i=0;i<num;i++){cout<<nearcntr[i]<<endl;}
}






/*void brownsys::shftcntr()
{
	main.newpos(cpos, "center"); // maybe work in a better check like check |cpos-> - main->|
	reaction=main.chkreac();
	for(int i=0;i<Ncr;i++){crowders[i].newpos(cpos, "center");}
	for(int i=0;i<3;i++){cpos[i]=0;}
}*/
	
