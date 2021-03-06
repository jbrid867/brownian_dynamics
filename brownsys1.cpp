// implemention files for brownsys class

#include "headers/brownsys.h"
#include "dem_funcs.cpp"


using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
///// CONSTRUCTORS  CONSTRUCTORS  CONSTRUCTORS CONSTRUCTORS........
////////////////////////////////////////////////////////////////////////////////////////

//default, one protein system constructor
brownsys::brownsys()
{
	crowds=false;
	reaction=false;
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
	main = protein(xx,yy,zz,M,r,false);
	for(int i=0;i<dim;i++){cvel.push_back(0);cpos.push_back(0);}
	crad=r;
	cmass=M;
}

//N protein system constructor
brownsys::brownsys(int num)
{
	
	Ncr=num;
	crowds=true;	
	reaction=false;
	escape=false;
	double xx,yy,zz,mag,mag2; // dummies
	int u,v,w,ran;
	double n=ceil(pow(num,1.0/3.0));
	int Np=n*n*n; // number of lattice points
	int dN=Np-num; // number of lattice points which should eventually not contain crowders
	double l=2*L/n; // cubic lattice spacing
	vector< int > rem(dN-1);

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
	main = protein(xx,yy,zz,M,r,false);

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
						crowders.push_back(protein(xx,yy,zz,Mc,rc,true));
						numb++;}// end k if
						else if(j!=v){// j if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,Mc,rc,true));
						numb++;}// end j if
						else if(i!=u){// i if
						xx=l/2.0 + i*l-L;
						yy=l/2.0 + j*l-L;
						zz=l/2.0 + k*l-L;
						crowders.push_back(protein(xx,yy,zz,Mc,rc,true));
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

void brownsys::moveall(mt19937& gen, normal_distribution<> distro)
{
	
	del_t=h; // del_t is the time remaining, t1 is the time until an event
	//bool mainmove, cmove, col, resolved, secondary = false;
	double v,tp,tm,t, t1, mag2, mag;
	int index, index1, index2, nnNum;
	double rad2, rad, vdr, crvdr, v2, vc2, r2, R, reactime=0;
	vector<double> move(3), mvVel(3), pos1(3), pos2(3), npos(3), npos2(3), dr(3);
	vector<int> nearest, nearest2;
	vector<double> PBCvec;
	bool collision, resolved = false;
	t1=del_t;

	/////////////// COLLISION RESOLUTION STUFF
	vector<double> rhat(3), vperp(3), vpar(3), v1f(3), v2f(3);
	double m1,m2,R1,R2;
	bool ismain;


	/*
	for now:
	 - Move main
	 - move center
	 - move crowders
	later:
	 - move them all randomly with checks
	 */
	 /////////////////////////////////////////////////////////////////////////
	 /////////////////////////////////////////////////////////////////////////
	 //Move main
	 /////////////////////////////////////////////////////////////////////////
	 /////////////////////////////////////////////////////////////////////////

	 nearest=main.getNNs(true);
	 pos1=main.getv("coords");
	 mag2=0;
	 for(int k=0;k<3;k++)
	 {
	 	move[k]=distro(gen); mvVel[k]=move[k]/h;npos[k]=move[k]+pos1[k]; 
	 	mag2+=npos[k]*npos[k]; v2+=mvVel[k]*mvVel[k]; vdr+=mvVel[k]*pos1[k];
	 	r2+=pos1[k]*pos1[k];
	 }
	 mag=pow(mag2,0.5);

	 // REACTION
	 if(mag<crad+main.getp("radius"))
	 {
	 	rad2=vdr*vdr - v2*(r2-(crad+main.getp("radius")));
	 	if(rad2>0)
	 	{
	 		rad=pow(rad2,0.5);
	 	}
	 	else
	 	{
	 		cout<<"err: negative in radical"<<endl;
	 	}
	 	tm=(-1.0*vdr - rad)/v2;
	 	tp=(-1.0*vdr + rad)/v2;
	 	if(tm<0)
	 	{
	 		t=tp;
	 	} 
	 	else
	 	{
	 		t=tm;
	 	}
	 	if(t<del_t)
	 	{
	 		reaction=true;
	 		t1=t;
	 		cout<<"potential reaction"<<endl;
	 	}
	 	else
	 	{
	 		cout<<"err: time is greater than time step"<<endl;
	 	}



	 } // closes reaction loop

	 // ESCAPE
	 else if(mag > q)
	 {
	 	rad2=vdr*vdr - v2*(r2-q*q);
	 	if(rad2>0)
	 	{
	 		rad=pow(rad2,0.5);
	 	}
	 	else
	 	{
	 		cout<<"err: negative in radical"<<endl;
	 	}
	 	tm=(-1.0*vdr - rad)/v2;
	 	tp=(-1.0*vdr + rad)/v2;
	 	if(tm<0)
	 	{
	 		t=tp;
	 	} 
	 	else
	 	{
	 		t=tm;
	 	}
	 	if(t<del_t)
	 	{
	 		escape=true;
	 		t1=t;
	 		cout<<"potential escape"<<endl;
	 	}
	 	else
	 	{
	 		cout<<"err: time is greater than time step"<<endl;
	 	}
	 } // closes escape loop

	 // collisions with crowders
	 R = main.getp("radius")+crowders[nearest[k]%Ncr].getp("radius");
	 for(int k=0;k<nearest.size();k++)
	 {
	 	vdr, mag2=0;
	 	pos2=crowders[nearest[k]%Ncr].getv("coords");
	 	PBCvec=PBCswitch(Ncr, nearest[k]);

	 	for(int i=0;i<3;i++) //
	 	{
	 		pos2[i]+=PBCvec[i]*L; // position of potential collidee
	 		dr[i]=pos2[i]-pos1[i];
	 		mag2+=dr[i]*dr[i]; // dist^2 between particles
	 		vdr+=mvVel[i]*dr[i];

	 	}
	 	mag=pow(mag2,0.5);
	 	if(mag<R)
	 	{
	 		
	 		rad2 = vdr*vdr - v2*(mag2-R*R);
	 		if(rad2>0)
	 		{
	 			rad=pow(rad2,0.5);
	 		}
	 		else
	 		{
	 			cout<<"err: negative in radical"<<endl;
	 		}
	 		tm = (vdr - rad)/v2;
	 		tp = (vdr + rad)/v2;
	 		if(tm<0)
	 		{
	 			t=tp;
	 		}
	 		else
	 		{
	 			t=tm;
	 		}
	 		if(t<t1)
	 		{
	 			t1 = t;
	 			reaction, escape = false;
	 			collision=true;
	 			index=nearest[k];
	 		}

	 	} // closes crowder collision loop
	 } // closes NN loop

	 if(escape)
	 {
	 	cout<<"No reaction"<<endl;
	 }
	 else if(reaction)
	 {
	 	cout<<"Reaction"<<endl;
	 }
	 else if(collision)
	 {
	 	index2=index;
	 	index1=-1; // use -1 to indicate that its main
	 	ismain=true;
	 	R1=main.getp("radius");
	 	R2=crowders[index%Ncr].getp("radius");
	 	m1=main.getp("mass");
		m2=crowders[index%Ncr].getp("mass");
		pos2=crowders[index%Ncr].getv("coords");
		PBCvec=PBCswitch(Ncr, index);
		del_t=del_t-t1; // remaining time, move forward to t1 seconds
	 	while(!resolved)
	 	{
		 	vdr=0;
		 	crvdr=0;
		 	v2,vc2,mag2=0;
		 	collision=false;
		 	
		 	for(int k=0;k<3;k++)
		 	{
		 		pos1[k]+=mvVel[k]*t1; // update new position
		 		rhat[k]=((pos2[k]+PBCvec[k]*L)-npos[k])/R; // unit vector from 1 to 2
		 		vperp[k]=mvVel[k]*rhat[k]; // normal component
		 		vpar[k]=mvVel[k]-vperp[k]; // tangent component
		 		v1f[k]=vpar[k] + vperp[k]*(m1-m2)/(m1+m2); // final vel of main
		 		v2f[k]=vperp[k]*2*m1/(m1+m2); // final vel of crowder

		 		// update some stuff to check for secondaries
		 		npos[k]=pos1[k]+v1f[k]*del_t;
		 		npos2[k]=pos2[k] + v2f[k]*del_t;
		 		v2+=v1f[k]*v1f[k];
		 		vc2+=v2f[k]*v2f[k];
		 		mag2+=npos[k]*npos[k];
		 		vdr+=pos1[k]*v1f[k];
		 		crvdr+=pos2[k]+v2f[k];
		 	}

		 	t1=del_t;
		 	if(ismain)
		 	{
			 	//check reaction
			 	mag=pow(mag2,0.5);
			 	if(mag<crad+R1)
				 {
				 	rad2=vdr*vdr - v2*(r2-(crad+main.getp("radius")));
				 	if(rad2>0)
				 	{
				 		rad=pow(rad2,0.5);
				 	}
				 	else
				 	{
				 		cout<<"err: negative in radical"<<endl;
				 	}
				 	tm=(-1.0*vdr - rad)/v2;
				 	tp=(-1.0*vdr + rad)/v2;
				 	if(tm<0)
				 	{
				 		t=tp;
				 	} 
				 	else
				 	{
				 		t=tm;
				 	}
				 	if(t<del_t)
				 	{
				 		reaction=true;
				 		t1=t;
				 		cout<<"potential reaction"<<endl;
				 	}
				 	else
				 	{
				 		cout<<"err: time is greater than time step"<<endl;
				 	}
				} // closes reaction loop

				// ESCAPE
		 		else if(mag > q)
		 		{
		 			rad2=vdr*vdr - v2*(r2-q*q);
				 	if(rad2>0)
				 	{
				 		rad=pow(rad2,0.5);
				 	}
				 	else
				 	{
				 		cout<<"err: negative in radical"<<endl;
				 	}
				 	tm=(-1.0*vdr - rad)/v2;
				 	tp=(-1.0*vdr + rad)/v2;
				 	if(tm<0)
				 	{
				 		t=tp;
				 	} 
				 	else if(tm>0)
				 	{
				 		t=tm;
				 	}
				 	else
				 	{
				 		cout<<"err: negative times"<<endl;
				 	}
				 	if(t<del_t)
				 	{
				 		escape=true;
				 		t1=t;
				 		cout<<"potential escape"<<endl;
				 	}
				 	else
				 	{
				 		cout<<"err: time is greater than remaining time"<<endl;
				 	}
		 		} // closes escape loop
		 	} // closes ismain if

	 		//check secondary collision for particle 1, use pos2 as dumby
	 		if(index1==-1){nearest=main.getNNs(true);}
	 		else{nearest=crowders[index1%Ncr]}
	 		for(int i=0;i<nearest.size();i++)
	 		{
	 			pos2=crowders[nearest[i]%Ncr].getv("coords");
	 			PBCvec=PBCswitch(Ncr, nearest[i]);
	 			R=main.getp("radius")+crowders[nearest[i]%Ncr].getp("radius");
	 			for(int i=0;i<3;i++) //
	 			{
			 		pos2[i]+=PBCvec[i]*L; // position of potential collidee
			 		dr[i]=pos2[i]-pos1[i];
			 		mag2+=dr[i]*dr[i]; // dist^2 between particles
			 		vdr+=mvVel[i]*dr[i];
			 	}
			 	mag=pow(mag2,0.5);
			 	if(mag<R)
			 	{
			 		rad2=(vdr*vdr - v2*(mag2-R*R))/v2;
			 		if(rad2>0)
			 		{
			 			rad=pow(rad2,0.5);
			 			tm=(vdr-rad)/v2;
			 			tp=(vdr+rad)/v2;
			 			if(tm<0)
			 			{
			 				t=tp;
			 			}
			 			else
			 			{
			 				t=tm;
			 			}
			 			if(t<t1)
			 			{
			 				t1=t;
			 				reaction, escape=false;
			 				secondary=true;

			 			}
			 		}
			 	} // close time calculation if
			} // close main nn loop
			if(secondary)
			{
				// set up to continue while
			}
			else if(reaction || escape)
			{
				resolved=true;
			}
 		

		}//closes secondary resolution
	}
			
				 	



		 








	 
	 else
	 {
	 	main.setpos(npos);
	 }



	 ////////////////////////////////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////////////////////
	 /////////////// MOVE CENTER ////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////////////////////
	 ////////////////////////////////////////////////////////////////////////////

	 for(int k=0;k<3;k++)
	 {
	 	move[k]=distro(gen); mvVel[k]=move[k]/h; npos[k]=move[k];
	 }

	 // reaction
	 pos1=main.getv("coords");
	 mag2=0;
	 for(int k=0;k<3;k++)
	 {
	 	mag2+=(npos[k]-pos1[k])*(npos[k]-pos1[k]);
	 }
	 mag=pow(mag2,0.5);
	 if(mag<main.getp("radius")+crad)
	 {
	 	// store time and compare with crowder collisions
	 }

	 // check for collision with crowders
	 for(int k=0; k<nearcntr.size(); k++)
	 {
	 	mag2=0;
	 	pos2=crowders[k].getv("coords");
	 	for(int i=0;i<3;i++)
	 	{
	 		mag2+=(pos2[k]-npos[k])*(pos2[k]-npos[k]);
	 	}
	 	mag=pow(mag2,0.5);
	 	if(mag<crad+crowders[k].getp("radius"))
	 	{
	 		// RESOLVE
	 	}
	 }

	 // update ALL positions, check for escape THEN check PBCs


	 //////////////////////////////////////////////////////////////////////////////////
	 //////////////////////////////////////////////////////////////////////////////////
	 // MOVE CROWDERS /////////////////////////////////////////////////////////////////
	 //////////////////////////////////////////////////////////////////////////////////
	 //////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////////////
//////// OLD MOVE SHIT ////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

/*	
	// move central particle. DONT NECESSARILY DO THIS FIRST
	for(int i=0;i<dim;i++){cpos[i]=distro(gen); cvel[i]=cpos[i]/del_t; v2+=cvel[i]*cvel[i];}
	
	
	// move main and check for collision/escape. definitely do this first
	//main.newpos(cpos,"center");
	
	
	
	
	// move crowders and resolve collisions
	while(!resolved){
	col=false;
	
	t=del_t;
	reactime=del_t;

	if(main.chkreac(cpos))
	{
		reaction=true;
		cout<<"Potential reaction"<<endl;
		t=coltime(main);
		reactime=t;
	}
	
	if(nearcntr.size()!=0){
	for(int i=0;i<nearcntr.size();i++)
	{
		
		//crowders[i].newpos(cpos, "center");
		if(crowders[nearcntr[i]].chkcntr(cpos))
		{
			col=true;
			t1=coltime(crowders[nearcntr[i]]);
			if(t1<t){t=t1; index=nearcntr[i]; reaction=false;}
		}
			
			
		
	}} // end find collision for and near center if

	if(del_t<pow(10.0,-18.0)){col=false; reaction=false;}
	if(col && !reaction)
	{
		//main.newpos(cvel, t,"center"); 
		
		for(int i=0;i<3;i++){cpos[i]=cvel[i]*t;}
		
		/*for(int i=0;i<Ncr;i++)
		{
			crowders[i].newpos(cvel, t, "center");
			crowders[i].update();
		}*//*

		//for(int i=0;i<3;i++){cpos[i]=cvel[i]*t;}
		
		resolvec(crowders[index],t); // check secondary co
		
		/*main.newpos(cpos, "center");
		for(int i=0;i<Ncr;i++)
		{
			crowders[i].newpos(cpos, "center");
		}*//*
		cout<<"collision resolved "<<index<<endl;
		
		
	}
	else if(reaction)
	{
		cout<<"reaction"<<endl;
		// throw end simulation flag
		resolved=true;
	}
	else
	{
		main.newpos(cpos, "center");
		for(int i=0; i<Ncr; i++)
		{
			crowders[i].newpos(cpos, "center");
		}
		
		
		cout<<"remaining time = "<<del_t<<endl;
		resolved=true;
	}//end else
	
	}//end while

*/

///////////////////////////////////////////////////////////////////////////
/////////MOVE MAIN/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/*
vector<double> disp(3), newpos(3), newvel(3);
vector<double> posm;
reactime=0;
mag2=0;
posm=main.getv("coords");
For(int i=0;i<3;i++){disp[i]=distro(gen); newpos[i]=posm[i]+disp[i]; newvel[i]=disp[i]/h; mag2+=newpos[i]*newpos[i];}
mag=pow(mag2,0.5);

if(mag<main.getp("radius")+crowders[1].getp("radius")){reaction=true;}
*/


	
	
}// end moveall



void brownsys::equilibrate(mt19937& gen, normal_distribution<> distro, int eqsteps)
{
	vector<double> pos1(3), pos2(3), posm(3), disp(3), PBCvec(3);
	vector<int> nns;
	int numnns, count, clashcount=0;
	bool clash;
	double mag, mag2, r1, r2;

	posm=main.getv("coords");

	vector<double> PDFunc=PDF1(crowders, Ncr);
	ofstream pdf;
	pdf.open("pdfbeg.txt");
	for(int i=0;i<PDFunc.size();i++)
	{
		pdf << PDFunc[i] << endl;
	}
	pdf.close();



	for(int i=0;i<eqsteps;i++) // should be while !equilibrated
	{cout<<"equilibration step # "<<i<<endl;
	for(int k=0;k<Ncr;k++){
		
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
			pos2=crowders[nns[count]%Ncr].getv("coords");
			PBCvec = PBCswitch(Ncr, nns[count]);

			for(int dumb=0;dumb<3;dumb++){mag2+=(pos1[dumb]-(pos2[dumb]+PBCvec[dumb]))*(pos1[dumb]-(pos2[dumb]+PBCvec[dumb])); } //
			mag=pow(mag2,0.5);
			if(mag<crowders[k].getp("radius")+crowders[nns[count]%Ncr].getp("radius"))
			{clash=true; clashcount++;}
			else{count++;}
		}
		if(!clash){crowders[k].setpos(pos1);}


	}}//ends eqsteps and Ncr for loops
cout<<"clash count = "<<clashcount<<endl;

PDFunc=PDF1(crowders, Ncr);
ofstream pdfe;
pdfe.open("./outs/pdfend10.txt");
for(int i=0;i<PDFunc.size();i++)
{
	pdfe << PDFunc[i] << endl;
}
pdfe.close();




}// ends equilibrate function

//////////////// Useful functions





//////////////////////////////////////////////////
	





void brownsys::resolvec(protein& two, double dt)
{
	// Need to add secondary collisions of crowder somehow. 
	// Give crowders velocity, check outside of the resolve function. dont update coordinates in resolvec. Check for another collision, 
	

	vector<double> dr0(3), newvel(3), rhat(3);
	vector<int> nns;
	vector<double> pos2=two.getv("coords");
	double mag, mag2, v1x, v2x, v0x=0;
	double r=two.getp("radius");
	double m=two.getp("mass");
	double R2=(r+crad)*(r+crad);
	for(int i=0;i<3;i++){dr0[i]=pos2[i]-cpos[i]; mag2+=dr0[i]*dr0[i];}
	mag=pow(mag2,0.5);
	for(int i=0;i<3;i++){rhat[i]=dr0[i]/mag; v0x+=rhat[i]*cvel[i];}
	v2x=(2*cmass)/(cmass+m) * v0x;
	v1x=(cmass-m)/(cmass+m) * v0x;
	
	for(int i=0;i<3;i++)
	{	
		cvel[i]+=(v1x-v0x)*rhat[i]; 
		newvel[i]=v2x*rhat[i];
		cpos[i]=cvel[i]*(del_t-dt);
		
		
	}
	two.newvel(newvel);
	two.newpos(newvel, del_t-dt, "new");
	//CHECK SECONDARY COLLISIONS
	two.sec_col(crowders, Ncr);
	del_t-=dt;	

}


double brownsys::coltime(protein two)
{
	//two.resetpos();
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




void brownsys::upall()
{
	main.update();
	for(int i=0;i<Ncr;i++){crowders[i].update();}
}


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
	
