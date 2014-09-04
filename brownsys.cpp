// implemention files for brownsys class

#include "headers/brownsys.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
///// CONSTRUCTORS  CONSTRUCTORS  CONSTRUCTORS CONSTRUCTORS........
////////////////////////////////////////////////////////////////////////////////////////

//default, one protein system constructor
brownsys::brownsys()
{
	crowds=false;
	reaction=2;
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
}

//N protein system constructor
brownsys::brownsys(int num)
{
	
	Ncr=num;
	crowds=true;	
	reaction=2;
	double xx,yy,zz,mag,mag2; // dummies
	int u,v,w,ran;
	double n=ceil(pow(num,1.0/3.0));
	int Np=n*n*n; // number of lattice points
	int dN=Np-num; // number of lattice points which should eventually not contain crowders
	double l=2*L/n; // cubic lattice spacing
	vector< int > rem(dN-1);

	
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
{if(crowds){
	vector<int> nns;
	vector<int> nnns;
	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2;
	dumby=main.getinfo("coords");
	xx=dumby[0];yy=dumby[1];zz=dumby[2];
	bool NN;	

	for(int u=0; u<Ncr; u++)
	{	
		NN=false;
		dumby=crowders[u].getinfo("coords");
		x=dumby[0];y=dumby[1];z=dumby[2];
		mag2=(xx-x)*(xx-x)+(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u);NN=true;}
		x-=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+Ncr);NN=true;}
		x+=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+2*Ncr);NN=true;}
		y-=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+3*Ncr);NN=true;}
		y+=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+4*Ncr);NN=true;}
		z-=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+5*Ncr);NN=true;}
		z+=2*L;
		mag2=(xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z);
		mag=pow(mag2,0.5);
		if(mag<cut){nns.push_back(u+6*Ncr);NN=true;}
		if(!NN){nnns.push_back(u);}
	}
	
	main.setNNs(nns, nnns);
	
	
	










}}// ends function and big if


void brownsys::updateNNs(double cut)
{
	vector<double> pos1, pos2;
	double dx,dy,dz,mag2,mag;
	vector<int> NNs;
	vector<int> nNNs;
	int index;
	int index2;
	//update main
	pos1=main.getinfo("coords");
	NNs=main.getNNs(true);
	nNNs=main.getNNs(false);
	int NNnum=NNs.size();
	int nNNnum=nNNs.size();
	
	
	for(int i=0; i<NNnum; i++)
	{
		
		index=NNs[i];
		
		if(index<Ncr)
		{
			pos2=crowders[index].getinfo("coords");
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=Ncr && index<2*Ncr) //-x border
		{
			index2=index-Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[0]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=2*Ncr && index<3*Ncr) //+x border
		{
			index2=index-2*Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[0]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=3*Ncr && index<4*Ncr) //-y border
		{
			index2=index-3*Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[1]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=4*Ncr && index<5*Ncr) //+y border
		{
			index2=index-4*Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[1]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=5*Ncr && index<6*Ncr) //-z border
		{
			index2=index-5*Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[2]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=6*Ncr && index<7*Ncr) //+z border
		{
			index2=index-6*Ncr;
			pos2=crowders[index2].getinfo("coords");
			pos2[2]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
	
	
	} // ends the NN removal for loop		

	for(int i=0; i<nNNnum; i++) //CAN MAKE THIS BETTER BY LINKING PARTICLES FOR WHICH IVE ALREADY FOUND NNS
	{
		index=nNNs[i];
		pos2=crowders[index].getinfo("coords");
		dx=pos1[0]-pos2[0];
		dy=pos1[1]-pos2[1];
		dz=pos1[2]-pos2[2];
		if(dx<cut)
		{if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//all ifs ended
		else if(dx<cut-2*L) // -x boundary
		{dx+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//-x ifs ended
		else if(dx>2*L-cut) // +x boundary
		{dx-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+2*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//+x ifs ended
		else if(dy<cut-2*L) // -y boundary
		{dy+=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+3*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//-y ifs ended
		else if(dy>2*L-cut) // +y boundary
		{dy-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+4*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//+y ifs ended
		else if(dz<cut-2*L) // -z boundary
		{dz+=2*L; if(pow(dz*dz+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+5*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//-z ifs ended
		else if(dz>2*L-cut) // +z boundary
		{dz-=2*L; if(pow(dx*dx+dy*dy,0.5)<cut)
		{if(pow(dx*dx+dy*dy+dz*dz,0.5)<cut){NNs.push_back(index+6*Ncr);nNNs.erase(nNNs.begin()+i);nNNnum-=1;}
		}}//+z ifs ended

	} // ends NN addition loop	
}

void brownsys::moveall(mt19937& gen, normal_distribution<> distro)
{
	double xx,yy,zz;
	vector<double> center(dim);
	// move central particle. DONT NECESSARILY DO THIS FIRST
	center[0]=-1.0*distro(gen);
	center[1]=-1.0*distro(gen);
	center[2]=-1.0*distro(gen);
	// move main and check for collision/escape
	main.newpos(center);
	
	if(main.chkcntr()){//check if collision happens first, if it doesnt
		reaction=1;}
	else if(main.chkesc()){reaction=0;}
	// move crowders and resolve collisions
	for(int i=0;i<Ncr;i++)
	{
		crowders[i].newpos(center);
		if(crowders[i].chkcntr()){/*resolve*/}
	}
	
	
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
	 
	test1=crowders[1].getinfo("coords");
	cout<<test1[0]<<" "<<test1[1]<<" "<<test1[2]<<endl;
	

}	






	
