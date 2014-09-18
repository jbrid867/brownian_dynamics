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
	for(int i=0;i<dim;i++){cvel.push_back(0);cpos.push_back(0);}
	crad=r;
	cmass=M;
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
{if(crowds){
	vector<int> nns;
	vector<int> nnns;
	vector<double> dumby;
	double xx, yy, zz, x, y, z,mag,mag2;
	dumby=main.getv("coords");
	xx=dumby[0];yy=dumby[1];zz=dumby[2];
	bool NN;	

	for(int u=0; u<Ncr; u++)
	{	
		NN=false;
		dumby=crowders[u].getv("coords");
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
	pos1=main.getv("coords");
	NNs=main.getNNs(true);
	nNNs=main.getNNs(false);
	int NNnum=NNs.size();
	int nNNnum=nNNs.size();
	
	
	for(int i=0; i<NNnum; i++)
	{
		
		index=NNs[i];
		
		if(index<Ncr)
		{
			pos2=crowders[index].getv("coords");
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=Ncr && index<2*Ncr) //-x border
		{
			index2=index-Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=2*Ncr && index<3*Ncr) //+x border
		{
			index2=index-2*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[0]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=3*Ncr && index<4*Ncr) //-y border
		{
			index2=index-3*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=4*Ncr && index<5*Ncr) //+y border
		{
			index2=index-4*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[1]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=5*Ncr && index<6*Ncr) //-z border
		{
			index2=index-5*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]-=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
		}
		else if(index>=6*Ncr && index<7*Ncr) //+z border
		{
			index2=index-6*Ncr;
			pos2=crowders[index2].getv("coords");
			pos2[2]+=2*L;
			mag2=(pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1])+(pos1[2]-pos2[2])*(pos1[2]-pos2[2]);
			mag=pow(mag2,0.5);
			if(mag>cut){NNs.erase(NNs.begin()+i);nNNs.push_back(index);NNnum-=1;}
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
	cout<<"move"<<endl;
	del_t=h;
	bool col, resolved = false;
	double xx,yy,zz,v,t1,t;
	int index, index1;
	double v2, reactime=0;
	
	// move central particle. DONT NECESSARILY DO THIS FIRST
	for(int i=0;i<dim;i++){cpos[i]=distro(gen); cvel[i]=cpos[i]/del_t; v2+=cvel[i]*cvel[i];}
	
	
	// move main and check for collision/escape. definitely do this first
	//main.newpos(cpos,"center");
	
	
	t=del_t+2;
	reactime=del_t+2;
	// move crowders and resolve collisions
	while(!resolved){
	col=false;

	if(main.chkreac(cpos))
	{
		reaction=true;
		cout<<"Potential reaction"<<endl;
		t=coltime(main);
		reactime=t;
	}
	
	
	for(int i=0;i<Ncr;i++)
	{
		
		//crowders[i].newpos(cpos, "center");
		if(crowders[i].chkcntr(cpos))
		{
			col=true;
			cout<<"COLLISION "<<i<<endl;
			index1=i;
			t1=coltime(crowders[i]);
			if(t1<t){t=t1; index=index1; reaction=false;}
		}
			
			
		
	} // end find collision for
	if(col && t<reactime)
	{
		main.newpos(cvel, t,"center"); 
		
		for(int i=0;i<Ncr;i++)
		{
			crowders[i].newpos(cvel, t, "center");
			crowders[i].update();
		}
		resolvec(crowders[index],t); // check secondary co
		main.newpos(cpos, "center");
		for(int i=0;i<Ncr;i++)
		{
			crowders[i].newpos(cpos, "center");
		}
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
		cout<<"all collisions resolved"<<endl;
		resolved=true;
	}//end else
	
	}//end while
	
	
}

void brownsys::resolvec(protein& two, double dt)
{
	// Need to add secondary collisions of crowder somehow. 
	// Give crowders velocity, check outside of the resolve function. dont update coordinates in resolvec. Check for another collision, 
	

	vector<double> dr0(3), newvel(3), rhat(3);
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
		two.newpos(newvel, del_t-dt, "normal");
	}
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
	

	//not doesnt work for dim~=3
	xy=2*pos[0]*pos[1]*cvel[0]*cvel[1];
	yz=2*pos[2]*pos[1]*cvel[2]*cvel[1];
	xz=2*pos[0]*pos[2]*cvel[0]*cvel[2];
	xy2=cvel[0]*cvel[0]*pos[1]*pos[1];
	xz2=cvel[0]*cvel[0]*pos[2]*pos[2];
	yz2=cvel[1]*cvel[1]*pos[2]*pos[2];
	yx2=cvel[1]*cvel[1]*pos[0]*pos[0];
	zx2=cvel[2]*cvel[2]*pos[0]*pos[0];
	zy2=cvel[2]*cvel[2]*pos[1]*pos[1];
	for(int i=0;i<3;i++)
	{
		cout<<i<<"coord="<<pos[i]<<", vel="<<cvel[i]<<endl;
	}	
	
	rdcl2=xy+xz+yz-xy2-xz2-yx2-yz2-zx2-zy2+R2*v2;
	rdcl=pow(rdcl2,0.5);
	if(rdcl2>0)
	{
		tm=(vdr-rdcl)/v2;
		cout<<"t minus = "<<tm<<endl;
		tp=(vdr+rdcl)/v2;
		cout<<"t plus = "<<tp<<endl;
		if(tm>0){t=tm;}
		else if(tp>0){t=tp;}
		else{cout<<"negative time"<<endl;}
		if(t>del_t){cout<<"too long"<<endl;}	
		
	}
	else{cout<<"complex time"<<endl;}
	cout<<"t = "<<t<<endl;
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


/*void brownsys::shftcntr()
{
	main.newpos(cpos, "center"); // maybe work in a better check like check |cpos-> - main->|
	reaction=main.chkreac();
	for(int i=0;i<Ncr;i++){crowders[i].newpos(cpos, "center");}
	for(int i=0;i<3;i++){cpos[i]=0;}
}*/
	
