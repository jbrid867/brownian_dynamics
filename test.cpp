#include "headers/brownsys.h"

int main(){
	bool equil = true;
	bool go = false;
	int count=0;
	int i;

	//brownsys baller(N);
	//double cutoff=6*r;
	cout<<"numcrowd = "<<N<<endl;
	random_device rd; // Seed RNG
	mt19937 gen(rd()); // start RNG
	normal_distribution<> distro(0,pow(2*D*h,0.5));

	
	for(i=0; i<1; i++)
	{
		brownsys baller(N);
		cout<<"1";
		baller.equilibrate(gen, distro, 10000);
		
	}

	/*baller.startNNs(cutoff); // NOT FINISHED. ONLY FINDS IF ITS NEAR MAIN
	baller.NCout();
	baller.equilibrate(gen,distro,1000);
 	baller.updateNNs(cutoff);baller.NCout();}*/

	
	
	 
	
return 0;}





