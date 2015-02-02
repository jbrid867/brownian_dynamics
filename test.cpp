#include "headers/brownsys.h"

int main(){
	bool equil = true;
	bool go = false;
	int count=0;

	brownsys baller;
	double cutoff=6*r;
	cout<<"numcrowd = "<<N<<endl;
	random_device rd; // Seed RNG
	mt19937 gen(rd()); // start RNG
	normal_distribution<> distro(0,pow(2*D*h,0.5));

	for(int i=0; i<1000; i++)
	{
		baller=brownsys(N);
		baller.startNNs(cutoff);
		baller.equilibrate(gen,distro,1000);
		baller.moveall(gen, distro, count);
	}
	cout<<"beta = "<<count/1000.0<<endl;

	/*baller.startNNs(cutoff); // NOT FINISHED. ONLY FINDS IF ITS NEAR MAIN
	baller.NCout();
	baller.equilibrate(gen,distro,1000);
 	baller.updateNNs(cutoff);baller.NCout();}*/

	
	
	 
	
return 0;}





