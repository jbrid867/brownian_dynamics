#include "headers/brownsys.h"

int main(){
	bool equil = true;
	bool go = false;


	double cutoff=4*r;
	cout<<"numcrowd = "<<N<<endl;
	random_device rd; // Seed RNG
	mt19937 gen(rd()); // start RNG
	normal_distribution<> distro(0,pow(2*D*h,0.5));
	brownsys baller(N);
	baller.moveall(gen, distro);


	/*baller.startNNs(cutoff); // NOT FINISHED. ONLY FINDS IF ITS NEAR MAIN
	baller.NCout();
	baller.equilibrate(gen,distro,1000);
 	baller.updateNNs(cutoff);baller.NCout();}*/

	
	
	 
	
return 0;}





