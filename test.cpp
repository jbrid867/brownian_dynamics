#include "headers/brownsys.h"

int main(){
	bool equil = true;
	bool go = false;
	int count=0;
	int i;

	//brownsys baller(N);
	//double cutoff=6*r;
	//cout<<"numcrowd = "<<N<<endl;
	random_device rd; // Seed RNG
	mt19937 gen(rd()); // start RNG
	normal_distribution<> distro(0,pow(2*D*h,0.5));
	cout<<"number of proteins: "<<N<<endl;
	brownsys baller(N,true);
	cout<<"1"<<endl;
	baller.diff_NNs(true);
	cout<<"2"<<endl;
	// baller.equilibrate(gen, distro,100);




	// #pragma omp parallel for schedule(dynamic) 
	// for(i=0; i<100; i++)
	// {
	// 	brownsys baller(N);
	// 	baller.startNNs();
	// 	baller.equilibrate(gen,distro,1000);
	// 	baller.moveall(gen, distro, count);
	// 	cout<<"beta = "<<(double)(count/(i+1))<<endl;

	// }
	// cout<<"final beta = "<<count/1000.0<<endl;

	/*baller.startNNs(cutoff); // NOT FINISHED. ONLY FINDS IF ITS NEAR MAIN
	baller.NCout();
	baller.equilibrate(gen,distro,1000);
 	baller.updateNNs(cutoff);baller.NCout();}*/

	
	
	 
	
return 0;}





