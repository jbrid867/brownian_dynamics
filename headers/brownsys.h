#include "proteins.h"
#ifndef _BROWNSYS_H_
#define _BROWNSYS_H_

using namespace std;

class brownsys
{
private:
	protein main;
	vector<protein> crowders;
	bool crowds; 
	int Ncr; //number of crowders
	int steps; //number of steps taken
	int reaction; // 0 for escape, one for reaction
	vector<double> cvel; //center velocity
	vector<double> cpos; //center coords before shift
	double crad;
	double cmass;
public:
	brownsys(); // creates a brownian system with one protein
	brownsys(int num); // creates brownian system with num crowders
	void something();

	// NNs
	void startNNs(double cut); // at phi=0.2, 25<NNs<50
	void updateNNs(double cut);

	// Dynamics
	void moveall(mt19937& gen, normal_distribution<> distro);
	void resolve(protein& one); // resolves collisions with center
	void shftcntr(); 

	// Access functions
	void upall();
	 
};
#endif


/*NEED (in order)

function to retrieve protein location






*/
