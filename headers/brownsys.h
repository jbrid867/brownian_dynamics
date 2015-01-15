#include "proteins.h"
#ifndef _BROWNSYS_H_
#define _BROWNSYS_H_

using namespace std;

class brownsys
{
private:
	protein main;
	vector<protein> crowders;
	vector<int> nearcntr; //stores indices of proteins near the center
	bool crowds; 
	int Ncr; //number of crowders
	int steps; //number of steps taken
	bool reaction, escape; // true if there is a reaction
	vector<double> cvel; //center velocity
	vector<double> cpos; //center coords before shift
	double crad;
	double cmass;
	double del_t; //goes from h -> 0 
public:
	brownsys(); // creates a brownian system with one protein
	brownsys(int num); // creates brownian system with num crowders
	void something();

	// NNs
	void startNNs(double cut); // at phi=0.2, 25<NNs<50
	void updateNNs(double cut);

	// Dynamics
	void moveall(mt19937& gen, normal_distribution<> distro); //String for equilibration or actual step
	void resolvec(protein& two, double dt); // resolves collisions with center
	void resolvec(protein& one, protein& two, double dt); // resolves collisions between two proteins
	void shftcntr(); 
	void printcoords(protein test);
	double coltime(protein one, protein two);
	double coltime(protein two);
	void equilibrate(mt19937& gen, normal_distribution<> distro, int eqsteps);

	void NCout();
	

	// Access functions
	void upall();
	 
};
#endif


/*NEED (in order)

function to retrieve protein location






*/
