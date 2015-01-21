#include "proteins.h"
#ifndef _BROWNSYS_H_
#define _BROWNSYS_H_

using namespace std;

class brownsys
{
private:
	mainP main;
	vector<protein> crowders;
	vector<int> nearcntr; //stores indices of proteins near the center
	bool crowds, centered; 
	int Ncr; //number of crowders
	int steps; //number of steps taken
	vector<double> cvel; //center velocity
	vector<double> cpos; //center coords before shift
	double crad;
	double cmass;
	double del_t; //goes from h -> 0 
	vector<bool> events; // contains bools for collision [0], escape [1], and reaction [2]
public:
	brownsys(); // creates a brownian system with one protein
	brownsys(int num); // creates brownian system with num crowders
	void something();

	// NNs
	void startNNs(double cut); // at phi=0.2, 25<NNs<50
	void updateNNs(double cut);

	// Dynamics
	void moveall(mt19937& gen, normal_distribution<> distro, int steps); //String for equilibration or actual step
	void shftcntr(); 
	void printcoords(protein test);
	double coltime(protein one, protein two);
	double coltime(protein two);
	void equilibrate(mt19937& gen, normal_distribution<> distro, int eqsteps);

	void mainRes(double& t_el, int index); // resolves the movement of main (esc/reac/collisions)
	bool mainColCheck(int& index, double& t_el);

	bool ColCheckAll(int& index1, int& index2, vector<int> moving_particles, double& t_el, vector<bool>& condition);

	void NCout();
	

	// Access functions
	void upall();
	 
};
#endif


/*NEED (in order)

function to retrieve protein location






*/
