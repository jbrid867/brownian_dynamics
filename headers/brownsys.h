#include "proteins.h"
#ifndef _BROWNSYS_H_
#define _BROWNSYS_H_

using namespace std;

class brownsys
{
private:
	// proteins
	mainP main;
	vector<protein> crowders;

	// central protein ish
	vector<int> nearcntr; //stores indices of proteins near the center
	bool crowds, centered; //eventually polymorph the crowds part
	vector<double> cvel; //center velocity
	vector<double> cpos; //center coords before shift
	vector<double> ncpos; // center newcoords
	vector<int> boxlist;

	// parameters
	int Ncr; //number of crowders
	double crad;
	double cmass;

	//double L;
	double l; //lattice spacing
	//int N;
	int n;

	// esc/reac stuff
	vector<bool> events; // contains bools for reaction [0] and escape [1]
public:
	// constructors
	brownsys(); // creates a brownian system with one protein
	// for reaction
	brownsys(int num); // creates brownian system with num crowders
	// for diffusion constant
	brownsys(int num, bool dumby);
	~brownsys();
	
	// NNs
	void startNNs(); // builds NN lists
	void diff_NNs(bool init); // utilizes grid
	void updateNNs(); // updates NN lists, consider using references/pointers

	// Dynamics
	void moveall(mt19937& gen, normal_distribution<> distro, int& betacount); //String for equilibration or actual step
	void equilibrate(mt19937& gen, normal_distribution<> distro, int eqsteps);
	void equil_diff(mt19937& gen, normal_distribution<> distro, int eqsteps);
	void upall(vector<int> moving_particles); // coords=newcoords for all
	
	// Collisions
	void mainRes(double& t_el, int index, int colSwitch, vector<int> moving_particles); // resolves the movement of main (esc/reac/collisions)
	bool mainColCheck(int& index, double& t_el); // initial collision check for main
	bool ColCheckAll(int& switcher, int& index1, int& index2, vector<int> moving_particles, double& t_el, vector<bool>& condition);
	
	// Crowder movement
	void mainCheck(int index, double& t_el);
	void moveCrowd(int j, mt19937& gen, normal_distribution<> distro);
	void crColCheck(int index);

	//center movement
	void cColCheck();
	void shftcntr(); // shifts reference frame to center, updates PBCs	
	// debugging
	void printcoords(protein test);	 
};

class Reaction_sys : public brownsys
{
private:
	mainP main;
};
#endif

