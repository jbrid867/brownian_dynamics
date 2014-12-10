//Defines the protein and crowder classes

#include "constants1.h"
#ifndef _PROTEINS_H_
#define _PROTEINS_H_

using namespace std;

class protein
{
private:
	vector<double> coordinates;
	vector<double> newcoords;
	vector<double> vel;
	vector<int> NNs;
	vector<int> notNNs;
	double mass;
	double radius;
	bool crowder; // tells if particle is a crowder.

public:
	
	//CONSTRUCTOR
	protein(); //default constructor. only contains mass and radius
	protein(double x,double y, double z, double ms, double rad, bool crow); 
	//~protein(); //destructor
	
	//Dynamics
	
	void move(mt19937& gen, normal_distribution<> distro);
	void resetpos(); //sets newcoords=coords
	void setpos(vector<double> pos); // sets coords = pos
	void newpos(vector<double> pos, string a); // adds newcoords + pos
	void newpos(vector<double> nvel, double dt, string a); // adds newcoords + vel*dt
	void newvel(vector<double> v);
	void update();
	bool chkreac(vector<double> center); //checks for collision with center
	bool chkcntr(vector<double> center);	
	bool chkcol(vector<double> pos2);
	bool sec_col(vector<protein> &prots, int cnum);

	void remove_neighbor(int index);
	void add_neighbor(int index);
	
	void NNout(); //Print NNs
	
	
	//vector<double> resolve(vector<double> p2); //resolves collision with another protein
	//vector<double> resolve(); //resolves collision with center
	
	
	//Access 
	vector<double> getv(string a); // coords, new coords, vel
	double getp(string a); // mass, radius
	
	 
	vector<int> getNNs(bool nn); //true=NNs, false=notNNs
	void setNNs(vector<int> nns, vector<int> nnns); // sets NNs and notNNs

	
	
	
	
	
	
	
	
	
};

#endif
