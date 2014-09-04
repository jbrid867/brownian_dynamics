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
	
	//Dynamic 
	vector<double> getinfo(const char* a); // a = coords, new coords or velocity
	void newpos(vector<double> pos); // sets newcoords
	void update();
	void PBC();
	
	//Parameters 
	double  getmass();
	double getrad();
	 
	vector<int> getNNs(bool nn); //true=NNs, false=notNNs
	void setNNs(vector<int> nns, vector<int> nnns); // sets NNs and notNNs

	//incomplete
	
	
	void move(mt19937& gen, normal_distribution<> distro); // updates coordinates and generates velocity
	
	
	
	
	
	
};

#endif
