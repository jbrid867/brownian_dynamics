//Defines the protein and crowder classes

#include "constants1.h"
#ifndef _PROTEINS_H_
#define _PROTEINS_H_

using namespace std;

class protein // for all intents and purposes, this is a crowder
{
protected:
	vector<double> coordinates;
	vector<double> newcoords;
	vector<double> colCoords;
	vector<double> vel;
	vector<int> NNs;
	vector<int> notNNs;
	vector<bool> Boundaries; //8 true/false elements
	double mass;
	double radius;
	bool nearcenter;
	double cdist; // distance from center	

public:	
	//CONSTRUCTOR
	protein(); //default constructor. only contains mass and radius
	protein(double x,double y, double z, double ms, double rad); 
	//~protein(); //destructor
	
	//Dynamics
	
	void move(mt19937& gen, normal_distribution<> distro); // brownian displacement for protein
	void setpos(vector<double> pos); // sets coords = pos
	void newvel(vector<double> v); // sets velocity
	void nudge(double t); //moves coords by v*t
	void update(); // sets coords=newcoords
	bool nearcntr(); // checks if protein is near the central particle
	void shift(vector<double> cpos, vector<protein>& crowders);
	
	

	vector<double> mvVel(double t); // returns coords of v*t_remaining
	vector<double> PBCswitch(int crowds, int index); // returns periodic coords

	//debugging and outputs
	void energy(); // prints the energy
	
	
	//Access 
	vector<double> getv(string a); // coords, new coords, vel
	double getp(string a); // mass, radius
	bool NearBound(int i); // returns if near the boundary
	//vector<int> NNindex(protein& other);
	bool IsNN(protein& other, int index, bool remover, int rindex);
	
	 
	vector<int> getNNs(bool nn); //true=NNs, false=notNNs
	void setNNs(vector<int> nnns); // sets NNs and notNNs
	void setNNs(vector<int> nns, vector<int> nnns);
	
	
	
};

class mainP : public protein
{
private:
	bool escape; // if main is near q, true
	double crad;

public:
	mainP();
	mainP(double x,double y, double z, double ms, double rad, double centerrad);

	bool nearEsc();
	void upmain();

};

class central : public protein
{
private:
	bool centered;


};


#endif
