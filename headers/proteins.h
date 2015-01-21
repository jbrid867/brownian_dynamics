//Defines the protein and crowder classes

#include "constants1.h"
#ifndef _PROTEINS_H_
#define _PROTEINS_H_

using namespace std;

class protein
{
protected:
	vector<double> coordinates;
	vector<double> newcoords;
	vector<double> colCoords;
	vector<double> vel;
	vector<int> NNs;
	vector<int> notNNs;
	double mass;
	double radius;
	bool collision;
	//bool crowder; // tells if particle is a crowder.
	//bool escape, reaction;

public:
	
	//CONSTRUCTOR
	protein(); //default constructor. only contains mass and radius
	protein(double x,double y, double z, double ms, double rad); 
	//~protein(); //destructor
	
	//Dynamics
	
	void move(mt19937& gen, normal_distribution<> distro);
	void setpos(vector<double> pos); // sets coords = pos
	void newpos(vector<double> pos); // adds newcoords + pos
	void newvel(vector<double> v);
	void nudge(double t);
	void update();
	bool colCheck(int& index, double& t, vector<protein> crowds);
	void resolve(protein& two, double t_el, int index);

	vector<double> mvVel(double t);
	vector<double> PBCswitch(int crowds, int index);

	void energy();

	// event handling

	// double escTime();
	// void reacTime
	// void colTime
	// void resVel
	

	
	void NNout(); //Print NNs
	
	
	//vector<double> resolve(vector<double> p2); //resolves collision with another protein
	//vector<double> resolve(); //resolves collision with center
	
	
	//Access 
	vector<double> getv(string a); // coords, new coords, vel
	double getp(string a); // mass, radius
	
	 
	vector<int> getNNs(bool nn); //true=NNs, false=notNNs
	void setNNs(vector<int> nns, vector<int> nnns); // sets NNs and notNNs

	
	
	
	
	
	
	
	
	
};

class mainP : public protein
{
private:
	bool escape, reaction;
	double crad;

public:
	mainP();
	mainP(double x,double y, double z, double ms, double rad, double centerrad);

	
	void EscReac(double& t_el, vector<bool>& events);
	void collisions(double& t_el, vector<protein>& crowds, vector<bool>& events);
};

class central : public protein
{
private:
	bool centered;
};


#endif
