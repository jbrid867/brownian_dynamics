// Constants.h
#if !defined(MYLIB_CONSTANTS1_H)
#define MYLIB_CONSTANTS1_H 1
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <cstdio>
#include <fstream>
#include <limits>
#include <vector>
#include <omp.h>
//#include <constants.h>
#include <algorithm>
#include <random> //NEED TO COMPILE WITH -std=c++11 flag



const double k =1; //boltzman constant
const double e = 8.8542*pow(10,-12); //vacuum permittivity

const double r =1; //protein radius
const double rc =1; //crowder radius all that matters is scale


const double eta =0.66*pow(10,-3); //dynamic viscosity of water
const double pi = 3.1415926535897;
const double ke=pow(4*pi*e,-1);
const double T =310; //physiological temperature

const double Vp = (4.0/3.0)*pi*r*r*r;
const double Vc = (4.0/3.0)*pi*rc*rc*rc;

const double D = k*T/(6*pi*eta*r); //protein diffusion constant
const double Dc = k*T/(6*pi*rc*eta);//crowder diffusion constant


const double M =30000*1.67*pow(10,-27); //protein mass
const double Mc =30000*1.67*pow(10,-27); //crowder mass
const double b =3*pow(10,-10); //inner region radius
const double q =8*pow(10,-10); //outer region radius
const double phi =0.1; //volume fraction of crowders
//const double L =q; // simulation box size for reaction 
//const int N=floor(phi*pow(2.0*L,3.0)/((4.0/3.0)*pi*pow(rc,3.0)));
const int N=300; // for diffusion
const double L=pow(N*Vc/phi,0.333333333); // for diffusion

const double cut=6*rc;


const double vrms = pow(3*k*T/M,0.5); //rms velocity
const double h = .000001; //time step
const int dim=3; //DIMENSION
#endif

