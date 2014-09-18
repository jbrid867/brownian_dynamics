#include "headers/brownsys.h"

int main(){
	random_device rd; // Seed RNG
	mt19937 gen(rd()); // start RNG
	normal_distribution<> distro(0,pow(2*D*h,0.5));
	brownsys baller(N);
	baller.startNNs(4*r);
	for(int i=0;i<5;i++){baller.moveall(gen,distro);}
	//shftcntr(baller); // seems clunky
	 // integrate with moveall eventually
	//baller.something();
	//baller.updateNNs(4*r);
	
	 
	
return 0;}

/* SOME PSUEDO

for the move-resolve-update stage

- move serially (possibly randomize the order)
- individually:
-- check for collisions, resolve
--- use NN list -> access other proteins in brownsys -> return their coords
--- check if there are collisions and which happened first
--- resolve earliest
--- if main, check for exit or reaction
--- check both particles involved (priority to main) for secondary collisions
--- resolve and repeat until dt -> h
-- check for periodic BCs, resolve
-- update coordinates for all
-- move on to next


*/
