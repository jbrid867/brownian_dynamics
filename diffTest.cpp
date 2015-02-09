#include "headers/brownsys.h"

int main()
{
	diffusion_sys test(10, 2);
	test.startNNs();
	test.updateNNs();
}