#include <stdio.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath.h>

using namespace std;

#define BLOCK_SIZE 256

__global__ void NN_make(float * NNmat, float coords_list, int numC) 
{
	// for one crowder, do all numC-1 NNs
}

void gpu_NN_start(float coords_list, float * NNlist, int numC) // wrapper to initiate kernels
{
	// coords_list needs to be 1D list of coords
	// allocate host and device memory
	// copy coords list to gpu
	// set up grid and blocks
	// execute kernel
	// copy nnlist back
	// free memory
	// exit
}