#include <stdio.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <time.h>

using namespace std;

#define BLOCK_SIZE 256

/*__global__ void NN_make(float * NNmat, float coords_list, int numC) 
{
	// for one crowder, do all numC-1 NNs
}

__global__ void start_sys(float * coords, int numC)
{

}*/

__global__ void init_RNG(curandState *rngStates,
							  const unsigned int seed, int xLen)
{
	int Row = blockIdx.y*blockDim.y + threadIdx.y;
	int Col = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = Row*xLen + Col;

	curand_init(seed+tid, tid, 0, &rngStates[tid]);
}

__global__ void rand_gen(curandState *state, float * rands, int xLen, int yLen)
{
	int Row = blockIdx.y*blockDim.y + threadIdx.y;
	int Col = blockIdx.x*blockDim.x + threadIdx.x;
	int idx=Row*xLen + Col;
	if(idx<yLen*xLen)
		rands[idx] = curand_uniform(&state[idx]);
}

void wrapper(int switcher)
{
	int seed=(int)time(NULL), xLen=100, yLen=10;
	curandState *d_state;
	curandState *state;
	float *rands;
	float *d_rands;

	rands=(float *)malloc(xLen*yLen*sizeof(float));
	state=(curandState *)malloc(xLen*yLen*sizeof(curandState));


	cudaMalloc((void **)&d_state, xLen*yLen*sizeof(curandState));
	cudaMalloc((void **)&d_rands, xLen*yLen*sizeof(float));

	cudaMemcpy(d_rands, rands, xLen*yLen*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_state, state, xLen*yLen*sizeof(curandState), cudaMemcpyHostToDevice);



	dim3 Grid((xLen-1)/32 + 1, (yLen-1)/32 +1);
	dim3 Blocks(32,32) ;
	init_RNG<<<Grid,Blocks>>>(d_state, seed, xLen);
	cudaMemcpy(state, d_state, xLen*yLen*sizeof(curandState), cudaMemcpyDeviceToHost);
	rand_gen<<<Grid,Blocks>>>(d_state, d_rands, xLen, yLen);
	cudaMemcpy(rands, d_rands, xLen*yLen*sizeof(float), cudaMemcpyDeviceToHost);
	for(int i=0;i<xLen*yLen;i+=10)
		printf("random number %d is %f.\n",i,rands[i]);

	cudaFree(d_state);
	cudaFree(d_rands);
	free(state);
	free(rands);
}

int main()
{
	wrapper(0);
	return 0;
}
/*void gpu_NN_start(float coords_list, float * NNlist, int numC) // wrapper to initiate kernels
{
	// coords_list needs to be 1D list of coords
	// allocate host and device memory
	// copy coords list to gpu
	// set up grid and blocks
	// execute kernel
	// copy nnlist back
	// free memory
	// exit
}*/