#include <stdio.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <time.h>
//#include <vector.h>

using namespace std;

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define BLOCK_SIZE 256

/*__global__ void NN_make(float * NNmat, float coords_list, int numC) 
{
	// for one crowder, do all numC-1 NNs
}

__global__ void start_sys(float * coords, int numC)
{

}*/

/*__global__ void init_RNG(curandState *rngStates,
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
		rands[idx] = curand_normal(&state[idx]);
}*/

/*__global__ void coord_ICs(float *coords, float params[])
{
	float space=params[0], Len=params[1];
	int N=params[2], n=params[3];
	int end=params[4]; // gives the index of the last element of params
	int indices[3];
	indices[0]=blockIdx.x*blockDim.x + threadIdx.x ;
	indices[1]=blockIdx.y*blockDim.y + threadIdx.y ;
	indices[2]=blockIdx.z*blockDim.z + threadIdx.z ;
	int arrIdx=(indices[0] + indices[1]*n + indices[2]*n*n)*3;
	if(indices[0]<n && indices[1]<n && indices[2]<n)
		for(int i=0;i<3;i++)
		{
			coords[arrIdx+i]=(indices[i]+0.5)*space - Len;
		}

}*/

__global__ void NN_lists(float *coords, int *NNs, int N, int n, float L)
{
	int xdx, ydx, count;
	xdx=blockIdx.x*blockDim.x + threadIdx.x;
	ydx=blockIdx.y*blockDim.y + threadIdx.y;
	float x1=0, y1=0, z1=0, x2=0, y2=0, z2=0,dx,dy,dz;
	float mag2;
	int index=xdx + n*ydx; // n should be sqrt(N)
	//float cut=2*L/n; // the lattice spacing
	float mag2s[10];
	float difference, difference2;
	int remove_index;
	bool remove_bool=false, first;
	

	x1=coords[3*index];
	y1=coords[3*index+1];
	z1=coords[3*index+2];

	count=0;
	for(int i=0; i<N; i++)
	{	
		first=true;
		if(i!=index)
		{
			x2=coords[3*i];
			y2=coords[3*i+1];
			z2=coords[3*i+2];
			dx=x2-x1;dy=y2-y1;dz=z2-z1;
			mag2=dx*dx+dy*dy+dz*dz;
			if(count<10)
			{
				NNs[10*index+count]=i;
				mag2s[count]=mag2;
				count++;
			}
			else
			{	
				
				remove_bool=false;
				first=true;
				for(int j=0;j<10;j++)
				{
					if(mag2<mag2s[j])
					{	
						difference=mag2s[j]-mag2;
						remove_bool=true;
						if(first)
						{
							remove_index=j;
							difference2=difference;
							first=false;
						}
						else if(difference<difference2)
						{
							remove_index=j;
							difference2=difference;
						}
					}
				}
				if(remove_bool)
				{
					NNs[10*index+remove_index]=i;
					mag2s[remove_index]=mag2;
				}
			}
		}
	}
}



void make_NNs(float *coords, float params[])
{
	printf("how about here????");
	float space=params[0];
	float length=params[1];
	int N=params[2], n=params[3];

	float *device_coords;
	int *device_NNs; // store X neighbors per crowder?
	int *system_NNs;
	int gridsize=ceil(pow(N,0.5));
	cudaMalloc((void **) &device_coords, N*3*sizeof(float));
	cudaMalloc((void **) &device_NNs, N*10*sizeof(int));

	cudaMemcpy(device_coords, coords, N*3*sizeof(float), cudaMemcpyHostToDevice);
	system_NNs=(int *)malloc(10*N*sizeof(int));

	dim3 Grid(gridsize, gridsize, 1);
	dim3 Blocks(16,16);

	NN_lists<<<Grid, Blocks>>>(device_coords, device_NNs, N, n, length);
	cudaMemcpy(system_NNs, device_NNs, N*10*sizeof(int), cudaMemcpyDeviceToHost);

	printf("Crowder #21's nearest neighbors are\n");
	for(int i=0;i<10;i++)
	{
		printf("%i\n", system_NNs[210+i]);
	}

	// cudaFree(device_coords);
	// cudaFree(device_NNs);
}



void crowd_build_wrap(float *coords, float params[])
{/*
	// params has spacing, total length, N, n = (N)^1/3
	//N needs to not include the subtracted points yet
	float *d_coords;
	float space=params[0], Len=params[1];
	int N=params[2], n=params[3];
	int end=params[4]; // gives the index of the last element of params
	//int N_red=0; // reduced number of proteins
	

	cudaMalloc((void **) &d_coords, N*3*sizeof(float));

	dim3 Grid((n-1)/4 +1, (n-1)/4 +1, (n-1)/4 + 1);
	dim3 Blocks(8,8,8);

	coord_ICs<<<Grid, Blocks>>>(d_coords, n, space, Len);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaMemcpy(coords, d_coords, N*3*sizeof(float), cudaMemcpyDeviceToHost));
	C
	printf("c=%f\n", coords[1]);*/
}

void rand_wrapper(int switcher)
{
	/*int seed=(int)time(NULL), xLen=100, yLen=10;
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
	free(rands);*/
}



/*int main()
{
	int N=10;
	int n=ceil(pow(N,0.33333));
	double L=10.0;
	double l= (double)2*L/n;
	float in[4]; // 0: space, 1: length, 2: N, 3: n
	in[0]=l;
	in[1]=L;
	in[2]=n*n*n;
	in[3]=n;
	float coords[n*n*n*3];
	for(int i=0;i<n*n*n*3;i++)
		coords[i]=0;
	crowd_build_wrap(coords, in);
	for(int i=0;i<n*n*n*3;i++)
		printf("c = %f\n",coords[i] );
	//rand_wrapper(0);
	//rand_wrapper(0);
	return 0;
}*/

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