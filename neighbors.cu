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
	
	__shared__ int coords_mz[300];
	__shared__ int coords_my[300];
	__shared__ int coords_mx[300];
	__shared__ int coords_pz[300];
	__shared__ int coords_py[300];
	__shared__ int coords_px[300];
	if(index<300){
		coords_mx[index]=coords[index*3]-2*L;
		coords_my[index]=coords[index*3+1]-2*L;
		coords_mz[index]=coords[index*3+2]-2*L;
		coords_px[index]=coords[index*3]+2*L;
		coords_py[index]=coords[index*3+1]+2*L;
		coords_pz[index]=coords[index*3+2]+2*L;
	

		__syncthreads();
		

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

		for(int i=0;i<N;i++) // check for neigbors through the boundaries. 
		{
			//-x
			x2=coords_mx[i];
			y2=coords[3*i+1];
			z2=coords[3*i+2];
			dx=x2-x1;dy=y2-y1;dz=z2-z1;
			mag2=dx*dx+dy*dy+dz*dz;
			if(i!=index)
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
					NNs[10*index+remove_index]=i+N;
					mag2s[remove_index]=mag2;
				}
				//-y////////////////////////////////////////////////////////////////
				////////////////////////////////////////////////////////////////////
				y2=coords_my[i];
				x2=coords[3*i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+2*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				y2=coords[3*i+1];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+3*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x
				////////////////////////////////////////////////////////////////
				z2=coords[3*i+2];
				x2=coords_px[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+4*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+y
				////////////////////////////////////////////////////////////////
				y2=coords_py[i];
				x2=coords[3*i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+5*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+z
				////////////////////////////////////////////////////////////////
				y2=coords[3*i+1];
				z2=coords_pz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+6*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x + y
				////////////////////////////////////////////////////////////////
				y2=coords_py[i];
				x2=coords_px[i];
				z2=coords[3*i+2];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+7*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x -y
				////////////////////////////////////////////////////////////////
				y2=coords_my[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+8*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x +y
				////////////////////////////////////////////////////////////////
				y2=coords_py[i];
				x2=coords_mx[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+9*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x -y
				////////////////////////////////////////////////////////////////
				y2=coords_my[i];
				x2=coords_mx[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+10*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x +z
				////////////////////////////////////////////////////////////////
				y2=coords[3*i+1];
				x2=coords_px[i];
				z2=coords_pz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+11*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x -z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+12*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x +z
				////////////////////////////////////////////////////////////////
				x2=coords_mx[i];
				z2=coords_pz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+13*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x -z
				////////////////////////////////////////////////////////////////
				x2=coords_mx[i];
				z2=coords_mz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+14*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+y +z
				////////////////////////////////////////////////////////////////
				y2=coords_py[i];
				x2=coords[3*i];
				z2=coords_pz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+15*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+y -z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+16*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-y +z
				////////////////////////////////////////////////////////////////
				y2=coords_my[i];
				z2=coords_pz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+17*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-y -z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+18*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x+y+z
				////////////////////////////////////////////////////////////////
				z2=coords_pz[i];
				y2=coords_py[i];
				x2=coords_px[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+19*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x-y-z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				y2=coords_my[i];
				x2=coords_mx[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+20*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x-y-z
				////////////////////////////////////////////////////////////////
				x2=coords_px[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+21*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x+y-z
				////////////////////////////////////////////////////////////////
				x2=coords_mx[i];
				y2=coords_py[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+22*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x-y+z
				////////////////////////////////////////////////////////////////
				z2=coords_pz[i];
				y2=coords_my[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+23*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x+y-z
				////////////////////////////////////////////////////////////////
				z2=coords_mz[i];
				y2=coords_py[i];
				x2=coords_px[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+24*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//+x-y+z
				////////////////////////////////////////////////////////////////
				z2=coords_pz[i];
				y2=coords_my[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+25*N;
					mag2s[remove_index]=mag2;
				}
				////////////////////////////////////////////////////////////////
				//-x+y+z
				////////////////////////////////////////////////////////////////
				x2=coords_mx[i];
				y2=coords_py[i];
				dx=x2-x1;dy=y2-y1;dz=z2-z1;
				mag2=dx*dx+dy*dy+dz*dz;		
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
					NNs[10*index+remove_index]=i+26*N;
					mag2s[remove_index]=mag2;
				}
			}			
		}
	}
}


void make_NNs(float *coords, int *NN, float params[])
{
	printf("how about here????");
	float space=params[0];
	float length=params[1];
	int N=params[2], n=params[3];

	float *device_coords;
	int *device_NNs; // store X neighbors per crowder?
	//int *system_NNs;
	int gridsize=ceil(pow(N,0.5));
	cudaMalloc((void **) &device_coords, N*3*sizeof(float));
	cudaMalloc((void **) &device_NNs, N*10*sizeof(int));

	cudaMemcpy(device_coords, coords, N*3*sizeof(float), cudaMemcpyHostToDevice);
	//system_NNs=(int *)malloc(10*N*sizeof(int));

	dim3 Grid(gridsize, gridsize, 1);
	dim3 Blocks(16,16);

	NN_lists<<<1, N>>>(device_coords, device_NNs, N, n, length);
	cudaMemcpy(NN, device_NNs, N*10*sizeof(int), cudaMemcpyDeviceToHost);

	printf("Crowder #1's nearest neighbors are\n");
	
	for(int i=0;i<10;i++)
	{
		printf("%i\n", NN[i]);
	}	
	cudaFree(device_coords);
	cudaFree(device_NNs);
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

__global__ void init_RNG(curandState *rngStates,
							  const unsigned int seed, int n)
{
	int Row = blockIdx.y*blockDim.y + threadIdx.y;
	int Col = blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int tid = Row*n + Col;

	curand_init(seed+tid, tid, 0, &rngStates[tid]);
}

__global__ void rand_gen(curandState *state, float * rands, int n, float sigma)
{
	int Row = blockIdx.y*blockDim.y + threadIdx.y;
	int Col = blockIdx.x*blockDim.x + threadIdx.x;
	int idx=Row*n + Col;
	float x,y,z,mag2=0,mag, random;
	if(idx<n*n)
	{
		random = curand_normal(&state[idx]);
		random=random*sigma;
		x=curand_uniform(&state[idx]);
		y=curand_uniform(&state[idx]);
		z=curand_uniform(&state[idx]);
		mag2=x*x+y*y+z*z;
		mag=sqrt(mag2);
		x=x/mag;y=y/mag;z=z/mag;
		rands[3*idx]=x*random;
		rands[3*idx+1]=y*random;
		rands[3*idx+2]=z*random;
	}
}

void rand_wrapper(int num, float *steparr, float sigma)
{
	int seed=(int)time(NULL);
	curandState *d_state;
	curandState *state;
	int n = ceil(pow(num,0.5));
	//float *rands;
	float *d_steps;

	//rands=(float *)malloc(xLen*yLen*sizeof(float));
	state=(curandState *)malloc(num*sizeof(curandState));


	cudaMalloc((void **)&d_state, num*sizeof(curandState));
	cudaMalloc((void **)&d_steps, 3*num*sizeof(float));

	cudaMemcpy(d_steps, steparr, 3*num*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_state, state, num*sizeof(curandState), cudaMemcpyHostToDevice);



	dim3 Grid((n/32) + 1, (n/32) + 1);
	dim3 Blocks(32,32) ;
	init_RNG<<<Grid,Blocks>>>(d_state, seed, n);
	cudaMemcpy(state, d_state, num*sizeof(curandState), cudaMemcpyDeviceToHost);
	rand_gen<<<Grid,Blocks>>>(d_state, d_steps, n, sigma);
	cudaMemcpy(steparr, d_steps, 3*num*sizeof(float), cudaMemcpyDeviceToHost);
	for(int i=0;i<3*num;i+=10)
		printf("random number %d is %f.\n",i,steparr[i]);

	cudaFree(d_state);
	cudaFree(d_steps);
	free(state);
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