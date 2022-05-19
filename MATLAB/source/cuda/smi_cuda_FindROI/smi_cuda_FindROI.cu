
//#include "cuda_runtime.h"
#include "definitions.h"
//#include "kernel.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>

        
        
__global__ void kernel_gaussMajor(float * d, int MajorSize, float b0, float b1, float b2, float b3, float B)
{
	//Gaussian filter along the Major dimension.  
	int MinorSize = blockDim.x;
	int idx = threadIdx.x;
	int idz = blockIdx.x;

	float w0, w1, w2, w3;
	float temp;
	int ii = 0;
	const int base = idz*MinorSize*MajorSize + idx*MajorSize;

	//forward
	w1 = w2 = w3 = d[base];
	for (ii = 0; ii<MajorSize; ii++)
	{
		w0 = d[base + ii];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}

	//backward
	w1 = w2 = w3 = d[base + MajorSize - 1];
	for (ii = MajorSize - 1; ii >= 0; ii--)
	{
		w0 = d[base + ii];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}
}


__device__ float kernel_norm(int x, float s1)
{
    return 1/sqrt(2*pi)/s1*exp( - x*x/(2*s1*s1));
}

__global__ void kernel_gaussMajor_sCMOS(const float * d, const float * v, float * d_out, int MajorSize, float Sigma)
{
	//Gaussian filter along the Major dimension.  
    // Noise-Weighted convolution by Gaussians, seperable
    // *d is the data
    // *v is a variance image
    // Sigma is the small kernel
         
    int MinorSize = blockDim.x;
	int idx = threadIdx.x;
	int idz = blockIdx.x;
    
    int st,en; 
	float weight, var, weightsum, varsum;
	int ii = 0;
	const int base = idz*MinorSize*MajorSize + idx*MajorSize;
    const int basev = idx*MajorSize;
    float winsize = 3* Sigma;

    // Each thread does one row.  
    // variance weighted  x =         
    for (ii = 0; ii< MajorSize; ii++){
        st = max(0,floor(ii-winsize));
        en = min(MajorSize-1,ii+winsize) ;
        varsum = 0;
        weightsum = 0;
        for (int jj = st; jj < en; jj++){
            var=kernel_norm(ii-jj,Sigma)/v[basev+jj];
            weight=var*d[base+jj];
            varsum += var;
            weightsum += weight;      
        }
        d_out[base+ii]=weightsum/varsum;
    }
}

__global__ void kernel_gaussMinor(float * d, int MinorSize, float b0, float b1, float b2, float b3, float B)
{
	//this kernel does gaussian filter along the Minor dimension.  
    //Note: faster to permute and use kernel_gaussMajor
	int MajorSize = blockDim.x;
	int idy = threadIdx.x;
	int idz = blockIdx.x;

	float w0, w1, w2, w3;
	float temp;
	int ii = 0;
	const int base = idz*MinorSize*MajorSize + idy;

	//forward
	w1 = w2 = w3 = d[base];
	for (ii = 0; ii<MinorSize; ii++)
	{
		w0 = d[base + ii*MajorSize];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii*MajorSize] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}

	//backward
	w1 = w2 = w3 = d[base + MajorSize*(MinorSize - 1)];
	for (ii = MinorSize - 1; ii >= 0; ii--)
	{
		w0 = d[base + ii*MajorSize];
		temp = w0*B + (b1*w1 + b2*w2 + b3*w3) / b0;
		d[base + ii*MajorSize] = temp;
		w3 = w2;
		w2 = w1;
		w1 = temp;
	}
}

__global__ void kernel_subtract(float * d_A, float * d_B)
{
	//Subtract second array from first array
    //Note: faster to subtract gpuArrays in MATLAB       
	int Xsize = blockDim.x;
	int Ysize = gridDim.x;
	int idx = threadIdx.x + Xsize*blockIdx.x + Xsize*Ysize*blockIdx.y;
	d_A[idx] = d_A[idx] - d_B[idx];
}

__global__ void kernel_LocalMaxFirstPass(const float * d_A, float * d_B, const int kernelsz, const float minval)
{
	//this kernel does max finding along the Major dimension.  
            //y-major, x-minor
	int MajorSize = blockDim.x;
	int MinorSize = gridDim.x;
	
    int x = blockIdx.x;  //minor
	int y = threadIdx.x; //major
	int z = blockIdx.y;

	//this is the pixel that we are searching around
	int idx = y + MajorSize*x + MinorSize*MajorSize*z;
	
	//define search only up to edges   
    int start = fmaxf(0, y - kernelsz);
	int end = fminf(MajorSize - 1, y + kernelsz);

	float maxval = minval;
	float inpixel = d_A[idx];

	for (int ii = start; ii<=end; ii++) 
		maxval = fmaxf(maxval, d_A[ii+ MajorSize*x + MinorSize*MajorSize*z]);

	//if any other pixel is larger set pixel idx to negative of that value, otherwise keep
    if (maxval>inpixel) d_B[idx]=-maxval;else d_B[idx]=maxval;

}

__global__ void kernel_LocalMaxSecondPass(const float * d_A, float * d_B, const int kernelsz, const float minval)
{
	//Max finding in second dimension, but input array should be permuted
    //so that operation is done along major axis.         
	int MajorSize = blockDim.x; 
	int MinorSize = gridDim.x;

	int x = blockIdx.x;  //minor
	int y = threadIdx.x; //major
	int z = blockIdx.y;

	//this is the pixel that we are searching around
	int idx = y + MajorSize*x + MinorSize*MajorSize*z;

	//define search only up to edges   
    int start = fmaxf(0, y - kernelsz);
	int end = fminf(MajorSize - 1, y + kernelsz);

	float maxval = minval;
	float inpixel = d_A[idx];

	//find the maximum absolute value in the filter window
	for (int ii = start; ii<=end; ii++)
		maxval = fmaxf(maxval, fabsf(d_A[ii+ MajorSize*x + MinorSize*MajorSize*z]));

	//if our pixel under test is equal to maximum, then flag that with '1', otherwise '0'
    if (fabsf(maxval-inpixel)<1e-6)d_B[idx] = (float)1;else d_B[idx] = (float)0;
	
}

