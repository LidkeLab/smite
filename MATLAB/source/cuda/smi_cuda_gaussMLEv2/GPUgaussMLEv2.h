/*!
 * \file GPUgaussMLEv2.h
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief Prototypes for the actual Cuda kernels.  
 */

#ifndef GPUGAUSSMLEV2_H
#define GPUGAUSSMLEV2_H

__global__ void kernel_MLEFit_noshared(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

__global__ void kernel_MLEFit(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

__global__ void kernel_MLEFit_sigma(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

__global__ void kernel_MLEFit_z(const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

__global__ void kernel_MLEFit_sigmaxy(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

#endif