/*!
 * \file GPUgaussMLEv2.cu
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief This file contains all of the Cuda kernels.  The helper functions
 * are defined in GPUgaussLib.cuh
 */

#include "definitions.h"
#include "MatInvLib.h"
#include "GPUgaussLib.cuh"
#include "GPUgaussMLEv2.h"

//*******************************************************************************************
//theta is: {x,y,N,bg}
__global__ void kernel_MLEFit_XYNB_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[NV_P*NV_P], Diag[NV_P], Minv[NV_P*NV_P];
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    const int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_P;
    float dudt[NV_P];
    float d2udt2[NV_P];
    float NR_Numerator[NV_P], NR_Denominator[NV_P];
    float theta[NV_P];
    float maxjump[NV_P]={1e0f, 1e0f, 1e2f, 2e0f};
    float Nmax;

    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_P*NV_P*sizeof(float));
	memset(Minv,0,NV_P*NV_P*sizeof(float));
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma); //Added 2* on 8.9.16 to account for smoothing filter.
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
		memset(NR_Numerator,0,NV_P*sizeof(float));
		memset(NR_Denominator,0,NV_P*sizeof(float));

        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
            PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], &d2udt2[1]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
         
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
        PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBS_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    
    //__shared__ float s_data[MEM];
    float M[NV_PS*NV_PS], Diag[NV_PS], Minv[NV_PS*NV_PS];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_PS;
    float dudt[NV_PS];
    float d2udt2[NV_PS];
    float NR_Numerator[NV_PS], NR_Denominator[NV_PS];
    float theta[NV_PS];
    float maxjump[NV_PS]={1e0f, 1e0f, 1e2f, 2e0f, 5e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_PS*NV_PS*sizeof(float));
	memset(Minv,0,NV_PS*NV_PS*sizeof(float));      
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV_PS*sizeof(float));
		memset(NR_Denominator,0,NV_PS*sizeof(float));
      
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], &d2udt2[4]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain Sigma
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]); //bug fix 8.9.16 
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]); //bug fix 8.9.16 
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBZ_(const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[5*5], Diag[5], Minv[5*5];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=5;
    float dudt[5];
    float d2udt2[5];
    float NR_Numerator[5], NR_Denominator[5];
    float theta[5];
    float maxjump[5]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
    theta[4]=0;
   
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, dudt, d2udt2);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating remaining derivatives
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -maxjump[4]), maxjump[4]);
        
        // Any other constraints
        theta[2]=max(theta[2], 1.0f);
        theta[3]=max(theta[3], 0.01f);
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay, Bx,By, gamma, d, &PSFx, &PSFy, dudt, NULL);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating remaining derivatives
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
       
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) 
    kernel_MatInvN(M, Minv, Diag, NV);
  
   //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_XYNBSXSY_(const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
 
    //__shared__ float s_data[MEM];
    float M[6*6], Diag[6], Minv[6*6];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=6;
    float dudt[6];
    float d2udt2[6];
    float NR_Numerator[6], NR_Denominator[6];
    float theta[6];
    float maxjump[6]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f,1e-1f};
    float Nmax;
    
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    
	//initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    theta[5]=PSFSigma;
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            
            //calculating derivatives
   
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], &d2udt2[4]);
            kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], &d2udt2[5]);
            
            
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=data/model-1;
            if (model>10e-3f) df=data/pow(model, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);
        theta[5]-=min(max(NR_Numerator[5]/NR_Denominator[5], -theta[5]), theta[5]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain SigmaX
        theta[5]=max(theta[5], 0.5f); //Constrain Sigma
        theta[5]=min(theta[5], sz/2.0f); //Constrain SigmaX
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], NULL);
        kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
   
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

// SCMOS Versions---------------------------------------

__global__ void kernel_MLEFit_SCMOSXYNB_(const float *d_data, const float *d_Coords, const float *d_GainRatio, 
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[NV_P*NV_P], Diag[NV_P], Minv[NV_P*NV_P];
    const int tx = threadIdx.x;
    const int bx = blockIdx.x;
    const int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_P;
    float dudt[NV_P];
    float d2udt2[NV_P];
    float NR_Numerator[NV_P], NR_Denominator[NV_P];
    float theta[NV_P];
    float maxjump[NV_P]={1e0f, 1e0f, 1e2f, 2e0f};
    float Nmax;
	float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_P*NV_P*sizeof(float));
	memset(Minv,0,NV_P*NV_P*sizeof(float));
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
	const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
		memset(NR_Numerator,0,NV_P*sizeof(float));
		memset(NR_Denominator,0,NV_P*sizeof(float));

        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
            PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			//GRind=(int)s_Coords[0];
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], &d2udt2[0]);//x
            kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], &d2udt2[1]);//y
            dudt[2] = PSFx*PSFy;// I
            d2udt2[2] = 0.0f;// I
            dudt[3] = 1.0f;// bg
            d2udt2[3] = 0.0f;// bg
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR); // add variance-gain ratio: v/g^2
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2); // add variance-gain ratio: v/g^2
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], PSFSigma);
        PSFy=kernel_IntGauss1D(jj, theta[1], PSFSigma);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], PSFSigma, theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], PSFSigma, theta[2], PSFx, &dudt[1], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/(model+gainR);// add gain ratio
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;// add gain ratio
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}


//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBS_(const float *d_data, const float *d_Coords, const float *d_GainRatio,
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view.
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    
    //__shared__ float s_data[MEM];
    float M[NV_PS*NV_PS], Diag[NV_PS], Minv[NV_PS*NV_PS];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=NV_PS;
    float dudt[NV_PS];
    float d2udt2[NV_PS];
    float NR_Numerator[NV_PS], NR_Denominator[NV_PS];
    float theta[NV_PS];
    float maxjump[NV_PS]={1e0f, 1e0f, 1e2f, 2e0f, 5e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV_PS*NV_PS*sizeof(float));
	memset(Minv,0,NV_PS*NV_PS*sizeof(float));      
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV_PS*sizeof(float));
		memset(NR_Denominator,0,NV_PS*sizeof(float));
      
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], &d2udt2[4]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
        // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain Sigma
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[4]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[4], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss2DSigma(ii, jj, theta[0], theta[1], theta[4], theta[2], PSFx, PSFy, &dudt[4], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/(model+gainR);
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBZ_(const float *d_data, const float *d_Coords, const float *d_GainRatio, const float *d_x0,
		const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
    //__shared__ float s_data[MEM];
    float M[5*5], Diag[5], Minv[5*5];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=5;
    float dudt[5];
    float d2udt2[5];
    float NR_Numerator[5], NR_Denominator[5];
    float theta[5];
    float maxjump[5]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
	const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	const float *z_initial = d_x0+(bx*BlockSize+tx);
    //initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*pi*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
    theta[4]=z_initial[0];
   
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, dudt, d2udt2);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating remaining derivatives
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -maxjump[4]), maxjump[4]);
        
        // Any other constraints
        theta[2]=max(theta[2], 1.0f);
        theta[3]=max(theta[3], 0.01f);
        
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        kernel_DerivativeIntGauss2Dz(ii, jj, theta, PSFSigma_x,PSFSigma_y, Ax,Ay, Bx,By, gamma, d, &PSFx, &PSFy, dudt, NULL);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating remaining derivatives
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
       
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) 
    kernel_MatInvN(M, Minv, Diag, NV);
  
   //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

//*******************************************************************************************
__global__ void kernel_MLEFit_SCMOSXYNBSXSY_(const float *d_data, const float *d_Coords, const float *d_GainRatio, 
	    const float PSFSigma, const int sz, const int Mapsz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
	/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_Coords array of subregions's pixel coordinates in original field of view.
	 * \param d_GainRatio calibration result of variance/gain^2 in each pixel of original field of view. 
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param Mapsz size of original field of view.
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */
 
    //__shared__ float s_data[MEM];
    float M[6*6], Diag[6], Minv[6*6];
    int tx = threadIdx.x;
    int bx = blockIdx.x;
    int BlockSize = blockDim.x;
    int ii, jj, kk, ll;
    float model, cf, df, data;
    float Div;
    float PSFy, PSFx;
    int NV=6;
    float dudt[6];
    float d2udt2[6];
    float NR_Numerator[6], NR_Denominator[6];
    float theta[6];
    float maxjump[6]={1e0f, 1e0f, 1e2f, 2e0f, 1e-1f,1e-1f};
    float Nmax;
    float gainR;
	int GRind;
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;
    
	memset(M,0,NV*NV*sizeof(float));
	memset(Minv,0,NV*NV*sizeof(float));      
    
    //load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);
    const float *s_Coords = d_Coords+(2*bx*BlockSize+2*tx);
	//initial values
    kernel_CenterofMass2D(sz, s_data, &theta[0], &theta[1]);
    kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &theta[3]);
    theta[2]=max(0.0f, (Nmax-theta[3])*2*2*pi*PSFSigma*PSFSigma);
    theta[4]=PSFSigma;
    theta[5]=PSFSigma;
    for (kk=0;kk<iterations;kk++) {//main iterative loop
        
        //initialize
  		memset(NR_Numerator,0,NV*sizeof(float));
		memset(NR_Denominator,0,NV*sizeof(float));
        
        for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
            PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
            PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
            
            model=theta[3]+theta[2]*PSFx*PSFy;
            data=s_data[sz*jj+ii];
            GRind=(s_Coords[1]+jj)*Mapsz+s_Coords[0]+ii;
			gainR=d_GainRatio[GRind];
            //calculating derivatives
            kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], &d2udt2[0]);
            kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], &d2udt2[1]);
            kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], &d2udt2[4]);
            kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], &d2udt2[5]);
            dudt[2] = PSFx*PSFy;
            d2udt2[2] = 0.0f;
            dudt[3] = 1.0f;
            d2udt2[3] = 0.0f;
            
            cf=0.0f;
            df=0.0f;
            if (model>10e-3f) cf=(data-model)/(model+gainR);
            if (model>10e-3f) df=(data+gainR)/pow(model+gainR, 2);
            cf=min(cf, 10e4f);
            df=min(df, 10e4f);
            
            for (ll=0;ll<NV;ll++){
                NR_Numerator[ll]+=dudt[ll]*cf;
                NR_Denominator[ll]+=d2udt2[ll]*cf-pow(dudt[ll], 2)*df;
            }
        }
        
         // The update: Editted 8.9.16 to allow for larger jumps in Photons/BG
        theta[0]-=min(max(NR_Numerator[0]/NR_Denominator[0], -maxjump[0]), maxjump[0]);
        theta[1]-=min(max(NR_Numerator[1]/NR_Denominator[1], -maxjump[1]), maxjump[1]);
        theta[2]-=min(max(NR_Numerator[2]/NR_Denominator[2], -theta[2]), theta[2]);
        theta[3]-=min(max(NR_Numerator[3]/NR_Denominator[3], -theta[3]), theta[3]);
        theta[4]-=min(max(NR_Numerator[4]/NR_Denominator[4], -theta[4]), theta[4]);
        theta[5]-=min(max(NR_Numerator[5]/NR_Denominator[5], -theta[5]), theta[5]);

        // Any other constraints
        theta[2]=max(theta[2], 1.0f); //Make sure Photons is postitve
        theta[3]=max(theta[3], 0.01f); //Make sure Background is postitve
        theta[4]=max(theta[4], 0.5f); //Constrain Sigma
        theta[4]=min(theta[4], sz/2.0f); //Constrain SigmaX
        theta[5]=max(theta[5], 0.5f); //Constrain Sigma
        theta[5]=min(theta[5], sz/2.0f); //Constrain SigmaX
    }
    
    // Calculating the CRLB and LogLikelihood
    Div=0.0f;
    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
        
        PSFx=kernel_IntGauss1D(ii, theta[0], theta[4]);
        PSFy=kernel_IntGauss1D(jj, theta[1], theta[5]);
        
        model=theta[3]+theta[2]*PSFx*PSFy;
        data=s_data[sz*jj+ii];
        
        //calculating derivatives
        kernel_DerivativeIntGauss1D(ii, theta[0], theta[4], theta[2], PSFy, &dudt[0], NULL);
        kernel_DerivativeIntGauss1D(jj, theta[1], theta[5], theta[2], PSFx, &dudt[1], NULL);
        kernel_DerivativeIntGauss1DSigma(ii, theta[0], theta[4], theta[2], PSFy, &dudt[4], NULL);
        kernel_DerivativeIntGauss1DSigma(jj, theta[1], theta[5], theta[2], PSFx, &dudt[5], NULL);
        dudt[2] = PSFx*PSFy;
        dudt[3] = 1.0f;
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= dudt[ll]*dudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if ((model+gainR)>0)
            if ((data+gainR)>0)Div+=(data+gainR)*log(model+gainR)-model-(data+gainR)*log(data+gainR)+data;
            else
                Div+=-model-gainR;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
   
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+BlockSize*bx+tx]=theta[kk];
    for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+BlockSize*bx+tx]=Diag[kk];
    d_LogLikelihood[BlockSize*bx+tx] = Div;
    return;
}

#define NUM_FIT 8
#define a1 -0.9911f
#define a2 -0.6763f
#define b1 0.8055f
#define b2 -1.2451f
	//********************************************************************************************************************************************
__global__ void kernel_gaussMFA(const float *d_data, const float PSFSigma, 
        const int sz, const int iterations,const int num, const float Nave, 
                 const float llThreshold, float *d_X, 
                        float *d_Y, float *d_b, float *d_PValue, int Nfits) 
{
	
	float Nmax;

	int tx = threadIdx.x; //fits
	int bx = blockIdx.x;
    int BlockSize = blockDim.x;

	int ii, jj, kk, qq, mm, nn, nnM, goodness;
	int nnfit=0;
	int breaktmp, nncount, rimsigny, rimsignx;
	float dLLx, dLLy, dLLb;
	float deflamodel;
	float tmp, signx, signy, meanx, meany;
	float b, x, y, xmin, xmax, ymin, ymax;
	float model, cf, df, data,bini;
	float Div;
	float bfit=0;
	float imdy, imdx,  imddx, imddy;
	float PSFy, PSFx, pval, zval,minDiv,maxpval;
	float stepxtot;
	float stepytot;
	float stepbgtot;
	float xarray[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};
	float yarray[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};
	float xarrayfit[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};
	float yarrayfit[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};
	float x_current[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};
	float y_current[NUM_FIT]={0, 0, 0, 0, 0, 0, 0, 0};

	//load data
    const float *s_data = d_data+(sz*sz*bx*BlockSize+sz*sz*tx);

	//Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	//initial values
	kernel_CenterofMass2D(sz, s_data, &x, &y);
	kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &b);
	bini=max( b, 0.0001f);

	PSFy = 0.0f; PSFx=0.0f;

	maxpval=-1.0f; //this keeps track of models' pvalue
	
    //loop over number of emitters
	for (mm=1;mm<=num;mm++) {
		breaktmp=0;
		goodness=0;
		b=bini;

		//INITIAL GUESSES

		//if fitting for a single fluorophore, use center of mass as initial guess
		if (mm==1){
			xarray[0]=x;
			yarray[0]=y;
		}

		//if fitting more than one fluorophore, use deflation method and find max
		else {
			tmp=Nave/2/pi/PSFSigma/PSFSigma/10;

			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
				model=b;
				nncount=0;
				meanx=0;
				meany=0;
				for (nn=0;nn<mm-1;nn++) {
					x=xarray[nn];
					y=yarray[nn];
					meanx+=xarray[nn];
					meany+=yarray[nn];
					nncount=nncount+1;
					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);					
					model+=Nave*PSFx*PSFy;
				}
				
				deflamodel=s_data[sz*jj+ii]-model;
				signx=0;
				signy=0;
				rimsignx=1;
				rimsigny=1;
				if (deflamodel>tmp) {
					if (ii==0) rimsignx=-1;
					if (ii==sz-1) rimsignx=-1;
					if (jj==0) rimsigny=-1;
					if (jj==sz-1) rimsigny=-1;

					signx=((ii-meanx/nncount) > 0)? 1 : -1;
					signy=((jj-meany/nncount) > 0)? 1 : -1;
					xarray[mm-1]=ii-0.5f*signx*rimsignx;
					yarray[mm-1]=jj-0.5f*signy*rimsigny;
					tmp=deflamodel;
				}
			}
			if (tmp==Nave/2/pi/PSFSigma/PSFSigma/10)
				breaktmp=1;
		}

		if (breaktmp==1)
			break;

		//MAIN ITERATIVE LOOP
		for (kk=0;kk<iterations;kk++){

			dLLb=0.0f;
			stepbgtot = 0;

			//generating the model and calc b update
			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

				model=b;
				for (nn=0;nn<mm;nn++) {
					x_current[nn]=xarray[nn];
					y_current[nn]=yarray[nn];

					x=xarray[nn];
					y=yarray[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					model+=Nave*PSFx*PSFy;
				}
				data=s_data[sz*jj+ii];
				cf=0.0f;
				df=0.0f;
				if (model>10e-3f) cf=data/model-1;
				if (model>10e-3f) df=data/pow(model, 2);
				cf=min( cf, 10e4f);
				df=min( df, 10e4f);
				stepbgtot += -df;
				dLLb+=cf;
			}
			b-=min( max( dLLb/stepbgtot, -1e0f), 1e0f)/mm/2;
			b=max( b, 0.001f);

			//This starts iterative routine for theta_i other than theta_bg which is calculated above.
			for (nn=0;nn<mm;nn++){

				dLLx=0.0f;
				dLLy=0.0f;
				stepxtot = 0.0f;
				stepytot = 0.0f;

				for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
					//generate model using initial value or iteration value.
					data=s_data[sz*jj+ii];
					//data=s_data[sz*sz*tx+sz*jj+ii];
					model=b;
					for (nnM=0;nnM<mm;nnM++) {
						x=x_current[nnM];
						y=y_current[nnM];
						PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
						PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
						model+=Nave*PSFx*PSFy;
					}

					cf=0.0f;
					df=0.0f;
					if (model>10e-3f) cf=data/model-1;
					if (model>10e-3f) df=data/pow( model, 2);
					cf=min( cf, 10e4f);
					df=min( df, 10e4f);

					x=x_current[nn];
					y=y_current[nn];

					PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
					PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
					
					kernel_DerivativeIntGauss1D(ii, x, PSFSigma, Nave, PSFy, &imdx, &imddx);
					kernel_DerivativeIntGauss1D(jj, y, PSFSigma, Nave, PSFx, &imdy, &imddy);

					//denominator
					stepxtot  += imddx*cf-pow( imdx, 2)*df;
					stepytot  += imddy*cf-pow( imdy, 2)*df;

					//numerator
					dLLx+=imdx*cf;
					dLLy+=imdy*cf;
				}

				x-=min( max( dLLx/stepxtot, -1e0f), 1e0f)/mm/2.0f;
				y-=min( max( dLLy/stepytot, -1e-0f), 1e-0f)/mm/2.0f;
				xarray[nn]=x;
				yarray[nn]=y;
			}
		}

		// calculating loglikelihood value
		Div=0.0f;
		xmin=1000;
		xmax=-1000;
		ymin=1000;
		ymax=-1000;

		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			// generating the model
			model=b;
			for (nn=0;nn<mm;nn++) {
				x=xarray[nn];
				y=yarray[nn];
                
                //reject this fit model if too far outside
				if (x>xmax) xmax=x;
				if (x<xmin) xmin=x;
				if (y>ymax) ymax=y;
				if (y<ymin) ymin=y;

				PSFx=kernel_IntGauss1D(ii, x, PSFSigma);
				PSFy=kernel_IntGauss1D(jj, y, PSFSigma);
				model+=Nave*PSFx*PSFy;
			}

			data=s_data[sz*jj+ii];

			if (data>0){
				Div+=-2*(data*log(model)-model-data*log(data)+data-0.5f*log(2*pi*data)-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}

			else{
				Div+=-2*(-model-model*log(model)+model+model*log(model)-model+0.5f*log(2*pi*model));}
		}

		zval=sqrt((float) Div)-sqrt((float) (sz*sz-2*mm-1));
		pval=(zval<0)*(1-0.5f*exp(b1*zval+a1*pow(zval,2)))+(zval>0)*(0.5f*exp(b2*zval+a2*pow(zval,2)));

        //if this is true, model is better than last
		if ((pval>maxpval) && (xmin>-round(1.5f*PSFSigma)) && (xmax <(sz-1+round(1.5f*PSFSigma))) && (ymin>-round(1.5f*PSFSigma)) && (ymax<(sz-1+round(1.5f*PSFSigma)))) 
        {
			for (qq=0;qq<mm;qq++)
			{
				xarrayfit[qq]=xarray[qq];
				yarrayfit[qq]=yarray[qq];
			}
			maxpval=pval;
			minDiv=Div;
			bfit=b;
			nnfit=mm;
		}
	}

	if ((xarrayfit[0]!=0)&&(maxpval>llThreshold*0.01)){
		b=bfit;
		for (qq=0;qq<nnfit;qq++)
		{
			xarray[qq]=xarrayfit[qq];
			yarray[qq]=yarrayfit[qq];
		}	
	}

	// output
	if (maxpval>llThreshold)
		goodness=1;

	int fitnum=num;

	for (nn=0;nn<fitnum;nn++)
		d_X[BlockSize*fitnum*bx+fitnum*tx+nn]=xarrayfit[nn]*goodness;

	for (nn=0;nn<fitnum;nn++)
		d_Y[BlockSize*fitnum*bx+fitnum*tx+nn]=yarrayfit[nn]*goodness;

	d_b[BlockSize*bx+tx]=bfit*goodness;
	d_PValue[BlockSize*bx+tx]=maxpval;
	return;
}
#define NMAX 5
#define MEM_CRLB (2*NMAX+1)
#define MATMEM (2*NMAX+1)*(2*NMAX+1)

__global__ void kernel_CRLB(const int sz, const int fitnum,const float Nave, 
							const float PSFSigma, const float *d_xarray, const float *d_yarray, 
							const float *d_barray, float *d_CRLBarray, 
							float *d_covariance, const int Nfits,  float *d_crlbx, float *d_crlby) 
{

	int tx = threadIdx.x; //fits
	int bx = blockIdx.x;
	int BlockSize = blockDim.x;

	int kk, ll, ii, jj,nn;
	float x, y;
	float PSFx, PSFy, imdx,imdy;
	float imddx, imddy;
	float imdbg;
	int fitstate, matsz;
	int matsize=0,numco;

	float model,disx,disy;//, *s_Diag=0, *s_Diag2=0;
	float s_Diag[MEM_CRLB]; //x,y,bg
	float s_temp[MATMEM];
	float s_fishermatrix[MATMEM];
	float zerofisher[MATMEM];

	for (kk=0;kk<MATMEM;kk++){
		s_fishermatrix[kk]=0;
		s_temp[kk]=0;
	}

	// data array locations
	const float *s_xarray = d_xarray+(fitnum*bx*BlockSize+fitnum*tx);
	const float *s_yarray = d_yarray+(fitnum*bx*BlockSize+fitnum*tx);
	const float *s_barray = d_barray+(bx*BlockSize+tx);

	// adding variables for CRLB_x and CRLB_y
	float *s_crlbx = d_crlbx+(fitnum*bx*BlockSize+fitnum*tx);
	float *s_crlby = d_crlby+(fitnum*bx*BlockSize+fitnum*tx);

	//Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nfits) return;

	//calculation starts.

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		// generate model for pixel ii jj

		model=s_barray[tx];

		for (nn=0;nn<fitnum;nn++){
			fitstate=0;
			x=s_xarray[nn];
			y=s_yarray[nn];
			if(x!=0) fitstate=1;
			PSFx=kernel_IntGauss1D(ii, x, PSFSigma)*fitstate;
			PSFy=kernel_IntGauss1D(jj, y, PSFSigma)*fitstate;
			model+=Nave*PSFx*PSFy;
		}

		matsize=0;
		for (nn=0;nn<fitnum;nn++){
			fitstate=0;
			x=s_xarray[nn];
			y=s_yarray[nn];
			if(x!=0) fitstate=1;
			PSFx=kernel_IntGauss1D(ii, x, PSFSigma)*fitstate;
			PSFy=kernel_IntGauss1D(jj, y, PSFSigma)*fitstate;
			
			kernel_DerivativeIntGauss1D(ii, x, PSFSigma, Nave, PSFy, &imdx, &imddx );
			kernel_DerivativeIntGauss1D(jj, y, PSFSigma, Nave, PSFx, &imdy, &imddy );
			imdx = imdx*fitstate;
			imdy = imdy*fitstate;

			//store in temp memory (which's allocate for LUDC method)
			s_temp[2*nn]=imdx;
			s_temp[2*nn+1]=imdy;
			matsize+=fitstate;
		}
		imdbg = 1.0f;
		// put derivitive of bg right after fitted parameters
		if (matsize>0)
			s_temp[2*matsize]=imdbg;

		for (kk=0;kk<(2*matsize+1);kk++)
			for (ll=0;ll<(2*matsize+1);ll++)
				s_fishermatrix[kk*(2*matsize+1)+ll]+=s_temp[kk]*s_temp[ll]/model;
	}

	for (kk=0;kk<MATMEM;kk++){

		s_temp[kk]=0;
	}
	/*-----------------------------------------------------------------------------------------------------------------------
	-----------------------------------------------------------------------------------------------------------------------*/
	//copy fisher information matrix into device memory for output.

	//modification of fisher information matrix
	matsz=matsize*2+1;
	if(matsize>0)
	{
		//copy fisher matrix
		for (ii=0;ii<matsz;ii++)
			for(kk=0;kk<matsz;kk++)
				zerofisher[kk*(2*matsize+1)+ii]=s_fishermatrix[kk*(2*matsize+1)+ii];

		//set covariance to zero
		for(kk=0;kk<matsize;kk++){
			for(jj=0;jj<matsize;jj++)
			{
				if(kk!=jj)
				{
					zerofisher[(2*kk+1)*(2*matsize+1)+(2*jj+1)]=0;
					zerofisher[(2*kk)*(2*matsize+1)+(2*jj)]=0;
				}
			}
		}

		// ********************************************************
		//inverse 0 fisher
		kernel_MatInvN(zerofisher, s_temp, s_Diag, matsz);

		//Calculate new fisher
		for(kk=0;kk<matsize;kk++)
		{
			s_crlbx[kk]=sqrt(fabs( s_temp[2*kk*matsz+2*kk]));        //crlbx
			s_crlby[kk]=sqrt(fabs( s_temp[(2*kk+1)*matsz+(2*kk+1)])); //crlby
		}


		for(kk=0;kk<matsize;kk++)
		{
			for(jj=0;jj<matsize;jj++)
			{
				if(kk!=jj)
				{
					disx=pow((fabs( (s_xarray[kk]-s_xarray[jj]))),2)/s_crlbx[kk]/s_crlbx[jj];
					disy=pow((fabs( (s_yarray[kk]-s_yarray[jj]))),2)/s_crlby[kk]/s_crlby[jj];

					s_fishermatrix[(2*kk+1)*(2*matsize+1)+(2*jj+1)]=disy/(disy+1)*s_fishermatrix[(2*kk+1)*(2*matsize+1)+(2*jj+1)];
					s_fishermatrix[(2*kk)*(2*matsize+1)+(2*jj)]=disx/(disx+1)*s_fishermatrix[(2*kk)*(2*matsize+1)+(2*jj)];

				}
			}

		}

		//inverse new fisher
		// matsize will be a indicator saying how big the matrix is for each thread. matsize*2+1 will be the scale of the matrix.
		// start calculating matrix inverse
		/*-----------------------------------------------------------------------------------------------------------------------
		-----------------------------------------------------------------------------------------------------------------------*/
		for (kk=0;kk<MATMEM;kk++){
			s_temp[kk]=0;
		}

		kernel_MatInvN(s_fishermatrix, s_temp, s_Diag, matsz);

	}
	// finished    
	// copy back to device memory
	numco=0;
	for (kk=0;kk<matsz;kk++)
	{
		d_CRLBarray[(2*fitnum+1)*bx*BlockSize+tx*(2*fitnum+1)+kk]=sqrt(fabs(s_temp[kk+kk*matsz]));
		if ((kk+1)%2==0)
		{
			d_covariance[fitnum*bx*BlockSize+tx*fitnum+numco]=s_temp[kk+(kk-1)*matsz];
			numco++;
		}
	}
	return;
}

