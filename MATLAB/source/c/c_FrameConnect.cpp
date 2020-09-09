/*Compile the mex file:
 *
 *
 *		[X,SE,NConnected,StartFrame,Photons,BG,ConnectID]=c_frameConnect(LOS,X_in,SE_in,FrameNum,Photons_in,BG_in)
 *
 *
 *		Inputs must come in sorted by time
 */

#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include <math.h>

#define pi 3.141592
#define DEBUG

float lognormpdf(float X, float U, float S)
{
    //extra terms not needed in LogLikelihood ratio
    return powf((X - U) / S, 2) / 2; // -logf(S); // -0.5*logf(2 * pi);
    
}

//The integral of the Chi-square PDF is the regularized gamma function.  I'm hard-coding for P(k/2,x/2) for k = 1,2,3 (up to 3 dimensions)
float X2_CDF(int k, float x){
    float out;
    switch (k){
        case 1:
            out = erf(sqrtf(x / 2)); break;
        case 2:
            out = 1 - expf(-x / 2); break;
        case 3:  //from wolfram alpha
            out = (2 * (-1 / 2 * sqrtf(pi)* (1 - erff(sqrtf(x) / sqrtf(2))) - (expf(-x / 2) *sqrtf(x)) / sqrtf(2) + sqrtf(pi) / 2)) / sqrtf(pi); break;
    }
    return out;
}

float pvalue(float *X1, float *SE1, float *X2, float *SE2, int N, float *Xout, float *SEout){
    
    //MLE of combined value
    for (int nn = 0; nn < N; nn++)
    {
        Xout[nn] = (X1[nn] / powf(SE1[nn], 2) + X2[nn] / powf(SE2[nn], 2)) / (1 / powf(SE1[nn], 2) + 1/ powf(SE2[nn], 2));
        SEout[nn] = sqrt(1 / (1 / powf(SE1[nn], 2) + 1/ powf(SE2[nn], 2)));
    }
    //Likelihood at MLE
    float logL = 0;
    
    
    for (int nn = 0; nn < N; nn++)
        logL += lognormpdf(X1[nn], Xout[nn], SE1[nn]) + lognormpdf(X2[nn], Xout[nn], SE2[nn]);
    
    
    //Likelihood of independent points- this is 0
    /*
     * float logL0 = 0;
     * for (int nn = 0; nn < N; nn++)
     * logL0 += lognormpdf(X1[nn], X1[nn], SE1[nn]) + lognormpdf(X2[nn], X2[nn], SE2[nn]);
     */
    
    //Calculate likelihood ratio :
    float R = 2 * logL;
    
    //mexPrintf("%g %g %g \n",logL,logL0,R);
    //Degrees of Freedom
    int DOF = N;
    return 1 - X2_CDF(DOF, R);
   }

float Weighted_ave(float *P_ave1, float *P_aveSE1, float *P_ave2, float *P_aveSE2, int N,float *P_Aveout, float *P_aveSEout ){
    float out;
    
    for (int nn = 0; nn < N; nn++)
        
    {
        P_Aveout[nn] = (P_ave1[nn] / powf(P_aveSE1[nn], 2) + P_ave2[nn] / powf(P_aveSE2[nn], 2)) / (1 / powf(P_aveSE1[nn], 2) + 1/ powf(P_aveSE2[nn], 2));
        P_aveSEout[nn] = sqrt(1 / (1 / powf(P_aveSE1[nn], 2) + 1/ powf(P_aveSE2[nn], 2)));
    }
    out=1;
    return out;
}


void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const	mxArray	*prhs[])
{
    float LOS;			//Level of Significance
    float *X;		//NxNDim input locations
    float *X_SE;		//NxNDim input Standard Errors
    float *P_Ave;		//NxNDim_Ave input like PSFSigma
    float *P_Ave_SE;     //NxNDim_Ave input PSFSigma Standard Errors
    float *P_Add;       //Nx2 input like Photons and Background
    int N;				//Number of input localizations
    int NDim;			//Number of dimensions for spatial coordinates
    int NDim_Ave;        //Number of dimensions for parameter (weighted average)
    int NDim_Add;        //Number of dimensions for parameter (Added)
    float MaxDistance;	//Maximum Distance between accepted connections
    int MaxFrameGap;	//Maximum Frame distance between accepted fits. 1 means consecutive
    
    unsigned int *FN_new;		//Frame number of last connected point
    int	CurrentFrame;			//Frame number of particle under test
    float *X_new;				//updated position from connection
    float *P_Ave_new;		    //parameters (weighted ave)
    float *P_Ave_SE_new;				//parameters (weighted ave) Standard Error
    float *P_Add_new;			//parameters (added)
    float PValue, out; 
    
    float tmp;
    float *SE_new;
    float *Xout, *SEout, *P_Aveout, *P_Ave_SEout, *P_Addout;
    unsigned int *FrameNum, *N_combined, *Nout,  *Fout;
    unsigned int *ConnectID, MaxConID, *CombinedID;
    int cnt, kk, ii, jj, nn, cont;
    const mwSize *datasize, *errorsize, *framesize,*P_avesize,*P_ave_SEsize,*P_addsize;
    mwSize outsize[3];
    
    int MaxDim = 3;
    float Xtmp1[3], SEtmp1[3], Xtmp2[3], SEtmp2[3], Xtmpout[3], SEtmpout[3];
    
    // So far we only have two parameter_ave: PSFSigma or PSFSigmaX & PSFSigmaY
    float P_Avetmp1[2], P_AveSEtmp1[2], P_Avetmp2[2], P_AveSEtmp2[2], P_Avetmpout[2], P_AveSEtmpout[2];
    
    //input checks
    if (nrhs != 10)
        mexErrMsgTxt("Inputs are: LOS,Coord_in,Coord_SE_in,FrameNum,P_Ave_in,P_Ave_SE_in,P_Add_in,MaxDistance,MaxFrameGap,MaxConID\n");
    
    
    datasize = mxGetDimensions(prhs[1]);
    errorsize = mxGetDimensions(prhs[2]);
    framesize = mxGetDimensions(prhs[3]);

    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        P_avesize = mxGetDimensions(prhs[4]);
        P_ave_SEsize = mxGetDimensions(prhs[5]);
    }
    P_addsize = mxGetDimensions(prhs[6]);
    //input checks
    if (datasize[0] != errorsize[0])
        mexErrMsgTxt("Size of X must match size of X_Variance.\n");
    
    if (datasize[1] != errorsize[1])
        mexErrMsgTxt("Size of X must match size of X_Variance.\n");
    
    if (framesize[0] != datasize[0])
        mexErrMsgTxt("Size of FrameNum must match size of X.\n");
    
    if (framesize[1] != 1)
        mexErrMsgTxt("FrameNum must be N x 1.\n");
    
    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        if (P_avesize[0] != datasize[0])
            mexErrMsgTxt("Size of P_Ave must match size of X.\n");
        if (P_avesize[0] != P_ave_SEsize[0])
            mexErrMsgTxt("Size of P_Ave must match size of P_Ave_Variance.\n");
        if (P_avesize[1] != P_ave_SEsize[1])
            mexErrMsgTxt("Size of P_Ave must match size of P_Ave_Variance.\n");
    }
    if (P_addsize[0] != datasize[0])
        mexErrMsgTxt("Size of P_Add must match size of X.\n");
    
    //Check input data types
    if (mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
        mexErrMsgTxt("X must be comprised of single floats\n");
    if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
        mexErrMsgTxt("X_Variance must be comprised of single floats\n");
    if (mxGetClassID(prhs[3]) != mxUINT32_CLASS)
        mexErrMsgTxt("FrameNum must be comprised of uint32\n");
    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        if (mxGetClassID(prhs[4]) != mxSINGLE_CLASS)
            mexErrMsgTxt("Parameter (weigthed ave) must be comprised of single floats\n");
        if (mxGetClassID(prhs[5]) != mxSINGLE_CLASS)
            mexErrMsgTxt("Parameter_Variance (weigthed ave) must be comprised of single floats\n");
    }
    if (mxGetClassID(prhs[6]) != mxSINGLE_CLASS)
        mexErrMsgTxt("Parameter (Added) must be comprised of single floats\n");
    
    //Get Sizes
    N = (int)datasize[0];
    NDim = (int)datasize[1];
if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0){
    NDim_Ave = (int) P_avesize[1];
}
    NDim_Add = (int) P_addsize[1];
    if (NDim > MaxDim)
        mexErrMsgTxt("Maximum dimension is 3!\n");
// mexPrintf("NDIM_ADD==%d and N==%d",NDim_Add,N);
    
    //Get input data
    LOS = (float)mxGetScalar(prhs[0]);
    X = (float *)mxGetData(prhs[1]);
    X_SE = (float *)mxGetData(prhs[2]);
    FrameNum = (unsigned int *)mxGetData(prhs[3]);
    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        P_Ave = (float *)mxGetData(prhs[4]);
        P_Ave_SE = (float *)mxGetData(prhs[5]);
    }
    P_Add =(float *)mxGetData(prhs[6]);
    MaxDistance = (float)mxGetScalar(prhs[7]);
    MaxFrameGap = (int)mxGetScalar(prhs[8]);
    MaxConID = (int)mxGetScalar(prhs[9]);

    //Copy input values to new arrays for processing
    X_new = (float*)mxCalloc(N*NDim, sizeof(float));
    for (ii = 0; ii<N*NDim; ii++)X_new[ii] = X[ii];
    
    SE_new = (float*)mxCalloc(N*NDim, sizeof(float));
    for (ii = 0; ii<N*NDim; ii++)SE_new[ii] = X_SE[ii];
    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        P_Ave_new = (float*)mxCalloc(N*NDim_Ave, sizeof(float));
        for (ii = 0; ii<N*NDim_Ave; ii++)P_Ave_new[ii] = P_Ave[ii];
        
        P_Ave_SE_new = (float*)mxCalloc(N*NDim_Ave, sizeof(float));
        for (ii = 0; ii<N*NDim_Ave; ii++)P_Ave_SE_new[ii] = P_Ave_SE[ii];
    }
    P_Add_new = (float*)mxCalloc(N*NDim_Add, sizeof(float));
    for (ii = 0; ii<N*NDim_Add; ii++) {
        P_Add_new[ii] = P_Add[ii];
//     mexPrintf("ii=%d and P_ADd_new{ii}==%g \n",ii,P_Add_new[ii]);
    }
    //keep track of number of points used
    N_combined = (unsigned int*)mxCalloc(N, sizeof(unsigned int));
    for (ii = 0; ii<N; ii++)N_combined[ii] = 1;
    
    //keep record of how orginal points were combined.  This is returned directly
    outsize[0] = N;
    outsize[1] = 1;
    plhs[7] = mxCreateNumericArray(2, outsize, mxUINT32_CLASS, mxREAL);
    ConnectID = (unsigned int *)mxGetData(plhs[7]);

    for (ii = 0; ii<N; ii++)ConnectID[ii] = 0;
            
    //keep track of intial frame number
    FN_new = (unsigned int*)mxCalloc(N, sizeof(unsigned int));
    for (ii = 0; ii<N; ii++)FN_new[ii] = FrameNum[ii];
    
    cnt = 0;
    //Logic in loop requires monotonically increasing frame number
    
    int Clustercnt = MaxConID+1;
//     mexPrintf("Max Connect ID %d and CurrentID %d\n",MaxConID,Clustercnt);

    //Main loop for combination
    for (ii = 0; ii<N - 1; ii++){
        CurrentFrame = (int)FN_new[ii];
//         mexPrintf("current localization is %d\n",ii);
//         mexPrintf("current connetID is %d\n",Clustercnt);
//         mexPrintf("CurrentFrame is %d\n",CurrentFrame);
//       
//Clustercnt++;
        for (jj = ii + 1; jj<N; jj++){
//             mexPrintf("Considering  jj %d and X_jj %g Y_jj %g and Connect ID of jj %d\n",jj,X_new[jj],X_new[N+jj],ConnectID[jj]);
//            mexPrintf("Considering  jj %d and SE_X_jj %g and SE_Y_jj %g\n",jj,SE_new[jj],SE_new[N+jj]);

            //fast reject on frame gap
            if ((int)FN_new[jj]>(CurrentFrame + MaxFrameGap))
            {
//                 mexPrintf("Considering %d and %d Frame Gap \n",ii,jj);
                break;
            }
            if ((int)FN_new[jj] == CurrentFrame)
            {
//                 mexPrintf("Considering %d and %d same Frame \n",ii,jj);
                continue;
            }
            //fast rejection on distance.  This is linear distance in each dimension.
            cont = 0;
            for (nn = 0; nn<NDim; nn++)
                if (fabs(X_new[nn*N + jj] - X_new[nn*N + ii])>MaxDistance) cont++;
            if (cont)
            {
//                 mexPrintf("Considering %d and %d Max Dis \n",ii,jj);
                
                continue;
            }
            //fast rejection on euclidian distance
            tmp = 0;
            for (nn = 0; nn<NDim; nn++) tmp += powf(X_new[nn*N + jj] - X_new[nn*N + ii], 2);
            if (sqrt(tmp) > MaxDistance)
            {
//                 mexPrintf("Considering %d and %d Max r \n",ii,jj);
                
                continue;
            }
            //mexPrintf("Considering %d and %d\n",ii,jj);
            
            //make temp arrays for position and standard error
            for (nn = 0; nn < NDim; nn++) Xtmp1[nn] = X_new[nn*N + ii];
            for (nn = 0; nn < NDim; nn++) Xtmp2[nn] = X_new[nn*N + jj];
            for (nn = 0; nn < NDim; nn++) SEtmp1[nn] = SE_new[nn*N + ii];
            for (nn = 0; nn < NDim; nn++) SEtmp2[nn] = SE_new[nn*N + jj];
           
//             mexPrintf("Considering %d and %d , X of ii %g, X of jj %g Y of ii %g, Y of jj %g\n ",ii,jj,Xtmp1[0],Xtmp2[0],Xtmp1[1],Xtmp2[1]);
//             mexPrintf("Considering %d and %d , SE_X of ii %g, SE_X of jj %g SE_Y of ii %g, SE_Y of jj %g\n ",ii,jj,SEtmp1[0],SEtmp2[0],SEtmp1[1],SEtmp2[1]);

            //calculate p-value
            PValue = pvalue(Xtmp1, SEtmp1, Xtmp2, SEtmp2, NDim, Xtmpout, SEtmpout);
//             mexPrintf("Considering %d and %d , p_value %g, \n",ii,jj,PValue);
            //reject on PValue
            if (PValue<LOS) continue;
            
//             mexPrintf("Considering %d and %d , p_value %g, LOS %g\n",ii,jj,PValue, LOS);
            
            //We have passed all tests, so now combine points into later coordinate and null previous
            for (nn = 0; nn < NDim; nn++){
                X_new[nn*N + jj] = Xtmpout[nn];
                SE_new[nn*N + jj] = SEtmpout[nn];
                X_new[nn*N + ii] = 0;
                SE_new[nn*N + ii] = 0;
            }
//             mexPrintf("considering ii=%d and Photon==%g and bg=%g\n",ii,P_Add_new[0*N+ii],P_Add_new[1*N+ii]);
            //sum parameter_Add (in this case photons and background) and null previous
            for (nn = 0; nn < NDim_Add; nn++){
                P_Add_new[nn*N+jj] += P_Add_new[nn*N+ii];
                    P_Add_new[nn*N+ii]=0;
            }
            
            //weighted average parameter_Ave (in this case PSFSigma) and null previous
            if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
            {
                //make temp arrays for parameter_Ave and standard error
                for (nn = 0; nn < NDim_Ave; nn++) P_Avetmp1[nn] = P_Ave_new[nn*N + ii];
                for (nn = 0; nn < NDim_Ave; nn++) P_Avetmp2[nn] = P_Ave_new[nn*N + jj];
                for (nn = 0; nn < NDim_Ave; nn++) P_AveSEtmp1[nn] = P_Ave_SE_new[nn*N + ii];
                for (nn = 0; nn < NDim_Ave; nn++) P_AveSEtmp2[nn] = P_Ave_SE_new[nn*N + jj];
                
                //calculate weighted average
                out = Weighted_ave(P_Avetmp1, P_AveSEtmp1, P_Avetmp2, P_AveSEtmp2, NDim_Ave,P_Avetmpout, P_AveSEtmpout);
                
                
                //We have passed all tests, so now combine points into later coordinate and null previous
                for (nn = 0; nn < NDim_Ave; nn++){
                    P_Ave_new[nn*N + jj] = P_Avetmpout[nn];
                    P_Ave_SE_new[nn*N + jj] = P_AveSEtmpout[nn];
                    P_Ave_new[nn*N + ii] = 0;
                    P_Ave_SE_new[nn*N + ii] = 0;
                }
            }
            
            //Count number of frames combined and null previous
            N_combined[jj] += N_combined[ii];
            N_combined[ii] = 0;
            
            //Propogate cluster ID
            if ((ConnectID[ii] == 0) && (ConnectID[jj] == 0)) {//neither are in cluster already
                ConnectID[ii] = Clustercnt;   ConnectID[jj] = Clustercnt; Clustercnt++;
            }
            
            if ((ConnectID[ii] == 0) && (ConnectID[jj]>0))  //jj is already in cluster - this can happen sometimes
                ConnectID[ii] = ConnectID[jj];
            
            if ((ConnectID[jj] == 0) && (ConnectID[ii]>0))
                ConnectID[jj] = ConnectID[ii];
            
            if ((ConnectID[jj] > 0) && (ConnectID[ii] > 0))// this is rare, but want to combine all with same cluster ID
            {// we'll use cluster number from ii.
                for (kk = 0; kk < jj; kk++) if (ConnectID[kk] == ConnectID[jj]) ConnectID[kk] = ConnectID[ii];
                ConnectID[jj] = ConnectID[ii];
            }
            
            
//             mexPrintf("ii: %d jj: %d ConnectID: %d Clustercnt: %d NCombined: %d\n", ii,jj,ConnectID[ii], Clustercnt,N_combined[jj]);
            
            FN_new[ii] = 0;
            
            cnt++;
            break;// I guess this is wrong!
        }
    }

//     mexPrintf("Currecnt Clustercnt: %d \n", Clustercnt);
    //assign cluster id to singletons
    for (ii = 0; ii<N; ii++)
    {
//         mexPrintf("Number of Localization %d \n",N);
        if (ConnectID[ii] == 0){ ConnectID[ii] = Clustercnt; Clustercnt++; }
//            mexPrintf("Localization %d and Currecnt Connect ID: %d \n",ii, ConnectID[ii]);

    }
//             mexPrintf("last Clustercnt: %d \n", Clustercnt);
            
    //Create ouput arrays
    int NumOut = N - cnt;
    outsize[0] = NumOut;
    outsize[1] = NDim;
    plhs[0] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    Xout = (float *)mxGetData(plhs[0]);
    
    outsize[0] = NumOut;
    outsize[1] = NDim;
    plhs[1] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    SEout = (float *)mxGetData(plhs[1]);
    
    outsize[0] = NumOut;
    outsize[1] = 1;
    plhs[2] = mxCreateNumericArray(2, outsize, mxUINT32_CLASS, mxREAL);
    Nout = (unsigned int *)mxGetData(plhs[2]);
    
    outsize[0] = NumOut;
    outsize[1] = 1;
    plhs[3] = mxCreateNumericArray(2, outsize, mxUINT32_CLASS, mxREAL);
    Fout = (unsigned int *)mxGetData(plhs[3]);
    
            if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
            {
    outsize[0] = NumOut;
    outsize[1] =NDim_Ave;
    plhs[4] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    P_Aveout = (float *)mxGetData(plhs[4]);

    outsize[0] = NumOut;
    outsize[1] =NDim_Ave;
    plhs[5] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    P_Ave_SEout = (float *)mxGetData(plhs[5]);
            }
    
    if (mxGetM(prhs[4]) == 0 && mxGetN(prhs[4]) == 0)
            {
    outsize[0] = 1;
    outsize[1] =1;
    plhs[4] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    P_Aveout = (float *)mxGetData(plhs[4]);
    
    outsize[0] = 1;
    outsize[1] =1;
    plhs[5] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    P_Ave_SEout = (float *)mxGetData(plhs[5]);
            }
    
    outsize[0] = NumOut;
    outsize[1] = NDim_Add;
    plhs[6] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
    P_Addout = (float *)mxGetData(plhs[6]);
    
    outsize[0] = NumOut;
    outsize[1] = 1;
    plhs[8] = mxCreateNumericArray(2, outsize, mxUINT32_CLASS, mxREAL);
    CombinedID = (unsigned int *)mxGetData(plhs[8]);

  
//     
    //write into compressed coords array
    //mexPrintf("Writing into output array: %d.\n",N-cnt);
//     
    kk = 0;
    for (ii = 0; ii<N; ii++) if (N_combined[ii]){//this means this point was not erased
        for (nn = 0; nn<NDim; nn++)
        {
            Xout[nn*NumOut +kk] = X_new[nn*N + ii];
            SEout[nn*NumOut +kk] = SE_new[nn*N + ii];
        }
        
        if (mxGetM(prhs[4]) == 0 && mxGetN(prhs[4]) == 0)
       {
           P_Aveout= 0;
            P_Ave_SEout=0;
        }
        
        if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
        {
            for (nn = 0; nn<NDim_Ave; nn++)
            {
                P_Aveout[nn*NumOut +kk] = P_Ave_new[nn*N + ii];
                P_Ave_SEout[nn*NumOut +kk] = P_Ave_SE_new[nn*N + ii];
            }
        }
        for (nn = 0; nn<NDim_Add; nn++)
        {
            P_Addout[nn*NumOut +kk] = P_Add_new[nn*N + ii];
        }
        
        Nout[kk] = N_combined[ii];
        Fout[kk] = FN_new[ii];
        CombinedID[kk]=ConnectID[ii];
        kk++;
    }

    mxFree(X_new);
    mxFree(SE_new);
    if (mxGetM(prhs[4]) != 0 && mxGetN(prhs[4]) != 0)
    {
        mxFree(P_Ave_new);
        mxFree(P_Ave_SE_new);
    }
    mxFree(P_Add_new);
    mxFree(FN_new);
    mxFree(N_combined);

    return;
};

