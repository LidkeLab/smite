/* c_HistRecon.cpp
use: HistIm = c_HistRecon(Ysize,Xsize,Y,X,histEqual);
 * INPUT
 *   Ysize - y size of output image 
 *   Xsize - x size of output image
 *   Y     - y coordinates of points (single floats)  
 *   X     - x coordinates of points (single floats)  
 *   HistEqual - flag for performing histogram equalization
 * OUTPUT
 *   HistIm = histogram image of points
*/

//

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

// Thread block size
#define CBLK 4000
#define pi 3.141592f

#define max(a,b) ( (a) >= (b) ? (a) : (b) )  
#define min(a,b) ( (a) < (b) ? (a) : (b) )  


void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const	mxArray	*prhs[])
{
	float *X, *Y;
	uint16_T *out;
	int i;
	int N, sizeX, sizeY;
	const mwSize *datasize;
	int Xtmp, Ytmp, useHE;
	mwSize outsize[2];
	int *Hist;
	float *CDF;
	int sum;
	int cdf_min;
	int Npixels;
	float *out2;
	const int NBIN = 1023;

	//check Inputs
	if (nrhs != 5)
		mexErrMsgIdAndTxt("MATLAB:WrongNumberOfInputs",
		"Must input sizeX,sizeY,X,Y,UseHistEqual!\n$Rev: 56 $\n$Date: 2011-09-28 12:29:18 -0700 (Wed, 28 Sep 2011) $\n$Author: jbyars $\n$HeadURL: https://rayleigh.phys.unm.edu/svn/cCode/cHistRecon.c $");
	if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cHistRecon:WrongPrecision", "X array inputs must be comprised of single floats!\n");
	if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cHistRecon:WrongPrecision", "Y array inputs must be comprised of single floats!\n");


	sizeY = (int)mxGetScalar(prhs[0]);
	sizeX = (int)mxGetScalar(prhs[1]);
	useHE = (int)mxGetScalar(prhs[4]);

	if (sizeX <= 0)
		mexErrMsgIdAndTxt("cHistRecon:InputOutofRange", "sizeX must be larger than 0\n");
	if (sizeY <= 0)
		mexErrMsgIdAndTxt("cHistRecon:InputOutofRange", "sizeY must be larger than 0\n");


	//get variables
	Y = (float *)mxGetData(prhs[2]);
	X = (float *)mxGetData(prhs[3]);

	datasize = mxGetDimensions(prhs[2]);
	//create output

	outsize[0] = sizeY;
	outsize[1] = sizeX;
	plhs[0] = mxCreateNumericArray(2, outsize, mxUINT16_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);

	out = (uint16_T *)mxGetData(plhs[0]);
	out2 = (float *)mxGetData(plhs[1]);

	sum = 0;
	for (i = 0; i<datasize[0]; i++){
		Xtmp = max(0, min(sizeX - 1, (int)floor(X[i])));
		Ytmp = max(0, min(sizeY - 1, (int)floor(Y[i])));
		out[Xtmp*sizeY + Ytmp] += 1;
	}

	if (useHE) //Histogram equalization
	{
		Hist = (int*)mxCalloc(NBIN, sizeof(int));
		CDF = (float*)mxCalloc(NBIN, sizeof(float));
		for (i = 0; i<NBIN; i++) Hist[i] = 0;

		Npixels = sizeX*sizeY;
		for (i = 0; i<Npixels; i++) if (out[i]) Hist[out[i]] += 1; //calc histogram

		i = 1; //ignoring zero pixels
		CDF[0] = 0;
		sum = 0;
		while (sum<datasize[0]) //calc CDF
		{
			CDF[i] = CDF[i - 1] + (float)Hist[i];
			sum += Hist[i] * i;
			i++;
		}
		N = (int)CDF[i - 1];
		//printf("%f %d %d\n",CDF[i-1],i-1,N);
		//printf("%f %f %f %f %f %f\n",CDF[0],CDF[1],CDF[2],CDF[3],CDF[4],CDF[5]);
		//convert pixels
		for (i = 0; i<Npixels; i++) if (out[i]) {
			out2[i] = NBIN*(CDF[out[i]] - 1) / (N - 1);
			out[i] = (uint16_T)NBIN*(CDF[out[i]] - 1) / (N - 1);
		}
		mxFree(Hist);
		mxFree(CDF);
	}
}














