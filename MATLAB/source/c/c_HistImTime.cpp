/* c_HistImTime

use: driftIm = c_HistImTime(ysize,xsize,y,x,t);
 * INPUT
 *   ysize - y size of output image 
 *   xsize - x size of output image
 *   y     - y coordinates of points (single floats)  
 *   x     - x coordinates of points (single floats)  
 *   t     - time stamp of points (single floats)
 * OUTPUT
 *   driftIm = image containing points with gray scale value based on timestamp
 
 */ 


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
	float *X, *Y, *T;
	float *out;
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
		"Must input sizeX,sizeY,X,Y,UseHistEqual!\n$Rev: 56 $\n$Date: 2011-09-28 13:29:18 -0600 (Wed, 28 Sep 2011) $\n$Author: jbyars $\n$HeadURL: https://abbe.phys.unm.edu/svn/cCode/cHistRecon.c $");
	if (mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cHistRecon:WrongPrecision", "Y array inputs must be comprised of single floats!\n");
	if (mxGetClassID(prhs[3]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cHistRecon:WrongPrecision", "X array inputs must be comprised of single floats!\n");
	if (mxGetClassID(prhs[4]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cHistRecon:WrongPrecision", "T array inputs must be comprised of single floats!\n");


	sizeY = (int)mxGetScalar(prhs[0]);
	sizeX = (int)mxGetScalar(prhs[1]);


	if (sizeX <= 0)
		mexErrMsgIdAndTxt("cHistRecon:InputOutofRange", "sizeX must be larger than 0\n");
	if (sizeY <= 0)
		mexErrMsgIdAndTxt("cHistRecon:InputOutofRange", "sizeY must be larger than 0\n");


	//get variables
	Y = (float *)mxGetData(prhs[2]);
	X = (float *)mxGetData(prhs[3]);
	T = (float *)mxGetData(prhs[4]);

	datasize = mxGetDimensions(prhs[2]);
	//create output

	outsize[0] = sizeY;
	outsize[1] = sizeX;
	plhs[0] = mxCreateNumericArray(2, outsize, mxSINGLE_CLASS, mxREAL);
	out = (float *)mxGetData(plhs[0]);



	for (i = 0; i<datasize[0]; i++){
		Xtmp = max(0, min(sizeX - 1, (int)floor(X[i])));
		Ytmp = max(0, min(sizeY - 1, (int)floor(Y[i])));
		out[Xtmp*sizeY + Ytmp] = T[i];
	}


}














