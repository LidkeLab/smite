// file c_GenColorChannels.c

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include "mex.h"
#include "matrix.h"

#define pi 3.141592f

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif



void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const	mxArray	*prhs[])
{
	float *in;
	double *cm;
	uint8_t *r, *g, *b;
	int i, j;
	int Nim;
	const mwSize *cmsize;
	int cmlength;
	float maxvalue, minvalue;
	float mintmp, maxtmp;
	float tmp, diff;

	//check Inputs
	if (nrhs<2)
		mexErrMsgIdAndTxt("cGenColorChannels:WrongNumberOfInputs",
		"Must input image,colormap, optional: minvalue,maxvalue");
	if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("cGenColorChannels:WrongPrecision", "input image must be comprised of single floats!\n");
	if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
		mexErrMsgIdAndTxt("cGenColorChannels:WrongPrecision", "colormap must be comprised of double floats!\n");

	//get variables
	in = (float *)mxGetData(prhs[0]);
	cm = (double *)mxGetData(prhs[1]);
	if (nrhs>2) minvalue = (float)mxGetScalar(prhs[2]);
	if (nrhs>3) maxvalue = (float)mxGetScalar(prhs[3]);

	Nim = (int)mxGetNumberOfElements(prhs[0]);

	//get and check color map size
	cmsize = mxGetDimensions(prhs[1]);
	if (cmsize[1] != 3) mexErrMsgIdAndTxt("cGenColorChannels:Wrong Input Size", "color map must be N x 3\n");
	cmlength = cmsize[0];

	//create output
	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxUINT8_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxUINT8_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]), mxGetDimensions(prhs[0]), mxUINT8_CLASS, mxREAL);

	r = (uint8_t *)mxGetData(plhs[0]);
	g = (uint8_t *)mxGetData(plhs[1]);
	b = (uint8_t *)mxGetData(plhs[2]);

	//find max,min value if nessesary
	if (nrhs<3)//find min value from data
	{
		minvalue = in[0];
		for (i = 1; i<Nim; i++)minvalue = min(minvalue, in[i]);
	}
	if (nrhs<4)//find min value from data
	{
		maxvalue = in[0];
		for (i = 1; i<Nim; i++)maxvalue = max(maxvalue, in[i]);
	}

	//calculate r,g,b values
	diff = 1 / (maxvalue - minvalue);
	for (i = 0; i<Nim; i++){
		tmp = (in[i] - minvalue)*diff;
		tmp = min(max(tmp, 0), 1);
		j = (int)floor((cmlength - 1)*tmp);
		r[i] = (uint8_t)(cm[j] * 255.0);
		g[i] = (uint8_t)(cm[j + cmlength] * 255.0);
		b[i] = (uint8_t)(cm[j + cmlength * 2] * 255.0);
	}
}














