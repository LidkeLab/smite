#include "device_launch_parameters.h"
#include "math.h"

__global__ void cuda_PSFSample3DBlob(const double *PX_in, const double *PY_in, const double *PZ_in, const double *PIntensity_in, const double *PBackground_in, const double XYSamPerPix, const double ZSamPerUnit, const int Row, const int Cent, const int *NN, const double *NotPSF, const int NFrames, float *M)
{

	// This function makes the model. It gets the position of the particle and then returns the image. This function is being called inside the loop, several times per iteration.
	// PX, PY, PZ, PIntensity, PBackground are respectively the X-position, Y-position, Z-position and intensity of the particle. PBackground is the background noise.
	// XYSamPerPix is the number of samples (samples of PFS) per pixel along the X and Y axes. ZSamPerUnit is the number of samples (PSF samples) along the Z-axis.
	// Row is the number of pixels along the X and Y axes, assuming that the image is a square image. Cent is the idex of the central column and central row of the NotPSF supposing that we have a square image.
	// NN is a 3-vector which has the 3 dimensions of the NotPSf array. The first element is the number of rows, the second element is the number of columns and the last one is the number of planes along the Z-axis.
	// NotPSF is a 3D array that contains the result of integrating the PSF samples. M is the output of the function, which is the model of the image. 

	double PX, PY, PZ, PIntensity, PBackground;
	double IntX, IntY, DecX, DecY, Axx, Ayy, XYL, YLowInd, XLowInd, PIndX1Y1, PIndX1Y2, PIndX2Y1, PIndX2Y2, CentZ, CeilZ;
	double IntZ, DecZ, Az, DecDecZ, DecDecX, DecDecY, PInterpolX, PInterpolY, PInterpolZ;
	int IndX1Y1, IndX1Y2, IndX2Y1, IndX2Y2, IndZ1CentZ, IndZ2CentZ, IndZ1, IndZ2;
	int Ind1, Ind2, Ind3, Ind4, Ind5, Ind6, Ind7, Ind8;
	int PixNum = Row*Row;

	// The three components of NN are given to the three following parameters.
	int NS1 = (int)NN[0];
	int NS2 = (int)NN[1];
	int NS3 = (int)NN[2];

	// Nn gives the number of pixels inside one of the planes along the Z-axis in the NotPSF variable.
	int Nn = NS1*NS2;

	int ThreadNumPerBlock = blockDim.x;
	int ThreadNumInBlock = threadIdx.x;
	int BlockNum = blockIdx.x;
	int Start = (BlockNum*ThreadNumPerBlock + ThreadNumInBlock)*PixNum;
	int PStart = BlockNum*ThreadNumPerBlock + ThreadNumInBlock;

    if (Start>((NFrames*PixNum)-1)) return;
        
	PX = PX_in[PStart]+1.0f/(2.0f*XYSamPerPix);
	PY = PY_in[PStart]+1.0f/(2.0f*XYSamPerPix);
	PZ = PZ_in[PStart];
	PIntensity = PIntensity_in[PStart];
	PBackground = PBackground_in[PStart];

	
    IntZ = floor(PZ);
	DecZ = PZ - IntZ;
	Az = floor(ZSamPerUnit*DecZ);
	DecDecZ = DecZ - Az / ZSamPerUnit;
	// I used the following three lines rather using round() function.
	CentZ = NS3 / 2.0f;
	CeilZ = ceil(CentZ);
	CentZ = CeilZ;
	// IndZ1 and IndZ2 are the indices of the first components of the PsfZ1 and PsfZ2.
	IndZ1 = (int)round(ZSamPerUnit*(PZ - DecDecZ));
	IndZ2 = IndZ1 + 1;
	IndZ1CentZ = IndZ1 + (int)CentZ;
	IndZ2CentZ = IndZ2 + (int)CentZ;

	IntX = floor(PX);
	IntY = floor(PY);
	DecX = PX - IntX;
	DecY = PY - IntY;
	Axx = floor(XYSamPerPix*DecX);
	Ayy = floor(XYSamPerPix*DecY);
	DecDecX = DecX - Axx / XYSamPerPix;
	DecDecY = DecY - Ayy / XYSamPerPix;

	//XLowInd and YLowInd are the indices of the first components the we pick from PsfZ along the X and Y axes.
	XYL = ceil(XYSamPerPix / 2);
	YLowInd = Cent + XYL - XYSamPerPix*(PY - DecDecY);
	XLowInd = Cent + XYL - XYSamPerPix*(PX - DecDecX);

	// The following for parameters are the linear indices of the first componenets for the PsfX1Y1, PsfX1Y2, PsfX2Y1, PsfX2Y2.
	PIndX1Y1 = (XLowInd)*(double)NS1 + YLowInd;
	PIndX1Y2 = (XLowInd)*(double)NS1 + YLowInd - 1.0f;
	PIndX2Y1 = (XLowInd - 1.0f)*(double)NS1 + YLowInd;
	PIndX2Y2 = (XLowInd - 1.0f)*(double)NS1 + YLowInd - 1.0f;

    PInterpolX = XYSamPerPix*(1.0f / XYSamPerPix - DecDecX);
    PInterpolY = XYSamPerPix*(1.0f / XYSamPerPix - DecDecY);
    PInterpolZ = ZSamPerUnit*(1.0f / ZSamPerUnit - DecDecZ);

for (int nn = 0; nn < PixNum; nn++)
	{
		// the interpolation along X and Y axes and gives the final model.
		double w = nn / Row;
		int xx = (int)floor(w);
		int yy = nn - Row*xx;

		IndX1Y1 = (int)PIndX1Y1 + (int)round(XYSamPerPix*(double)(NS1*xx)) + (int)round(XYSamPerPix*(double)yy);
		IndX1Y2 = (int)PIndX1Y2 + (int)round(XYSamPerPix*(double)(NS1*xx)) + (int)round(XYSamPerPix*(double)yy);
		IndX2Y1 = (int)PIndX2Y1 + (int)round(XYSamPerPix*(double)(NS1*xx)) + (int)round(XYSamPerPix*(double)yy);
		IndX2Y2 = (int)PIndX2Y2 + (int)round(XYSamPerPix*(double)(NS1*xx)) + (int)round(XYSamPerPix*(double)yy);
		Ind1 = Nn*(IndZ1CentZ-1)+IndX1Y1;
        Ind2 = Nn*(IndZ2CentZ-1)+IndX1Y1;
        Ind3 = Nn*(IndZ1CentZ-1)+IndX2Y1;
        Ind4 = Nn*(IndZ2CentZ-1)+IndX2Y1;
        Ind5 = Nn*(IndZ1CentZ-1)+IndX1Y2;
        Ind6 = Nn*(IndZ2CentZ-1)+IndX1Y2;
        Ind7 = Nn*(IndZ1CentZ-1)+IndX2Y2;
        Ind8 = Nn*(IndZ2CentZ-1)+IndX2Y2;

		M[Start + nn] = PIntensity*(PInterpolY * (PInterpolX*(PInterpolZ*NotPSF[Ind1] + (1.0f-PInterpolZ)*NotPSF[Ind2])
			+ (1.0f-PInterpolX)*(PInterpolZ*NotPSF[Ind3] + (1.0f-PInterpolZ)*NotPSF[Ind4])) 
            + (1.0f-PInterpolY)*(PInterpolX*(PInterpolZ*NotPSF[Ind5] + (1.0f-PInterpolZ)*NotPSF[Ind6])
			+ (1.0f-PInterpolX)*(PInterpolZ*NotPSF[Ind7] + (1.0f-PInterpolZ)*NotPSF[Ind8]))) + PBackground;

	}
}
