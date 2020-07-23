#define pi 3.141592
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

//kernel_guassiansampleblobs(sz,Nframes,d_xarray,d_yarray,d_Narray,d_xsigma,d_ysigma,d_covariance,d_im);  

__global__ void kernel_guassiansampleblobs( const int sz, const int Nframes, const float *d_xarray, const float *d_yarray, const float *d_Narray, const float *d_Barray, const float *d_xsigma, const float *d_ysigma, const float *d_covariance,      float *d_im  ) {
	int tx = threadIdx.x; 
	int bx = blockIdx.x;
    int BlockSize = blockDim.x;
	float x,y,xsigma,ysigma,covariance,N,Bg,A;
	int ii,jj,pixelx,pixely,Idx,Idxii;
	
    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nframes) return;
    
	//import datas from device to shared memory
    x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
    Bg=d_Barray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];
	covariance=d_covariance[bx*BlockSize+tx];
    
    // precalculate for speed
    A = N/(2*pi*xsigma*ysigma*sqrt(1-pow(covariance,2)));
    Idx = bx*BlockSize*sz*sz+tx*sz*sz;
	
    for (ii=0;ii<sz;ii++) {
        pixelx=ii; // sample at pixel center
        Idxii = ii*sz;
        for(jj=0;jj<sz;jj++) {
            pixely=jj; // sample at pixel center
            // generate model for pixel ii jj
            d_im[Idx+Idxii+jj] = Bg + (A* exp( -1/(2*(1-pow(covariance,2))) * ( pow(x-pixelx-0.5,2)/pow(xsigma,2) + pow(y-pixely-0.5,2)/pow(ysigma,2) - 2*covariance*(x-pixelx-0.5)*(y-pixely-0.5)/(xsigma*ysigma) ) ) );
        }
    }
	return;
}


//kernel_guassianintegrateblobs(sz,Nframes,d_xarray,d_yarray,d_Narray,d_xsigma,d_ysigma,d_im);

__global__ void kernel_guassianintegrateblobs( const int sz, const int Nframes, const float *d_xarray, const float *d_yarray, const float *d_Narray, const float *d_Barray, const float *d_xsigma, const float *d_ysigma,      float *d_im  ) {

	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
    int BlockSize = blockDim.x;
	float x,y,xsigma,ysigma,N,Bg,A,B,SX,SY;
	int ii,jj,pixelx,pixely,Idx,Idxii;

    //Prevent read/write past end of array
    if ((bx*BlockSize+tx)>=Nframes) return;

	//import datas from device to shared memory
	x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
    Bg=d_Barray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];

    // precalculate these for speed
    A = N/4;
    Idx = bx*BlockSize*sz*sz+tx*sz*sz;
    SX = 1/sqrt(2*pow(xsigma,2));
    SY = 1/sqrt(2*pow(ysigma,2));
    
	for (ii=0;ii<sz;ii++) {
        pixelx=ii;
        Idxii = ii*sz;
        B = erf((x-pixelx)*SX)-erf((x-pixelx-1)*SX);
        for(jj=0;jj<sz;jj++) {
            pixely=jj;
            // generate model for pixel ii jj
            d_im[Idx+Idxii+jj] = Bg + (A*B*(erf((y-pixely)*SY)-erf((y-pixely-1)*SY)));  
        }
    }
	return;
}



