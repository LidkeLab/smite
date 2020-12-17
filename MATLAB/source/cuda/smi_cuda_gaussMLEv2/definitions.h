/*!
 * \file definitions.h
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief The constants used throughout the project. 
 */

#define BSZ 128			//!< max number of threads per block 
#define MEM 3872        //3872		11616//!< shared 
#define IMSZ 11			//!< not used
#define IMSZBIG 21		//!< maximum fitting window size
#define NK 128			//!< number of blocks to run in each kernel
#define pi 3.141592f	//!< ensure a consistent value for pi
#define NV_P 4			//!< number of fitting parameters for MLEfit (x,y,bg,I)
#define NV_PS 5			//!< number of fitting parameters for MLEFit_sigma (x,y,bg,I,Sigma)
#define NV_PZ 5			//!< not used (x,y,bg,I,z)
#define NV_PS2 6		//!< number of fitting parameters for MLEFit_sigmaxy (x,y,bg,I,Sx,Sy)
