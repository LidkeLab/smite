/*!
 * \file definitions.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  If you are defining any constants they better be in this file and only this file.
 */

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#define pi 3.141592f	//!< ensure a consistent value for pi

//pi is defined in <mex.h>
//#ifndef pi
//#define pi 3.141592f;

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif


#endif
