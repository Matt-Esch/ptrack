#include "pitchutils.h"
#include <stdlib.h>
#include <math.h>

//------------------------------------------------------------------------------
// Performs quadratic interpolation on 3 points spaced 1 apart. The return value
// is the second derivative, indicating whether or not the point is a local
// maxima (<0) or local minima (>0). Unreliable interpolations return 0.
// Adpated from: http://www.ebyte.it/library/codesnippets/P3Interpolation.html
//------------------------------------------------------------------------------
REAL interpolate(REAL xc, REAL yl, REAL yc, REAL yu, REAL* xe, REAL* ye)
{
    REAL d1,d2;
    d2 = yu-yc+yl-yc;
    d1 = (REAL)0.5*(yu-yl);
    if (d2) {
        *xe = xc - d1/d2;
        *ye = yc + (REAL)0.5*d1*(*xe-xc);
    } else {
        *xe = xc;
        *ye = yc;
        return 0.0;
    }
    if (fabs(xc-*xe)>1)
    {
        return 0.0;
    }
    else 
    {
        return d2;
    }
}

//------------------------------------------------------------------------------
// Returns the absolute (positive) value of a given sample.
//------------------------------------------------------------------------------
REAL getabs( REAL val )
{
    if( val < 0 )
    {
        return -1*val;
    }
    else
    {
        return val;
    }
}

//------------------------------------------------------------------------------
// Generates a sine wave with a specified frequency in a given buffer length.
//------------------------------------------------------------------------------
void sinewave(SIGNAL* buffer,int bufferSize, int sampleRate, int frequency) {
    int i;
    REAL ratio = (REAL)frequency/(REAL)sampleRate;
    for (i = 0; i < bufferSize; i++) {
        buffer[i] = (SIGNAL)MAX_SIGNAL*(SIGNAL)sin(2 * PI * i * ratio);
    }
}

//------------------------------------------------------------------------------
// Indicates the sign of the sample. -1 = -ve, +1 = +ve
//------------------------------------------------------------------------------
REAL sgn( REAL sample )
{
    if( sample < 0 )
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

//------------------------------------------------------------------------------
// Computes the absolute difference between two given samples. Used to determine
// a change in sign in the zero crossing rate estimator.
//------------------------------------------------------------------------------
REAL absdif(REAL x, REAL y)
{
    REAL r = x-y;
	if(r<0)
	{
		return -1*r;
	}
	else
	{
		return r;
	}
}

//------------------------------------------------------------------------------
// Sets the trailing values of a buffer to 0
//------------------------------------------------------------------------------
void zeroPad(SIGNAL* buff, long buffLen, int startPos)
{
    long i;
    for(i=startPos; i<buffLen; i++)
    {
        buff[i]=0;
    }
}

REAL Square(REAL x)
{
    return x * x;
}

REAL* GetGaussWindow(int bufferLen)
{
    REAL sigma = 0.4f;

    REAL* wbuffer = (REAL*)malloc(bufferLen*sizeof(REAL));

	int k1 = (bufferLen - 1)/2;
	REAL k2 = sigma * (REAL)(bufferLen-1)/2;
	int i;
	for( i = 0; i < bufferLen; i++ )
	{
		wbuffer[i] = (REAL)exp(-0.5*pow(((i-k1)/(k2)),2));
	}

    return wbuffer;
}

REAL* GetHannWindow(int bufferLen)
{
    REAL* wbuffer = (REAL*)malloc(bufferLen*sizeof(REAL));
    int i;
	for( i = 0; i < bufferLen; i++ )
	{
		wbuffer[i] = (REAL)0.5*(1-cos((REAL)(2*PI*i)/(REAL)(bufferLen-1)));
	}

    return wbuffer;
}