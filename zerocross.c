#include <stdlib.h>
#include "pitchutils.h"
#include "zerocross.h"
//------------------------------------------------------------------------------
// Computes the zero crossing rate estimator of a given input buffer. This
// input can be pre-processed with a lowpass filter to increase robustness to
// noise.
//------------------------------------------------------------------------------
REAL zeroCross(SIGNAL* input, int length, int sampleRate)
{
    int i;
    
    // interpolate crossings
    int firstCross = -1;
    int lastCross = -1;
    float iLeft;
    float iRight;

    int count = 0;
    REAL sgn0 = sgn(input[0]);
    REAL sgn1 = 0;

    for(i=1; i<length; i++)
    {
        sgn1 = sgn(input[i]);
        if(absdif(sgn0,sgn1)>1)
        {
            if(firstCross<0)
            {
                firstCross = i;
                lastCross = firstCross;
                iLeft = (i-1) + (input[i-1] / (REAL)(input[i-1]-input[i]));
            }
            else
            {
                lastCross = i;
            }
            count++;
        }
        sgn0 = sgn1;
    }


    if( lastCross - firstCross <= 0)
    {
        return 0;
    }
    else
    {
        iRight = (lastCross-1) + (input[lastCross-1] / (REAL)(input[lastCross-1]-input[lastCross]));
        return sampleRate*(count-1)/((REAL)2.0*(iRight-iLeft));
    }
}

ZeroCross* NewZeroCross(int windowSize, int sampleRate)
{
    ZeroCross* z = (ZeroCross*)malloc(sizeof(ZeroCross));
    z->windowSize = windowSize;
    z->sampleRate = sampleRate;
    return z;
}

void FreeZeroCross(ZeroCross* zeroCross)
{
    free(zeroCross);
}
