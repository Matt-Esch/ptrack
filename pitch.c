#include <stdlib.h>
#include <math.h>
#include "pitch.h"

REAL zerocross( SIGNAL* input, ZeroCross* z )
{
    return zeroCross(input,z->windowSize,z->sampleRate);
}

//------------------------------------------------------------------------------
// Performs a correlation (autocorrelation, ASDF, AMDF or YIN) on a given input.
// A threshold for maxima/minima is specified to prevent location of suboptimal
// peaks. 
//------------------------------------------------------------------------------
REAL correlation( SIGNAL* input, Correlation* c )
{
    REAL timePeriod;

    switch(c->correlationType) 
    {
        case AUTO:
            autocorrelation(input,c);
            timePeriod = firstMaxima(c);
            break;
        case AMDF:
            amdf(input,c);
            timePeriod = firstMinima(c);
            break;
        case ASDF:
            asdf(input,c);
            timePeriod = firstMinima(c);
            break;
        case YIN:
            yin(input,c);
            timePeriod = firstMinima(c);
            break;
        case FAST_AUTO:
            fast_autocorrelation(input,c);
            timePeriod = firstMaxima(c);
            break;
        case FAST_ASDF:
            fast_asdf(input,c);
            timePeriod = firstMinima(c);
            break;
        case FAST_YIN:
            fast_yin(input,c);
            timePeriod = firstMinima(c);
            break;
        default:
            return 0;
    }

    if(timePeriod < 1)
    {
        // No frequency detected;
        return 0;
    }
    return (REAL)c->sampleRate/timePeriod;
}

//------------------------------------------------------------------------------
// Computes the zero crossing rate estimator of a given input buffer. This
// input can be pre-processed with a lowpass filter to increase robustness to
// noise.
//------------------------------------------------------------------------------
REAL interpolateFFT( SIGNAL* input, QIFFT* qifft )
{
    int i;
    int bin;
    REAL i_bin = 0;
    REAL i_val = 0;
    REAL sqrtnfft = sqrt((REAL)qifft->nfft);

    // Apply window function to input and copy to FFT input
    if( qifft->windowType == Rectangular )
    {
        for( i=0; i<qifft->nfft; i++ )
        {
            qifft->fft->realSignal[i] = (kiss_fft_scalar)input[i];
        }
    }
    else
    {
        for( i=0; i<qifft->nfft; i++ )
        {
            qifft->fft->realSignal[i] = input[i]*qifft->window[i];
        }
    }


    // Transform the signal
    FFT_Transform(qifft->fft);

    // Compute the real spectrum
    for( i=0; i<qifft->nfft/2 + 1; i++ )
    {
            qifft->realSpectrum[i] = cpx_abs(qifft->fft->spectrum[i])/sqrtnfft;
    }


    // Identify the most signficant bin above the threshold
    bin = getBin( qifft );

    if( bin == 0 )
    {
        // No frequency detected
        return 0;
    }
    else if( bin > qifft->nfft/2 )
    {
        i_bin = (REAL)bin;
    }
    else
    {
       interpolate( (REAL)bin,
                    logf(qifft->realSpectrum[bin-1]+1),
                    logf(qifft->realSpectrum[bin]+1),
                    logf(qifft->realSpectrum[bin+1]+1),
                    &i_bin,
                    &i_val );
    }

    if( i_bin >= 1 && i_bin <= qifft->nfft/2 )
    {
        if(i_bin>90)
        {
            return ((REAL)qifft->sampleRate/(REAL)qifft->nfft) * i_bin;
        }
        return ((REAL)qifft->sampleRate/(REAL)qifft->nfft) * i_bin;
    }

    return 0;
}


