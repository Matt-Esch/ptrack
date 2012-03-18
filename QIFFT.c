#include <stdlib.h>
#include "QIFFT.h"

QIFFT* NewQIFFT(int nfft, T_WINDOW windowType, T_BINSELECT binSelection, REAL threshold, int sampleRate)
{
    QIFFT* qifft = (QIFFT*)malloc(sizeof(QIFFT));
    qifft->fft = NewFFT(nfft);
    qifft->nfft = nfft;
    qifft->fft->realSignal = (kiss_fft_scalar*)malloc(nfft*sizeof(kiss_fft_scalar));

    qifft->realSpectrum = (REAL*)malloc(((nfft/2)+1)*sizeof(REAL));
    qifft->binSelection = binSelection;
    qifft->windowType = windowType;
    
    switch(windowType)
    {
    case Gauss:
        qifft->window = GetGaussWindow(nfft);
        break;
    case Hann:
        qifft->window = GetHannWindow(nfft);
        break;
    default:
        qifft->window = NULL;
        break;
    }

    switch(binSelection)
    {
    case MaxBin:
        qifft->hps = NULL;
        break;
    case HPS_2:
        qifft->hps = (REAL*)malloc(((nfft/2)-1)*sizeof(REAL));
        break;
    case HPS_3:
        qifft->hps = (REAL*)malloc(((nfft/3)-1)*sizeof(REAL));
    }

    qifft->threshold = threshold;
    qifft->sampleRate = sampleRate;

    return qifft;
}

void FreeQIFFT( QIFFT* qifft )
{
    if( qifft->realSpectrum != NULL )
    {
        free(qifft->realSpectrum);
    }
    if( qifft->hps != NULL )
    {
        free(qifft->hps);
    }
    if( qifft->window != NULL )
    {
        free(qifft->window);
    }
    if( qifft->fft != NULL )
    {
        FreeFFT(qifft->fft);
    }
    free(qifft);
}

int getBin( QIFFT* qifft )
{
    switch( qifft->binSelection )
    {
    case MaxBin:
        return getMaxBin( qifft );
        break;
    case HPS_2:
        return getHPSBin(qifft, 2);
        break;
    case HPS_3:
        return getHPSBin(qifft, 3);
        break;
    default:
        return 0;
        break;
    }
}

int getMaxBin( QIFFT* qifft )
{
    int i;
    int max = 1;
    
    for( i = 2; i<=qifft->nfft/2; i++ )
    {
        if( qifft->realSpectrum[i] > qifft->realSpectrum[max] )
        {
            max = i;
        }
    }

    if( qifft->realSpectrum[max] >= qifft->threshold )
    {
        return max;
    }

    return 0;
}

int getHPSBin( QIFFT* qifft, int levels )
{
    int i, j;
    int minBin = 1;
    int topBin = qifft->nfft/(2*levels);
    int maxBin = minBin;

    for(i=minBin; i<topBin; i++)
    {
        for(j=1; j<=levels; j++)
        {
            if(j==1)
            {
                qifft->hps[i-minBin] = qifft->realSpectrum[i];
            }

            qifft->hps[i-minBin] *= qifft->realSpectrum[i*j];

        }
        if(qifft->hps[i-minBin] > qifft->hps[maxBin-minBin])
        {
            maxBin = i;
        }
    }


    if( qifft->realSpectrum[maxBin] >= qifft->threshold )
    {
        return maxBin;
    }
    return 0;
}
