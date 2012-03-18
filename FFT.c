#include <stdlib.h>
#include "FFT.h"

FFT* NewFFT( int numSamples )
{
	FFT* fft = (FFT*)malloc(sizeof(FFT));

    fft->direction = FORWARD;
    
	fft->config = kiss_fftr_alloc(numSamples, 0, NULL, NULL);
	fft->spectrum = (kiss_fft_cpx*) malloc(sizeof(kiss_fft_cpx) * numSamples);
    fft->realSignal = NULL;
	fft->numSamples = numSamples;
    fft->usescale = USESCALE_NO;
    fft->scale = 1;
	return fft;
}

FFT* NewIFFT( int numSamples )
{
	FFT* fft = (FFT*)malloc(sizeof(FFT));

    fft->direction = INVERSE;

	fft->config = kiss_fftr_alloc(numSamples, 1, NULL, NULL);
	fft->spectrum = NULL;
    fft->realSignal = (kiss_fft_scalar*)malloc(sizeof(kiss_fft_scalar)*numSamples);
	fft->numSamples = numSamples;
    fft->usescale = USESCALE_YES;
    fft->scale = 1.0f/numSamples;
	return fft;
}

void FreeFFT( FFT* fft )
{
    switch(fft->direction)
    {
        case FORWARD:
            free(fft->config);
	        free(fft->spectrum);
	        free(fft);
            break;
        case INVERSE:
            free(fft->config);
	        free(fft->realSignal);
	        free(fft);
            break;
    }
}

REAL cpx_abs(kiss_fft_cpx val)
{
    return (REAL)sqrtf(val.r*val.r + val.i*val.i);
}

void scale( FFT* fft )
{
    int i;
    if( fft->usescale == USESCALE_YES )
    {
        switch( fft->direction )
        {
        case FORWARD:
            for( i = 0; i<fft->numSamples; i++ )
            {
                fft->spectrum[i].i *= fft->scale;
                fft->spectrum[i].r *= fft->scale;
            }
            break;
        case INVERSE:
            for( i = 0; i<fft->numSamples; i++ )
            {
                fft->realSignal[i] *= fft->scale;
            }
            break;
        }
    }       
}

void FFT_Transform( FFT* fft )
{
    switch( fft->direction )
    {
        case FORWARD:
            kiss_fftr(fft->config,fft->realSignal,fft->spectrum);
            scale(fft);
            break;
        case INVERSE:
            kiss_fftri(fft->config,fft->spectrum,fft->realSignal);
            scale(fft);
            break;
    }
}