#include <kiss_fftr.h>
#include <kiss_fft.h>
#include "pitchutils.h"

#ifndef FFT_H_
#define FFT_H_
#ifdef __cplusplus
extern "C" {
#endif

#define USESCALE_YES 'y'
#define USESCALE_NO  'n'

typedef enum fft_direction{
    FORWARD,
    INVERSE
} FFTDIRECT;

typedef struct{
    FFTDIRECT direction;
	kiss_fftr_cfg config;
    kiss_fft_scalar* realSignal;
	kiss_fft_cpx* spectrum;
	int numSamples;
    char usescale;
    float scale;
} FFT;

FFT* NewFFT( int numSamples );
FFT* NewIFFT( int numSamples );
void FreeFFT( FFT* fft );

REAL cpx_abs(kiss_fft_cpx val);

void FFT_Transform( FFT* fft );

#ifdef __cplusplus
}
#endif
#endif