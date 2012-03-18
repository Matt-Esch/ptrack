#include "FFT.h"

#ifndef QIFFT_H_
#define QIFFT_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    Rectangular,
    Gauss,
    Hann
} T_WINDOW;

typedef enum {
    MaxBin,
    HPS_2,
    HPS_3
} T_BINSELECT;

typedef struct {
    int nfft;
    FFT* fft;
    T_WINDOW windowType;
    REAL* window;
    REAL* realSpectrum;
    REAL threshold;
    T_BINSELECT binSelection;   // Indicates which bins selection to use
    REAL* hps;
    int sampleRate;
} QIFFT;

QIFFT* NewQIFFT(int nfft, T_WINDOW windowType, T_BINSELECT binSelection, REAL threshold, int sampleRate );
void FreeQIFFT( QIFFT* qifft );
int getBin( QIFFT* qifft );
int getMaxBin( QIFFT* qifft );
int getHPSBin( QIFFT* qifft, int levels );
#ifdef __cplusplus
}
#endif
#endif



