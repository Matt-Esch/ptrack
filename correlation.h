#include "pitchutils.h"
#include "FFT.h"

#ifndef CORRELATION_H_
#define CORRELATION_H_
#ifdef __cplusplus
extern "C" {
#endif

/* Standard Correlation */
typedef enum correlation_type{
    AUTO,
    AMDF,
    ASDF,
    YIN,
    FAST_AUTO,
    FAST_ASDF,
    FAST_YIN
} CORR;

typedef struct
{
    CORR correlationType;
    REAL* correlation;
    int length;
    int searchStart;
    int searchEnd;
    int sampleRate;
    SIGNAL threshold;
    FFT* fft_x1;
    FFT* fft_x2;
    FFT* ifft;
    SIGNAL* x2;
    REAL* buffer;   // General purpose buffer for efficiency;
} Correlation;


Correlation* NewCorrelation( CORR correlationType, 
                             int length, 
                             int searchStart, 
                             int searchEnd, 
                             int sampleRate, 
                             REAL threshold );
Correlation* NewFastCorrelation( CORR correlationType, 
                                 int length, 
                                 int searchStart, 
                                 int searchEnd, 
                                 int sampleRate, 
                                 REAL threshold );
void FreeCorrelation( Correlation* c );



REAL firstMinima( Correlation* c );
REAL firstMaxima( Correlation* c );

void autocorrelation(SIGNAL* input, Correlation* c);
void amdf(SIGNAL* input, Correlation* c);
void asdf(SIGNAL* input, Correlation* c);
void yin(SIGNAL* input, Correlation* c);
void fast_autocorrelation(SIGNAL* input, Correlation* c);
void fast_asdf(SIGNAL*, Correlation* c);
void fast_yin(SIGNAL*, Correlation* c);

#ifdef __cplusplus
}
#endif
#endif