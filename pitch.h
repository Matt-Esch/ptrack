#include "pitchutils.h"
#include "zerocross.h"
#include "correlation.h"
#include "FFT.h"
#include "QIFFT.h"

#ifndef PITCH_H_
#define PITCH_H_
#ifdef __cplusplus
extern "C" {
#endif

REAL zeroCross(REAL* input, int length, int sampleRate);
REAL correlation( REAL* input, Correlation* c );
REAL interpolateFFT( REAL* input, QIFFT* c );

#ifdef __cplusplus
}
#endif
#endif

