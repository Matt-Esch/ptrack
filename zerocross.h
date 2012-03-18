#ifndef ZEROCROSS_H_
#define ZEROCROSS_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    int windowSize;
    int sampleRate;
} ZeroCross;

ZeroCross* NewZeroCross(int windowSize, int sampleRate);
void FreeZeroCross(ZeroCross* zeroCross);
REAL zeroCross(SIGNAL* input, int length, int sampleRate);

#ifdef __cplusplus
}
#endif
#endif