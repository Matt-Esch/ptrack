#ifndef PITCHUTILS_H_
#define PITCHUTILS_H_
#ifdef __cplusplus
extern "C" {
#endif

#define PI (3.141592653589793)
    
typedef float REAL;

#ifdef FIXED_POINT
typedef int32 SIGNAL;
#define MAX_SIGNAL (32768)
#else
typedef float SIGNAL;
#define MAX_SIGNAL (1.0f)
#endif

REAL interpolate(REAL xc, REAL yl, REAL yc, REAL yu, REAL* xe, REAL* ye);
REAL getabs( REAL val );
void sinewave(REAL* buffer, int bufferSize, int sampleRate, int frequency);
REAL sgn( REAL sample );
REAL absdif(REAL x, REAL y);
void zeroPad(REAL* buff, long buffLen, int startPos);
REAL Square(REAL x);
REAL* GetGaussWindow(int bufferLen);
REAL* GetHannWindow(int bufferLen);

#ifdef __cplusplus
}
#endif
#endif
