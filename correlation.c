#include "correlation.h"
#include "pitchutils.h"


Correlation* NewCorrelation( CORR correlationType, 
                             int length, 
                             int searchStart, 
                             int searchEnd, 
                             int sampleRate, 
                             REAL threshold )
{
    Correlation* c = (Correlation*)malloc(sizeof(Correlation));
    c->correlationType = correlationType;
    c->length = length;
    c->searchStart = searchStart;
    c->searchEnd = searchEnd;
    c->sampleRate = sampleRate;
    c->threshold = (SIGNAL)threshold*MAX_SIGNAL;
    c->correlation = (REAL*)malloc(sizeof(REAL)*length);
    c->x2 = NULL;
    c->fft_x1 = NULL;
    c->fft_x2 = NULL;
    c->ifft = NULL;
    c->buffer = NULL;
    return c;
}

Correlation* NewFastCorrelation( CORR correlationType, 
                                 int length, 
                                 int searchStart, 
                                 int searchEnd, 
                                 int sampleRate, 
                                 REAL threshold )
{
    Correlation* c = (Correlation*)malloc(sizeof(Correlation));
    c->correlationType = correlationType;
    c->length = length;
    c->searchStart = searchStart;
    c->searchEnd = searchEnd;
    c->sampleRate = sampleRate;
    c->threshold = threshold;    
    
    /* Initialize new FFTs and IFFT */
    c->fft_x1 = NewFFT( 2*length );
    c->fft_x2 = NewFFT( 2*length );
    c->ifft = NewIFFT( 2*length );
    c->ifft->spectrum = c->fft_x1->spectrum;

    /* Calloc important for 0 padded signal */
    c->x2 = (REAL*)calloc(2*length,sizeof(REAL));
    c->fft_x2->realSignal = c->x2;

    c->correlation = c->ifft->realSignal;

    if( correlationType == FAST_AUTO )
    {
        c->buffer = NULL;
    }
    else
    {
        c->buffer = (REAL*)malloc(2*length*sizeof(REAL));
    }
    return c;
}

void FreeCorrelation( Correlation* c )
{
    if( c == NULL )
    {
        return;
    }

    if( c->correlationType == FORWARD && c->correlation != NULL )
    {
        free(c->correlation);
    }

    if( c->fft_x1 != NULL )
    {
        FreeFFT( c->fft_x1 );
    }
    if( c->fft_x2 != NULL )
    {
        FreeFFT( c->fft_x2 );
    }
    if( c->ifft != NULL )
    {
        FreeFFT( c->ifft );
    }
    if( c->buffer != NULL )
    {
        free( c->buffer );
    }
    free(c);
}

/*
 *------------------------------------------------------------------------------
 * Returns the first minima from a given correlation function. This would      
 * normally come from the AMDF or YIN functions which generate minima at       
 * multiples of the time period. The threshold determines the value below which
 * the minima must be, ignoring any minima below this value. The interpolated  
 * location (quadratic) is returned. If no local minima is found the global 
 * minimum is returned.
 *------------------------------------------------------------------------------
 */
REAL firstMinima( Correlation* c )
{
    int i;
    REAL xe = 0;
    REAL ye = 0;

    // Compute direction of 1st derivative
    REAL d1_prev = 1;
    REAL d1_next = -1;


    for( i = c->searchStart+1; i<c->searchEnd; i++ )
    {
        if(c->correlation[i-1] <= c->threshold || 
           c->correlation[i] <= c->threshold ||
           c->correlation[i+1] <= c->threshold )
        {
            d1_next = c->correlation[i+1]-c->correlation[i];

            if(d1_prev < 0 && d1_next >= 0)
            {

                interpolate((REAL)i,
                            c->correlation[i-1],
                            c->correlation[i],
                            c->correlation[i+1],
                            &xe,
                            &ye);
                
                return xe;
            }

            if( d1_next < 0 )
            {
                d1_prev = d1_next;
            }

        }
    }

    // not found, return first maxima

    // Compute 1st derivative
    d1_prev = 1;

    for( i = c->searchStart+1; i<c->searchEnd; i++ )
    {
        d1_next = c->correlation[i+1]-c->correlation[i];

            if(d1_prev < 0 && d1_next >= 0)
            {
                interpolate((REAL)i,
                            c->correlation[i-1],
                            c->correlation[i],
                            c->correlation[i+1],
                            &xe,
                            &ye);

                return xe;
            }

            if( d1_next < 0 )
            {
                d1_prev = d1_next;
            }
    }


    // still not found, return 0
    return 0;
}

//------------------------------------------------------------------------------
// Returns the first maxima from a given correlation function. This would 
// normally come from the autocorrelation function which generates maxima at 
// multiples of the time period. The threshold determines the value above which 
// the maxima must be, ignoring any maxima below this value. The interpolated
// location (quadratic) is returned. If no local maxima is found the global 
// maximum is returned.
//------------------------------------------------------------------------------
REAL firstMaxima( Correlation* c )
{
    int i;
    REAL xe = 0;
    REAL ye = 0;

    // Compute 1st derivative
    REAL d1_prev = -1;
    REAL d1_next = 1;


    for( i = c->searchStart+1; i<c->searchEnd; i++ )
    {
        if(c->correlation[i-1] >= c->threshold || 
           c->correlation[i] >= c->threshold ||
           c->correlation[i+1] >= c->threshold )
        {
            d1_next = c->correlation[i+1]-c->correlation[i];

            if(d1_prev > 0 && d1_next <= 0)
            {

                interpolate((REAL)i,
                    c->correlation[i-1],
                    c->correlation[i],
                    c->correlation[i+1],
                    &xe,
                    &ye);
                    return xe;
            }

            if( d1_next > 0 )
            {
                d1_prev = d1_next;
            }

        }
    }

    // not found, return first maxima

    // Compute 1st derivative
    d1_prev = -1;

    for( i = c->searchStart+1; i<c->searchEnd; i++ )
    {
        d1_next = c->correlation[i+1]-c->correlation[i];

        if(d1_prev > 0 && d1_next <= 0)
        {
            interpolate((REAL)i,
                c->correlation[i-1],
                c->correlation[i],
                c->correlation[i+1],
                &xe,
                &ye);
            return xe;
        }

        if( d1_next > 0 )
        {
            d1_prev = d1_next;
        }
    }

    // still not found, return 0
    return 0;

}

//------------------------------------------------------------------------------
// Computes the autocorrelation function directly from the definition.
//------------------------------------------------------------------------------
void autocorrelation(SIGNAL* input, Correlation* c)
{
    int t;
    int i;
    static int k = 0;
    FILE* debug = fopen("correlationslow.txt","a");
    for (t = c->searchStart; t <= c->searchEnd; t++) {
        c->correlation[t] = 0;
        for (i = 0; i <= c->searchEnd; i++) {
            c->correlation[t] += (REAL)(input[i]*input[i + t]);
        }        
        c->correlation[t] = c->correlation[t]/(REAL)c->length;
        if(k==0)
        {
            fprintf(debug,"%f\n", c->correlation[t]);
        }
    }
    fclose(debug);
    k++;

}

//------------------------------------------------------------------------------
// Computes the AMDF function directly from the definition.
//------------------------------------------------------------------------------
void amdf(SIGNAL* input, Correlation* c) 
{
    int t;
    int i;    
    //FILE* debug = fopen("amdfslow.txt","w");
    c->correlation[0]=0;
    for (t = c->searchStart; t <= c->searchEnd; t++) {
        c->correlation[t] = 0;
        for (i = 0; i <= c->length; i++) {
            c->correlation[t] += getabs(input[i] - input[i + t]);
        }        
        c->correlation[t] = c->correlation[t]/(REAL)c->length;
        //fprintf(debug,"%f\n", c->correlation[t]);
    }
    //fclose(debug);
}

//------------------------------------------------------------------------------
// Computes the ASDF function directly from the definition.
//------------------------------------------------------------------------------
void asdf(SIGNAL* input, Correlation* c) 
{
    int t;
    int i;
    //FILE* debug = fopen("asdfslow.txt","w");
    c->correlation[0]=0;
    for (t = c->searchStart; t <= c->searchEnd; t++) {
        c->correlation[t] = 0;
        for (i = 0; i <= c->length; i++) {
            c->correlation[t] += Square(input[i] - input[i + t]);
        }        
        c->correlation[t] = c->correlation[t]/(REAL)c->length;
        //fprintf(debug,"%f\n",c->correlation[t]);
    }
    //fclose(debug);
}

//------------------------------------------------------------------------------
// Computes the YIN function directly from the definition.
//------------------------------------------------------------------------------
void yin(SIGNAL* input, Correlation* c)
{
    int t;
    int i;
    REAL totalSum = 0;
    //FILE* debug = fopen("yinslow.txt","w");

    c->correlation[0]=0;
    for (t = 1; t <= c->searchEnd; t++) {
        c->correlation[t]=0;
        for(i=0; i <= c->length; i++)
        {
            // Compute the average square difference function
            c->correlation[t]+=Square(input[i] - input[i + t]);
        }
        
        // Compute the cumulative mean normalized difference function
        totalSum += c->correlation[t];
        if(totalSum == 0)
        {
            c->correlation[t] = 0;
        }
        else
        {
            c->correlation[t] *= t/totalSum;
        }
        //fprintf(debug,"%f\n",c->correlation[t]);
    }
    //fclose(debug);
}

//------------------------------------------------------------------------------
// Computes the autocorrelation function with the FFT.
//------------------------------------------------------------------------------
void fast_autocorrelation(SIGNAL* input, Correlation* c)
{
    int i;
    int N = 2*c->length;
    kiss_fft_cpx temp;
    //FILE* debug = fopen("correlationfast.txt","w");

    /* compute the FFT of the input */
    c->fft_x1->realSignal = input;    
    FFT_Transform( c->fft_x1 );

    /* compute the FFT of half the input with zero padding */
    for(i=0; i<c->length; i++)
    {
        c->x2[i] = input[i];
    }
    FFT_Transform( c->fft_x2 );

    /* Compute the correlation in the frequency domain */
    for( i=0; i<N; i++ )
    {
        temp.r = (c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].r + 
                  c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].i);
        temp.i = (c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].r - 
                  c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].i);
        c->ifft->spectrum[i].r = temp.r;
        c->ifft->spectrum[i].i = temp.i;
    }    
    
    FFT_Transform( c->ifft );

    /* Half of the correlation is unique */
    for( i=0; i<c->length; i++ )
    {        
        c->correlation[i] = c->correlation[i]/(REAL)c->length;
        //fprintf(debug,"%f\n",c->correlation[i]);
    }
    //fclose(debug);
}

//------------------------------------------------------------------------------
// Computes the ASDF function with the FFT.
//------------------------------------------------------------------------------
void fast_asdf(SIGNAL* input, Correlation* c)
{
    int i;
    int t;
    int N = 2*c->length;
    kiss_fft_cpx temp;
    REAL sumx1squared = 0;
    REAL sumx2squared = 0;
    //FILE* debug = fopen("asdffast.txt","w");

    /* compute the FFT of the input */
    c->fft_x1->realSignal = input;    
    FFT_Transform( c->fft_x1 );

    /* compute the FFT of half the input with zero padding */
    for(i=0; i<c->length; i++)
    {
        c->x2[i] = input[i];
    }
    FFT_Transform( c->fft_x2 );

    /* Compute the sum of x_1 and x_2 squared */
    for( i=0; i<N; i++ )
    {        
        c->buffer[i] = input[i]*input[i];
        if(i<c->length)
        {
            sumx1squared += c->buffer[i];
        }
    }
    sumx2squared = sumx1squared;

    /* Compute the correlation in the frequency domain */
    for( i=0; i<N; i++ )
    {
        temp.r = (c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].r + 
                  c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].i);
        temp.i = (c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].r - 
                  c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].i);
        c->ifft->spectrum[i].r = temp.r;
        c->ifft->spectrum[i].i = temp.i;
    }    
    
    FFT_Transform( c->ifft );

    /* Half of the correlation is unique */
    for( i=0; i<c->length; i++ )
    {
        #ifdef FIXED_POINT
        c->correlation[i] = (-2.0f*c->correlation[i] + sumx1x2_squared)/(REAL)N;
        #else
        // Compute the sum of the square of the current window
        if( i > 0 )
        {
            sumx2squared -= c->buffer[i-1];
            sumx2squared += c->buffer[i-1+c->length];
        }
        c->correlation[i] = (-2.0*c->correlation[i] + sumx1squared + sumx2squared)/c->length;
        //fprintf(debug,"%f\n",c->correlation[i]);
        #endif
    }

    //fclose(debug);
}

//------------------------------------------------------------------------------
// Computes the YIN function with the FFT.
//------------------------------------------------------------------------------
void fast_yin(SIGNAL* input, Correlation* c)
{
    int i;
    int t;
    int N = 2*c->length;
    kiss_fft_cpx temp;
    REAL totalSum = 0;
    REAL sumx1squared = 0;
    REAL sumx2squared = 0;
    //FILE* debug = fopen("yinfast.txt","w");

    /* compute the FFT of the input */
    c->fft_x1->realSignal = input;    
    FFT_Transform( c->fft_x1 );

    /* compute the FFT of half the input with zero padding */
    for(i=0; i<c->length; i++)
    {
        c->x2[i] = input[i];
    }
    FFT_Transform( c->fft_x2 );

    /* Compute the sum of x_1 and x_2 squared */
    for( i=0; i<N; i++ )
    {        
        c->buffer[i] = input[i]*input[i];
        if(i<c->length)
        {
            sumx1squared += c->buffer[i];
        }
    }
    sumx2squared = sumx1squared;

    /* Compute the correlation in the frequency domain */
    for( i=0; i<N; i++ )
    {
        temp.r = (c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].r + 
                  c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].i);
        temp.i = (c->fft_x1->spectrum[i].i*c->fft_x2->spectrum[i].r - 
                  c->fft_x1->spectrum[i].r*c->fft_x2->spectrum[i].i);
        c->ifft->spectrum[i].r = temp.r;
        c->ifft->spectrum[i].i = temp.i;
    }    
    
    FFT_Transform( c->ifft );

    /* Half of the correlation is unique */
    for( i=0; i<c->length; i++ )
    {
        // Compute the sum of the square of the current window
        if( i > 0 )
        {
            sumx2squared -= c->buffer[i-1];
            sumx2squared += c->buffer[i-1+c->length];
        }
        c->correlation[i] = (-2.0*c->correlation[i] + sumx1squared + sumx2squared);        
    }

    // Compute the cumulative mean normalized difference function
    c->correlation[0]=0;
    for (t = 1; t < c->length; t++) {
        // Compute the cumulative mean normalized difference function
        totalSum += c->correlation[t];
        if(totalSum == 0)
        {
            c->correlation[t] = 0;
        }
        else
        {
            c->correlation[t] *= t/totalSum;
        }
        //fprintf(debug,"%f\n",c->correlation[t]);
    }
    //fclose(debug);
}
