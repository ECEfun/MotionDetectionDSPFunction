//
//  DSPFunction
//
//  Created by Linh Nguyen on 10/23/15.
//  Copyright (c) 2015 Linh Nguyen. All rights reserved.
//

#ifndef __DSPFunction__functions__
#define __DSPFunction__functions__
// #include <application.h>
#include <stdio.h>
#include <string>
#include <math.h>
#include <stdlib.h>
#endif /* defined(__DSPFunction__functions__) */

class Detection
{
public:
    void killall();
    
    void conv(float *A, float *B, int lenA, int lenB,float *C, int *lenC);
    
    float Mean();
    
    float Variane();
    
    float SampleVariane();
    
    int SetValues(volatile float *p, int count);
    
    float GetStandardDeviation();
    
    float GetSampleStandardDeviation();
    
    void WMA(volatile float *signal,float *wsignal, int lensignal, int window);
    
    int Detector(volatile float *XBsignal, volatile float *PIRsignal, int lenXBsignal, int lenPIRsignal, float noise_Var, float Vt, float signal_Var, int Fs,int *PIRsum);
    
private:
    
};