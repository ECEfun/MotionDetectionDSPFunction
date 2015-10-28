//
//  DSPFunction.cpp
//
//  Created by Linh Nguyen on 10/23/15.
//  Copyright (c) 2015 Linh Nguyen. All rights reserved.
//

#include "DSPFunction.h"

// #include <sstream>
// #include <math.h>
#include <array>

int max;
float *value = NULL;
float mean;

//convolution algorithm
void Detection::killall()
{
    if (value != NULL)
    {
        free(value);
        value = NULL;
    }
}
void Detection::conv(float *A, float *B, int lenA, int lenB,float *C, int *lenC)
{
    int nconv;
    int i, j, i1;
    float tmp;
    // float *C;
    
    //allocated convolution array
    nconv = lenA+lenB-1;
    // C = (float*) calloc(nconv, sizeof(float));
    
    //convolution process
    for (i=0; i<nconv; i++)
    {
        i1 = i;
        tmp = 0.0;
        for (j=0; j<lenB; j++)
        {
            if(i1>=0 && i1<lenA)
                tmp = tmp + (A[i1]*B[j]);
            
            i1 = i1-1;
            C[i] = tmp;
        }
    }
    
    //get length of convolution array
    (*lenC) = nconv;
    
    //return convolution array
    // return(C);
}

float Detection::Mean()
{
    float sum = 0;
    for(int i = 0; i < max; i++)
        sum += value[i];
    return (sum / max);
}

float Detection::Variane()
{
    mean = Mean();
    
    float temp = 0;
    for(int i = 0; i < max; i++)
    {
        temp += (value[i] - mean) * (value[i] - mean) ;
    }
    return temp / max;
}

float Detection::SampleVariane()
{
    mean = Mean();
    
    double temp = 0;
    for(int i = 0; i < max; i++)
    {
        temp += (value[i] - mean) * (value[i] - mean) ;
    }
    return temp / (max - 1);
}

int Detection::SetValues(volatile float *p, int lenp)
{
    if (value != NULL)
    {
        free(value);
        value = NULL;
        value = (float*) malloc(lenp*sizeof(float));
    }
    else
    {
        value = (float*) malloc(lenp*sizeof(float));
    }
    // if(count > lenp)
    //     return -1;
    max = lenp;
    for(int i = 0; i < lenp; i++)
        value[i] = p[i];
    return 0;
}

float Detection::GetStandardDeviation()
{
    return sqrt(Variane());
}

float Detection::GetSampleStandardDeviation()
{
    return sqrt(SampleVariane());
}

/*Weighted Moving Average
 Inputs:
 signal: signal array
 wsignal: output filtered array
 lensignal: length of signal
 window: moving average window
 */

void Detection::WMA(volatile float *signal,float *wsignal,int lensignal, int window)
{
    //Set WMA polynomial
    // float *wsignal;
    // wsignal = (float*) calloc(lensignal, sizeof(float));
    
    float *binomialCoeff = NULL;
    int lenbinomialCoeff;
    float h[2] = {0.5,0.5};
    binomialCoeff = (float*) malloc(window*sizeof(float));
    
    
    conv(h,h,2,2,binomialCoeff,&lenbinomialCoeff);
    
    for (int n = 1;n<=window-3;n++)
    {
        conv(binomialCoeff,h,lenbinomialCoeff,2,binomialCoeff,&lenbinomialCoeff);
    }
    
    // initialize sum. In theory, the first filter output element are sum from 1 to number of window - window/2
    
    for (int i = 0;i<(lensignal);i++)
    {
        float sum = 0;
        if (i < window - int(window/2))
        {
            for (int n = (window/2)-i;n < window;n++)
            {
                sum = sum+binomialCoeff[n]*(signal[n-window/2 + i]-mean);
            }
            wsignal[i] = sum;  //should be devide by sum of weight, but sum of weight = 1 in this case
        }
        else if (i < lensignal - int(round(window/2)))
        {
            for (int n = 0;n<window;n++)
            {
                sum = sum+binomialCoeff[n]*(signal[n+i-int(window/2)]-mean);
            }
            wsignal[i] = sum;  //should be devide by sum of weight, but sum of weight = 1 in this case
        }
        else
        {
            wsignal[i] = signal[i] - mean;
        }
        
    }
    free(binomialCoeff);
    binomialCoeff = NULL;
    // return wsignal;
    
}

/*
 Detector function
 Input:
 XBsignal: XBand signal array
 PIRsignal: PIR signal array
 lenXBsignal: length of XBand array
 lenPIRsignal: length of PIR array
 noise_Var: Noise variance
 Vt: Threshold using log likelihood ratio
 signal_Var: Signal variance
 Fs: Sampling frequency
 meanXBVal : mean value of XBand, can also use Mean function but have to setValue first.
 PIRsum: sum of binary integration of PIR sensor
 
 Output: return sum of PIR and Xband binary integration for decision making (EX: 4 out of 5 -> target)
 */
int Detection::Detector(volatile float *XBsignal,volatile float *PIRsignal, int lenXBsignal, int lenPIRsignal, float noise_Var, float Vt, float signal_Var,int Fs,float meanXBVal,int *PIRsum)
{
    int Fo = 100;  //default for stage 2, 10 observation
    int Fo1 = 10;  // default stage 1, 10 observation
    int *stage1Detection = NULL;
    stage1Detection = (int*) malloc(Fo*sizeof(int));
    int window = 7;  //moving average window, using Matlab to calculate window
    float *binomialCoeff = NULL;
    int lenbinomialCoeff;
    float h[2] = {0.5,0.5};
    binomialCoeff = (float*) malloc(window*sizeof(float));
    
    
    conv(h,h,2,2,binomialCoeff,&lenbinomialCoeff);
    
    for (int n = 1;n<=window-3;n++)
    {
        conv(binomialCoeff,h,lenbinomialCoeff,2,binomialCoeff,&lenbinomialCoeff);
    }
    
    float lrt =(2*pow(noise_Var,2)*pow(signal_Var,2)/(pow(signal_Var,2) - pow(noise_Var,2)))*log(Vt) - (Fs/Fo)*log(noise_Var/signal_Var);
    
    // Detection::SetValues(XBsignal, lenXBsignal);
    mean = meanXBVal;
    // XFiltered = (float*) malloc(Fs*sizeof(float));
    // Detection::WMA(XBsignal,XFiltered,lenXBsignal,7);
    // do not need mean value anymore free memory
    //120 usecond up to here
    int sumXB = 0;
    int sumPIR = 0;
    int stage2Observation = 0;
    int sumXBandstage2 = 0;
    int sumPIRstage2 = 0;
    int k;
    for (int observation = 0;observation < Fo; observation++)
    {
        float sumxobssq = 0;
        // calculate weighted moving average for each element
        for (int i = 0; i < Fs/Fo; i++)
        {
            float WMA_k = 0;
            k  = i + (Fs/Fo)*observation;
            
            if (k < window - int(window/2))
            {
                for (int n = (window/2)-i;n < window;n++)
                {
                    WMA_k = WMA_k+binomialCoeff[n]*(XBsignal[n-window/2 + k]-mean);
                }
                // wsignal[k] = WMA_k;  //should be devide by sum of weight, but sum of weight = 1 in this case
            }
            else if (k < Fs - int(round(window/2)))
            {
                for (int n = 0;n<window;n++)
                {
                    WMA_k = WMA_k+binomialCoeff[n]*(XBsignal[n+k-int(window/2)]-mean);
                }
                // wsignal[k] = WMA_k;  //should be devide by sum of weight, but sum of weight = 1 in this case
            }
            else
            {
                WMA_k = XBsignal[k] - mean;
            }
            //end weighted moving average
            sumxobssq = sumxobssq + pow(WMA_k,2);  // calculate N*R^2 for log likelihood ratio test
        }
        if (sumxobssq > lrt)
        {
            sumXB = sumXB + 1;
        }
        sumPIR = sumPIR + PIRsignal[observation];
        stage2Observation++;
        if (stage2Observation == Fo1)
        {
            if (sumXB > 7)
                sumXBandstage2 =  sumXBandstage2 + 1;
            if (sumPIR > 7)
                sumPIRstage2 = sumPIRstage2 + 1;
            
            stage2Observation = 0;
            sumPIR = 0;
            sumXB = 0;
        }
    }
    
    (*PIRsum) = sumPIRstage2;
    free(stage1Detection);
    //  free(XFiltered);
    stage1Detection = NULL;
    //  XFiltered = NULL;
    return (sumPIRstage2 + sumXBandstage2);
}



