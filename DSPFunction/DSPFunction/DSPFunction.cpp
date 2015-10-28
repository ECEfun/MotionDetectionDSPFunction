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
 PIRsum: sum of binary integration of PIR sensor
 
 Output: return sum of PIR and Xband binary integration for decision making (EX: 4 out of 5 -> target)
 */
int Detection::Detector(volatile float *XBsignal,volatile float *PIRsignal, int lenXBsignal, int lenPIRsignal, float noise_Var, float Vt, float signal_Var,int Fs,int *PIRsum)
{
    int Fo = 100;
    int Fo1 = 10;
    int *BinaryIntegrationXB = NULL;
    BinaryIntegrationXB = (int*) malloc(Fo1*sizeof(int));
    int *BinaryIntegrationPIR = NULL;
    BinaryIntegrationPIR = (int*) malloc(Fo1*sizeof(int));
    int lendetection = 20;
    int *stage1Detection = NULL;
    //    int *detection;
    //    detection = (int*) calloc(lendetection,sizeof(int));
    stage1Detection = (int*) malloc(Fo*sizeof(int));
    float *XFiltered = NULL;
    int sum = 0;
    
    float lrt =(2*pow(noise_Var,2)*pow(signal_Var,2)/(pow(signal_Var,2) - pow(noise_Var,2)))*log(Vt) - (Fs/Fo)*log(noise_Var/signal_Var);
    
    Detection::SetValues(XBsignal, lenXBsignal);
    mean = Detection::Mean();
    XFiltered = (float*) malloc(Fs*sizeof(float));
    Detection::WMA(XBsignal,XFiltered,lenXBsignal,9);
    
    
    for (int observation = 0;observation < Fo; observation++)
    {
        float sumxobssq = 0;
        for (int i = 0; i < Fs/Fo; i++)
        {
            sumxobssq = sumxobssq + pow(XFiltered[i + (Fs/Fo)*observation],2);
        }
        if (sumxobssq > lrt)
        {
            stage1Detection[observation] = 1;
        }
        else
            stage1Detection[observation] = 0;
    }
    
    
    for (int observation = 0; observation < Fo1;observation++)
    {
        float sumXB = 0;
        float sumPIR = 0;
        for (int i = 0; i < Fo1; i++)
        {
            sumXB = sumXB + stage1Detection[i+observation*Fo1];
            sumPIR = sumPIR + PIRsignal[i+observation*Fo1];
            
        }
        if (sumXB > 7)
            BinaryIntegrationXB[observation] = 1;
        else
            BinaryIntegrationXB[observation] = 0;
        if (sumPIR > 4)
            BinaryIntegrationPIR[observation] = 1;
        else
            BinaryIntegrationPIR[observation] = 0;
    }
    // for testing array of decition
    //    for (int i = 0; i < lendetection/2;i++)
    //    {
    //        detection[i] = BinaryIntegrationXB[i];
    //        detection[lendetection-i-1] = BinaryIntegrationPIR[i];
    //    }
    //
    //    return detection;
    int sumPIR = 0;
    for (int i = 0; i < lendetection/2; i++)
    {
        sum = sum + BinaryIntegrationXB[i] + BinaryIntegrationPIR[i];
        sumPIR = sumPIR + BinaryIntegrationPIR[i];
    }
    (*PIRsum) = sumPIR;
    free(BinaryIntegrationXB);
    free(BinaryIntegrationPIR);
    free(stage1Detection);
    free(XFiltered);
    free(value);
    BinaryIntegrationXB = NULL;
    BinaryIntegrationPIR = NULL;
    stage1Detection = NULL;
    XFiltered = NULL;
    value = NULL;
    return sum;
    // if (sum > 10)
    //     return 1;
    // else
    //     return 0;
}