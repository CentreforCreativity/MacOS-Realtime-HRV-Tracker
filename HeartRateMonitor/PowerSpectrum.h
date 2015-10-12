//
//  PowerSpectrum.h
//  HRV_Project
//
//  Created by Gareth Loudon on 07/09/2015.
//  Copyright (c) 2015 Gareth Loudon. All rights reserved.
//

#ifndef __HRV_Project__PowerSpectrum__
#define __HRV_Project__PowerSpectrum__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* macros */
#define TWO_PI (6.2831853071795864769252867665590057683943L)
#define NFFT 128              // number of samples in FFT
#define SAMPLING_FREQ 2       // resampling frequency of the heart rate data
#define HP_FILTER_CUTOFF 0.005

#define LF_MIN_FREQ 0.07
#define LF_MAX_FREQ 0.15
#define HF_MIN_FREQ 0.15
#define HF_MAX_FREQ 0.40



#define WINDOW_INCREMENT 8    // number of seconds to increment window
#define LOCAL_WINDOW (SAMPLING_FREQ * 3)
#define ARTEFACT_THRESHOLD 350


/* data structures */
struct HRVData
{
    double heartRate = 0.0;
    double averagePower[2];
    double peakPower = 0.0;
    double peakPowerRatio = 0.0;
    int attentionLevel = 0;
    int calmingLevel = 0;
};



//void analyseData( int thedata);
int analyseData(float timeStampreceived, NSInteger heartBeatreceived, int *calmingLevel);

void sendData(char dataString, void *sender);

/* function prototypes */

void updateInterface();

void PowerSpectrum(double hrCurrent, double timeCurrent);
void ReadHeartRateFileData(int N, double *heartrate, double *timeStamp);

int ResampleHeartRateData(int N, double *heartrate, double *timeStamp, double *interval);
void HighPassFilterHeartRateData(double *interval, double *filteredArray);
void ArtefactCorrection(double *interval);

void CalcAverageHeartRate(double *interval);
double ApplyHanningWindow(double *interval, double (*windowedData)[2]);
void FastFourierTransform(double (*x)[2], double (*X)[2]);
void fft_rec(int N, int offset, int delta, double (*x)[2], double (*X)[2], double (*XX)[2]);
void CalcPowerSpectrum(double (*fftdata)[2], double s2, double *powerspectrum);

int CalcAttentionLevel(double *averagePower);
int CalcCalmingLevel(double peakPowerRatio);

void SaveHeartRateData(double timeStampreceived, int heartBeatreceived);
void SavePowerSpectrumData(double *powerspectrum);



#endif /* defined(__HRV_Project__PowerSpectrum__) */
