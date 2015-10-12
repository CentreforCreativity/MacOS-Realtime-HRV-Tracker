//
//  PowerSpectrum.c
//  HRV_Project
//
//  Created by Gareth Loudon on 07/09/2015.
//  Copyright (c) 2015 Gareth Loudon. All rights reserved.
//


#include "PowerSpectrum.h"
#import "HeartRateMonitorAppDelegate.h"


struct HRVData hrvValues;
FILE *fpPS, *fpHR;

/*---------------------- Testing communication with the Main App Delegate ---------------------------*/
 int analyseData(float timeStampreceived, NSInteger heartBeatreceived, int *calmingLevel)
 {
     
     SaveHeartRateData((double)timeStampreceived, (int)heartBeatreceived);
     
     PowerSpectrum((double)heartBeatreceived, (double)timeStampreceived);
     
     *calmingLevel = hrvValues.calmingLevel;
     
     return hrvValues.attentionLevel;
     
     
 }



//update the interface after performing the calculation
void sendData(char dataString, void *sender)
{
    
    
 //   [sender updateInterface: dataString];
    
    
}




/*---------------------- Testing communication with the Main App Delegate ---------------------------*/

void PowerSpectrum(double hrCurrent, double timeCurrent)
{
    static int numSamples = 0;
    static double heartrate[NFFT*2];
    static double timeStamp[NFFT*2];
    double interval[NFFT*2];
    int numIntervals, i;
    double windowedData[NFFT][2];
    double fftdata[NFFT][2];
    double powerspectrum[NFFT];
    double s2;
    
    
    heartrate[numSamples] = hrCurrent;
    timeStamp[numSamples] = timeCurrent;
    
    numSamples++;
    
    if(numSamples > (NFFT/SAMPLING_FREQ))
    {
        numIntervals = ResampleHeartRateData(numSamples, heartrate, timeStamp, interval);
        
        if(numIntervals > NFFT)
        {
            CalcAverageHeartRate(interval);
            
            ArtefactCorrection(interval);
            
            s2 = ApplyHanningWindow(interval, windowedData);
            
            FastFourierTransform(windowedData, fftdata);
            
            CalcPowerSpectrum(fftdata, s2, powerspectrum);
            
            hrvValues.attentionLevel = CalcAttentionLevel(hrvValues.averagePower);
            
            hrvValues.calmingLevel = CalcCalmingLevel(hrvValues.peakPowerRatio);
            
            SavePowerSpectrumData(powerspectrum);
            
            for(i=0; i<(numSamples - WINDOW_INCREMENT); i++)
            {
                heartrate[i] = heartrate[i+WINDOW_INCREMENT];
                timeStamp[i] = timeStamp[i+WINDOW_INCREMENT];
            }
            numSamples = numSamples - WINDOW_INCREMENT;
        }
        
    }
    
}



/* -------------------------- Resample heart rate data at the SAMPLING_FREQ ------------------------------*/

int ResampleHeartRateData(int N, double *heartrate, double *timeStamp, double *interval)
{
    int i, k;
    double sampleTime, timeDiff;
    double newTotalTime, newTime[500*SAMPLING_FREQ];
    
    interval[0] = ((double)60000.0/heartrate[0]);
    newTime[0] = timeStamp[0];
    newTotalTime = newTime[0] + (1.0 / (double)SAMPLING_FREQ);
    
    k = 1;
    
    for(i=1; i<N; i++)
    {
        timeDiff = timeStamp[i] - timeStamp[i-1];
        sampleTime = newTotalTime - timeStamp[i-1];
        
        while(timeDiff > sampleTime)
        {
            interval[k] = heartrate[i-1] + (((heartrate[i] - heartrate[i-1])*sampleTime)/timeDiff);
            interval[k] = ((double)60000.0/interval[k]);
            
            sampleTime = sampleTime + (1.0 / (double)SAMPLING_FREQ);
            newTime[k] = newTime[k-1] + sampleTime;
            newTotalTime = newTotalTime + (1.0 / (double)SAMPLING_FREQ);
            k++;
        }
    }
    
    return k;
}





/* ----------------------- Higher Pass Filter the RR intervals with HP_FILTER_CUTOFF setting ----------*/

void HighPassFilterHeartRateData(double *interval, double *filteredArray)
{
    double RC, dt, alpha;
    int i;
    
    RC = 1.0/(HP_FILTER_CUTOFF*TWO_PI);
    dt = 1.0/SAMPLING_FREQ;
    alpha = RC/(RC + dt);
    
    filteredArray[0] = 0.0;
    
    //   filteredArray[0] = interval[0];
    
    for (i = 1; i<NFFT; i++)
    {
        //    filteredArray[i] = alpha * (filteredArray[i-1] + interval[i] - interval[i-1]);
        filteredArray[i] = interval[i] - interval[0];
        
    }
    
}


/* ----------------------- Correct for artefacts in interval data ----------*/

void ArtefactCorrection(double filteredData[NFFT])
{
    int i, j;
    double localMeanInterval;
    
    
    for (i=LOCAL_WINDOW; i<(NFFT-LOCAL_WINDOW); i++)
    {
        localMeanInterval = 0.0;
        
        for(j = i; j < (i+LOCAL_WINDOW); j++)
        {
            localMeanInterval = localMeanInterval + filteredData[j-(LOCAL_WINDOW/2)];
        }
        
        localMeanInterval = localMeanInterval / (double)LOCAL_WINDOW;
        
        
        if(filteredData[i] > (ARTEFACT_THRESHOLD+localMeanInterval))
        {
            filteredData[i] = (filteredData[i-1] + filteredData[i+1])/2;
        }
        else if(filteredData[i] < (localMeanInterval-ARTEFACT_THRESHOLD))
        {
            filteredData[i] = (filteredData[i-1] + filteredData[i+1])/2;
        }
    }
    
}


/* --------------------------------- Calculate Average Heart Rate ---------------------------------------*/


void CalcAverageHeartRate(double *interval)
{
    int i;
    
    hrvValues.heartRate = 0.0;
    
    for(i = 0; i < NFFT; i++)
    {
        hrvValues.heartRate = hrvValues.heartRate + ((double)60000.0/interval[i]);
    }
    
    hrvValues.heartRate = hrvValues.heartRate / (double)NFFT;
}




/* ---------------------------------Hanning Window ------------------------------------------------*/

double ApplyHanningWindow(double *interval, double (*windowedData)[2])
{
    int i;
    double windowValue, s1, s2;
    
    s1 = s2 = 0.0;
    
    for(i = 0; i < NFFT; i++)
    {
        windowValue = (double)0.54 - ((double)0.46 * cos(  (  (double)TWO_PI * i)   /(NFFT-1))  );
        windowedData[i][0] = interval[i] * windowValue;
        windowedData[i][1] = 0.0;
        
        s1 = windowValue + s1;
        s2 = (windowValue*windowValue) + s2;
        
    }
    
    return s2;
}

/* ---------------Calculate power spectrum and total power in LF and HF ranges ----------------------*/


void CalcPowerSpectrum(double (*fftdata)[2], double s2, double *powerspectrum)
{
    int i, lowFreqBinCount, highFreqBinCount, localPeakBin;
    double freqValue, freqRangeBin, timePeriod, deltaTime;
    double totalPower[2], localPeakPower;
    static double averagePS[NFFT];
    static int count = 0;
    int minBinValue, maxBinValue;
    
    totalPower[0] = 0.0;
    totalPower[1] = 0.0;
    
    lowFreqBinCount = 0;
    highFreqBinCount = 0;
    localPeakBin = 0;
    minBinValue = 0;
    maxBinValue = 0;
    
    freqRangeBin = (double)SAMPLING_FREQ / (double)NFFT;
    
    freqValue = 0.0;
    
    timePeriod = (double)NFFT / (double)SAMPLING_FREQ;
    deltaTime = (double)1.0 / (double)SAMPLING_FREQ;
    
    localPeakPower = 0.0;
    
    for(i=0; i < (NFFT/2); i++)
    {
        powerspectrum[i] = (double)2.0*((fftdata[i][0] * fftdata[i][0]) + (fftdata[i][1] * fftdata[i][1])) / (s2 * (double)SAMPLING_FREQ);
        
        if(count == 0)
        {
            averagePS[i] = powerspectrum[i];
        }
        else
        {
            averagePS[i] = (1.00*powerspectrum[i]) + (0.00*averagePS[i]);
        }
        
        if(freqValue > LF_MIN_FREQ && freqValue < LF_MAX_FREQ)
        {
            totalPower[0] = totalPower[0] + averagePS[i];
            
            if(lowFreqBinCount == 0)
                minBinValue = i;
            maxBinValue = i;

            lowFreqBinCount++;
            
            if(powerspectrum[i] > localPeakPower)
            {
                localPeakPower = powerspectrum[i];
                localPeakBin = i;
            }
        }
        
        if(freqValue > HF_MIN_FREQ && freqValue < HF_MAX_FREQ)
        {
            totalPower[1] = totalPower[1] + averagePS[i];
            highFreqBinCount++;
        }
        
        freqValue = freqValue + freqRangeBin;
        
    }
    
    if(count == 0)
    {
        count++;
    }
    
    hrvValues.averagePower[0] = totalPower[0] / timePeriod;
    hrvValues.averagePower[1] = totalPower[1] / timePeriod;
    
    hrvValues.peakPower = powerspectrum[localPeakBin] / timePeriod;
 
    
    if((NFFT/SAMPLING_FREQ) > 60)
    {
        if(localPeakBin > minBinValue && localPeakBin < maxBinValue)
        {
            hrvValues.peakPowerRatio = (powerspectrum[localPeakBin-1] + powerspectrum[localPeakBin] + powerspectrum[localPeakBin+1])/(totalPower[0] - powerspectrum[localPeakBin-1] - powerspectrum[localPeakBin] - powerspectrum[localPeakBin+1] + totalPower[1]);
        }
        else if(localPeakBin == minBinValue && localPeakBin < maxBinValue)
        {
            hrvValues.peakPowerRatio = (powerspectrum[localPeakBin] + powerspectrum[localPeakBin+1])/(totalPower[0] - powerspectrum[localPeakBin] - powerspectrum[localPeakBin+1] + totalPower[1]);
        }
        else if(localPeakBin > minBinValue && localPeakBin == maxBinValue)
        {
            hrvValues.peakPowerRatio = (powerspectrum[localPeakBin-1] + powerspectrum[localPeakBin])/(totalPower[0] - powerspectrum[localPeakBin-1] - powerspectrum[localPeakBin] + totalPower[1]);
        }
        else
        {
            hrvValues.peakPowerRatio = (powerspectrum[localPeakBin])/(totalPower[0] - powerspectrum[localPeakBin] + totalPower[1]);
        }
    }
    else
    {
        hrvValues.peakPowerRatio = (powerspectrum[localPeakBin])/(totalPower[0] - powerspectrum[localPeakBin] + totalPower[1]);
    }
    
}


/* ------------------------ Calculate the attention level from range 1 - 5  ---------------------------------*/

int CalcAttentionLevel(double *averagePower)
{
    int attentionLevel;
    
    attentionLevel = 0;
    
    if(averagePower[0] > 200)
    {
        attentionLevel = 1;
    }
    else if (averagePower[0] > 150)
    {
        attentionLevel = 2;
    }
    else if (averagePower[0] > 100)
    {
        attentionLevel = 3;
    }
    else if (averagePower[0] > 50)
    {
        attentionLevel = 4;
    }
    else
    {
        attentionLevel = 5;
    }
    
    
    return attentionLevel;
}




/* ------------------------ Calculate the calming level from range 1 - 5  ---------------------------------*/


int CalcCalmingLevel(double peakPowerRatio)
{
    int calmingLevel;
    
    calmingLevel = 0;
    
    if(peakPowerRatio < 0.50)
    {
        calmingLevel = 1;
    }
    else if (peakPowerRatio < 1.25)
    {
        calmingLevel = 2;
    }
    else if (peakPowerRatio < 2.50)
    {
        calmingLevel = 3;
    }
    else if (peakPowerRatio < 7.00)
    {
        calmingLevel = 4;
    }
    else
    {
        calmingLevel = 5;
    }
    
    
    return calmingLevel;
}

/* ----------------------------------------  FFT --------------------------------------------------------------------*/

void FastFourierTransform(double (*intervalsTime)[2], double (*IntervalsFreq)[2])
{
    double freqValues[2*NFFT][2];
    
    /* Calculate FFT by a recursion. */
    fft_rec(NFFT, 0, 1, intervalsTime, IntervalsFreq, freqValues);
    
}

/* FFT recursion */
void fft_rec(int N, int offset, int delta,
             double (*intervalsTime)[2], double (*intervalsFreq)[2], double (*freqValues)[2])
{
    int N2 = N/2;            /* half the number of points in FFT */
    int k;                   /* generic index */
    double cs, sn;           /* cosine and sine */
    int k00, k01, k10, k11;  /* indices for butterflies */
    double tmp0, tmp1;       /* temporary storage */
    
    if(N != 2)  /* Perform recursive step. */
    {
        /* Calculate two (N/2)-point DFT's. */
        fft_rec(N2, offset, 2*delta, intervalsTime, freqValues, intervalsFreq);
        fft_rec(N2, offset+delta, 2*delta, intervalsTime, freqValues, intervalsFreq);
        
        /* Combine the two (N/2)-point DFT's into one N-point DFT. */
        for(k=0; k<N2; k++)
        {
            k00 = offset + k*delta;
            k01 = k00 + N2*delta;
            
            k10 = offset + 2*k*delta;
            k11 = k10 + delta;
            
            cs = cos(TWO_PI*k/(double)N);
            sn = sin(TWO_PI*k/(double)N);
            
            tmp0 = cs * freqValues[k11][0] + sn * freqValues[k11][1];
            tmp1 = cs * freqValues[k11][1] - sn * freqValues[k11][0];
            
            intervalsFreq[k01][0] = freqValues[k10][0] - tmp0;
            intervalsFreq[k01][1] = freqValues[k10][1] - tmp1;
            intervalsFreq[k00][0] = freqValues[k10][0] + tmp0;
            intervalsFreq[k00][1] = freqValues[k10][1] + tmp1;
        }
    }
    else  /* Perform 2-point DFT. */
    {
        k00 = offset;
        k01 = k00 + delta;
        
        intervalsFreq[k01][0] = intervalsTime[k00][0] - intervalsTime[k01][0];
        intervalsFreq[k01][1] = intervalsTime[k00][1] - intervalsTime[k01][1];
        intervalsFreq[k00][0] = intervalsTime[k00][0] + intervalsTime[k01][0];
        intervalsFreq[k00][1] = intervalsTime[k00][1] + intervalsTime[k01][1];
    }
}


/*---------------------- Save Heart Rate Data ---------------------------*/

void SaveHeartRateData(double timeStampreceived, int heartBeatreceived)
{
    static int count = 0;
    char *timeString, filename[255];
    time_t rawtime;
    struct tm *timeinfo;
    int len;
    
    if(count == 0)
    {
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        timeString = asctime(timeinfo);
        len = strlen(timeString);
        timeString[len-1] = '\0';
        strcpy(filename,"/Users/garethloudon/Dropbox/Innovate UK Project Code/Logged_HeartRate_Data/HR ");
        strcat(filename,timeString);
        strcat(filename, ".txt");
        fpHR = fopen(filename, "w");
        count++;
    }
    
    fprintf(fpHR,"%.3lf %d\n", timeStampreceived, heartBeatreceived);
    
}
/*---------------------- Save Power Spectrum Data ---------------------------*/

void SavePowerSpectrumData(double *powerspectrum)
{
    static int count = 0;
    char *timeString, filename[255];
    time_t rawtime;
    struct tm *timeinfo;
    int len, i;
    
    if(count == 0)
    {
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        timeString = asctime(timeinfo);
        len = strlen(timeString);
        timeString[len-1] = '\0';
        strcpy(filename,"/Users/garethloudon/Dropbox/Innovate UK Project Code/Logged_PowerSpec_Data/PS ");
        strcat(filename,timeString);
        strcat(filename, ".txt");
        fpPS = fopen(filename, "w");
        count++;
    }
    
    
    fprintf(fpPS,"%.3lf, ", hrvValues.heartRate);
    
    fprintf(fpPS,"%.3lf, ", hrvValues.averagePower[0]);
    fprintf(fpPS,"%.3lf, ", hrvValues.averagePower[1]);
    
    fprintf(fpPS,"%.3lf, ", hrvValues.peakPower);
    
    fprintf(fpPS,"%.3lf, ", hrvValues.peakPowerRatio);
    
    fprintf(fpPS,"%d, ", hrvValues.attentionLevel);
    fprintf(fpPS,"%d\n", hrvValues.calmingLevel);
    
    for(i=0; i<(NFFT/2); i++)
    {
        fprintf(fpPS,"%.3f, ", powerspectrum[i]);

    }
    fprintf(fpPS,"\n");


    
}


