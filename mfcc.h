#pragma once
#ifndef _MFCC_H_H
#define _MFCC_H_H


#include"hmath.hpp"
#include<cmath>

//����MFCC����
//MFCC�˲�������Ϣ�������˳�ʼ��
typedef struct {
	int frameSize;       /* speech frameSize */
	int numChans;        /* number of channels */
	double sampPeriod;     /* sample period */
	int fftN;            /* fft size */
	int klo, khi;         /* lopass to hipass cut-off fft indices */
	int usePower;    /* use power rather than magnitude */
	int takeLogs;    /* log filterbank channels */
	double fres;          /* scaled fft resolution */
	hmath::Vector cf;           /* array[1..pOrder+1] of centre freqs */
	hmath::ShortVec loChan;     /* array[1..fftN/2] of loChan index */
	hmath::Vector loWt;         /* array[1..fftN/2] of loChan weighting */
	hmath::Vector x;            /* array[1..fftN] of fftchans */
}FBankInfo;

double Mel(int k, double fres);
/*
return mel-frequency corresponding to given FFT index k.
Resolution is normally determined by fres field of FBankInfo
record.
*/

FBankInfo InitFBank(int frameSize, double sampPeriod, int numChans,
	double lopass, double hipass, int usePower, int takeLogs,
	int doubleFFT,
	double alpha, double warpLowCut, double warpUpCut);
/*
Initialise an FBankInfo record prior to calling Wave2FBank.
*/

void Wave2FBank(hmath::Vector s, hmath::Vector fbank, double *te,double* te2, FBankInfo info);
/*
Convert given speech frame in s into mel-frequency filterbank
coefficients.  The total frame energy is stored in te.  The
info record contains precomputed filter weights and should be set
prior to using Wave2FBank by calling InitFBank.
*/

void FBank2MFCC(hmath::Vector fbank, hmath::Vector c, int n);
/*
Apply the DCT to fbank and store first n cepstral coeff in c.
Note that the resulting coef are normalised by sqrt(2/numChans)
*/

double FBank2C0(hmath::Vector fbank);
/* EXPORT->FBank2C0: return zero'th cepstral coefficient */

void FFT(hmath::Vector s, int invert);
/*
When called s holds nn complex values stored in the
sequence   [ r1 , i1 , r2 , i2 , .. .. , rn , in ] where
n = VectorSize(s) DIV 2, n must be a power of 2. On exit s
holds the fft (or the inverse fft if invert == 1)
*/

void Realft(hmath::Vector s);
/*
When called s holds 2*n real values, on exit s holds the
first  n complex points of the spectrum stored in
the same format as for fft
*/



#endif