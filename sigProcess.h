#pragma once
#define _CRT_SECURE_NO_WARNINGS
#ifndef _SIGPROCESS_H_
#define _SIGPROCESS_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"hmath.h"


void circshift(Vector v, int shift);
int find(Vector v, double thre, int FrontOrEndFlag);
void pad_signal(Vector* yP, Vector x, int Npad);
void unpad_signal(Vector* yP, Vector x, int res, int target_sz );

/*frame the origin int signal according to frame size and frame shift size*/
/*return a double matrix,each row contains a frame*/
Matrix frameRawSignal(IntVec v, int wlen, int inc,double preEmphasiseCoefft,int enableHamWindow);




/*---------预处理----------------*/
/*zero mean a complete speech waveform nSamples long*/
void ZeroMean(IntVec data);
/*预加重,k一般取0.9-1，如果取k=0,不做任何预加重*/
void PreEmphasise(Vector s, double k);

/*-----------------计算谱中心以及子带能量---------------*/
/*计算谱中心,相对位置*/
double calBrightness(Vector fftx);
/*计算子带能量,输出子带能量占总能量的比值*/
void calSubBankE(Vector fftx, Vector subBankEnergy);
/*计算过零率，返回过零次数除以一帧的采样点的个数*/
double zeroCrossingRate(Vector s, int frameSize);

/*计算差分系数*/
/*参数分别为向量，需要差分的每帧的采样点数，帧数，每次移动的步长，差分系数和原信号的距离，0，0，是否简单差分*/
void Regress(double* data, int vSize, int n, int step, int offset, int delwin, int head, int tail, int simpleDiffs);

void RegressMat(Matrix* m,int delwin, int regressOrder);

void NormaliseLogEnergy(double *data, int n, int step, double silFloor, double escale);

void ZNormalize(double *data, int vSize, int n, int step);

/* GenHamWindow: generate precomputed Hamming window function */
Vector GenHamWindow(int frameSize);
/*Apply Hamming Window to Speech frame s*/
void Ham(Vector s);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // !_SIGPROCESS_H_
