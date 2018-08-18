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




/*---------Ԥ����----------------*/
/*zero mean a complete speech waveform nSamples long*/
void ZeroMean(IntVec data);
/*Ԥ����,kһ��ȡ0.9-1�����ȡk=0,�����κ�Ԥ����*/
void PreEmphasise(Vector s, double k);

/*-----------------�����������Լ��Ӵ�����---------------*/
/*����������,���λ��*/
double calBrightness(Vector fftx);
/*�����Ӵ�����,����Ӵ�����ռ�������ı�ֵ*/
void calSubBankE(Vector fftx, Vector subBankEnergy);
/*��������ʣ����ع����������һ֡�Ĳ�����ĸ���*/
double zeroCrossingRate(Vector s, int frameSize);

/*������ϵ��*/
/*�����ֱ�Ϊ��������Ҫ��ֵ�ÿ֡�Ĳ���������֡����ÿ���ƶ��Ĳ��������ϵ����ԭ�źŵľ��룬0��0���Ƿ�򵥲��*/
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
