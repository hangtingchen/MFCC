#pragma once

#ifndef _SIGPROCESS_H_
#define _SIGPROCESS_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include"hmath.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>



/*---------预处理----------------*/
/*htk模块*/
/*zero mean a complete speech waveform nSamples long*/
void ZeroMean(short *data, long nSamples);
/*Apply Hamming Window to Speech frame s*/
void PreEmphasise(Vector s, double k);
/* GenHamWindow: generate precomputed Hamming window function */
void GenHamWindow(int frameSize);
/*Apply Hamming Window to Speech frame s*/
void Ham(Vector s);

/*自己的模块*/
/*读入二进制文件，生成十进制文件*/
Vector pcmToData(FILE* f, long fileLength);
Vector* pcmToData2(FILE* f, long fileLength, int biteNum, int perSample, int vecNum);
/*计算过零率*/
double zeroCrossingRate(Vector s, int frameSize);

/*----------------后期处理-------------*/
/*htk模块*/
/* EXPORT->WeightCepstrum: Apply cepstral weighting to c */
void WeightCepstrum(Vector c, int start, int count, int cepLiftering);
/* EXPORT->UnWeightCepstrum: Undo cepstral weighting of c */
void UnWeightCepstrum(Vector c, int start, int count, int cepLiftering);
/* EXPORT->FZeroMean: Zero mean the given data sequence */
void FZeroMean(double *data, int vSize, int n, int step);
/* EXPORT->AddRegression: add regression vector at +offset from source vector */
void AddRegression(double *data, int vSize, int n, int step, int offset,int delwin, int head, int tail, int simpleDiffs);
/* EXPORT->AddHeadRegress: add regression at start of data */
void AddHeadRegress(double *data, int vSize, int n, int step, int offset,int delwin, int simpleDiffs);
/* EXPORT->AddTailRegress: add regression at end of data */
void AddTailRegress(double *data, int vSize, int n, int step, int offset,int delwin, int simpleDiffs);
/* EXPORT->NormaliseLogEnergy: normalise log energy to range -X .. 1.0 */
void NormaliseLogEnergy(double *data, int n, int step, double silFloor, double escale);
/* Regression: add regression vector at +offset from source vector.  If head
or tail is less than delwin then duplicate first/last vector to compensate */
void Regress(double *data, int vSize, int n, int step, int offset, int delwin, int head, int tail, int simpleDiffs);

/*自己的模块*/
/*使用FZroMean后方差归一*/
void FNormalize(double *data, int vSize, int n, int step);
/*EXPORT->NormaliseLogEnergy: normalise log energy to range 0 .. 1.0*/
void NormaliseLogEnergy2(double *data, int n, int step );


/*-----------------计算谱中心以及子带能量---------------*/
/*自己的模块*/
/*计算谱中心*/
void calBrightness(Vector fftx, double *b,double te);
/*计算子带能量*/
void calSubBankE(Vector fftx, Vector subBankEnergy,double te);

#ifdef __cplusplus
}
#endif // __cplusplus

#endif