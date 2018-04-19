#pragma once
#ifndef _WAVE_H_
#define _WAVE_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifdef _MSC_VER  
typedef __int32 int32_t; 
typedef unsigned __int32 uint32_t; 
typedef __int64 int64_t; 
typedef unsigned __int64 uint64_t;  
typedef unsigned __int16 uint16_t;
#else 
#include <stdint.h> 
#endif


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include "hmath.h"

/*RIFF chunk*/
typedef struct {
	char ckID[4];
	uint32_t cksize;
	char WAVEID[4];
}RIFF_t;

/*fmt chunk*/
typedef struct {
	char ckID[4];
	uint32_t cksize;
	uint16_t wFormatTag;
	uint16_t nChannels;
	uint32_t nSamplesPerSec;
	uint32_t nAvgBytesPerSec;
	uint16_t nBlockAlign;
	uint16_t wBitsPerSample;
}fmt_t;

/*data chunk�����е�data��һ���ŵ���*����������Int�;���*/
typedef struct {
	char ckID[4];
	uint32_t cksize;
	int** data;
}DATA_t;

/*WAVEParams_t�Ǹ���������Ƶ���洢����Ƶ����Ҫ����*/
/*��Щ�������Ǹ���WAVE��ͷֱ�Ӷ�ȡ�����߼�Ӽ����*/
/*��Щ�����У�����������������������Ƶ�ʣ�ÿ����������ֽڴ�С��ÿ���������ֽڴ�С��һ���ǰ����ͬ��*/
typedef struct {
	int numChannels;
	int numSamples;
	int sampleRate;
	int sampleLengthInByte;
	int containerLengthInByte;
}WAVEParams_t;

/*WAVE����Ҫ�ṹ*/
typedef struct {
	RIFF_t RIFF;
	fmt_t fmt;
	DATA_t DATA;
	WAVEParams_t WAVEParams;
	int myEndian;
}WAVE_t;

/*��ʼ��WAVE�ṹ*/
WAVE_t initWAVE_t();

/*���ڶ�ȡWAVEͷ��һЩ�Ǹ�������Ϣ*/
long int readWAVE(FILE * f, size_t sizeInByte, int EndianFlag);

/*����WAVE�ļ���WAVE�ṹ*/
void loadWAVEFile(WAVE_t* w, FILE* f);

/*�ж�WAVE�ļ��Ƿ����Ҫ��*/
int WAVEParamsCheck(WAVEParams_t w1, WAVEParams_t w2);

/*��ӡWAVE�Ļ�����Ϣ*/
void print_WAVE(WAVE_t w);

/*�ͷ�WAVE*/
void free_WAVE(WAVE_t* w);

/*д��WAVE�ļ�*/
void writeWaveFile(FILE* f, WAVEParams_t params, IntMat m);

#ifdef _cplusplus
}
#endif // _cplusplus

#endif // !_WAVE_H_
