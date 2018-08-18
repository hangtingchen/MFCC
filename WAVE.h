#pragma once
#ifndef _WAVE_H_
#define _WAVE_H_

#ifdef __cplusplus
extern "C" {
#endif // _cplusplus

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

/*data chunk，其中的data是一个信道数*采样点数的Int型矩阵*/
typedef struct {
	char ckID[4];
	uint32_t cksize;
	int** data;
}DATA_t;

/*WAVEParams_t是根据输入音频，存储的音频的主要参数*/
/*这些参数都是根据WAVE的头直接读取，或者间接计算的*/
/*这些参数有：声道数，采样点数，采样频率，每个采样点的字节大小，每个容器的字节大小（一般和前者相同）*/
typedef struct {
	int numChannels;
	int numSamples;
	int sampleRate;
	int sampleLengthInByte;
	int containerLengthInByte;
}WAVEParams_t;

/*WAVE的主要结构*/
typedef struct {
	RIFF_t RIFF;
	fmt_t fmt;
	DATA_t DATA;
	WAVEParams_t WAVEParams;
	int myEndian;
}WAVE_t;

/*初始化WAVE结构*/
WAVE_t initWAVE_t();

/*用于读取WAVE头中一些非负整型信息*/
long int readWAVE(FILE * f, size_t sizeInByte, int EndianFlag);

/*读入WAVE文件至WAVE结构*/
void loadWAVEFile(WAVE_t* w, FILE* f);

/*判断WAVE文件是否符合要求*/
int WAVEParamsCheck(WAVEParams_t w1, WAVEParams_t w2);

/*打印WAVE的基本信息*/
void print_WAVE(WAVE_t w);

/*释放WAVE*/
void free_WAVE(WAVE_t* w);

/*写入WAVE文件*/
void writeWaveFile(FILE* f, WAVEParams_t params, IntMat m);

#ifdef __cplusplus
}
#endif // _cplusplus

#endif // !_WAVE_H_
