#pragma once
#ifndef _MFCC_H_H
#define _MFCC_H_H


#include"hmath.hpp"
#include"sigProcess.hpp"
#include"WAVE.hpp"
#include"fileIO.hpp"
#include"cnpy.hpp"
#include<cmath>

#define maxBuffLength 600


namespace hMFCC {
	/*MFCC提取的配置
	sampleRate 采样率
	samplePeriod 采样周期，再checkConfig中初始化
	hipassfre 高频截断
	lowpassfre 低频截断
	preemphasise 预加重系数，=0不预加重，=0.97一般设置对于语音
	zeroMeanSigFlag 对音频信号零均值
	wlen 窗长度
	inc 位移量
	fbankFlag 是否仅计算fbank
	bankNum mel滤波器组的个数
	MFCCNum MFCC参量的个数
	delwin 计算加速系数时的窗口大小
	energyFlag 是否计算能量
	zeroCrossingFlag 是否计算过零率
	brightFlag 是否计算谱中心
	subBankEFlag 是否计算子带能量以及个数
	regreOrder 加速系数的阶数，1为静态参量，2为静态参量和一阶加速系数，3为静态参量和一阶加速系数和二阶加速系数，以此类推
	znormFlag 对本段音频Z标准化
	MFCC0thFlag 是否包含DCT的0阶MFCC，当且仅当fbankFlag=0时有效
	fftLength 输出原始的fft系数，仅用作调试
	saveType 保存格式，=0 数字格式，=1 科学计数，=2 numpy格式，=3 二进制文件
	vecNum 输出声道数，=1 单声道，=2 双声道，=4 四声道
	channels 音频的声道
	sampleNum 采样点数，可以不初始化
	numThreads 线程数
	fileList 文件列表位置
	*/
	typedef struct {
		int sampleRate;
		double samplePeriod;
		double hipassfre;
		double lowpassfre;
		double preemphasise;
		int zeroMeanSigFlag;
		int wlen;
		int inc;
		int fbankFlag;
		int bankNum;
		int MFCCNum;
		int delwin;
		int energyFlag;
		int zeroCrossingFlag;
		int brightFlag;
		int subBandEFlag;
		int regreOrder;
		int znormFlag;
		int MFCC0thFlag;
		int fftLength;
		int saveType;
		int vecNum;
		int channels;
		int sampleNum;
		int numThreads;
		char fileList[maxBuffLength];
	}Config;

	/*检查配置*/
	void configCheck(Config* config);

	/*读取配置，使用inih*/
	static int handler(void* user, const char* section, const char* name, const char* value) {
		Config* pconfig = (Config*)user;
#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
		if (MATCH("Frame", "sampleRate")) { pconfig->sampleRate = atoi(value); }
		else if (MATCH("Frame", "hipassfre")) { pconfig->hipassfre = atof(value); }
		else if (MATCH("Frame", "lowpassfre")) { pconfig->lowpassfre = atof(value); }
		else if (MATCH("Frame", "preemphasise")) { pconfig->preemphasise = atof(value); }
		else if (MATCH("Frame", "zeroMeanSigFlag")) { pconfig->zeroMeanSigFlag = atoi(value); }
		else if (MATCH("Frame", "wlen")) { pconfig->wlen = atoi(value); }
		else if (MATCH("Frame", "inc")) { pconfig->inc = atoi(value); }
		else if (MATCH("Frame", "vecNum")) { pconfig->vecNum = atoi(value); }
		else if (MATCH("MFCC", "fbankFlag")) { pconfig->fbankFlag = atoi(value); }
		else if (MATCH("MFCC", "bankNum")) { pconfig->bankNum = atoi(value); }
		else if (MATCH("MFCC", "MFCCNum")) { pconfig->MFCCNum = atoi(value); }
		else if (MATCH("MFCC", "MFCC0thFlag")) { pconfig->MFCC0thFlag = atoi(value); }
		else if (MATCH("Others", "energyFlag")) { pconfig->energyFlag = atoi(value); }
		else if (MATCH("Others", "zeroCrossingFlag")) { pconfig->zeroCrossingFlag = atoi(value); }
		else if (MATCH("Others", "brightFlag")) { pconfig->brightFlag = atoi(value); }
		else if (MATCH("Others", "subBandEFlag")) { pconfig->subBandEFlag = atoi(value); }
		else if (MATCH("Others", "fftLength")) { pconfig->fftLength = atoi(value); }
		else if (MATCH("Regression", "znormFlag")) { pconfig->znormFlag = atoi(value); }
		else if (MATCH("Regression", "regreOrder")) { pconfig->regreOrder = atoi(value); }
		else if (MATCH("Regression", "delwin")) { pconfig->delwin = atoi(value); }
		else if (MATCH("IO", "saveType")) {
			if (strcmp(value, "f\0") == 0)pconfig->saveType = 0;
			else if (strcmp(value, "n\0") == 0) pconfig->saveType = 2;
			else if (strcmp(value, "b\0") == 0) pconfig->saveType = 3;
			else pconfig->saveType = 1;
		}
		else if (MATCH("IO", "fileList")) { strcpy(pconfig->fileList, value); }
		else if (MATCH("IO", "numThreads")) { pconfig->numThreads = atoi(value); }
		else return 0;
		return 1;
	}

	/*MFCC滤波器的信息，依靠此初始化
	仿照HTK的设置
	*/
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

	/*为多线程准备的暂存变量*/
	typedef struct {
		FBankInfo* fbVec;
		hmath::Vector hamWin;
		hmath::Matrix subBankEnergy;
		hmath::Matrix fbank;
		hmath::Matrix d2;
		hmath::Matrix d3;
	}MFCCWapperTempStruct;

	/*Mel-Realft函数参考HTK的源码，有部分修改*/
	/*
	return mel-frequency corresponding to given FFT index k.
	Resolution is normally determined by fres field of FBankInfo
	record.
	*/
	double Mel(int k, double fres);

	double WarpFreq(double fcl, double fcu, double freq, double minFreq, double maxFreq, double alpha);

	/*
	Initialise an FBankInfo record prior to calling Wave2FBank.
	*/
	FBankInfo InitFBank(int frameSize, double sampPeriod, int numChans,
		double lopass, double hipass, int usePower, int takeLogs,
		int doubleFFT,
		double alpha, double warpLowCut, double warpUpCut);

	/*
	Convert given speech frame in s into mel-frequency filterbank
	coefficients.  The total frame energy is stored in te.  The
	info record contains precomputed filter weights and should be set
	prior to using Wave2FBank by calling InitFBank.
	*/
	void Wave2FBank(hmath::Vector s, hmath::Vector fbank, double *te, double* te2, FBankInfo info);

	/*
	Apply the DCT to fbank and store first n cepstral coeff in c.
	Note that the resulting coef are normalised by sqrt(2/numChans)
	*/
	void FBank2MFCC(hmath::Vector fbank, hmath::Vector c, int n);

	/* EXPORT->FBank2C0: return zero'th cepstral coefficient */
	double FBank2C0(hmath::Vector fbank);

	/*
	When called s holds nn complex values stored in the
	sequence   [ r1 , i1 , r2 , i2 , .. .. , rn , in ] where
	n = VectorSize(s) DIV 2, n must be a power of 2. On exit s
	holds the fft (or the inverse fft if invert == 1)
	*/
	void FFT(hmath::Vector s, int invert);

	/*
	When called s holds 2*n real values, on exit s holds the
	first  n complex points of the spectrum stored in
	the same format as for fft
	*/
	void Realft(hmath::Vector s);

	/*MFCC暂存变量的初始化*/
	MFCCWapperTempStruct MFCCWapperTempInit(Config config);

	/*释放暂存变量*/
	void MFCCWapperTempFree(MFCCWapperTempStruct* mwts, Config config);

	/*MFCC提取的Wapper，
	输入音频位置，输出位置，暂存变量，配置，线程ID，wav检查变量
	1.读取WAV文件,并转化为十进制
	2.计算MFCC系数以及其他参量
	3.能量归一化
	4.计算加速系数
	5.MFCC参量归一归整
	6.写入目标文件
	最后参量的排列分别为 ( MFCCNum+energyFlag+zeroCrossingFlag+brightFlag+subBankEFlag )*regreOrder
	最后参量的排列分别为 ( MFCC+能量+过零率+谱中心+子带能量 )*阶数
	*/
	void MFCCWapper(const char* inputWAV, const char* outputFile,
		MFCCWapperTempStruct mwts, Config config, int threadNum,
		const hWAVE::WAVEParams_t * const wParamStd = NULL);
}
#endif