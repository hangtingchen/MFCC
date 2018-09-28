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
	/*MFCC��ȡ������
	sampleRate ������
	samplePeriod �������ڣ���checkConfig�г�ʼ��
	hipassfre ��Ƶ�ض�
	lowpassfre ��Ƶ�ض�
	preemphasise Ԥ����ϵ����=0��Ԥ���أ�=0.97һ�����ö�������
	zeroMeanSigFlag ����Ƶ�ź����ֵ
	wlen ������
	inc λ����
	fbankFlag �Ƿ������fbank
	bankNum mel�˲�����ĸ���
	MFCCNum MFCC�����ĸ���
	delwin �������ϵ��ʱ�Ĵ��ڴ�С
	energyFlag �Ƿ��������
	zeroCrossingFlag �Ƿ���������
	brightFlag �Ƿ����������
	subBankEFlag �Ƿ�����Ӵ������Լ�����
	regreOrder ����ϵ���Ľ�����1Ϊ��̬������2Ϊ��̬������һ�׼���ϵ����3Ϊ��̬������һ�׼���ϵ���Ͷ��׼���ϵ�����Դ�����
	znormFlag �Ա�����ƵZ��׼��
	MFCC0thFlag �Ƿ����DCT��0��MFCC�����ҽ���fbankFlag=0ʱ��Ч
	fftLength ���ԭʼ��fftϵ��������������
	saveType �����ʽ��=0 ���ָ�ʽ��=1 ��ѧ������=2 numpy��ʽ��=3 �������ļ�
	vecNum �����������=1 ��������=2 ˫������=4 ������
	channels ��Ƶ������
	sampleNum �������������Բ���ʼ��
	numThreads �߳���
	fileList �ļ��б�λ��
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

	/*�������*/
	void configCheck(Config* config);

	/*��ȡ���ã�ʹ��inih*/
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

	/*MFCC�˲�������Ϣ�������˳�ʼ��
	����HTK������
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

	/*Ϊ���߳�׼�����ݴ����*/
	typedef struct {
		FBankInfo* fbVec;
		hmath::Vector hamWin;
		hmath::Matrix subBankEnergy;
		hmath::Matrix fbank;
		hmath::Matrix d2;
		hmath::Matrix d3;
	}MFCCWapperTempStruct;

	/*Mel-Realft�����ο�HTK��Դ�룬�в����޸�*/
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

	/*MFCC�ݴ�����ĳ�ʼ��*/
	MFCCWapperTempStruct MFCCWapperTempInit(Config config);

	/*�ͷ��ݴ����*/
	void MFCCWapperTempFree(MFCCWapperTempStruct* mwts, Config config);

	/*MFCC��ȡ��Wapper��
	������Ƶλ�ã����λ�ã��ݴ���������ã��߳�ID��wav������
	1.��ȡWAV�ļ�,��ת��Ϊʮ����
	2.����MFCCϵ���Լ���������
	3.������һ��
	4.�������ϵ��
	5.MFCC������һ����
	6.д��Ŀ���ļ�
	�����������зֱ�Ϊ ( MFCCNum+energyFlag+zeroCrossingFlag+brightFlag+subBankEFlag )*regreOrder
	�����������зֱ�Ϊ ( MFCC+����+������+������+�Ӵ����� )*����
	*/
	void MFCCWapper(const char* inputWAV, const char* outputFile,
		MFCCWapperTempStruct mwts, Config config, int threadNum,
		const hWAVE::WAVEParams_t * const wParamStd = NULL);
}
#endif