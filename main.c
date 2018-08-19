#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"ini.h"
#include "hmath.h"
#include "mfcc.h"
#include"sigProcess.h"
#include"WAVE.h"


#define maxBuffLength 600

/*含义
pi
wlen 窗长度
inc 位移量
bankNum mel滤波器组的个数
MFCCNum MFCC参量的个数
delwin 计算加速系数时的窗口大小
energyFlag 是否计算能量
zeroCrossingFlag 是否计算过零率
brightFlag 是否计算谱中心
subBankEFlag 是否计算子带能量以及个数
regreOrder 加速系数的阶数，1为静态参量，2为静态参量和一阶加速系数，3为静态参量和一阶加速系数和二阶加速系数
*/

/*  最后参量的排列分别为 ( MFCCNum+energyFlag+zeroCrossingFlag+brightFlag+subBankEFlag )*regreOrder  */
/*  最后参量的排列分别为 ( MFCC+能量+过零率+谱中心+子带能量 )*阶数  */

/*程序基本结构
1.读取pcm文件,并转化为十进制
2.计算MFCC系数以及其他参量
3.能量归一化
4.计算加速系数
5.MFCC参量归一归整
6.写入目标文件
*/

typedef struct {
	int sampleRate;
	double samplePeriod;
	double hipassfre;
	double lowpassfre;
	double preemphasise;
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
	int MFCC0thFlag;
	int fftLength;
	int saveType;
	int vecNum;
	int channels;
	int sampleNum;
	char fileList[maxBuffLength];
}Config;

static int handler(void* user, const char* section, const char* name, const char* value) {
	Config* pconfig = (Config*)user;
	#define MATCH(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
	if (MATCH("Frame", "sampleRate")) { pconfig->sampleRate = atoi(value); }
	else if (MATCH("Frame", "hipassfre")) { pconfig->hipassfre = atof(value);   }
	else if (MATCH("Frame", "lowpassfre")) { pconfig->lowpassfre = atof(value); }
	else if (MATCH("Frame", "preemphasise")) { pconfig->preemphasise = atof(value); }
	else if (MATCH("Frame", "wlen")) { pconfig->wlen = atoi(value);}
	else if (MATCH("Frame", "inc")) { pconfig->inc = atoi(value); }
	else if (MATCH("Frame", "saveType")) { 
		if(strcmp(value,"f\0")==0)pconfig->saveType=0 ;
		else if(strcmp(value, "n\0") == 0) pconfig->saveType=2;
		else pconfig->saveType = 1;
	}
	else if (MATCH("Frame", "vecNum")) { pconfig->vecNum = atoi(value); }
	else if (MATCH("Frame", "fileList")) { strcpy(pconfig->fileList, value);}
	else if (MATCH("MFCC", "fbankFlag")) { pconfig->fbankFlag = atoi(value); }
	else if (MATCH("MFCC", "bankNum")) { pconfig->bankNum = atoi(value); }
	else if (MATCH("MFCC", "MFCCNum")) { pconfig->MFCCNum = atoi(value);  }
	else if (MATCH("MFCC", "MFCC0thFlag")) { pconfig->MFCC0thFlag = atoi(value);}
	else if (MATCH("Others", "energyFlag")) { pconfig->energyFlag = atoi(value); }
	else if (MATCH("Others", "zeroCrossingFlag")) { pconfig->zeroCrossingFlag = atoi(value); }
	else if (MATCH("Others", "brightFlag")) { pconfig->brightFlag = atoi(value); }
	else if (MATCH("Others", "subBandEFlag")) { pconfig->subBandEFlag = atoi(value); }
	else if (MATCH("Others", "fftLength")) { pconfig->fftLength = atoi(value);  }
	else if (MATCH("Regression", "regreOrder")) { pconfig->regreOrder = atoi(value); }
	else if (MATCH("Regression", "delwin")) { pconfig->delwin = atoi(value); }
	else return 0;
	return 1;
}

int main(int argc, char** argv) {

	Config config;
	double samplePeriod;
	double hipassfre;
	double lowpassfre;
	int wlen;
	int inc;
	int bankNum;
	int MFCCNum;
	int delwin;
	int energyFlag;
	int zeroCrossingFlag;
	int brightFlag;
	int subBandEFlag;
	int regreOrder;
	int MFCC0thFlag;
	int fftLength;
	FILE *f ;FILE* f_result;//分别是原文件,处理后的十进制文件
	FILE *fList;
	char* pcmFile1;char* pcmFile2;char fileNameBuf[maxBuffLength];
	Vector d2=NULL , d3=NULL ,fbank=NULL,subBankEnergy=NULL,hamWin=NULL;
	Vector* d1=NULL;
	FBankInfo info;
	Vector dpostProc;
	int otherFeatureNum;
	int curr_pos = 1;
	int vSize,rowNum,step;
	double zc,te,te2,brightness;
	int i,j,i0,j0;

	i0 = ini_parse(argv[1], handler, &config);
	if (i0<0) {
		printf("Can't load '.ini'\n");
		system("pause");
		return 1;
	}
	else if (i0 > 0) {
		printf("Unknown input config settings in line %d",i0);
		system("pause");
		return 1;
	}
	else printf("finish reading the config.ini\n");
	config.samplePeriod = (double)1e7 / (double)config.sampleRate;
	if (config.fbankFlag) {
		config.MFCCNum = config.bankNum; config.MFCC0thFlag = 0;
	}

	samplePeriod = config.samplePeriod;
	hipassfre = config.hipassfre; lowpassfre = config.lowpassfre;
	wlen = config.wlen;
	inc = config.inc;
	bankNum = config.bankNum;
	MFCCNum = config.MFCCNum;
	delwin = config.delwin;
	energyFlag = config.energyFlag;
	zeroCrossingFlag = config.zeroCrossingFlag;
	brightFlag = config.brightFlag;
	subBandEFlag = config.subBandEFlag;
	regreOrder = config.regreOrder;
	MFCC0thFlag = config.MFCC0thFlag;
	fftLength = config.fftLength;
	fList=fopen(config.fileList,"r");
	otherFeatureNum = MFCC0thFlag + energyFlag + zeroCrossingFlag + brightFlag + subBandEFlag + fftLength;//计算其他参量的个数

	subBankEnergy = CreateVector(subBandEFlag);
	d2 = CreateVector(wlen);
	fbank = CreateVector(bankNum);
	d3 = CreateVector(MFCCNum);

	info = InitFBank(wlen, samplePeriod, bankNum, lowpassfre, hipassfre, 1, 1, 0, 1.0, 0, 0);
	//初始化，其中wlen，采样间隔=(10E7/16000),bankNum,最低频率,最高频率,是否计算能量，是否加log，是否doublefft，是否resacle(1.0不进行)
	//产生ham窗
	hamWin=GenHamWindow(wlen);

	while(fgets(fileNameBuf,maxBuffLength,fList)!=NULL){
		curr_pos = 1;
		if(fileNameBuf[0]=='\n')break;
		pcmFile1=strtok(fileNameBuf,"\t");
		pcmFile2=strtok(NULL,"\n");
		/*1.读取wav文件, 并转化为十进制*/
		//打开原文件
		printf("Convert %s to %s\n", pcmFile1, pcmFile2);
		f = fopen(pcmFile1, "rb");
		if (!f) { printf("open .wav failed\n"); system("pause");  return -1; }

		printf("including : \nMFCCNum	%d\nenergyFlag %d\nzeroCrossingFlag %d\nbrightFlag %d\nsubBandEFlag %d\n", MFCCNum, energyFlag, zeroCrossingFlag, brightFlag, subBandEFlag);
		printf("the frame feature dimension is %d\n",(MFCCNum + otherFeatureNum) * config.vecNum* regreOrder);
		printf("order %d\n", regreOrder);
		printf("start...\n");

		//读取wave文件
		WAVE_t wavfile = initWAVE_t();
		loadWAVEFile(&wavfile, f); print_WAVE(wavfile);
		if (config.sampleRate != wavfile.WAVEParams.sampleRate) {
			printf("Sample rate not right");
			system("pause");
			return 1;
		}
		config.sampleNum=wavfile.WAVEParams.numSamples; config.channels = wavfile.WAVEParams.numChannels;
		if(config.channels==2)d1 = CreateMatrix(wavfile.WAVEParams.numChannels*2, wavfile.WAVEParams.numSamples);
		else d1 = CreateMatrix(wavfile.WAVEParams.numChannels, wavfile.WAVEParams.numSamples);
		for (i0 = 1; i0 <= wavfile.WAVEParams.numChannels; i0++)
			for (i = 1; i <= wavfile.WAVEParams.numSamples; i++) {
				d1[i0][i] = (double)wavfile.DATA.data[i0][i];
			}
		if (config.channels == 2)for (i = 1; i <= wavfile.WAVEParams.numSamples; i++) {
			d1[3][i] = 0.5*(d1[1][i] + d1[2][i]);
			d1[4][i] = d1[1][i] - d1[2][i];
		}
		for(i0 = 1; i0 <= NumRows(d1); i0++) PreEmphasise(d1[i0],config.preemphasise);
		fclose(f); free_WAVE(&wavfile);

		dpostProc = CreateVector(((MFCCNum + otherFeatureNum) * config.vecNum* regreOrder *(int)((config.sampleNum - (wlen - inc)) / inc)));
		printf("total coef size: %d\n", VectorSize(dpostProc));

		/*2.计算MFCC系数以及其他参量*/
		//d2,d3分别为经过mel滤波器组的信号，MFCC参数
		//最终的数据
		for (j = 0; j <= (config.sampleNum - wlen); j += inc) {
			if (config.vecNum > NumRows(d1)) {
				printf("vecnum %d > num_channels %d ",config.vecNum,config.channels);
				system("pause");
				return 1;
			}
			for (i0 = 1; i0 <= config.vecNum; i0++) {
				for (i = 1; i <= wlen; i++)d2[i] = d1[i0][i + j];
//				if (j == 0)ShowVector(d2);
				//计算过零率
				zc = zeroCrossingRate(d2, wlen);
				//加窗
				Ham(d2,hamWin,wlen);
				//te和te2分别是根据信号和fft后的信号计算的能量
				//经过mel滤波器组
				Wave2FBank(d2, fbank, &te, &te2, info);
//				if (j == 0)ShowVector(fbank);
				//计算谱中心和子带能量，都是百分数
				brightness = calBrightness(info.x); 
				if (subBandEFlag)calSubBankE(info.x, subBankEnergy);
				//计算MFCC系数
				if (config.fbankFlag) CopyVector(fbank, d3);
				else FBank2MFCC(fbank, d3, MFCCNum);;
				//迁移数据，并且根据flag加上部分特征
				for (j0 = 1; j0 <= MFCCNum; j0++, curr_pos++)dpostProc[curr_pos] = d3[j0];
				if (MFCC0thFlag) { dpostProc[curr_pos] = FBank2C0(fbank); curr_pos++; }
				if (energyFlag) { dpostProc[curr_pos] = log(te); curr_pos++; }
				if (zeroCrossingFlag) { dpostProc[curr_pos] = zc; curr_pos++; } //printf("%f\n", zc);
				if (brightFlag) { dpostProc[curr_pos] = brightness; curr_pos++; };
				if (subBandEFlag) { for (j0 = 1; j0 <= subBandEFlag; j0++, curr_pos++)dpostProc[curr_pos] = subBankEnergy[j0]; }
				if (fftLength) { for (j0 = 1; j0 <= fftLength; j0++, curr_pos++)dpostProc[curr_pos] = sqrt(info.x[2 * j0 - 1] * info.x[2 * j0 - 1] + info.x[2 * j0] * info.x[2 * j0]); }
			}
			curr_pos += (MFCCNum + otherFeatureNum) * (regreOrder - 1)* config.vecNum ;
		}
		
		printf("post-processing...\n");
		//计算加速参量
		vSize = (MFCCNum + otherFeatureNum)*config.vecNum; rowNum = (config.sampleNum - (wlen - inc)) / inc; step = (MFCCNum + otherFeatureNum) * regreOrder*config.vecNum;
	//	NormaliseLogEnergy(&dpostProc[1 + MFCCNum], rowNum, step, 50.0, 0.1);
		/*3.能量归一化*/
	//	NormaliseLogEnergy2(&dpostProc[1 + MFCCNum], rowNum, step);
		/*4.计算加速系数*/
		if(regreOrder>=2)Regress(&dpostProc[1], vSize , rowNum, step, vSize, delwin, 0, 0,0);
		if(regreOrder>=3)Regress(&dpostProc[1+vSize], vSize, rowNum, step, vSize, delwin, 0, 0, 0);
		/*5.MFCC参量归一归整*/
		//MFCC系数均值方差归整，不对其他参数做此操作
	//	FZeroMean(&dpostProc[1], MFCCNum, rowNum, step); FNormalize(&dpostProc[1], MFCCNum, rowNum, step);
	//	FZeroMean(&dpostProc[1], step, rowNum, step);FNormalize(&dpostProc[1], step, rowNum, step);
	//	if (regreOrder >= 2){ FZeroMean(&dpostProc[1+vSize], MFCCNum, rowNum, step); FNormalize(&dpostProc[1+vSize], MFCCNum, rowNum, step); }
	//	if (regreOrder >= 3) { FZeroMean(&dpostProc[1+2*vSize], MFCCNum, rowNum, step); FNormalize(&dpostProc[1+2*vSize], MFCCNum, rowNum, step); }
		
	//	扩大
	//	for (int i = 1; i <= VectorSize(dpostProc); i++)dpostProc[i] *= 10.0;

		/*6.写入目标文件*/
		f_result= fopen(pcmFile2, "w");
		if (!f_result) { printf("open result.dat failed\n"); system("pause");  return -1; }
		printf("writing the doc...\n");
		if(config.saveType==0){
			for(i=1;i<=VectorSize(dpostProc);i++){
				fprintf(f_result, "%f\t", dpostProc[i]);
				if (i % step == 0)fprintf(f_result, "\n");
			}
		}
		else if(config.saveType==1){
			for(i=1;i<=VectorSize(dpostProc);i++){
				fprintf(f_result, "%e\t", dpostProc[i]);
				if (i % step == 0)fprintf(f_result, "\n");
			}
		}
		fclose(f_result);
		FreeVector(dpostProc);
		FreeMatrix(d1);
	}

	FreeVector(d2);
	FreeVector(fbank);
	FreeVector(d3);
	FreeVector(subBankEnergy);
	FreeVector(hamWin);
	fclose(fList);
	//system("pause");
	return 0;
}
