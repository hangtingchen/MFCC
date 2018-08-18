#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<string>
#include"ini.h"
#include "hmath.h"
#include "mfcc.h"
#include"sigProcess.h"
#include"WAVE.h"
#include"cnpy.hpp"

#define maxBuffLength 600

using namespace hmath;

/*����
pi
wlen ������
inc λ����
bankNum mel�˲�����ĸ���
MFCCNum MFCC�����ĸ���
delwin �������ϵ��ʱ�Ĵ��ڴ�С
energyFlag �Ƿ��������
zeroCrossingFlag �Ƿ���������
brightFlag �Ƿ����������
subBankEFlag �Ƿ�����Ӵ������Լ�����
regreOrder ����ϵ���Ľ�����1Ϊ��̬������2Ϊ��̬������һ�׼���ϵ����
3Ϊ��̬������һ�׼���ϵ���Ͷ��׼���ϵ��
*/

/*  �����������зֱ�Ϊ ( MFCCNum+energyFlag+zeroCrossingFlag+brightFlag+subBankEFlag )*regreOrder  */
/*  �����������зֱ�Ϊ ( MFCC+����+������+������+�Ӵ����� )*����  */

/*��������ṹ
1.��ȡwav�ļ�,��ת��Ϊʮ����
2.����MFCCϵ���Լ���������
3.������һ��
4.�������ϵ��
5.MFCC������һ����
6.д��Ŀ���ļ�
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
	FILE *f ;FILE* f_result;//�ֱ���ԭ�ļ�,������ʮ�����ļ�
	FILE *fList;
	char* pcmFile1;char* pcmFile2;char fileNameBuf[maxBuffLength];
	Vector d2 , d3 ,fbank,subBankEnergy=NULL;
	Vector* d1;
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
	otherFeatureNum = MFCC0thFlag + energyFlag + zeroCrossingFlag + brightFlag + subBandEFlag + fftLength;//�������������ĸ���

	subBankEnergy = CreateVector(subBandEFlag);
	d2 = CreateVector(wlen);
	fbank = CreateVector(bankNum);
	d3 = CreateVector(MFCCNum);

	info = InitFBank(wlen, samplePeriod, bankNum, lowpassfre, hipassfre, 1, 1, 0, 1.0, 0, 0);
	//��ʼ��������wlen���������=(10E7/16000),bankNum,���Ƶ��,���Ƶ��,�Ƿ����������
	//�Ƿ��log���Ƿ�doublefft���Ƿ�resacle(1.0������)
	//����ham��
	GenHamWindow(wlen);

	while(fgets(fileNameBuf,maxBuffLength,fList)!=NULL){
		curr_pos = 1;
		if(fileNameBuf[0]=='\n')break;
		pcmFile1=strtok(fileNameBuf,"\t");
		pcmFile2=strtok(NULL,"\n");
		/*1.��ȡwav�ļ�, ��ת��Ϊʮ����*/
		//��ԭ�ļ�
		printf("Convert %s to %s\n", pcmFile1, pcmFile2);
		f = fopen(pcmFile1, "rb");
		if (!f) { printf("open .wav failed\n"); system("pause");  return -1; }

		printf("including : \nMFCCNum	%d\nenergyFlag %d\nzeroCrossingFlag %d\nbrightFlag %d\nsubBandEFlag %d\n",
			MFCCNum, energyFlag, zeroCrossingFlag, brightFlag, subBandEFlag);
		printf("the frame feature dimension is %d\n",(MFCCNum + otherFeatureNum) 
			* config.vecNum* regreOrder);
		printf("order %d\n", regreOrder);
		printf("start...\n");

		//��ȡwave�ļ�
		WAVE_t wavfile = initWAVE_t();
		loadWAVEFile(&wavfile, f); print_WAVE(wavfile);
		if (config.sampleRate != wavfile.WAVEParams.sampleRate) {
			printf("Sample rate not right");
			system("pause");
			return 1;
		}
		config.sampleNum=wavfile.WAVEParams.numSamples; 
		config.channels = wavfile.WAVEParams.numChannels;
		if(config.channels==2)d1 = CreateMatrix(wavfile.WAVEParams.numChannels*2, 
			wavfile.WAVEParams.numSamples);
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

		dpostProc = CreateVector(((MFCCNum + otherFeatureNum) * config.vecNum* 
			regreOrder *(int)((config.sampleNum - (wlen - inc)) / inc)));
		printf("total coef size: %d\n", VectorSize(dpostProc));

		/*2.����MFCCϵ���Լ���������*/
		//d2,d3�ֱ�Ϊ����mel�˲�������źţ�MFCC����
		//���յ�����
		for (j = 0; j <= (config.sampleNum - wlen); j += inc) {
			if (config.vecNum > NumRows(d1)) {
				printf("vecnum %d > num_channels %d ",config.vecNum,config.channels);
				system("pause");
				return 1;
			}
			for (i0 = 1; i0 <= config.vecNum; i0++) {
				if (config.vecNum == 1 && config.channels == 2)i0 = 3;
				for (i = 1; i <= wlen; i++)d2[i] = d1[i0][i + j];
//				if (j == 0)ShowVector(d2);
				//���������
				zc = zeroCrossingRate(d2, wlen);
				//�Ӵ�
				Ham(d2);
				//te��te2�ֱ��Ǹ����źź�fft����źż��������
				//����mel�˲�����
				Wave2FBank(d2, fbank, &te, &te2, info);
//				if (j == 0)ShowVector(fbank);
				//���������ĺ��Ӵ����������ǰٷ���
				brightness = 0.0; 
				calBrightness(info.x, &brightness, te2); 
				if (subBandEFlag)calSubBankE(info.x, subBankEnergy, te2);
				//����MFCCϵ��
				if (config.fbankFlag) CopyVector(fbank, d3);
				else FBank2MFCC(fbank, d3, MFCCNum);;
				//Ǩ�����ݣ����Ҹ���flag���ϲ�������
				for (j0 = 1; j0 <= MFCCNum; j0++, curr_pos++)dpostProc[curr_pos] = d3[j0];
				if (MFCC0thFlag) { dpostProc[curr_pos] = FBank2C0(fbank); curr_pos++; }
				if (energyFlag) { dpostProc[curr_pos] = log(te); curr_pos++; }
				if (zeroCrossingFlag) { dpostProc[curr_pos] = zc; curr_pos++; } //printf("%f\n", zc);
				if (brightFlag) { dpostProc[curr_pos] = brightness; curr_pos++; };
				if (subBandEFlag) { for (j0 = 1; j0 <= subBandEFlag; j0++, curr_pos++)
					dpostProc[curr_pos] = subBankEnergy[j0]; }
				if (fftLength) { for (j0 = 1; j0 <= fftLength; j0++, curr_pos++)
					dpostProc[curr_pos] = sqrt(info.x[2 * j0 - 1] * info.x[2 * j0 - 1] + info.x[2 * j0] * info.x[2 * j0]); }
			}
			curr_pos += (MFCCNum + otherFeatureNum) * (regreOrder - 1)* config.vecNum ;
		}
		
		printf("post-processing...\n");
		//������ٲ���
		vSize = (MFCCNum + otherFeatureNum)*config.vecNum; rowNum = (config.sampleNum - (wlen - inc)) / inc; step = (MFCCNum + otherFeatureNum) * regreOrder*config.vecNum;
	//	NormaliseLogEnergy(&dpostProc[1 + MFCCNum], rowNum, step, 50.0, 0.1);
		/*3.������һ��*/
	//	NormaliseLogEnergy2(&dpostProc[1 + MFCCNum], rowNum, step);
		/*4.�������ϵ��*/
		if(regreOrder>=2)Regress(&dpostProc[1], vSize , rowNum, step, vSize, delwin, 0, 0,0);
		if(regreOrder>=3)Regress(&dpostProc[1+vSize], vSize, rowNum, step, vSize, delwin, 0, 0, 0);
		/*5.MFCC������һ����*/
		//MFCCϵ����ֵ������������������������˲���
	//	FZeroMean(&dpostProc[1], MFCCNum, rowNum, step); FNormalize(&dpostProc[1], MFCCNum, rowNum, step);
	//	FZeroMean(&dpostProc[1], step, rowNum, step);FNormalize(&dpostProc[1], step, rowNum, step);
	//	if (regreOrder >= 2){ FZeroMean(&dpostProc[1+vSize], MFCCNum, rowNum, step); FNormalize(&dpostProc[1+vSize], MFCCNum, rowNum, step); }
	//	if (regreOrder >= 3) { FZeroMean(&dpostProc[1+2*vSize], MFCCNum, rowNum, step); FNormalize(&dpostProc[1+2*vSize], MFCCNum, rowNum, step); }
		
	//	����
	//	for (int i = 1; i <= VectorSize(dpostProc); i++)dpostProc[i] *= 10.0;

		/*6.д��Ŀ���ļ�*/

		if(config.saveType==0){
			f_result = fopen(pcmFile2, "w");
			if (!f_result) { printf("open result.dat failed\n"); system("pause");  return -1; }
			printf("writing the doc...\n");
			for(i=1;i<=VectorSize(dpostProc);i++){
				fprintf(f_result, "%f\t", dpostProc[i]);
				if (i % step == 0)fprintf(f_result, "\n");
			}
			fclose(f_result);
		}
		else if(config.saveType==1){
			f_result = fopen(pcmFile2, "w");
			if (!f_result) { printf("open result.dat failed\n"); system("pause");  return -1; }
			printf("writing the doc...\n");
			for(i=1;i<=VectorSize(dpostProc);i++){
				fprintf(f_result, "%e\t", dpostProc[i]);
				if (i % step == 0)fprintf(f_result, "\n");
			}
			fclose(f_result);
		}
		else {
			cnpy::npy_save(pcmFile2, &dpostProc[1], { (size_t)rowNum ,(size_t)step },"w");
		}

		FreeVector(dpostProc);
		FreeMatrix(d1);
	}

	FreeVector(d2);
	FreeVector(fbank);
	FreeVector(d3);
	FreeVector(subBankEnergy);
	fclose(fList);
	//system("pause");
	return 0;
}
