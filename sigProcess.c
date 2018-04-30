#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "hmath.h"
#include "sigProcess.h"

using namespace hmath;

void ZeroMean(short *data, long nSamples)
{
	long i, hiClip = 0, loClip = 0;
	short *x;
	double sum = 0.0, y, mean;

	x = data;
	for (i = 0; i<nSamples; i++, x++)
		sum += *x;
	mean = sum / (double)nSamples;
	x = data;
	for (i = 0; i<nSamples; i++, x++) {
		y = (double)(*x) - mean;
		if (y<-32767.0) {
			y = -32767.0; ++loClip;
		}
		if (y>32767.0) {
			y = 32767.0; ++hiClip;
		}
		*x = (short)((y>0.0) ? y + 0.5 : y - 0.5);
	}
	if (loClip>0)
		printf("ZeroMean: %d samples too -ve\n", loClip);
	if (hiClip>0)
		printf("ZeroMean: %d samples too +ve\n", hiClip);
}

double zeroCrossingRate(Vector s, int frameSize) {
	int count = 0;int i;
	for (i = 1; i < frameSize; i++)  if ((s[i] * s[i + 1]) < 0.0)count++; 
	return ((double)count)/(double)(frameSize-1) ;
}

/* EXPORT->PreEmphasise: pre-emphasise signal in s */
void PreEmphasise(Vector s, double k)
{
	int i;
	double preE;

	preE = k;
	if (k == 0.0)return;
	for (i = VectorSize(s); i >= 2; i--)
		s[i] -= s[i - 1] * preE;
	s[1] *= 1.0 - preE;
}

//将PCM文件转化为十进制
Vector pcmToData(FILE* f, long fileLength) {
	//short int正好是两个字节
	short* data = (short *)malloc(sizeof(short)*fileLength / 2);
	int i = 0;
	Vector pcmVector = CreateVector(fileLength / 2);
	for (i = 0; i<fileLength / 2; i++) {
		//每次读取两个字节,存储到data中
		fread(&data[i], 2, 1, f);
	}
	ZeroMean(data, fileLength / 2);
	for (i = 0; i<fileLength / 2; i++)	pcmVector[i + 1] = data[i];
	free(data);
	PreEmphasise(pcmVector, 0);//不做预加重
	fclose(f);
	return pcmVector;
}

Vector* pcmToData2(FILE* f, long fileLength,int biteNum,int perSample,int vecNum) {
	//short int正好是两个字节
	short* data = (short *)malloc(sizeof(short)*fileLength / biteNum);
	int i = 0;
	Vector pcmVector1 = CreateVector(fileLength / perSample); Vector pcmVector2 = CreateVector(fileLength / perSample);
	Vector pcmVector3 = CreateVector(fileLength / perSample); Vector pcmVector4 = CreateVector(fileLength / perSample);
	Vector* pcmVector = (Vector*)malloc(sizeof(Vector) * 4);
	pcmVector[0] = pcmVector3; pcmVector[1] = pcmVector1; pcmVector[2] = pcmVector2; pcmVector[3] = pcmVector4;
	for (i = 0; i<fileLength / biteNum; i++) {
		//每次读取两个字节,存储到data中
		fread(&data[i], 2, 1, f);
	}
	ZeroMean(data, fileLength / biteNum);
	for (i = 0; i < fileLength / perSample; i++) { pcmVector1[i + 1] = data[2 * i]; pcmVector2[i + 1] = data[2 * i + 1]; }
	for (i = 1; i <= fileLength / perSample; i++) { pcmVector3[i] = (double)(int)(0.5*pcmVector1[i] + 0.5*pcmVector2[i] + 0.5); pcmVector4[i] = pcmVector1[i] - pcmVector2[i]; }
	free(data);
	fclose(f); 
	for(i=0;i<4;i++)PreEmphasise(pcmVector[i], 0);
	if (vecNum == 1||vecNum==4)return pcmVector;
	else if (vecNum == 2)return &(pcmVector[1]);
	else return NULL;
}


static int hamWinSize = 0;          /* Size of current Hamming window */
static Vector hamWin = NULL;        /* Current Hamming window */

/* GenHamWindow: generate precomputed Hamming window function */
void GenHamWindow(int frameSize)
{
	int i;
	double a;

	if (hamWin == NULL || VectorSize(hamWin) < frameSize)
		hamWin = CreateVector(frameSize);
	a = 2 * pi / (frameSize - 1);
	for (i = 1; i <= frameSize; i++)
		hamWin[i] = 0.54 - 0.46 * cos(a*(i - 1));
	hamWinSize = frameSize;
}

/* EXPORT->Ham: Apply Hamming Window to Speech frame s */
void Ham(Vector s)
{
	int i, frameSize;
	frameSize = VectorSize(s);
	if (hamWinSize != frameSize)
		GenHamWindow(frameSize);
	for (i = 1; i <= frameSize; i++) {
		s[i] *= hamWin[i];
		//		printf("%d %f\n", i,s[i]);
	}
}


static int cepWinSize = 0;            /* Size of current cepstral weight window */
static int cepWinL = 0;               /* Current liftering coeff */
static Vector cepWin = NULL;        /* Current cepstral weight window */

									/* GenCepWin: generate a new cep liftering vector */
void GenCepWin(int cepLiftering, int count)
{
	int i;
	double a, Lby2;

	if (cepWin == NULL || VectorSize(cepWin) < count)
		cepWin = CreateVector( count);
	a = pi / cepLiftering;
	Lby2 = cepLiftering / 2.0;
	for (i = 1; i <= count; i++)
		cepWin[i] = 1.0 + Lby2*sin(i * a);
	cepWinL = cepLiftering;
	cepWinSize = count;
}

/* EXPORT->WeightCepstrum: Apply cepstral weighting to c */
void WeightCepstrum(Vector c, int start, int count, int cepLiftering)
{
	int i, j;

	if (cepWinL != cepLiftering || count > cepWinSize)
		GenCepWin(cepLiftering, count);
	j = start;
	for (i = 1; i <= count; i++)
		c[j++] *= cepWin[i];
}

/* EXPORT->UnWeightCepstrum: Undo cepstral weighting of c */
void UnWeightCepstrum(Vector c, int start, int count, int cepLiftering)
{
	int i, j;

	if (cepWinL != cepLiftering || count > cepWinSize)
		GenCepWin(cepLiftering, count);
	j = start;
	for (i = 1; i <= count; i++)
		c[j++] /= cepWin[i];
}

/* The following operations apply to a sequence of n vectors step apart.
They are used to operate on the 'columns' of data files
containing a sequence of feature vectors packed together to form a
continguous block of doubles.  The logical size of each vector is
vSize (<=step) */

/* EXPORT->FZeroMean: Zero mean the given data sequence */
void FZeroMean(double *data, int vSize, int n, int step)
{
	double sum;
	double *fp, mean;
	int i, j;

	for (i = 0; i<vSize; i++) {
		/* find mean over i'th component */
		sum = 0.0;
		fp = data + i;
		for (j = 0; j<n; j++) {
			sum += *fp; fp += step;
		}
		mean = sum / (double)n;
		/* subtract mean from i'th components */
		fp = data + i;
		for (j = 0; j<n; j++) {
			*fp -= mean; fp += step;
		}
	}
}

void FNormalize(double * data, int vSize, int n, int step)
{
	double sum;
	double *fp, sd;
	int i,j;
	for (i = 0; i < vSize; i++) {
		sum = 0.0;
		fp = data + i;
		for (j = 0; j < n; j++) { sum += (*fp)*(*fp); fp += step; }
		sd = sqrt(sum / (double)n);
		fp = data + i;
		for (j = 0; j < n; j++) {
			*fp /= sd; fp += step;
		}
	}
}

/* Regression: add regression vector at +offset from source vector.  If head
or tail is less than delwin then duplicate first/last vector to compensate */
void Regress(double *data, int vSize, int n, int step, int offset,int delwin, int head, int tail, int simpleDiffs)
{
	double *fp, *fp1, *fp2, *back, *forw;
	double sum, sigmaT2;
	int i, t, j;

	sigmaT2 = 0.0;
	for (t = 1; t <= delwin; t++)
		sigmaT2 += t*t;
	sigmaT2 *= 2.0;
	fp = data;
	for (i = 1; i <= n; i++) {
		fp1 = fp; fp2 = fp + offset;
		for (j = 1; j <= vSize; j++) {
			back = forw = fp1; sum = 0.0;
			for (t = 1; t <= delwin; t++) {
				if (head + i - t > 0)     back -= step;
				if (tail + n - i + 1 - t > 0) forw += step;
				if (!simpleDiffs) sum += t * (*forw - *back);
			}
			if (simpleDiffs)
				*fp2 = (*forw - *back) / (2 * delwin);
			else
				*fp2 = sum / sigmaT2;
			++fp1; ++fp2;
		}
		fp += step;
	}
}

void calBrightness(Vector fftx, double * b,double te)
{
	int i;
	double sum = 0.0;
	const double energyfloor = 1.0;if (te <= energyfloor)te = energyfloor;
	if (((int)VectorSize(fftx)) % 2 != 0)printf("something wrong in cal brightness");
	for (i = 1; i <=((int)VectorSize(fftx)) / 2; i++) {
		sum += (fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i])*(double)i;
	}
	*b = sum / te;
	*b = (*b) / ((double)VectorSize(fftx) / 2.0);
}


void calSubBankE(Vector fftx, Vector subBankEnergy,double te)
{
	const double energyfloor = 1.0;
	int i;
	int numBank = VectorSize(subBankEnergy); int bankSize = (int)VectorSize(fftx) / (2*numBank);
	int bankNum = 1;
	double sum = 0.0; 
	for (i = 1; i <= (int)VectorSize(fftx) / 2; i++) {
		if (i <= bankNum*bankSize) {
			sum += fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i];
		}
		else {
			subBankEnergy[bankNum] = sum/te;
			//printf("sum: %f\n", sum/te);
			bankNum++; sum = 0.0; i--;
		}
	}
	if (te < energyfloor || sum < energyfloor)subBankEnergy[bankNum] = 0.0;
	else subBankEnergy[bankNum] = sum / te;

}


/* EXPORT->AddRegression: add regression vector at +offset from source vector */
void AddRegression(double *data, int vSize, int n, int step, int offset,int delwin, int head, int tail, int simpleDiffs)
{
	Regress(data, vSize, n, step, offset, delwin, head, tail, simpleDiffs);
}

/* EXPORT->AddHeadRegress: add regression at start of data */
void AddHeadRegress(double *data, int vSize, int n, int step, int offset,int delwin, int simpleDiffs)
{
	double *fp, *fp1, *fp2;
	int i, j;

	fp = data;
	if (delwin == 0) {
		for (i = 1; i <= n; i++) {
			fp1 = fp; fp2 = fp + offset;
			for (j = 1; j <= vSize; j++) {
				*fp2 = *(fp1 + step) - *fp1;
				++fp1; ++fp2;
			}
			fp += step;
		}
	}
	else {
		Regress(data, vSize, n, step, offset, delwin, 0, delwin, simpleDiffs);
	}
}

/* EXPORT->AddTailRegress: add regression at end of data */
void AddTailRegress(double *data, int vSize, int n, int step, int offset,int delwin, int simpleDiffs)
{
	double *fp, *fp1, *fp2;
	int i, j;

	fp = data;
	if (delwin == 0) {
		for (i = 1; i <= n; i++) {
			fp1 = fp; fp2 = fp + offset;
			for (j = 1; j <= vSize; j++) {
				*fp2 = *fp1 - *(fp1 - step);
				++fp1; ++fp2;
			}
			fp += step;
		}
	}
	else {
		Regress(data, vSize, n, step, offset, delwin, delwin, 0, simpleDiffs);
	}
}

/* EXPORT->NormaliseLogEnergy: normalise log energy to range -X .. 1.0 */
void NormaliseLogEnergy(double *data, int n, int step, double silFloor, double escale)
{
	double *p, max, min;
	int i;

	/* find max log energy */
	p = data; max = *p;
	for (i = 1; i<n; i++) {
		p += step;                   /* step p to next e val */
		if (*p > max) max = *p;
	}
	min = max - (silFloor*log(10.0)) / 10.0;  /* set the silence floor */
											  /* normalise */
	p = data;
	for (i = 0; i<n; i++) {
		if (*p < min) *p = min;          /* clamp to silence floor */
		*p = 1.0 - (max - *p) * escale;  /* normalise */
		p += step;
	}
}

void NormaliseLogEnergy2(double * data, int n, int step)
{
	double *p, max;
	int i;
	p = data; max = *p;
	for (i = 1; i<n; i++) {
		p += step;                   /* step p to next e val */
		if (*p > max) max = *p;
	}
	p = data;
	for (i = 0; i < n; i++) {
		*p /= max;
		p += step;
	}
}
