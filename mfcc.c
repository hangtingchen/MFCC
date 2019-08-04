/* -------------------- MFCC Related Operations -------------------- */

/* EXPORT->Mel: return mel-frequency corresponding to given FFT index */
#include"mfcc.h"
#include"hmath.hpp"
#include"sigProcess.hpp"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

using namespace hmath;

void hMFCC::configCheck(hMFCC::Config* config){
	if (!(config->sampleRate == 44100 || config->sampleRate == 8000
		|| config->sampleRate == 16000 || config->sampleRate == 48000)) {
		printf("Unusual sample rate %d\n", config->sampleRate);
	}
	config->samplePeriod = (double)1e7 / (double)config->sampleRate;
	if (config->hipassfre > config->sampleRate ||
		config->hipassfre < 100 || config->lowpassfre < 0 ||
		config->hipassfre < config->lowpassfre) {
		printf("Unexpected high bound and low bound %f %f\n", config->lowpassfre, config->hipassfre);
		exit(1);
	}
	if (config->preemphasise < 0) {
		printf("Unexpected preemphasise %f\n", config->preemphasise);
		exit(1);
	}
	if (config->fbankFlag) {
		config->MFCCNum = config->bankNum;
		config->MFCC0thFlag = 0;
	}
	if (config->fftLength) {
		printf("You want to extract fft coefft. Don't recommand to do this\n");
	}
	if (!(config->vecNum == 1 || config->vecNum == 2 || config->vecNum == 4)) {
		printf("Unexpected vecNum %d\n", config->vecNum);
		exit(1);
	}
	if (config->sampleRate <= 0 || config->wlen <= 0 || config->inc <= 0 || config->fftLength < 0 || config->numThreads <= 0) {
		printf("Check your sampleRate wlen inc fftLength numThreads\n");
		exit(1);
	}
	if (config->regreOrder <=0 ) {
		printf("RegreOrder should be >=1 and 1 means no diff coeffs.");
		exit(1);
	}
}

double hMFCC::Mel(int k, double fres)
{
	return 1127 * log(1 + (k - 1)*fres);
}

/* EXPORT->WarpFreq: return warped frequency */
double hMFCC::WarpFreq(double fcl, double fcu, double freq, double minFreq, double maxFreq, double alpha)
{
	if (alpha == 1.0)
		return freq;
	else {
		double scale = 1.0 / alpha;
		double cu = fcu * 2 / (1 + scale);
		double cl = fcl * 2 / (1 + scale);

		double au = (maxFreq - cu * scale) / (maxFreq - cu);
		double al = (cl * scale - minFreq) / (cl - minFreq);

		if (freq > cu)
			return  au * (freq - cu) + scale * cu;
		else if (freq < cl)
			return al * (freq - minFreq) + minFreq;
		else
			return scale * freq;
	}
}

/* EXPORT->InitFBank: Initialise an FBankInfo record */
hMFCC::FBankInfo hMFCC::InitFBank(int frameSize, double sampPeriod, int numChans,
	double lopass, double hipass, int usePower, int takeLogs,
	int doubleFFT,
	double alpha, double warpLowCut, double warpUpCut)
{
	hMFCC::FBankInfo fb;
	double mlo, mhi, ms, melk;
	int k, chan, maxChan, Nby2;

	/* Save sizes to cross-check subsequent usage */
	fb.frameSize = frameSize; fb.numChans = numChans;
	fb.sampPeriod = sampPeriod;
	fb.usePower = usePower; fb.takeLogs = takeLogs;
	/* Calculate required FFT size */
	fb.fftN = 2;
	while (frameSize>fb.fftN) fb.fftN *= 2;
	if (doubleFFT)
		fb.fftN *= 2;
	Nby2 = fb.fftN / 2;
	fb.fres = 1.0E7 / (sampPeriod * fb.fftN * 700.0);//700源于mel scale
	maxChan = numChans + 1;
	/* set lo and hi pass cut offs if any */
	fb.klo = 2; fb.khi = Nby2;       /* apply lo/hi pass filtering */
	mlo = 0; mhi = hMFCC::Mel(Nby2 + 1, fb.fres);
	if (lopass >= 0.0) {
		mlo = 1127 * log(1 + lopass / 700.0);
		fb.klo = (int)((lopass * sampPeriod * 1.0e-7 * fb.fftN) + 2.5);
		if (fb.klo<2) fb.klo = 2;
	}
	if (hipass >= 0.0) {
		mhi = 1127 * log(1 + hipass / 700.0);
		//printf("%f\n", 1127 * log(1 + hipass / 700.0));
		fb.khi = (int)((hipass * sampPeriod * 1.0e-7 * fb.fftN) + 0.5);
		if (fb.khi>Nby2) fb.khi = Nby2;
	}
		printf("FFT passband %d to %d out of 1 to %d\n", fb.klo, fb.khi, Nby2);
		printf("Mel passband %f to %f\n", mlo, mhi);
	/* Create vector of fbank centre frequencies */
	fb.cf = CreateVector(maxChan);
	ms = mhi - mlo;
	for (chan = 1; chan <= maxChan; chan++) {
		if (alpha == 1.0) {
			fb.cf[chan] = ((double)chan / (double)maxChan)*ms + mlo;
		}
		else {
			/* scale assuming scaling starts at lopass */
			double minFreq = 700.0 * (exp(mlo / 1127.0) - 1.0);
			double maxFreq = 700.0 * (exp(mhi / 1127.0) - 1.0);
			double cf = ((double)chan / (double)maxChan) * ms + mlo;

			cf = 700 * (exp(cf / 1127.0) - 1.0);

			fb.cf[chan] = 1127.0 * log(1.0 + hMFCC::WarpFreq(warpLowCut, warpUpCut, cf, minFreq, maxFreq, alpha) / 700.0);
		}
	}

	/* Create loChan map, loChan[fftindex] -> lower channel index */
	fb.loChan = CreateIntVec(Nby2);
	for (k = 1, chan = 1; k <= Nby2; k++) {
		melk = hMFCC::Mel(k, fb.fres);
		if (k<fb.klo || k>fb.khi) fb.loChan[k] = -1;
		else {
			while (fb.cf[chan] < melk  && chan <= maxChan) ++chan;
			fb.loChan[k] = chan - 1;
		}
	}

	/* Create vector of lower channel weights */
	fb.loWt = CreateVector(Nby2);
	for (k = 1; k <= Nby2; k++) {
		chan = fb.loChan[k];
		if (k<fb.klo || k>fb.khi) fb.loWt[k] = 0.0;
		else {
			if (chan>0)
				fb.loWt[k] = ((fb.cf[chan + 1] - hMFCC::Mel(k, fb.fres)) /
				(fb.cf[chan + 1] - fb.cf[chan]));
			else
				fb.loWt[k] = (fb.cf[1] - hMFCC::Mel(k, fb.fres)) / (fb.cf[1] - mlo);
		}
	}
	/* Create workspace for fft */
	fb.x = CreateVector(fb.fftN);
	return fb;
}

/* EXPORT->Wave2FBank:  Perform filterbank analysis on speech s */
void hMFCC::Wave2FBank(Vector s, Vector fbank, double *te,double* te2, hMFCC::FBankInfo info)
{
	const double melfloor = 1.0;
	const double energyfloor = 1.0;
	int k, bin;
	double t1, t2;   /* real and imag parts */
	double ek;      /* energy of k'th fft channel */

				   /* Check that info record is compatible */
	if (info.frameSize != VectorSize(s))
		printf("Wave2FBank: frame size mismatch");
	if (info.numChans != VectorSize(fbank))
		printf("Wave2FBank: num channels mismatch");
	/* Compute frame energy if needed */
	if (te != NULL) {
		*te = 0.0;
		for (k = 1; k <= info.frameSize; k++)
			*te += (s[k] * s[k]);
	}
	if (*te < energyfloor)*te = energyfloor;
	/* Apply FFT */
	for (k = 1; k <= info.frameSize; k++)
		info.x[k] = s[k];    /* copy to workspace */
	for (k = info.frameSize + 1; k <= info.fftN; k++)
		info.x[k] = 0.0;   /* pad with zeroes */
	Realft(info.x);                            /* take fft */

											   /* Fill filterbank channels */
	ZeroVector(fbank); 
	for (k = info.klo; k <= info.khi; k++) {             /* fill bins */
		t1 = info.x[2 * k - 1]; t2 = info.x[2 * k];
		if (info.usePower)
			ek = t1*t1 + t2*t2;
		else
			ek = sqrt(t1*t1 + t2*t2);

		bin = info.loChan[k];
		t1 = info.loWt[k] * ek;
		if (bin>0) fbank[bin] += t1;
		if (bin<info.numChans) fbank[bin + 1] += ek - t1;
	}
	*te2 = 0.0; //printf("te2 %f\n", *te2);
	for (k = 1; k <= info.fftN/2; k++) {
		*te2 += info.x[2 * k - 1] * info.x[2 * k - 1] + info.x[2 * k] * info.x[2 * k];
//		printf("te2 %d %f\n",k, *te2);
	}
	if (*te2 < energyfloor)*te2 = energyfloor;

	/* Take logs */
	if (info.takeLogs)
		for (bin = 1; bin <= info.numChans; bin++) {
			t1 = fbank[bin];
			if (t1<melfloor) t1 = melfloor;
			fbank[bin] = log(t1);
		}
}

/* EXPORT->FBank2MFCC: compute first n cepstral coeff */
void hMFCC::FBank2MFCC(Vector fbank, Vector c, int n)
{
	int j, k, numChan;
	double mfnorm, pi_factor, x;

	numChan = VectorSize(fbank);
	mfnorm = sqrt(2.0 / (double)numChan);
	pi_factor = pi / (double)numChan;
	for (j = 1; j <= n; j++) {
		c[j] = 0.0; x = (double)j * pi_factor;
		for (k = 1; k <= numChan; k++)
			//c[j] += log(fbank[k]) * cos(x*(k - 0.5));
			c[j] += fbank[k] * cos(x*(k - 0.5));
		c[j] *= mfnorm;
	}
}

double hMFCC::FBank2C0(Vector fbank)
{
	int k, numChan;
	double mfnorm, sum;

	numChan = VectorSize(fbank);
	mfnorm = sqrt(2.0 / (double)numChan);
	sum = 0.0;
	for (k = 1; k <= numChan; k++)
		sum += fbank[k];
	return sum * mfnorm;
}

/* EXPORT-> FFT: apply fft/invfft to complex s */
void hMFCC::FFT(Vector s, int invert)
{
	int ii, jj, n, nn, limit, m, j, inc, i;
	double wx, wr, wpr, wpi, wi, theta;
	double xre, xri, x;

	n = VectorSize(s);
	nn = n / 2; j = 1;
	for (ii = 1; ii <= nn; ii++) {
		i = 2 * ii - 1;
		if (j>i) {
			xre = s[j]; xri = s[j + 1];
			s[j] = s[i];  s[j + 1] = s[i + 1];
			s[i] = xre; s[i + 1] = xri;
		}
		m = n / 2;
		while (m >= 2 && j > m) {
			j -= m; m /= 2;
		}
		j += m;
	};
	limit = 2;
	while (limit < n) {
		inc = 2 * limit; theta = 2*pi / limit;
		if (invert) theta = -theta;
		x = sin(0.5 * theta);
		wpr = -2.0 * x * x; wpi = sin(theta);
		wr = 1.0; wi = 0.0;
		for (ii = 1; ii <= limit / 2; ii++) {
			m = 2 * ii - 1;
			for (jj = 0; jj <= (n - m) / inc; jj++) {
				i = m + jj * inc;
				j = i + limit;
				xre = wr * s[j] - wi * s[j + 1];
				xri = wr * s[j + 1] + wi * s[j];
				s[j] = s[i] - xre; s[j + 1] = s[i + 1] - xri;
				s[i] = s[i] + xre; s[i + 1] = s[i + 1] + xri;
			}
			wx = wr;
			wr = wr * wpr - wi * wpi + wr;
			wi = wi * wpr + wx * wpi + wi;
		}
		limit = inc;
	}
	if (invert)
		for (i = 1; i <= n; i++)
			s[i] = s[i] / nn;

}

/* EXPORT-> Realft: apply fft to real s */
void hMFCC::Realft(Vector s)
{
	int n, n2, i, i1, i2, i3, i4;
	double xr1, xi1, xr2, xi2, wrs, wis;
	double yr, yi, yr2, yi2, yr0, theta, x;

	n = VectorSize(s) / 2; n2 = n / 2;
	theta = pi / n;
	hMFCC::FFT(s, 0);
	x = sin(0.5 * theta);
	yr2 = -2.0 * x * x;
	yi2 = sin(theta); yr = 1.0 + yr2; yi = yi2;
	for (i = 2; i <= n2; i++) {
		i1 = i + i - 1;      i2 = i1 + 1;
		i3 = n + n + 3 - i2; i4 = i3 + 1;
		wrs = yr; wis = yi;
		xr1 = (s[i1] + s[i3]) / 2.0; xi1 = (s[i2] - s[i4]) / 2.0;
		xr2 = (s[i2] + s[i4]) / 2.0; xi2 = (s[i3] - s[i1]) / 2.0;
		s[i1] = xr1 + wrs * xr2 - wis * xi2;
		s[i2] = xi1 + wrs * xi2 + wis * xr2;
		s[i3] = xr1 - wrs * xr2 + wis * xi2;
		s[i4] = -xi1 + wrs * xi2 + wis * xr2;
		yr0 = yr;
		yr = yr * yr2 - yi  * yi2 + yr;
		yi = yi * yr2 + yr0 * yi2 + yi;
	}
	xr1 = s[1];
	s[1] = xr1 + s[2];
	s[2] = 0.0;
}

hMFCC::MFCCWapperTempStruct hMFCC::MFCCWapperTempInit(hMFCC::Config config) {
	hMFCC::FBankInfo* fbVec = (hMFCC::FBankInfo*)malloc(sizeof(hMFCC::FBankInfo)*config.numThreads);
	hmath::Vector hamWin = hsigProcess::GenHamWindow(config.wlen);
	hmath::Matrix subBankEnergy = NULL;
	hmath::Matrix fbank = CreateMatrix(config.numThreads, config.bankNum);
	hmath::Matrix d2 = CreateMatrix(config.numThreads, config.wlen);
	hmath::Matrix d3 = CreateMatrix(config.numThreads, config.MFCCNum);
	hMFCC::MFCCWapperTempStruct mfccWapperTempStruct = { NULL,NULL,NULL,NULL,NULL,NULL };
	if (config.subBandEFlag > 0)subBankEnergy = CreateMatrix(config.numThreads, config.subBandEFlag);

	for (int i = 0; i < config.numThreads; i++) {
		fbVec[i] = InitFBank(config.wlen, config.samplePeriod,
			config.bankNum, config.lowpassfre, config.hipassfre,
			1, 1, 0, 1.0, 0, 0);
	}
	mfccWapperTempStruct = { fbVec ,hamWin ,subBankEnergy ,fbank,d2,d3 };
	return mfccWapperTempStruct;
}


void hMFCC::MFCCWapperTempFree(hMFCC::MFCCWapperTempStruct* mwts, hMFCC::Config config)
{
	FreeMatrix(mwts->d2); mwts->d2 = NULL;
	FreeMatrix(mwts->d3); mwts->d3 = NULL;
	FreeMatrix(mwts->fbank); mwts->fbank = NULL;
	FreeVector(mwts->hamWin); mwts->hamWin = NULL;
	for (int i = 0; i < config.numThreads; i++) {
		FreeVector(mwts->fbVec[i].cf); mwts->fbVec[i].cf = NULL;
		FreeIntVec(mwts->fbVec[i].loChan); mwts->fbVec[i].loChan = NULL;
		FreeVector(mwts->fbVec[i].loWt); mwts->fbVec[i].loWt = NULL;
		FreeVector(mwts->fbVec[i].x); mwts->fbVec[i].x = NULL;
	}
	if (config.subBandEFlag) {
		FreeMatrix(mwts->subBankEnergy);
		mwts->subBankEnergy = NULL;
	}
}

void hMFCC::MFCCWapper(const char* inputWAV, const char* outputFile, hMFCC::MFCCWapperTempStruct mwts,
	hMFCC::Config config,int threadNum, const hWAVE::WAVEParams_t * const wParamStd) {
	FILE* fin = NULL,*fout=NULL;
	int otherFeatureNum = config.MFCC0thFlag + config.energyFlag +
		config.zeroCrossingFlag + config.brightFlag + config.subBandEFlag + config.fftLength;
	Vector* d1; Vector dpostProc = NULL;
	double te = 0.0, te2 = 0.0, brightness = 0.0;
	hfileIO::BinaryFileStruct bf = { -1,-1,-1,-1,nullptr };
	int curr_pos = 1;
	printf("Convert %s to %s\n", inputWAV, outputFile);
	fin = fopen(inputWAV, "rb");
	if (!fin) { printf("open .wav failed\n"); system("pause");  exit(1); }
	printf("including : \nMFCCNum	%d\nenergyFlag %d\nzeroCrossingFlag %d\nbrightFlag %d\nsubBandEFlag %d\n",
		config.MFCCNum, config.energyFlag,config.zeroCrossingFlag,config.brightFlag, config.subBandEFlag);
	printf("the frame feature dimension is %d\n", (config.MFCCNum + otherFeatureNum) * config.vecNum* config.regreOrder);
	hWAVE::WAVE_t wavfile = hWAVE::initWAVE_t();
	hWAVE::loadWAVEFile(&wavfile, fin); 
	hWAVE::print_WAVE(wavfile);
	if (wParamStd)
		if (!hWAVE::WAVEParamsCheck(wavfile.WAVEParams, *wParamStd)) {
			printf("input file has an unexpected format %s\n", inputWAV);
		}

	if(config.zeroMeanSigFlag || wavfile.WAVEParams.containerLengthInByte==1)
		for (int i = 1; i <= wavfile.WAVEParams.numChannels; i++) {
			hsigProcess::ZeroMean(wavfile.DATA.data[i]);
		}
	if (config.sampleRate != wavfile.WAVEParams.sampleRate) {
		printf("Sample rate not right");
		system("pause");
		exit(1);
	}
	config.sampleNum = wavfile.WAVEParams.numSamples; 
	config.channels = wavfile.WAVEParams.numChannels;

	int vSize = (config.MFCCNum + otherFeatureNum)*config.vecNum;
	int rowNum = (config.sampleNum - (config.wlen - config.inc)) / config.inc +
		((config.sampleNum - (config.wlen - config.inc)) / config.inc>0);
	int step = (config.MFCCNum + otherFeatureNum) * config.regreOrder*config.vecNum;
	if (config.channels == 2)d1 = CreateMatrix(wavfile.WAVEParams.numChannels * 2, wavfile.WAVEParams.numSamples);
	else d1 = CreateMatrix(wavfile.WAVEParams.numChannels, wavfile.WAVEParams.numSamples);
	for (int i0 = 1; i0 <= wavfile.WAVEParams.numChannels; i0++)
		for (int i = 1; i <= wavfile.WAVEParams.numSamples; i++) {
			d1[i0][i] = (double)wavfile.DATA.data[i0][i];
		}
	if (config.channels == 2)for (int i = 1; i <= wavfile.WAVEParams.numSamples; i++) {
		d1[3][i] = 0.5*(d1[1][i] + d1[2][i]);
		d1[4][i] = d1[1][i] - d1[2][i];
	}
	for (int i0 = 1; i0 <= NumRows(d1); i0++) hsigProcess::PreEmphasise(d1[i0], config.preemphasise);
	fclose(fin); hWAVE::free_WAVE(&wavfile);

	dpostProc = CreateVector(((config.MFCCNum + otherFeatureNum) * config.vecNum* config.regreOrder *rowNum));
	printf("The feature if of dimension %d * (%d * %d * %d)\n",rowNum,config.vecNum,config.regreOrder,(config.MFCCNum + otherFeatureNum));
	printf("total coef size: %d\n", VectorSize(dpostProc));

	/*2.计算MFCC系数以及其他参量*/
	//d2,d3分别为经过mel滤波器组的信号，MFCC参数
	//最终的数据
	for (int j = 0; j < rowNum; j ++) {
		if (config.vecNum > NumRows(d1)) {
			printf("vecnum %d > num_channels %d ", config.vecNum, config.channels);
			system("pause");
			exit(1);
		}
		for (int i0 = 1; i0 <= config.vecNum; i0++) {
			ZeroVector(mwts.d2[threadNum]);
			if (config.vecNum == 1 && config.channels == 2)i0 = 3;
			for (int i = 1; i <= config.wlen && j*config.inc+i<= config.sampleNum ; i++)
				mwts.d2[threadNum][i] = d1[i0][i + j*config.inc];
			//				if (j == 0)ShowVector(d2);
			//计算过零率
			double zc = hsigProcess::zeroCrossingRate(mwts.d2[threadNum], config.wlen);
			//加窗
			hsigProcess::Ham(mwts.d2[threadNum],mwts.hamWin,config. wlen);
			//te和te2分别是根据信号和fft后的信号计算的能量
			//经过mel滤波器组
			Wave2FBank(mwts.d2[threadNum], mwts.fbank[threadNum], &te, &te2, mwts.fbVec[threadNum-1]);
			//				if (j == 0)ShowVector(fbank);
			//计算谱中心和子带能量，都是百分数
			brightness = hsigProcess::calBrightness(mwts.fbVec[threadNum-1].x);
			if (config.subBandEFlag)hsigProcess::calSubBankE(mwts.fbVec[threadNum-1].x, mwts.subBankEnergy[threadNum]);
			//计算MFCC系数
			if (config.fbankFlag) CopyVector(mwts.fbank[threadNum],mwts.d3[threadNum]);
			else hMFCC::FBank2MFCC(mwts.fbank[threadNum], mwts.d3[threadNum], config.MFCCNum);
			//迁移数据，并且根据flag加上部分特征
			for (int j0 = 1; j0 <= config.MFCCNum; j0++, curr_pos++)dpostProc[curr_pos] = mwts.d3[threadNum][j0];
			if (config.MFCC0thFlag) { dpostProc[curr_pos] = hMFCC::FBank2C0(mwts.fbank[threadNum]); curr_pos++; }
			if (config.energyFlag) { dpostProc[curr_pos] = log(te); curr_pos++; }
			if (config.zeroCrossingFlag) { dpostProc[curr_pos] = zc; curr_pos++; } //printf("%f\n", zc);
			if (config.brightFlag) { dpostProc[curr_pos] = brightness; curr_pos++; };
			if (config.subBandEFlag) { for (int j0 = 1; j0 <= config.subBandEFlag; j0++, curr_pos++)dpostProc[curr_pos] = mwts.subBankEnergy[threadNum][j0]; }
			if (config.fftLength) { for (int j0 = 1; j0 <= config.fftLength; j0++, curr_pos++)dpostProc[curr_pos] = sqrt(mwts.fbVec[threadNum-1].x[2 * j0 - 1] * mwts.fbVec[threadNum-1].x[2 * j0 - 1] + mwts.fbVec[threadNum-1].x[2 * j0] * mwts.fbVec[threadNum-1].x[2 * j0]); }
		}
		curr_pos += (config.MFCCNum + otherFeatureNum) * (config.regreOrder - 1)* config.vecNum;
	}
	printf("post-processing...\n");
	//计算加速参量
	//	NormaliseLogEnergy(&dpostProc[1 + MFCCNum], rowNum, step, 50.0, 0.1);
	/*3.能量归一化*/
	//	NormaliseLogEnergy2(&dpostProc[1 + MFCCNum], rowNum, step);
	/*4.计算加速系数*/
	for (int j = 1; j < config.regreOrder; j++) {
		hsigProcess::Regress(&dpostProc[1+(j-1)*vSize], vSize, rowNum, step, vSize, config.delwin, 0, 0, 0);
	}
	/*5.MFCC参量归一归整*/
	//MFCC系数均值方差归整，不对其他参数做此操作
	if (config.znormFlag) {
		hsigProcess::ZNormalize(&dpostProc[1], step, rowNum, step);
	}
	//	扩大
	//	for (int i = 1; i <= VectorSize(dpostProc); i++)dpostProc[i] *= 10.0;
	/*6.写入目标文件*/
	if (config.saveType == 0) {
		fout = fopen(outputFile, "w");
		if (!fout) { printf("open result.dat failed\n"); system("pause");  exit(1); }
		printf("writing the doc...\n");
		for (int i = 1; i <= VectorSize(dpostProc); i++) {
			fprintf(fout, "%f\t", dpostProc[i]);
			if (i % step == 0)fprintf(fout, "\n");
		}
		fclose(fout);
	}
	else if (config.saveType == 1) {
		fout = fopen(outputFile, "w");
		if (!fout) { printf("open result.dat failed\n"); system("pause");  exit(1); }
		printf("writing the doc...\n");
		for (int i = 1; i <= VectorSize(dpostProc); i++) {
			fprintf(fout, "%e\t", dpostProc[i]);
			if (i % step == 0)fprintf(fout, "\n");
		}
		fclose(fout);
	}
	else if (config.saveType == 3) {
		bf.numFrames = rowNum; bf.lengthFrame = 10000;
		bf.sizeFrameInByte = step * 4; bf.typeFlag = 0x49;
		fout = fopen(outputFile, "wb");
		if (!fout) { printf("open %s failed\n",outputFile); system("pause");  exit(1); }
		printf("writing the doc...\n");
		writeBinaryFile(fout, bf, dpostProc);
		fclose(fout);
	}
	else if (config.saveType == 2) {
		cnpy::npy_save(outputFile, &dpostProc[1], { (size_t)rowNum ,(size_t)step }, "w");
	}
	FreeVector(dpostProc);
	FreeMatrix(d1);
}
