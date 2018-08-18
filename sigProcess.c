#include "sigProcess.h"

/*circle shift the signal*/
void circshift(Vector v, int shift)
{
	int i = 1; Vector v_temp = CreateVector(VectorSize(v));
	if (shift < 0)do { shift += VectorSize(v); } while (shift < 0);
	if (shift >= VectorSize(v))do { shift -= VectorSize(v); } while (shift >= VectorSize(v));
	for (i = 1; (i + shift) <= VectorSize(v); i++)v_temp[i + shift] = v[i];
	for (; i <= VectorSize(v); i++)v_temp[i + shift - VectorSize(v)] = v[i];
	CopyVector(v_temp, v);
	FreeVector(v_temp);
}

/*find the first index of abs|sample| exceeding the thre from the front or end*/
int find(Vector v, double thre, int FrontOrEndFlag)
{
	int i; int m = 0;
	if (FrontOrEndFlag == 1) {
		for (i = 1; i <= VectorSize(v); i++)if (fabs(v[i]) > thre) { m = i; break; }
		return m;
	}
	else if (FrontOrEndFlag == -1) {
		for (i = VectorSize(v); i >= 1; i--)if (fabs(v[i]) > thre) { m = i; break; }
		return m;
	}
	else
		return 0;
}

void pad_signal(Vector * yP, Vector x, int Npad)
{
	int i = 0; int j = 0;
	int orig_sz = VectorSize (x);
	int Norig = VectorSize(x);
	int end = Norig + (int)floor((double)(Npad - Norig) / 2.0);
	int end2 = (int)ceil((double)(Npad - Norig) / 2.0);
	int end3 = Norig + (int)floor((double)(Npad - Norig) / 2.0) + 1;
	IntVec ind0 = CreateIntVec(2 * Norig);
	IntVec conjugate0 = CreateIntVec(2 * Norig);
	IntVec conjugate = CreateIntVec(Npad);
	IntVec ind = CreateIntVec(Npad);
	IntVec src = CreateIntVec(end - Norig);
	IntVec dst = CreateIntVec(end - Norig);
	ZeroIntVec(ind); ZeroIntVec(conjugate);
	for (i = 1; i <= Norig; i++)ind0[i] = i;
	for (i=Norig; i >= 1; i--)ind0[2 * Norig - i + 1] = i;
	for (i = 1; i <= Norig; i++)conjugate0[i] = 0;
	for (i = Norig + 1; i <= 2 * Norig; i++)conjugate0[i] = 1;
	for (i = 1; i <= Norig; i++)ind[i] = i;
	for (i = 1; i <= VectorSize(src); i++)src[i] = (Norig + i-1) % (VectorSize(ind0)) + 1;
	for (i = 1; i <= VectorSize(dst); i++)dst[i] = Norig + i;
	for (i = 1; i <= VectorSize(src); i++)ind[dst[i]] = ind0[src[i]];
	for (i = 1; i <= VectorSize(src); i++)conjugate[dst[i]] = conjugate0[src[i]];
	FreeIntVec(src); FreeIntVec(dst);
	src = CreateIntVec(end2); dst = CreateIntVec(Npad - end3 + 1);
	for (i = 1; i <= VectorSize(src); i++) {
		if((VectorSize(ind0) - i)>=0)src[i] = ((VectorSize(ind0) - i) % (VectorSize(ind0))) + 1;
		else src[i] = ((VectorSize(ind0) - i+ VectorSize(ind0)) % (VectorSize(ind0))) + 1;
	}
	for (i = Npad, j = 1; i >= end3; i--, j++)dst[j] = i;
	for (i = 1; i <= VectorSize(src); i++)ind[dst[i]] = ind0[src[i]];
	for (i = 1; i <= VectorSize(src); i++)conjugate[dst[i]] = conjugate0[src[i]];
	*yP = CreateVector(VectorSize(ind));
	for (i = 1; i <= VectorSize(ind); i++)(*yP)[i] = x[ind[i]];
	FreeIntVec(ind0); FreeIntVec(conjugate0); FreeIntVec(conjugate); FreeIntVec(ind); FreeIntVec(src); FreeIntVec(dst);
}

void unpad_signal(Vector * yP, Vector x, int res, int target_sz)
{
	int i = 0;
	int padded_sz = VectorSize(x);
	double offset = 0;
	int offset_ds = 0;
	int target_sz_ds = 1 + floor((double)(target_sz - 1) / pow(2.0, (double)res));
	*yP = CreateVector(target_sz_ds);
	for (i = 1; i <= VectorSize(*yP); i++)(*yP)[i] = x[i];
}

Matrix frameRawSignal(IntVec v, int wlen, int inc, double preEmphasiseCoefft, int enableHamWindow)
{
	int numSamples = VectorSize(v);
	int numFrames = (numSamples - (wlen - inc)) / inc;
	Matrix m = NULL; 
	Vector v1 = NULL; Vector HamWindow = NULL;
	int i = 0, j = 0, pos = 1; double a = 0;

	v1 = CreateVector(numSamples);
	for (i = 1; i <= numSamples; i++)v1[i] = (double)v[i];
	PreEmphasise(v1, preEmphasiseCoefft);

	HamWindow = CreateVector(wlen);
	a = 2 * pi / (wlen - 1);
	if(enableHamWindow)for (i = 1; i <= wlen; i++)HamWindow[i] = 0.54 - 0.46 * cos(a*(i - 1));
	else for (i = 1; i <= wlen; i++)HamWindow[i] = 1.0;

	if ((numSamples - (inc - wlen)) % inc != 0)numFrames++;
	m = CreateMatrix(numFrames, wlen);
	for (i = 1; i <= numFrames; i++) {
		pos = (i - 1)*inc + 1;
		for (j = 1; j <= wlen; j++,pos++) {
			if (pos > numSamples)m[i][j] = 0.0;
			else m[i][j] = (double)v1[pos]*HamWindow[j];
		}
	}
	FreeVector(v1); FreeVector(HamWindow);
	return m;
}

/*FFT from HTK*/
/* EXPORT-> FFT: apply fft/invfft to complex s */
/*
void FFT(Vector s, int invert)
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
		inc = 2 * limit; theta = 2 * pi / limit;
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
	if (invert) {
		for (i = 1; i <= n; i++)
			s[i] = s[i] / nn;
	}
	for (i = 3; i <= n / 2; i++) {
		if (mod(i, 2) == 0)ii = n + 4 - i;
		else ii = n + 2 - i;
		x = s[ii]; s[ii] = s[i]; s[i] = x;
	}
}
*/



void ZeroMean(IntVec data)
{
	long i, hiClip = 0, loClip = 0;
	int *x;
	double sum = 0.0, y, mean;
	int nSamples = VectorSize(data);

	x = &data[1];
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
		*x = (int)((y>0.0) ? y + 0.5 : y - 0.5);
	}
	if (loClip>0)
		printf("ZeroMean: %d samples too -ve\n", loClip);
	if (hiClip>0)
		printf("ZeroMean: %d samples too +ve\n", hiClip);
}

double zeroCrossingRate(Vector s, int frameSize) {
	int count = 0; int i;
	for (i = 1; i < frameSize; i++)  if ((s[i] * s[i + 1]) < 0.0)count++;
	return (double)count / (double)(frameSize - 1);
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

double calBrightness(Vector fftx)
{
	int i;
	double sum = 0.0;
	double te = 0.0;
	double b = 0;
	if (((int)VectorSize(fftx)) % 2 != 0)printf("something wrong in cal brightness");
	for (i = 1; i <= ((int)VectorSize(fftx)) / 2; i++) {
		sum += (fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i])*(double)i;
		te += fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i];
	}
	b = sum / te;
	b = b / ((double)VectorSize(fftx) / 2.0);
	return b;
}


void calSubBankE(Vector fftx, Vector subBankEnergy)
{
	int i;
	int numBank = VectorSize(subBankEnergy); int bankSize = (int)VectorSize(fftx) / (2 * numBank);
	int bankNum = 1;
	double te = 0.0;
	double sum = 0.0;
	for (i = 1; i <= (int)VectorSize(fftx) / 2; i++)te+= fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i];
	for (i = 1; i <= (int)VectorSize(fftx) / 2; i++) {
		if (i <= bankNum*bankSize) {
			sum += fftx[2 * i - 1] * fftx[2 * i - 1] + fftx[2 * i] * fftx[2 * i];
		}
		else {
			subBankEnergy[bankNum] = sum / te;
			//printf("sum: %f\n", sum/te);
			bankNum++; sum = 0.0; i--;
		}
	}
	subBankEnergy[bankNum] = sum / te;

}

void Regress(double* data, int vSize, int n, int step, int offset, int delwin, int head, int tail, int simpleDiffs)
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

void RegressMat(Matrix* m, int delwin,int regressOrder)
{
	Vector v = NULL;
	int dimOrigin = NumCols(*m), numFrames = NumRows(*m); int dimAfter = dimOrigin*(1+regressOrder);
	int i = 0, j = 0;
	
	if (regressOrder == 0)return;

//	printf("%d\t%d\n", NumRows(*m), NumCols(*m));
	v = CreateVector(dimOrigin*numFrames*(regressOrder+1));
	for(i=1;i<=numFrames;i++)for (j = 1; j <= dimOrigin; j++) v[j + dimAfter*(i - 1)] = (*m)[i][j];

	for (i = 1; i <= regressOrder; i++) {
		Regress(&v[1+(i-1)*dimOrigin], dimOrigin, numFrames, dimAfter, dimOrigin, delwin, 0, 0, 0);
	}
	FreeMatrix(*m);

	*m = CreateMatrix(numFrames, dimAfter);
//	printf("%d\t%d\n", NumRows(*m), NumCols(*m));
	for (i = 1; i <= numFrames; i++)for (j = 1; j <= dimAfter; j++)  (*m)[i][j]= v[j + dimAfter*(i - 1)]  ;
	FreeVector(v);
}

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

/*Z-normalization, Not tested */
void ZNormalize(double * data, int vSize, int n, int step)
{
	double sum1,sum2;
	double *fp, sd,mean;
	int i,j;
	for (i = 0; i < vSize; i++) {
		sum1 = 0.0;sum2=0.0;
		fp = data + i;
		for (j = 0; j < n; j++) { sum1+=(*fp);sum2 += (*fp)*(*fp); fp += step; }
		mean=sum1/(double)n;
		sd = sqrt(sum2 / (double)n-mean*mean);
		fp = data + i;
		for (j = 0; j < n; j++) {
			*fp = ((*fp)-mean)/sd; fp += step;
		}
	}
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
