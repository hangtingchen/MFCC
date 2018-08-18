#include "hmath.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//XFFT

void InitXFFT(XFFT * xfftP, int N)
{
	(*xfftP).N = N;
	(*xfftP).real = CreateVector(N);
	ZeroVector((*xfftP).real);
	(*xfftP).imag = CreateVector(N);
	ZeroVector((*xfftP).imag);
}

void ShowXFFT(XFFT xf)
{
	int i;
	for (i = 1; i <= xf.N; i++)printf("%f+%fi\t", xf.real[i], xf.imag[i]);
	printf("\n");
}

void ShowXFFTE(XFFT xf)
{
	int i;
	for (i = 1; i <= xf.N; i++) { 
		if(xf.imag[i]>=0)printf("%e+%ei\t", xf.real[i], xf.imag[i]); 
		else printf("%e-%ei\t", xf.real[i], -xf.imag[i]);
	}
	printf("\n");
}

void FreeXFFT(XFFT* xfftP)
{
	FreeVector((*xfftP).real);
	FreeVector((*xfftP).imag);
	(*xfftP).N = -1;
}

int XFFTSize(XFFT x)
{
	return x.N;
}

void XFFTToVector(XFFT xf, Vector * vp,int power2Flag)
{
	int i = 0; int N = xf.N; int N2 = N;
	if (power2Flag)N2 = (int)pow(2.0,ceil(log((double)N)/log(2)));
	(*vp) = CreateVector(2 * N2);
	for (i = 1; i <= N; i++) {
		(*vp)[2 * i - 1] = xf.real[i];
		(*vp)[2 * i] = xf.imag[i];
	}
	if (power2Flag) {
		for (i = N + 1; i <= N2; i++) {
			(*vp)[2 * i - 1] = 0;
			(*vp)[2 * i] = 0;
		}
	}
}

void VectorToXFFT(XFFT * xfp, Vector v)
{
	int i = 0; int N = VectorSize(v) / 2;
	InitXFFT(xfp, N);
	for (i = 1; i <= N; i++) {
		(*xfp).real[i] = v[2 * i - 1];
		(*xfp).imag[i] = v[2 * i];
	}
}

//IntVec

IntVec CreateIntVec(int n)
{
	int t;
	IntVec v = (int*)malloc(sizeof(int)*(n + 1));
	*(int *)v= n; for (t = 1; t <= n; t++)v[t] = 0;
	return v;
}

void FreeIntVec(IntVec v) { free(v); }

void ShowIntVec(IntVec v)
{
	int i=0;
	for (i = 1; i <= VectorSize(v); i++)printf("%d\t", v[i]);
	printf("\n");
}

void WriteIntVec(FILE * f, IntVec v)
{
	int i=0;
	for (i = 1; i <= VectorSize(v); i++)fprintf(f, "%d\t", v[i]);
	fprintf(f, "\n");
}

void ZeroIntVec(IntVec v)
{
	int i, n;
	n = VectorSize(v);
	for (i = 1; i <= n; i++) v[i] = 0;
}

void CopyIntVec(IntVec v1, IntVec v2)
{
	int i = 0;
	if (VectorSize(v1) > VectorSize(v2)) { printf("v1 and v2 size dismatch when copying vector and size(v1)>size(v2)\n"); 	system("pause"); exit(-1); }
	if (VectorSize(v1) != VectorSize(v2))printf("v1 and v2 size dismatch. If you do it on purpose, ignore this warning\n");
	for (i = 1; i <= VectorSize(v1); i++)v2[i] = v1[i];
}

//IntMat

IntMat CreateIntMat(int nrows, int ncols)
{
	size_t vsize;
	int *i, j;
	IntVec *m;
	char *p;//使用char是因为其是单个字节
	if ((nrows==0)||(ncols == 0))return NULL;
	p = (char*)malloc(sizeof(int)*((ncols + 1)*nrows) + (nrows + 1) * sizeof(IntVec));
	//其中第一个储存行数,其后作为索引,其后矩阵的第一列储存列数
	i = (int *)p; *i = nrows;
	vsize = (ncols + 1) * sizeof(int);
	m = (IntVec *)p;
	p += (nrows + 1) * sizeof(IntVec);
	for (j = 1; j <= nrows; j++, p += vsize) {
		i = (int *)p; *i = ncols;
		m[j] = (IntVec)p;
	}
	return m;
}

void FreeIntMat(IntMat m) { free(m); }

void ShowIntMat(IntMat m)
{
	int r = NumRows(m), c = NumCols(m);int i=0,j=0;
	for (i = 1; i <= r; i++) {
		for (j = 1; j <= c; j++)printf("%d\t", m[i][j]);
		printf("\n");
	}
}

void ZeroIntMat(IntMat m)
{
	int r = NumRows(m); int c = NumCols(m); int r0 = 0, c0 = 0;
	for (r0 = 1; r0 <= r; r0++)for (c0 = 1; c0 <= c; c0++)m[r0][c0] = 0;
}

void WriteIntMat(FILE * f, IntMat m)
{
	int i = 0, j = 0;
	for (i = 1; i <= NumRows(m); i++) {
		for (j = 1; j <= NumCols(m); j++)fprintf(f, "%d\t", m[i][j]);
		fprintf(f, "\n");
	}
}

//Vector

Vector CreateVector(int n)
{
	int *i; int j;
	Vector v = (double*)malloc(sizeof(double)*(n + 1));
	i = (int *)v; *i = n; for (j = 1; j <= n; j++)v[j] = 0.0;
	return v;
}

int VectorSize(Vector v){
//	printf("%d\n", *((int*)v));
	return *((int*)v);
}

void ZeroVector(Vector v)
{
	int i, n;
	n = VectorSize(v);
	for (i = 1; i <= n; i++) v[i] = 0.0;
}

void ShowVector(Vector v)
{
	int i=0;
	for (i = 1; i <= VectorSize(v); i++)printf("%f\t", v[i]);
	printf("\n");
}

void ShowVectorE(Vector v)
{
	int i=0;
	for (i = 1; i <= VectorSize(v); i++)printf("%e\t", v[i]);
	printf("\n");
}

void FreeVector(Vector v) { free(v); }

double FindMax(Vector v)
{
	int i; double m=v[1];
	for (i = 1; i <= VectorSize(v); i++)if (v[i] > m)m = v[i];
	return m;
}

int FindMaxIndex(Vector v)
{
	int i; int m = 1;
	for (i = 1; i <= VectorSize(v); i++)if (v[i] > v[m])m = i;
	return m;
}

void CopyVector(Vector v1, Vector v2)
{
	int i = 0;
	if (VectorSize(v1) > VectorSize(v2)) { printf("v1 and v2 size dismatch when copying vector and size(v1)>size(v2)\n"); 	system("pause"); exit(-1); }
	if (VectorSize(v1) != VectorSize(v2))printf("v1 and v2 size dismatch. If you do it on purpose, ignore this warning\n");
	for (i = 1; i <= VectorSize(v1); i++)v2[i] = v1[i];
}

void CopyVector2(Vector v1, IntVec ind1, Vector v2, IntVec ind2)
{
	int reMalloc = 0, maxIndex2 = 0; int i = 0; Vector temp = NULL;
	if (VectorSize(ind1) != VectorSize(ind2)) { printf("index size dismatch when copying vector \n"); 	system("pause"); exit(-1); }
	for (i = 1; i <= VectorSize(ind1); i++)if (ind2[i] > VectorSize(v2)) { reMalloc = 1; if (ind2[i] > maxIndex2)maxIndex2 = ind2[i]; }
	if (reMalloc) {
		temp = CreateVector(maxIndex2); ZeroVector(temp);
		CopyVector(v2, temp);
		for (i = 1; i <= VectorSize(ind2); i++)temp[ind2[i]] = v1[ind1[i]];
		free(v2); v2 = temp;
	}
	else {
		for (i = 1; i <= VectorSize(ind2); i++)temp[ind2[i]] = v1[ind1[i]];
	}
}

void WriteVectorE(FILE * f, Vector v)
{
	int i = 0;
	for (i = 1; i <= VectorSize(v); i++)fprintf(f, "%e\t", v[i]);
	fprintf(f, "\n");
}

void WriteVector(FILE * f, Vector v)
{
	int i = 0;
	for (i = 1; i <= VectorSize(v); i++)fprintf(f, "%f\t", v[i]);
	fprintf(f, "\n");
}

void LoadVector(FILE * f, Vector v)
{
	int i = 0;
	int size = VectorSize(v);
	for (i = 1; i <= size; i++)fscanf(f, "%lf\t", &v[i]);
}

void LoadVectorE(FILE * f, Vector v)
{
	int size = VectorSize(v); int i = 0;
	for (i = 1; i <= size; i++)fscanf(f, "%le\t", &v[i]);
}

//SVector

Vector CreateSVector(int n)
{
	SVector v;
	Ptr *p;
	int *i;
	p = (Ptr *)malloc(sizeof(double)*(n+1)+sizeof(Ptr)*2);
	v = (SVector)(p + 2);
	i = (int *)v; *i = n;
	SetHook(v, NULL);
	SetUse(v, 0);
	return v;
}

void FreeSVector(SVector v)
{
	Ptr* ptr= (Ptr*)v - 2;
	free(ptr);
}

//Methods for Shared Vec or Mat

Ptr GetHook(Ptr m)
{
	Ptr *p = (Ptr*)m; p -= 2;
	return *p;//interesting
	//m-=2;return m;
}

void SetHook(Ptr m, Ptr ptr)
{
	Ptr *p = (Ptr*)m; p -= 2; *p = ptr;
}

void SetUse(Ptr m, int n)
{
	Ptr* p = (Ptr*)m; --p; *((int *)p) = n;
}

//Matrix

Matrix CreateMatrix(int nrows, int ncols)
{
	size_t vsize;
	int *i, j;
	Vector *m;
	char *p;//使用char是因为其是单个字节
	p = (char*)malloc(sizeof(double)*((ncols + 1)*nrows) + (nrows + 1) * sizeof(Vector));
	//其中第一个储存行数,其后作为索引,其后矩阵的第一列储存列数
	i = (int *)p; *i = nrows;
	vsize = (ncols+1)*sizeof(double);
	m = (Vector *)p;
	p += (nrows + 1) * sizeof(Vector);
	for (j = 1; j <= nrows; j++, p += vsize) {
		i = (int *)p; *i = ncols;
		m[j] = (Vector)p;
	}
	return m;
}

int NumRows(Matrix m)
{
	int *nrows;
	if (m == NULL)return 0;
	nrows = (int *)m;
	return *nrows;
}

int NumCols(Matrix m)
{
	int *ncols;
	if (m == NULL)return 0;
	ncols = (int *)m[1];
	return *ncols;
}

void ShowMatrix(Matrix m)
{
	int r = NumRows(m), c = NumCols(m);int i=0,j=0;
	for (i = 1; i <= r; i++) {
		for (j = 1; j <= c; j++)printf("%f\t", m[i][j]);
		printf("\n");
	}
}

void FreeMatrix(Matrix m) { free(m); }

void ZeroMatrix(Matrix m)
{
	int r = NumRows(m); int c = NumCols(m); int r0 = 0, c0 = 0;
	for (r0 = 1; r0 <= r; r0++)for (c0 = 1; c0 <= c; c0++)m[r0][c0] = 0.0;
}

void CopyMatrix(Matrix m1, Matrix m2)
{
	int r = NumRows(m1); int c = NumCols(m1); int r0 = 0, c0 = 0;
	if ((NumRows(m1) != NumRows(m2)) || (NumCols(m1) != NumCols(m2))) { printf("m1 and m2 size dismatch when copying matrix\n"); 	system("pause"); exit(-1); }
	for (r0 = 1; r0 <= r; r0++)for (c0 = 1; c0 <= c; c0++)m2[r0][c0] = m1[r0][c0];
}

void CopyMatToTri(Matrix m1, STriMat m2)
{
	int i = 0, j = 0;
	if (NumRows(m1) != STriMatSize(m2)) {
		printf("dismatch when copying mat to trimat\n"); exit(-1);
	}
	for (i = 1; i <= STriMatSize(m2); i++) {
		for (j = 1; j <= i; j++)m2[i][j] = m1[i][j];
	}
}

void WriteMatrix(FILE * f, Matrix m)
{
	int i = 0, j = 0;
	for (i = 1; i <= NumRows(m); i++) {
		for (j = 1; j <= NumCols(m); j++)fprintf(f, "%f\t", m[i][j]);
		fprintf(f, "\n");
	}
}

void LoadMatrix(FILE * f, Matrix m)
{
	int i = 0, j = 0;
	for (i = 1; i <= NumRows(m); i++) {
		for (j = 1; j <= NumCols(m); j++)fscanf(f, "%lf\t", &m[i][j]);
		fscanf(f, "\n");
	}
}

//SMatrix

SMatrix CreateSMatrix(int nrows, int ncols)
{
	size_t vsize;
	int *i, j;
	Vector *m;
	char *p;

	p = (char*)malloc(sizeof(double)*((ncols + 1)*nrows) + (nrows + 3) * sizeof(Vector))+2 * sizeof(Ptr *);
	i = (int *)p; *i = nrows;
	vsize = (ncols + 1) * sizeof(double);
	m = (Vector *)p;
	p += (nrows + 1) * sizeof(Vector);
	for (j = 1; j <= nrows; j++, p += vsize) {
		i = (int *)p; *i = ncols;
		m[j] = (Vector)p;
	}
	SetHook(m, NULL);
	SetUse(m, 0);
	return m;
}

void FreeSMatrix(SMatrix m)
{
	Ptr* ptr = (Ptr*)m - 2;
	free(ptr);
}

//STriMat

Matrix CreateSTriMat(int size)
{
	int *i, j;
	Vector *m;
	char *p;
	p = (char *)malloc(size*(sizeof(double) * 2 + (size + 1) * sizeof(double)) / 2+ (size + 1) * sizeof(Vector) + 2 * sizeof(Ptr)) + 2 * sizeof(Ptr *);
	i = (int *)p; *i = size;
	m = (Vector *)p;
	p += (size + 1) * sizeof(Vector);
	for (j = 1; j <= size; j++) {
		i = (int *)p; *i = j;
		m[j] = (Vector)p; p += sizeof(double)*(j+1);
	}
	SetHook(m, NULL);
	SetUse(m, 0);
	return m;
}

int STriMatSize(STriMat m)
{
	int *nrows;
	nrows = (int *)m;
	return *nrows;
}

void ShowSTriMat(STriMat m)
{
	int i=0,j=0;
	for (i = 1; i <= STriMatSize(m); i++) {
		for (j = 1; j <= i; j++)printf("%f\t", m[i][j]);
		printf("\n");
	}
}

void FreeSTriMat(STriMat m)
{
	Ptr* ptr = (Ptr*)m - 2;
	free(ptr);
}

void ZeroSTriMat(STriMat m)
{
	int i, j, size;
	size = STriMatSize(m);
	for (i = 1; i <= size; i++)
		for (j = 1; j <= i; j++) m[i][j] = 0.0;
}

void CopySTriMat(STriMat m1, STriMat m2)
{
	int i, size;
	size = STriMatSize(m1);
	if (size != STriMatSize(m2))printf("dim dismatch when copying stri");
	for (i = 1; i <= size; i++)
		CopyVector(m1[i], m2[i]);
}

void WriteSTriMat(FILE * f, STriMat m)
{
	int size = STriMatSize(m);int i=0,j=0;double t=0;
	for (i = 1; i <= size; i++) {
		for (j = 1; j <= i; j++) { t = m[i][j]; fprintf(f, "%f\t",t); }
		fprintf(f, "\n");
	}
}

void LoadStriMat(FILE * f, STriMat m)
{
	int size = STriMatSize(m);int i=0,j=0;
	for (i = 1; i <= size; i++) {
		for (j = 1; j <= i; j++) { fscanf(f, "%lf\t", &m[i][j]); }
		fscanf(f, "\n");
	}
}



int Choleski(STriMat A, Matrix L)
{
	int size, i, j, k;
	double sum;

	size = STriMatSize(A);
	for (i = 1; i <= size; i++)
		for (j = 1; j <= i; j++) {
			sum = A[i][j];
			for (k = 1; k < j; k++) {
				sum -= (L[i][k] * L[j][k]);
			}
			if ((i == j) && (sum <= 0.0))
				return 0;
			else if (i == j)
				sum = sqrt(sum);
			else if (L[j][j] == 0.0)
				return 0;
			else
				sum /= L[j][j];
			L[i][j] = sum;
		}
	for (i = 1; i <= size; i++)
		for (j = i + 1; j <= size; j++)
			L[i][j] = 0.0;
	return 1;
}

void MSolve(Matrix L, int i, Vector x, Vector y)
{
	int nr, j, k;
	double sum;

	nr = NumRows(L);
	for (j = 1; j<i; j++) y[j] = 0.0;  /* forward sub */
	y[i] = 1.0 / L[i][i];
	for (j = i + 1; j <= nr; j++) {
		sum = 0.0;
		for (k = i; k<j; k++)
			sum -= L[j][k] * y[k];
		y[j] = sum / L[j][j];
	}
	x[nr] = y[nr] / L[nr][nr];         /* backward sub */
	for (j = nr - 1; j >= 1; j--) {
		sum = y[j];
		for (k = j + 1; k <= nr; k++)
			sum -= L[k][j] * x[k];
		x[j] = sum / L[j][j];
	}
}

logdouble CovInvert(STriMat c, STriMat invc)
{
	Matrix l;     /* Lower Tri Choleski Matrix */
	Vector x, y;   /* for f/b substitution */
	logdouble ldet = 0.0;
	int i, j, n;
	int isTri=1;
//	for (int i = 1; i <= STriMatSize(c); i++)c[i][i] += 0.000001;

	n = STriMatSize(c); 
	x = CreateVector( n);
	y = CreateVector( n);
	l = CreateMatrix(n, n);
	if (Choleski(c, l)) {
		for (j = 1; j <= n; j++) {
			MSolve(l, j, x, y);
			for (i = isTri ? j : 1; i <= n; i++)
				invc[i][j] = x[i];
			ldet += log(l[j][j]);
		}
	}
	else {
		FILE* fcov = fopen("cov.dat", "w");
		WriteSTriMat(fcov, c);
		fclose(fcov);
		printf("the matrix is not inversible\n");
		system("pause");
		exit(-1);
	}
	FreeVector(x); FreeVector(y);
	FreeMatrix( l);    /* cut back stack to entry state */
	return 2.0*ldet;
}

logdouble CovDet(STriMat c)
{
	Matrix l;  /* Lower Tri Choleski Matrix */
	logdouble ldet = 0.0;
	int j, n;

	n = STriMatSize(c);
	l = CreateMatrix( n, n);
	if (Choleski(c, l)) {
		for (j = 1; j <= n; j++)
			ldet += log(l[j][j]);
	}
	else
		printf("the matrix is not inversible\n");
	FreeMatrix( l);
	return 2.0*ldet;
}

int mod(int a, int b)
{
	if (b <= 0)printf("Modulus calculation may be wrong because the dividend <=0 , please check\n");
	if (a < 0)do { a += b; } while (a < 0);
	return a%b;
}

void reshape(Matrix * mp, Vector v, int r, int c, int dim)
{
	int i = 0; int n = VectorSize(v); int index_r, index_c;
	if (n != (r*c))printf("reshape when numRow*numColumn!=VECTORSIZE !\n");
	(*mp) = CreateMatrix(r, c);
	if (dim == 1) {
		for (i = 1; i <= n; i++) {
			index_r = ((i - 1) % r) + 1;
			index_c = ((i - 1) / r) + 1;
			(*mp)[index_r][index_c] = v[i];
		}
	}
	else if (dim == 2) {
		for (i = 1; i <= n; i++) {
			index_r = ((i - 1) / c) + 1;
			index_c = ((i - 1) % c) + 1;
			(*mp)[index_r][index_c] = v[i];
		}
	}
	else { printf("reshape:dim wrong\n"); exit(-1); }
}

