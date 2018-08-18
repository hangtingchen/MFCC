#pragma once
#define _CRT_SECURE_NO_WARNINGS
#ifndef _HMATH_H_
#define _HMATH_H_

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#include<stdio.h>
#include<math.h>
typedef double logdouble;
typedef void* Ptr;
typedef double* Vector;
typedef double** Matrix;
typedef int* IntVec;
typedef int* ShortVec;
typedef int** IntMat;
//#define pi 3.1415926
#define pi 3.1415926535897932384626433832795028841971

typedef Vector SVector;
typedef Matrix STriMat;
typedef Matrix SMatrix;

typedef struct {	/*用于储存复数的数组结构，主要用于频域和快速傅里叶变换*/
	Vector real;	/*实数和复数部分被分别存储在两个向量中，大小都为N*/
	Vector imag;
	int N;
}XFFT;

/*------------------向量部分----------------------------------------*/
/*-----------------------------------------------------------------*/
/*可供所有向量通用的函数VectorSize*/
/*可供Vector，SVector向量通用的函数，除了Create和Free都可以*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/*------------------Vector------------------*/
/*------------------double型向量------------------*/
Vector CreateVector(int n);
int VectorSize(Vector v);	/*返回向量的大小*/
void ZeroVector(Vector v);	/*向量归零*/
void ShowVector(Vector v);
void ShowVectorE(Vector v);	/*用科学计数格式打印*/
void FreeVector(Vector v);
double FindMax(Vector v);	/*找到向量中的最大元素*/
int FindMaxIndex(Vector v);	/*找到向量中最大元素的下标*/
void CopyVector(Vector v1, Vector v2);	/*将v1复制到v2，v2要事先创建好*/
void CopyVector2(Vector v1, IntVec ind1, Vector v2, IntVec ind2);	/*根据下标ind1和ind2来复制v1到v2*/
void WriteVectorE(FILE* f, Vector v);	/*用科学计数将向量写入文件*/
void WriteVector(FILE* f, Vector v);
void LoadVector(FILE * f, Vector v);	/*从文件中载入向量*/
void LoadVectorE(FILE * f, Vector v);	/*用科学计数从文件中载入向量*/

/*------------------IntVec------------------*/
/*------------------整形向量------------------*/
//VectorSize(IntVec v) 来读取大小
IntVec CreateIntVec(int n);	/*创建一个整形向量*/
void FreeIntVec(IntVec v);
void ShowIntVec(IntVec v);
void WriteIntVec(FILE* f, IntVec v);	/*将一个向量储存到文件f中*/
void ZeroIntVec(IntVec v);
void CopyIntVec(IntVec v1, IntVec v2);

/*------------------SVector------------------*/
/*------------------double型共享向量------------------*/
//其他操作可以使用Vector的函数
SVector CreateSVector(int n);
void FreeSVector(SVector v);


/*------------------复数向量部分----------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/


/*------------------XFFT------------------*/
/*------------------存储复数向量------------------*/
void InitXFFT(XFFT* xfftP, int N);	/*输入N来初始化XFFT*/
void ShowXFFT(XFFT xf);	/*在屏幕上打印一组复数*/
void ShowXFFTE(XFFT xf);	/*在屏幕上打印复数，用科学计数*/
void FreeXFFT(XFFT* xfftP);	
int XFFTSize(XFFT x);	/*输出复数向量的大小*/
void XFFTToVector(XFFT xf, Vector* vp, int power2Flag);	/*将复数结构转化为数组形式，数组的长度为复数数组的长度的两倍，间隔存储实部和虚部*/
void VectorToXFFT(XFFT* xfp, Vector v);	/*将长度为2N的数组转化为XFFT形式，要求在v中间隔存储实部和虚部*/



/*------------------矩阵部分----------------------------------------*/
/*-----------------------------------------------------------------*/
/*可供IntMat，Matrix，SMatrix通用的函数NumRows，NumCols*/
/*可供Matrix，SMatrix通用的函数，除了Create和Free都可以*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/*------------------IntMat------------------*/
/*------------------整形矩阵------------------*/
//NumRows(IntMat m),NumCols(IntMat m)来读取矩阵的行数和列数
IntMat CreateIntMat(int nrows, int ncols);
void FreeIntMat(IntMat m);
void ZeroIntMat(IntMat m);
void WriteIntMat(FILE* f, IntMat m);
void ShowIntMat(IntMat m);


/*------------------Matrix------------------*/
/*------------------double型矩阵------------------*/
Matrix CreateMatrix(int nrows, int ncols);	/*创建矩阵，需要行数和列数的值的输入*/
int NumRows(Matrix m);	/*返回矩阵m的行数*/
int NumCols(Matrix m);	/*返回矩阵m的列数*/
void ShowMatrix(Matrix m);
void FreeMatrix(Matrix m);
void ZeroMatrix(Matrix m);
void CopyMatrix(Matrix m1, Matrix m2);	/*从矩阵m1复制到m2，m2要实现创建好*/
void CopyMatToTri(Matrix m1, STriMat m2);	/*将对称矩阵转换为下三角矩阵*/
void WriteMatrix(FILE* f, Matrix m);
void LoadMatrix(FILE* f, Matrix m);

/*------------------SMatrix------------------*/
/*------------------double型共享矩阵------------------*/
//其他操作可以使用Matrix的函数
SMatrix CreateSMatrix(int nrows, int ncols);
void FreeSMatrix(SMatrix m);

/*------------------STriMat------------------*/
/*------------------double型共享下三角矩阵------------------*/
Matrix CreateSTriMat(int size);
int STriMatSize(STriMat m);	/*下三角矩阵一定是方阵，行和列数目相同*/
void ShowSTriMat(STriMat m);
void FreeSTriMat(STriMat m);
void ZeroSTriMat(STriMat m);
void CopySTriMat(STriMat m1, STriMat m2);
void WriteSTriMat(FILE* f, STriMat m);
void LoadStriMat(FILE*f, STriMat m);

/*------------------Methods for Shared Vec or Mat------------------*/
Ptr GetHook(Ptr m);
void SetHook(Ptr m, Ptr ptr);
void SetUse(Ptr m, int n);


/*------------------求矩阵的逆和行列式的相关------------------*/
/*------------------算法限制为对称的方阵，因此只能用于下三角矩阵------------------*/
int Choleski(STriMat A, Matrix L);
void MSolve(Matrix L, int i, Vector x, Vector y);
logdouble CovInvert(STriMat c, STriMat invc);	/*计算矩阵c的逆矩阵invc，并返回log(det(c))*/
logdouble CovDet(STriMat c);	/*返回log(det(c))*/

/*------------------一些其他函数------------------*/
int mod(int a, int b);	/*计算a%b，主要考虑a为负值的情况*/
void reshape(Matrix* mp, Vector v, int r, int c,int dim);	/*将v按照dim重新排列，生成矩阵mp，mp不需要事先分配空间*/

#ifdef __cplusplus
}
#endif // __cplusplus

#endif
