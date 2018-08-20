#include<cstdio>
#include<cstdlib>
#include<cmath>
#include<string>
#include<fstream>
#include"omp.h"
#include"ini.h"
#include "hmath.hpp"
#include "mfcc.h"
#include"sigProcess.hpp"
#include"WAVE.hpp"
#include"cnpy.hpp"
#include"fileIO.hpp"

using namespace hmath;
using namespace hWAVE;
using namespace hsigProcess;
using namespace hfileIO;

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




int main(int argc, char** argv) {

	Config config;

	int i0 = ini_parse(argv[1], handler, &config);
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

	std::ifstream fileList(config.fileList);
	std::vector<std::string> fileListString;
	char fileNameBuf[maxBuffLength];
	
	if (!fileList.is_open()) {
		printf("Fail to open %s", config.fileList);
	}
	while (fileList.getline(fileNameBuf, maxBuffLength)) {
		fileListString.push_back( strtok(fileNameBuf, "\t"));
		fileListString.push_back( strtok(NULL, "\n"));
	}
	fileList.close();

	MFCCWapperTempStruct mwts = MFCCWapperTempInit(config);

	omp_set_num_threads(config.numThreads);
	#pragma omp parallel for
	for (int i = 0; i < fileListString.size(); i = i + 2) {
		int ID = omp_get_thread_num()+1;
		printf("Thread ID %d\n", ID);
		MFCCWapper(fileListString[i].c_str(), fileListString[i + 1].c_str(),mwts,config,ID);
	}
	MFCCWapperTempFree(&mwts,config);
	//system("pause");
	return 0;
}
