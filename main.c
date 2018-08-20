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
regreOrder ����ϵ���Ľ�����1Ϊ��̬������2Ϊ��̬������һ�׼���ϵ����3Ϊ��̬������һ�׼���ϵ���Ͷ��׼���ϵ��
*/

/*  �����������зֱ�Ϊ ( MFCCNum+energyFlag+zeroCrossingFlag+brightFlag+subBankEFlag )*regreOrder  */
/*  �����������зֱ�Ϊ ( MFCC+����+������+������+�Ӵ����� )*����  */

/*��������ṹ
1.��ȡpcm�ļ�,��ת��Ϊʮ����
2.����MFCCϵ���Լ���������
3.������һ��
4.�������ϵ��
5.MFCC������һ����
6.д��Ŀ���ļ�
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
