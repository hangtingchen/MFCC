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
using namespace hMFCC;


int main(int argc, char** argv) {

	Config config;
    	char * pchs;
    	char * pcht;
	
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

	configCheck(&config);

	std::ifstream fileList(config.fileList);
	std::vector<std::string> fileListString;
	char fileNameBuf[maxBuffLength];
	
	if (!fileList.is_open()) {
		printf("Fail to open %s", config.fileList);
	}
	while (fileList.getline(fileNameBuf, maxBuffLength)) {
        pchs=strtok(fileNameBuf,"\t");pcht=strtok(NULL,"\n");
        	if(pcht==NULL){pchs=strtok(fileNameBuf," ");pcht=strtok(NULL,"\n");}
        	if(pcht==NULL){printf("Unable to found delimeter tab or space in filelist\n");return 1;}
        	fileListString.push_back(pchs);
		fileListString.push_back(pcht);
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
