#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <math.h>
//#include"tool.h"


#define PI acos(-1)
#define ARRAY_SIZE 786432
#define MAX_LINE 1024

typedef struct complex {  /*定义复数结构体*/
	double re;
	double im;
}COMPLEX;


typedef struct calibResult1 {
	double AngleMateMat[16][12];
	double RangeMat[16][12];
	double PeakValMat[16][12];
	double RxMismatch[16];
	double TxMismatch[3][4];
	double Rx_fft[1280][16][12];
}CalibResult;


typedef struct binfilePath1 {
	char master[100];
	char masterIdxFile[100];
	char slave1[100];
	char slave1IdxFile[100];
	char slave2[100];
	char slave2IdxFile[100];
	char slave3[100];
	char slave3IdxFile[100];
	char dataFolderName[100];
}BinfilePath;


typedef struct calibrationObj1 {
	BinfilePath binfilePath;
	char calibrationfilePath[100];
	double frameIdx;
	double adcCalibrationOn;
	double numSamplePerChirp;
	double nchirp_loops;
	double numChirpsPerFrame;
	double TxToEnable[12];
	double Slope_calib;
	double Sampling_Rate_sps;
	double fs_calib;
	double chirpSlope;
	double calibrationInterp;
	double TI_Cascade_RX_ID[16];
	double RxForMIMOProcess[16];
	double IdTxForMIMOProcess[12];
	double numRxToEnable;
	double phaseCalibOnly;
	double ADVANCED_FRAME_CONFIG;
	double N_TXForMIMO;
	double NumAnglesToSweep;
	char dataPlatform[10];
	double RxOrder[16];
	double NumDevices;
	double debugMode;
	double enable;
	char name[100];
	char pfile[100];
	//double parameters; //结构体
}CalibrationObj;


typedef struct rangeFFTObj1 {
	double enable;
	double rangeFFTSize;
	double rangeWindowCoeffVec;
	double FFTOutScaleOn;
	double scaleFactorRange;
}RangeFFTObj;


auto readBinFile(const char* fileFullPath, int frameIdx, int numSamplePerChirp, int numChirpPerLoop, int numLoops, int numRXPerDevice, int numDevices) {
	int Expected_Num_SamplesPerFrame = numSamplePerChirp * numChirpPerLoop * numLoops * numRXPerDevice * 2;
	printf("Expected_Num_SamplesPerFrame: %d\n", Expected_Num_SamplesPerFrame);
	FILE* fidlist = NULL;
	unsigned short buf;  /*缓冲区*/
	char dir[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\空时杂波数据\\实测数据_20210330\\STAP_CROSSROAD_20\\";
	char new_dir[200];
	sprintf_s(new_dir, "%s%s", dir, fileFullPath);
	puts(new_dir);
	fidlist = fopen(new_dir, "rb");
	//unsigned short e;
	int i = 0, j = 0;
	int* adcData1 = 0;  /*初始数据*/
	complex* adcData2 = 0; /*转复数*/
	unsigned short* neg = 0;
	adcData1 = (int*)malloc(sizeof(int) * 2048 * 1024);  /*动态分配*/
	adcData2 = (complex*)malloc(sizeof(complex) * 2048 * 1024);  /*复数数组*/
	neg = (unsigned short*)malloc(sizeof(unsigned short) * 2048 * 1024);
	fseek(fidlist, (frameIdx - 1) * Expected_Num_SamplesPerFrame * 2, SEEK_SET); /*定位*/
	// 循环读取文件
	while (!feof(fidlist) && i < Expected_Num_SamplesPerFrame)
	{
		// 定义函数返回值
		unsigned short rc = fread(&buf, sizeof(unsigned short), 1, fidlist);
		// 读取到文件的结束，退出循环
		if (rc <= 0)
			break;
		// 输出读取的结果到屏幕
		buf >= 32768 ? neg[i] = 1 : neg[i] = 0; /*逻辑数组，大于2`15*/
		buf >= 32768 ? adcData1[i] = int(buf) - 65536 : adcData1[i] = int(buf);
		//printf("%d %u\n", adcData1[i], neg[i]);
		if (i % 2 != 0) {
			int re = adcData1[i - 1];  /*实部*/
			int im = adcData1[i];  /*虚部*/
			adcData2[j].re = re;
			adcData2[j].im = im;
			//printf("\nindex %d: %d + %di", j+1, re, im);
			j++;
		}
		i++;
	}
	//reshape(numRXPerDevice, numSamplePerChirp, numChirpPerLoop, numLoops)
	complex(*adcData1Complex)[256][12][64];
	adcData1Complex = (complex(*)[256][12][64])malloc(2 * 786432 * sizeof(complex));
	int k = 0;
	if (adcData1Complex == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		for (int p = 0; p < 64; p++) {
			for (int q = 0; q < 12; q++) {
				for (int j = 0; j < 4; j++) {
					for (int i = 0; i < 256; i++) {  /*matlab列优先*/
						adcData1Complex[j][i][q][p].re = adcData2[k].re;
						adcData1Complex[j][i][q][p].im = adcData2[k].im;
						k++;
						//printf("index %d: %d + %di\n", k, adcData1Complex[j][i][q][p].re, adcData1Complex[j][i][q][p].im);
					}
				}
			}
		}
		for (int i = 52; i < 64; i++) {
			printf("reshape last %d:  %.1lf + %.1lfi\n", i + 1, adcData1Complex[3][255][11][i].re, adcData1Complex[3][255][11][i].im);
		}
	}

	//permute(adcData1Complex, [2 4 1 3])
	k = 0;
	typedef complex(*mytype)[64][4][12];
	mytype adcData1Complex1;
	adcData1Complex1 = (complex(*)[64][4][12])malloc(786432 * sizeof(complex));
	if (adcData1Complex1 == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		for (int i = 0; i < 256; i++) {  /*matlab列优先*/
			for (int j = 0; j < 4; j++) {
				for (int q = 0; q < 12; q++) {
					for (int p = 0; p < 64; p++) {
						adcData1Complex1[i][p][j][q].re = adcData1Complex[j][i][q][p].re;
						adcData1Complex1[i][p][j][q].im = adcData1Complex[j][i][q][p].im;
						k++;
						//printf("index %d: %d + %di\n", k, adcData1Complex1[q][j][p][i].re, adcData1Complex1[q][j][p][i].im);
					}
				}
			}
		}
		for (int i = 0; i < 12; i++) {
			printf("last %d:  %.1lf + %.1lfi\n", i + 1, adcData1Complex1[255][63][3][i].re, adcData1Complex1[255][63][3][i].im);
		}
	}
	//scanf("~~");
	puts("\n");
	fclose(fidlist);
	free(adcData1);
	adcData1 = NULL;
	free(adcData1Complex);
	adcData1Complex = NULL;
	free(adcData2);
	adcData2 = NULL;
	//free(adcData1Complex1);
	return adcData1Complex1;
}



auto read_ADC_bin_TDA2_separateFiles() {
	int frameIdx = 5, numSamplePerChirp = 256, numChirpPerLoop = 12, numLoops = 64, numRXPerDevice = 4, numDevices = 1;
	auto radar_data_Rxchain_master = readBinFile("master_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices);
	auto radar_data_Rxchain_slave1 = readBinFile("slave1_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices);
	auto radar_data_Rxchain_slave2 = readBinFile("slave2_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices);
	auto radar_data_Rxchain_slave3 = readBinFile("slave3_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices);
	typedef complex(*mytype)[64][16][12];
	mytype radar_data_Rxchain;
	radar_data_Rxchain = (complex(*)[64][16][12])malloc(4 * 786432 * sizeof(complex));
	if (radar_data_Rxchain == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		int k = 0;
		// 合并
		for (int i = 0; i < 64; i++) {
			for (int j = 0; j < 256; j++) {  /*列优先*/
				for (int m = 0; m < 4; m++) {
					for (int n = 0; n < 12; n++) {
						//k++;
						radar_data_Rxchain[j][i][m][n] = radar_data_Rxchain_master[j][i][m][n];
						radar_data_Rxchain[j][i][m + 4][n] = radar_data_Rxchain_slave1[j][i][m][n];
						radar_data_Rxchain[j][i][m + 8][n] = radar_data_Rxchain_slave2[j][i][m][n];
						radar_data_Rxchain[j][i][m + 12][n] = radar_data_Rxchain_slave3[j][i][m][n];
						//printf("data%d: %d + %di\n", k, radar_data_Rxchain[j][i][m + 12][n].re, radar_data_Rxchain[j][i][m + 12][n].im);
					}
				}
			}
		}
		for (int n = 0; n < 12; n++) {
			k++;
			printf("%d: %.1lf + %.1lfi\n", k, radar_data_Rxchain[255][63][3][n].re, radar_data_Rxchain[255][63][3][n].im);
		}
	}
	return radar_data_Rxchain;
}


auto datapath(CalibrationObj* obj) {
	//参数区
	char dataPlatform[10] = "TDA2";
	//char buf[MAX_LINE];  /*缓冲区*/
	int TxToEnable[12] = { 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0 };
	const int numSamplePerChrip = 256, adcCalibrationOn = 1, phaseCalibOnly = 1;
	int numTX = sizeof(TxToEnable) / sizeof(int);  /*12*/
	const double fs_calib = 8000000, Sampling_Rate_sps = 22500000, chirpSlope = 8.999300000000000e+13, Slope_calib = 7.898600000000000e+13;
	//double PeakValMat[12][16];
	double(*PeakValMat)[16] = (double(*)[16])malloc(sizeof(double) * 12 * 16);
	double(*RangeMat)[16] = (double(*)[16])malloc(sizeof(double) * 12 * 16);
	double(*freq_calib) = (double(*))malloc(sizeof(double) * 16);
	double(*phase_calib) = (double(*))malloc(sizeof(double) * 16);
	double(*correction_vec)[256] = (double(*)[256])malloc(sizeof(double) * 16 * 256);
	double(*freq_correction_mat)[256][64] = (double(*)[256][64])malloc(sizeof(double) * 16 * 256 * 64);
	double(*phase_correction_mat)[256][64] = (double(*)[256][64])malloc(sizeof(double) * 16 * 256 * 64);
	double(*freq_correction_mat_p)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 16 * 256 * 64);
	double(*phase_correction_mat_p)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 16 * 256 * 64);
	const double numChirpsPerFrame = 768, nchirp_loops = 64, numSamplePerChirp = 256, calibrationInterp = 5;
	double numSamplePerChirp_m[256];  /*numSamplePerChirp=256*/
	for (int i = 0; i < numSamplePerChirp; i++) {
		numSamplePerChirp_m[i] = i;
	}

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 16; j++) {
			PeakValMat[i][j] = -1000000;
		}
	}
	int TX_ref = TxToEnable[0];
	for (int iTX = 0; iTX < numTX; iTX++) {
		int TXind = TxToEnable[iTX];
		double phase_calib[] = { 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1 };
	}
	//typedef COMPLEX(*mytype)[64][16][12];
	//mytype outData;
	COMPLEX(*outData)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*输出；256-64-16-12*/
	COMPLEX(*outData_re)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*输出；256-64-16-12*/

	//工作区
	if (outData == NULL || outData_re == NULL || freq_correction_mat_p == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		if (strcmp(dataPlatform, "TDA2") == 0) {
			double numChirpPerLoop = numChirpsPerFrame / nchirp_loops;
			double numLoops = nchirp_loops;
			double numRXPerDevice = 4;
			auto radar_data_Rxchain = read_ADC_bin_TDA2_separateFiles(); //读取初始数据
			if (adcCalibrationOn == 0) {
				auto outData = radar_data_Rxchain;
			}
			else
			{
				int TX_ref = TxToEnable[0];
				for (int iTX = 0; iTX < numTX; iTX++) {  /*12~1*/
					int TXind = TxToEnable[iTX];
					for (int i = 0; i < 16; i++) {
						freq_calib[i] = (RangeMat[TXind][i] - RangeMat[TX_ref][0]) * fs_calib / Sampling_Rate_sps * chirpSlope / Slope_calib;
						freq_calib[i] = 2 * PI * freq_calib[i] / (numSamplePerChirp * calibrationInterp);
					}
					int k = 0;
					for (int i = 0; i < 16; i++) {
						for (int j = 0; j < 256; j++) {
							//带一个虚部i
							correction_vec[i][j] = exp(freq_calib[i] * numSamplePerChirp_m[j]);
							//printf("c%d: %d\n", k, correction_vec[i][j]);
							k++;
							//repat(correction_vec, 1, 1, nchirp_loops), 16-256-64
							for (int q = 0; q < 64; q++) {
								//虚部复数矩阵, 实部为0
								freq_correction_mat[i][j][q] = correction_vec[i][j];
							}
						}
					}
					//permute(freq_correction_mat, [2 3 1]), 256-64-16
					for (int i = 0; i < 64; i++) {  /*列优先*/
						for (int j = 0; j < 256; j++) {
							for (int q = 0; q < 16; q++) {
								freq_correction_mat_p[j][i][q] = freq_correction_mat[q][j][i];
							}
						}
					}
					//scanf("~~");
					typedef complex(*mytype)[64][16];
					mytype outData1TX;
					outData1TX = (complex(*)[64][16])malloc(262144 * sizeof(complex));  /*256-64-16*/
					//construct the phase compensation matrix
					for (int i = 0; i < 16; i++) {
						phase_calib[i] = PeakValMat[TX_ref][0] / PeakValMat[TXind][i];
						if (phaseCalibOnly == 1) {
							phase_calib[i] = phase_calib[i] / double(abs(int(phase_calib[i])));
						}
						//repmat(phase_calib.', 1,numSamplePerChirp, nchirp_loops), 16-256-64
						for (int j = 0; j < 256; j++) {
							for (int q = 0; q < 64; q++) {
								phase_correction_mat[i][j][q] = phase_calib[i];
							}
						}
					}
					//permute(phase_correction_mat, [2 3 1]),256-64-16
					for (int i = 0; i < 64; i++) {  /*列优先*/
						for (int j = 0; j < 256; j++) {
							for (int q = 0; q < 16; q++) {
								phase_correction_mat_p[j][i][q] = phase_correction_mat[q][j][i];
							}
						}
					}
					//截取矩阵切片
					for (int i = 0; i < 64; i++) {  /*列优先*/
						for (int j = 0; j < 256; j++) {
							for (int m = 0; m < 16; m++) {
								//复数乘法, 实部=虚*虚*(-1) + 实*实，虚部=实*虚 + 虚*实
								//outData1TX = radar_data_Rxchain(:,:,:,iTX)*freq_correction_mat * phase_correction_mat;
								double ir = radar_data_Rxchain[j][i][m][iTX].re * freq_correction_mat_p[j][i][m];
								double ii = radar_data_Rxchain[j][i][m][iTX].im * freq_correction_mat_p[j][i][m];

								outData1TX[j][i][m].re = ((-1) * ii + 0) * phase_correction_mat_p[j][i][m];
								outData1TX[j][i][m].im = (ir + 0) * phase_correction_mat_p[j][i][m];
							}
						}
					}
					//矩阵合并
					for (int i = 0; i < 64; i++) {
						for (int j = 0; j < 256; j++) {
							for (int m = 0; m < 16; m++) {
								outData[j][i][m][iTX] = outData1TX[j][i][m];
							}
						}
					}
					free(outData1TX);
					outData1TX = NULL;
				}
			}
		}
		else {
			printf("error: Not supported data capture platform!");
		}
		// 天线重排
		int  RxForMIMOProcess[16] = { 12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7 }; /*重排顺序*/
		for (int n = 0; n < 16; n++) {
			int k = RxForMIMOProcess[n];
			for (int i = 0; i < 64; i++) {
				for (int j = 0; j < 256; j++) {
					for (int m = 0; m < 12; m++) {
						outData_re[j][i][n][m] = outData[j][i][k][m];
					}
				}
			}
		}
	}

	free(freq_correction_mat);
	free(freq_correction_mat_p);
	free(phase_correction_mat);
	free(phase_correction_mat_p);
	free(outData);
	//freq_correction_mat = NULL;
	//freq_correction_mat_p = NULL;
	//phase_correction_mat = NULL;
	//phase_correction_mat_p = NULL;
	//outData = NULL;
	return outData_re;
}


auto rangeProcCascade_datapath(COMPLEX input[256][64][16]) {
	int numLines = 64, numAnt = 16, enable = 1, rangeFFTSize = 256;
	if (enable == 1) {
		double(*out)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 256 * 64 * 16);
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 64; j++) {
				for (int m = 0; m < 16; m++) {
					out[i][j][m] = 0;
				}
			}
		}
		COMPLEX(*inputMat)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
		for (int i_an = 0; i_an < numAnt; i_an++) {
			//inputMat = squeeze(input(:,:,i_an))  256,64,16 -> 256,64
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j] = input[i][j][i_an];
				}
			}
			//inputMat    = bsxfun(@minus, inputMat, mean(inputMat))
			double mean_re = 0, mean_im = 0, total_re = 0, total_im = 0;
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					total_re += inputMat[i][j].re;
					total_im += inputMat[i][j].im;
				}
			}
			mean_re = total_re / (256 * 64);
			mean_im = total_im / (256 * 64);
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j].re -= mean_re;
					inputMat[i][j].im -= mean_im;
				}
			}
			//inputMat    = bsxfun(@times, inputMat, obj.rangeWindowCoeffVec)

		}
	}

	return 0;
}


auto DopplerProcClutterRemove_datapath() {

	return 0;
}


int main() {
	//参数区
	BinfilePath binfilePath = {
		"master_0000_data.bin", "master_0000_idx.bin",
		"slave1_0000_data.bin", "slave1_0000_idx.bin",
		"slave2_0000_data.bin", "slave2_0000_idx.bin",
		"slave3_0000_data.bin", "slave3_0000_idx.bin",
		"E:\\MIMO_Clutter_Space_Time_Distribution\\空时杂波数据\\实测数据_20210330\\STAP_CROSSROAD_20\\"
	};
	CalibrationObj calibrationObj = {
		binfilePath, "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\calibrateResults_dummy.mat",
		5, 1, 256, 64, 768, {11,10,9,8,7,6,5,4,3,2,1,0}, 7.8986e+13, 22500000, 8000000, 
		8.9993e+13, 5, {12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7}, {12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7},
		{11,10,9,8,7,6,5,4,3,2,1,0}, 16, 1, 0, 0, 0, "TDA2", {12,13,14,15,0,1,2,3,8,9,10,11,4,5,6,7}, 4, 0, 1,
		"calibrationCascade", "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\test1_param.m"
	};
	//printf("%s!!!!!!!", calibrationObj.dataPlatform);
	//char pro_path[] = "E:\\MIMO_Clutter_Space_Time_Distribution";
	char input_path[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\";
	const char test_list[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\testList.txt";
	char pathGenParaFile[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\test1_param.m";
	//char pfile[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\test1_param.m";
	//char calibrationfilePath[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\main\\input\\calibrateResults_dummy.mat";
	int PARAM_FILE_GEN_ON = 1;
	FILE* fidlist = NULL;
	fidlist = fopen(test_list, "r");
	char dataFolder_test[1000];
	char dataFolder_calib[1000];
	char module_param_file[1000];
	//char pathGenParaFile[1000];
	char testID[2] = "1";
	//工作区
	while (!feof(fidlist)) {
		fgets(dataFolder_test, 1000, fidlist);
		fgets(dataFolder_calib, 1000, fidlist);
		fgets(module_param_file, 1000, fidlist);
		printf("%s\n%s\n%s\n", dataFolder_test, dataFolder_calib, module_param_file);
		sprintf_s(pathGenParaFile, "%s%s%s%s", input_path, "test", testID, "_param.m");
		if (PARAM_FILE_GEN_ON == 1) {
			int parameter_file_gen_json;  /*不明函数*/
		}
		auto adcData = datapath(&calibrationObj);
		for (int i_tx = 0; i_tx < 12; i_tx++) {
			COMPLEX(*adcData_slice)[64][16] = (COMPLEX(*)[64][16])malloc(262144  * sizeof(COMPLEX));  /*输出；256-64-16*/
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					for (int m = 0; m < 16; m++) {
						adcData_slice[i][j][m] = adcData[i][j][m][i_tx];
					}
				}
			}
			auto rangeFFTOut = rangeProcCascade_datapath(adcData_slice);
			free(adcData_slice);
			auto DopplerFFTOut = DopplerProcClutterRemove_datapath();
		}

	}
	return 0;
}

