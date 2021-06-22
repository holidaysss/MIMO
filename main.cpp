#include<stdio.h>
#include<string.h>
#include <stdlib.h>
#include <math.h>
#include "tool.h"


#define PI acos(-1)
#define ARRAY_SIZE 786432
#define MAX_LINE 1024

typedef struct complex {  /*定义复数结构体*/
	double re;
	double im;
}COMPLEX, Complex;


typedef struct detectionResult {
	int rangeInd;
	double range;
	int dopplerInd_org;
	int dopplerInd;
	double doppler;
	double doppler_corr;
	double doppler_corr_overlap;
	double doppler_corr_FFT;
	double overlapTests;
	double overlapTestsVal;
	double noise_var;
	COMPLEX *bin_val;
	double estSNR;
}Detection, DR;

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


void readBinFile(const char* fileFullPath, int frameIdx, int numSamplePerChirp, int numChirpPerLoop, int numLoops, int numRXPerDevice, int numDevices, COMPLEX adcData1Complex[][64][4][12]) {
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
	//permute(adcData1Complex, [2 4 1 3])
	int k = 0;
	if (adcData1Complex == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		for (int j = 0; j < 64; j++) {
			for (int p = 0; p < 12; p++) {
				for (int i = 0; i < 256; i++) {
					for (int q = 0; q < 4; q++) {  /*matlab列优先*/
						adcData1Complex[i][j][q][p].re = adcData2[k].re;
						adcData1Complex[i][j][q][p].im = adcData2[k].im;
						k++;
						//printf("index %d: %d + %di\n", k, adcData1Complex[j][i][q][p].re, adcData1Complex[j][i][q][p].im);
					}
				}
			}
		}
		//for (int i = 0; i < 4; i++) {
		//	printf("reshape last %d:  %.1lf + %.1lfi\n", i + 1, adcData1Complex[255][63][i][0].re, adcData1Complex[255][63][i][0].im);
		//}
		//scanf("end");
	}
	puts("\n");
	fclose(fidlist);
	free(adcData1);
	free(adcData2);
	free(neg);
	adcData1 = NULL;
	adcData2 = NULL;
	neg = NULL;
}


void read_ADC_bin_TDA2_separateFiles(COMPLEX radar_data_Rxchain[][64][16][12]) {
	int frameIdx = 5, numSamplePerChirp = 256, numChirpPerLoop = 12, numLoops = 64, numRXPerDevice = 4, numDevices = 1;
	COMPLEX(*radar_data_Rxchain_master)[64][4][12] = (complex(*)[64][4][12])malloc(786432 * sizeof(complex));
	COMPLEX(*radar_data_Rxchain_slave1)[64][4][12] = (complex(*)[64][4][12])malloc(786432 * sizeof(complex));
	COMPLEX(*radar_data_Rxchain_slave2)[64][4][12] = (complex(*)[64][4][12])malloc(786432 * sizeof(complex));
	COMPLEX(*radar_data_Rxchain_slave3)[64][4][12] = (complex(*)[64][4][12])malloc(786432 * sizeof(complex));
	readBinFile("master_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices, radar_data_Rxchain_master);
	readBinFile("slave1_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices, radar_data_Rxchain_slave1);
	readBinFile("slave2_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices, radar_data_Rxchain_slave2);
	readBinFile("slave3_0000_data.bin", frameIdx, numSamplePerChirp, numChirpPerLoop, numLoops, numRXPerDevice, numDevices, radar_data_Rxchain_slave3);
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
	free(radar_data_Rxchain_master);
	free(radar_data_Rxchain_slave1);
	free(radar_data_Rxchain_slave2);
	free(radar_data_Rxchain_slave3);
	radar_data_Rxchain_master = NULL;
	radar_data_Rxchain_slave1 = NULL;
	radar_data_Rxchain_slave2 = NULL;
	radar_data_Rxchain_slave3 = NULL;

}


void calibration_datapath(CalibrationObj* obj, COMPLEX outData_re[][64][16][12]) {
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
	double(*phase_correction_mat)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 16 * 256 * 64);
	double(*freq_correction_mat_p)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 16 * 256 * 64);
	//double(*phase_correction_mat_p)[64][16] = (double(*)[64][16])malloc(sizeof(double) * 16 * 256 * 64);
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
	COMPLEX(*outData)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*输出；256-64-16-12*/
	//工作区
	if (outData == NULL || freq_correction_mat_p == NULL) {
		printf("内存分配不成功！\n");
	}
	else {
		if (strcmp(dataPlatform, "TDA2") == 0) {
			double numChirpPerLoop = numChirpsPerFrame / nchirp_loops;
			double numLoops = nchirp_loops;
			double numRXPerDevice = 4;
			COMPLEX(*radar_data_Rxchain)[64][16][12] = (complex(*)[64][16][12])malloc(4 * 786432 * sizeof(complex));
			read_ADC_bin_TDA2_separateFiles(radar_data_Rxchain); //读取初始数据
			//for (int i = 0; i < 16; i++) {
			//	printf("radar_data_Rxchain: re %.5lf, im %.5lf\n", radar_data_Rxchain[255][63][i][0].re, radar_data_Rxchain[255][63][i][0].im);
			//}
			//scanf("end");
			if (adcCalibrationOn == 0) {
				auto outData = radar_data_Rxchain;
			}
			else
			{
				int TX_ref = TxToEnable[0];
				typedef complex(*mytype)[64][16];
				mytype outData1TX;
				outData1TX = (complex(*)[64][16])malloc(262144 * sizeof(complex));  /*256-64-16*/
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
							k++;
							//repat(correction_vec, 1, 1, nchirp_loops), 16-256-64
							for (int q = 0; q < 64; q++) {
								//当虚部不为0时：虚部复数矩阵, 实部为0；当虚部为0时：只有实部
								freq_correction_mat[i][j][q] = correction_vec[i][j];
							}
						}
					}
					//permute(freq_correction_mat, [2 3 1]), 256-64-16, 存疑！
					for (int j = 0; j < 256; j++) {
						for (int q = 0; q < 16; q++) {
							for (int i = 0; i < 64; i++) {
								freq_correction_mat_p[j][i][q] = freq_correction_mat[q][j][i];
							}
						}
					}
					//outData1TX = (complex(*)[64][16])realloc(outData1TX, 262144 * sizeof(complex));  /*256-64-16*/

					//construct the phase compensation matrix
					//repmat(phase_calib.', 1,numSamplePerChirp, nchirp_loops), 16-256-64
					//permute(phase_correction_mat, [2 3 1]),256-64-16， 存疑！
					for (int j = 0; j < 256; j++) {
						for (int i = 0; i < 16; i++) {
							phase_calib[i] = PeakValMat[TX_ref][0] / PeakValMat[TXind][i];
							if (phaseCalibOnly == 1) {
								phase_calib[i] = phase_calib[i] / double(abs(int(phase_calib[i])));
							}
							for (int q = 0; q < 64; q++) {
								phase_correction_mat[j][q][i] = phase_calib[i];
							}
						}
					}
					//截取矩阵切片
					for (int i = 0; i < 64; i++) {  /*列优先*/
						for (int j = 0; j < 256; j++) {
							for (int m = 0; m < 16; m++) {
								//复数乘法, 实部=虚*虚*(-1) + 实*实，虚部=实*虚 + 虚*实
								//outData1TX = radar_data_Rxchain(:,:,:,iTX)*freq_correction_mat * phase_correction_mat;
								double rr = radar_data_Rxchain[j][i][m][iTX].re * freq_correction_mat_p[j][i][m];
								double ir = radar_data_Rxchain[j][i][m][iTX].im * freq_correction_mat_p[j][i][m];
								outData1TX[j][i][m].re = (rr + 0) * phase_correction_mat[j][i][m];
								outData1TX[j][i][m].im = (ir + 0) * phase_correction_mat[j][i][m];
								//printf("radar_data_rxchain: re %.5lf, im %.5lf\n", outData1TX[j][i][m].re, outData1TX[j][i][m].im);
							}
						}
					}
					//for (int i = 0; i < 16; i++) {
					//	printf("outData1TX: re %.5lf, im %.5lf\n", outData1TX[255][63][i].re, outData1TX[255][63][i].im);
					//}
					//scanf("end");

					//矩阵合并
					for (int i = 0; i < 64; i++) {
						for (int j = 0; j < 256; j++) {
							for (int m = 0; m < 16; m++) {
								outData[j][i][m][iTX] = outData1TX[j][i][m];
							}
						}
					}
					
				}
				free(outData1TX);
				outData1TX = NULL;
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

	free(PeakValMat);
	free(RangeMat);
	free(freq_calib);
	free(phase_calib);
	free(correction_vec);
	free(freq_correction_mat);
	free(phase_correction_mat);
	free(freq_correction_mat_p);
	//free(phase_correction_mat_p);
	free(outData);
	PeakValMat = NULL;
	RangeMat = NULL;
	freq_calib = NULL;
	phase_calib = NULL;
	correction_vec = NULL;
	freq_correction_mat = NULL;
	phase_correction_mat = NULL;
	freq_correction_mat_p = NULL;
	outData = NULL;
}


void rangeProcCascade_datapath(COMPLEX input[256][64][16], COMPLEX out[][64][16]) {
	double numLines = 64, numAnt = 16, enable = 1, rangeFFTSize = 256, FFTOutScaleOn = 0, scaleFactorRange = 0.0156;
	COMPLEX(*inputMat)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
	COMPLEX(*fftOutput)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
	if (enable == 1) {
		for (int i_an = 0; i_an < numAnt; i_an++) {
			//inputMat = squeeze(input(:,:,i_an))  256,64,16 -> 256,64
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j] = input[i][j][i_an];			
				}
			}
			//inputMat = bsxfun(@minus, inputMat, mean(inputMat))
			double mean_re[64], mean_im[64], total_re = 0, total_im = 0;
			for (int j = 0; j < 64; j++) {
				for (int i = 0; i < 256; i++) {
					total_re += inputMat[i][j].re;
					total_im += inputMat[i][j].im;
				}
				mean_re[j] = total_re / 256;
				mean_im[j] = total_im / 256;
				total_re = 0;
				total_im = 0;
			}
			double rangeWindowCoeffVec[] = { 0.0001, 0.0006, 0.0013, 0.0024, 0.0037, 0.0054, 0.0073, 0.0095, 0.0121, 0.0149, 0.0180, 0.0214, 0.0250, 0.0290, 0.0332, 0.0378, 0.0426, 0.0476, 0.0530, 0.0586, 0.0645, 0.0706, 0.0770, 0.0836, 0.0905, 0.0977, 0.1050, 0.1126, 0.1205, 0.1286, 0.1369, 0.1454, 0.1541, 0.1630, 0.1721, 0.1815, 0.1910, 0.2007, 0.2106, 0.2206, 0.2308, 0.2412, 0.2518, 0.2625, 0.2733, 0.2842, 0.2953, 0.3065, 0.3179, 0.3293, 0.3408, 0.3525, 0.3642, 0.3760, 0.3879, 0.3998, 0.4118, 0.4239, 0.4360, 0.4481, 0.4603, 0.4725, 0.4847, 0.4969, 0.5092, 0.5214, 0.5336, 0.5458, 0.5579, 0.5701, 0.5821, 0.5942, 0.6061, 0.6181, 0.6299, 0.6417, 0.6533, 0.6649, 0.6764, 0.6878, 0.6991, 0.7102, 0.7213, 0.7322, 0.7429, 0.7535, 0.7640, 0.7743, 0.7844, 0.7944, 0.8042, 0.8138, 0.8232, 0.8324, 0.8415, 0.8503, 0.8589, 0.8673, 0.8755, 0.8835, 0.8912, 0.8987, 0.9059, 0.9130, 0.9197, 0.9262, 0.9325, 0.9385, 0.9442, 0.9497, 0.9549, 0.9599, 0.9645, 0.9689, 0.9730, 0.9768, 0.9804, 0.9836, 0.9866, 0.9892, 0.9916, 0.9937, 0.9955, 0.9970, 0.9982, 0.9991, 0.9997, 1.0000, 1.0000, 0.9997, 0.9991, 0.9982, 0.9970, 0.9955, 0.9937, 0.9916, 0.9892, 0.9866, 0.9836, 0.9804, 0.9768, 0.9730, 0.9689, 0.9645, 0.9599, 0.9549, 0.9497, 0.9442, 0.9385, 0.9325, 0.9262, 0.9197, 0.9130, 0.9059, 0.8987, 0.8912, 0.8835, 0.8755, 0.8673, 0.8589, 0.8503, 0.8415, 0.8324, 0.8232, 0.8138, 0.8042, 0.7944, 0.7844, 0.7743, 0.7640, 0.7535, 0.7429, 0.7322, 0.7213, 0.7102, 0.6991, 0.6878, 0.6764, 0.6649, 0.6533, 0.6417, 0.6299, 0.6181, 0.6061, 0.5942, 0.5821, 0.5701, 0.5579, 0.5458, 0.5336, 0.5214, 0.5092, 0.4969, 0.4847, 0.4725, 0.4603, 0.4481, 0.4360, 0.4239, 0.4118, 0.3998, 0.3879, 0.3760, 0.3642, 0.3525, 0.3408, 0.3293, 0.3179, 0.3065, 0.2953, 0.2842, 0.2733, 0.2625, 0.2518, 0.2412, 0.2308, 0.2206, 0.2106, 0.2007, 0.1910, 0.1815, 0.1721, 0.1630, 0.1541, 0.1454, 0.1369, 0.1286, 0.1205, 0.1126, 0.1050, 0.0977, 0.0905, 0.0836, 0.0770, 0.0706, 0.0645, 0.0586, 0.0530, 0.0476, 0.0426, 0.0378, 0.0332, 0.0290, 0.0250, 0.0214, 0.0180, 0.0149, 0.0121, 0.0095, 0.0073, 0.0054, 0.0037, 0.0024, 0.0013, 0.0006, 0.0001 };
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j].re -= mean_re[j];
					inputMat[i][j].im -= mean_im[j];
					//inputMat = bsxfun(@times, inputMat, obj.rangeWindowCoeffVec)
					inputMat[i][j].re *= rangeWindowCoeffVec[i];
					inputMat[i][j].im *= rangeWindowCoeffVec[i];
					//printf("inputMat:%.5lf + %.5lf\n", inputMat[i][j].re, inputMat[i][j].im);
				}
			}
			//scanf("end");
			//fftOutput = fft(inputMat, obj.rangeFFTSize)
			double pr[256], pi[256], fr[256], fi[256];
			int k = 0;
			for (int j = 0; j < 64; j++) {
				k = 0;
				for (int i = 0; i < 256; i++) {
					pr[k] = inputMat[i][j].re;
					pi[k] = inputMat[i][j].im;
					//printf("p[k]:%.5lf + %.5lf\n", pr[k], pi[k]);
					k++;
				}
				kfft(pr, pi, 256, 8, fr, fi); //fix: 每次输入一行256个点
				k = 0;
				for (int i = 0; i < 256; i++) {
					fftOutput[i][j].re = fr[k];
					fftOutput[i][j].im = fi[k];
					//printf("f[k]:%.5lf + %.5lf\n", fr[k], fi[k]);
					k++;
				}
				//scanf(":e");
			}
			//for (int i = 0; i < 10; i++) {
			//	printf("fftOutput:%.5lf + %.5lf\n", fftOutput[2][i].re, fftOutput[2][i].im);
			//}
			//scanf("end");
			
			if (FFTOutScaleOn == 1) {
				for (int i = 0; i < 64; i++) {
					for (int j = 0; j < 256; j++) {
						fftOutput[j][i].re *= scaleFactorRange;
						fftOutput[j][i].im *= scaleFactorRange;
					}
				}
			}
			//out(:,:,i_an) = fftOutput
			for (int i = 0; i < 64; i++) {
				for (int j = 0; j < 256; j++){
					out[j][i][i_an] = fftOutput[j][i];
				}
			}
		}
		//for (int i = 0; i < 10; i++) {
		//	printf("out:%.5lf + %.5lf\n", out[255][i][0].re, out[255][i][0].im);
		//}
		//scanf("end");
	}
	else
	{
		out = input;
	}
	free(inputMat);
	free(fftOutput);
	inputMat = NULL;
	fftOutput = NULL;
}


void fftshift_2(COMPLEX input[][64], COMPLEX out[][64]) {  //2维 中心对调
	int k = 0;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 64; j++) {
			k = (j + 32) % 64;
			out[i][j] = input[i][k];
		}
	}
}


void fftshift_2_1(COMPLEX input[][7], COMPLEX out[][7]) {  //2维 上下对调
	int k = 0;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 7; j++) {
			out[i][j] = input[255-i][j];
		}
	}
}

void fftshift_2_2(COMPLEX input[][256], COMPLEX out[][256]) {  //2维 上下对调
	int k = 0;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 256; j++) {
			out[i][j] = input[i][255-j];
		}
	}
}


void dopplerProcClutterRemove_datapath(COMPLEX input[256][64][16], COMPLEX out[][64][16]) {
	int numLines = 256, numAnt = 16;
	double enable = 1, clutterRemove = 0, dopplerFFTSize = 64, FFTOutScaleOn = 0, scaleFactorDoppler = 0.0625; /*obj属性*/
	COMPLEX(*inputMat)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
	COMPLEX(*inputMat_t)[256] = (COMPLEX(*)[256])malloc(sizeof(COMPLEX) * 256 * 64);
	COMPLEX(*fftOutput)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
	COMPLEX(*fftOutput1)[64] = (COMPLEX(*)[64])malloc(sizeof(COMPLEX) * 256 * 64);
	double pr[64], pi[64], fr[64], fi[64];
	int k;
	if (enable == 1) {
		for (int i_an = 0; i_an < numAnt; i_an++) {
			//inputMat = squeeze(input(:,:,i_an))
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j] = input[i][j][i_an];
				}
			}
			//inputMat = bsxfun(@times, inputMat, obj.dopplerWindowCoeffVec.')
			double dopplerWindowCoeffVec[64];
			for (int i = 0; i < 64; i++) {
				dopplerWindowCoeffVec[i] = 1;
			}
			for (int i = 0; i < 256; i++) {
				for (int j = 0; j < 64; j++) {
					inputMat[i][j].re *= dopplerWindowCoeffVec[j];
					inputMat[i][j].im *= dopplerWindowCoeffVec[j];
					inputMat_t[j][i] = inputMat[i][j]; /*inputMat矩阵的转置:64-256*/
				}
			}
			if (clutterRemove == 1) {
				//inputMat = inputMat - (repmat(mean(inputMat'),size(inputMat,2),1))'
			}
			//fftOutput = fft(inputMat, obj.dopplerFFTSize, 2)
			for (int j = 0; j < 256; j++) {
				k = 0;
				for (int i = 0; i < 64; i++) {
					pr[k] = inputMat_t[i][j].re;
					pi[k] = inputMat_t[i][j].im;
					k++;
				}
				kfft(pr, pi, 64, 6, fr, fi);  //pr,pi -> input ; fr,fi -> output ; 采样点数 64=2~6
				k = 0;
				for (int i = 0; i < 64; i++) {
					fftOutput[j][i].re = fr[k];
					fftOutput[j][i].im = fi[k];
					k++;
				}
			}
			//fftshift
			fftshift_2(fftOutput, fftOutput1); /*左右行互换*/
			if (FFTOutScaleOn == 1) {
				for (int i = 0; i < 256; i++) {
					for (int j = 0; j < 64; j++) {
						fftOutput1[i][j].re = fftOutput1[i][j].re * scaleFactorDoppler;
						fftOutput1[i][j].im = fftOutput1[i][j].im * scaleFactorDoppler;
					}
				}
			}
			//out(:,:,i_an) = fftOutput
			for (int j = 0; j < 256; j++) {
				for (int i = 0; i < 64; i++) {
					out[j][i][i_an] = fftOutput1[j][i];
				}
			}
		}
	}
	else
	{
		out = input;
	}
	free(inputMat);
	free(inputMat_t);
	free(fftOutput);
	free(fftOutput1);
	inputMat = NULL;
	inputMat_t = NULL;
	fftOutput = NULL;
	fftOutput1 = NULL;
}


void writetxt(double input[], size_t size, const char name[])
{
	FILE* fp;
	unsigned int i = 0;
	fp = fopen(name, "w");
	//printf("longggggggggg: %d\n", size);
	while (i< size) {
		fprintf(fp, "%.4lf\n", input[i]);
		i++;
	}
	fclose(fp);
}


void writetxt_2d(double input[][64], size_t size, const char name[])
{
	FILE* fp;
	fp = fopen(name, "w");
	//printf("longggggggggg: %d", size);
	for (unsigned int i = 0; i < size; i++) {
		for (int j = 0; j < 64; j++) {
			fprintf(fp, "%.4lf\n", input[i][j]);
		}
	}
	fclose(fp);
}


void CFAR_CASO_Range_2D(double sig_integrate[256][64], int* N_obj, int Ind_obj[][2], double* noise_obj, double* CFAR_SNR) {
	//参数区
	char home[] = "E:\\MIMO_Clutter_Space_Time_Distribution\\modules\\detection";
	int obj_refWinSize[2] = { 8, 4 };
	int obj_guardWinSize[2] = { 8, 0 };
	int obj_KO[2] = { 5, 3 };
	int cellNum = obj_refWinSize[0];
	int gapNum = obj_guardWinSize[0];
	double KO = obj_KO[0];
	int M_samp = 256, N_pul = 64;
	int gaptot = gapNum + cellNum; //16
	//int N_obj = 0;     /*检测到的目标数*/
	int discardCellLeft = 10, discardCellRight = 10;
	int discardCellTop = 10, discardCellBottom = 20;
	double *sigv = (double(*))malloc(sizeof(double) * 64);
	double *v = (double(*))malloc(sizeof(double) * 44);
	double *vecLeft = (double(*))malloc(sizeof(double) * 16);
	double *vecRight = (double(*))malloc(sizeof(double) * 16);
	double(*vec)[76] = (double(*)[76])malloc(sizeof(double) * M_samp * 76); /*左右补数据*/


	//工作区
	for (int k = 0; k < M_samp; k++) {  /*行操作*/
		for (int i = 0; i < 64; i++) {
			sigv[i] = sig_integrate[k][i];
		}
		for (int i = 0; i < 44; i++) {
			v[i] = sigv[i + discardCellLeft];
			vec[k][i + 16] = v[i];  /*居中*/
		}
		for (int i = 0; i < 16; i++) {
			vecLeft[i] = v[i];
			vecRight[i] = v[i + 28];  /*44-16=28*/
			vec[k][i] = vecLeft[i];
			vec[k][i+16+44] = vecRight[i];
		}
	}

	int m = 256, n = 76; /*size(vec)*/
	double(*sigv_1) = (double(*))malloc(sizeof(double) * 256);
	double(*v_1) = (double(*))malloc(sizeof(double) * 226);

	double(*vector)[76] = (double(*)[76])malloc(sizeof(double) * 258 * 76); /*左右补数据, 226+16+16=258*/
	for (int k = 0; k < n; k++) {  /*列操作*/
		for (int i = 0; i < 256; i++) {
			sigv_1[i] = vec[i][k];
		}
		for (int i = 0; i < 226; i++) {
			v_1[i] = sigv_1[i + discardCellTop];
			vector[i+16][k] = v_1[i];  /*居中*/
		}
		for (int i = 0; i < 16; i++) {
			vecLeft[i] = v_1[i];
			vecRight[i] = v_1[i + 210];  /*226-16=210*/
			vector[i][k] = vecLeft[i];
			vector[i + 16 + 226][k] = vecRight[i];
		}
	}
	m = 258, n = 76;
	//operator = ones(2*gaptot + 1)/((gaptot*2+1)^2-(gapNum*2+1)^2)
	double(*operator_1)[33] = (double(*)[33])malloc(sizeof(double*) * 33 * 33*5); /* 2*gaptot+1=33 */
	for (int i = 0; i < 2 * gaptot + 1; i++) {
		for (int j = 0; j < 2 * gaptot + 1; j++) {
			operator_1[i][j] = 1 / (pow(gaptot * 2 + 1, 2) - pow(gapNum * 2 + 1, 2));
		}
	}
	//operator(cellNum+1:cellNum+2*gapNum+1,cellNum+1:cellNum+2*gapNum+1) = 0
	for (int i = cellNum; i < cellNum + 2 * gapNum+1; i++) {
		for (int j = cellNum; j < cellNum + 2 * gapNum+1; j++) {
			operator_1[i][j] = 0;
		}
	}
	//meanCell = conv2(vector,operator);
	/* step1: padding */
	double(*meanCell_pre)[140] = (double(*)[140])malloc(sizeof(double) * 322 * 140*5);
	for (int i = 0; i < 322; i++) {
		for (int j = 0; j < 140; j++) {
			meanCell_pre[i][j] = 0;
		}
	}
	for (int i = 32; i < 290; i++) {
		for (int j = 32; j < 108; j++) {
			meanCell_pre[i][j] = vector[i-32][j-32];
		}
	}
	/* step2: whirl operator 180` */
	double(*operator_2)[33] = (double(*)[33])malloc(sizeof(double) * 33 * 33 * 5);
	for (int i = 0; i < 33; i++) {
		for (int j = 0; j < 33; j++) {
			operator_2[i][j] = operator_1[32 - i][32 - j];
		}
	}
	/* step3: convolution */
	double(*meanCell)[108] = (double(*)[108])malloc(sizeof(double) * 290 * 108);
	for (int i = 0; i < 290; i++) {
		for (int j = 0; j < 108; j++) {
			meanCell[i][j] = 0;  //初始化
			for (int m = i; m < i + 33; m++) {
				for (int n = j; n < j + 33; n++) {
					meanCell[i][j] += meanCell_pre[m][n] * operator_2[m - i][n - j];
				}
			}
		}
	}
	//result = K0*meanCell(gaptot + 1:m - gaptot , gaptot + 1:n - gaptot)< vector(gaptot + 1:m - gaptot, gaptot + 1 : n - gaptot);
	int(*result)[44] = (int(*)[44])malloc(sizeof(int) * 226 * 44);
	int cout=0;  /*TrueNum*/
	for (int i = 0; i < 226; i++) {
		for (int j = 0; j < 44; j++) {
			if (KO*meanCell[i+gaptot][j+gaptot]<vector[i+gaptot][j+gaptot]) {
				result[i][j] = 1;
				cout++;
			}
			else
			{
				result[i][j] = 0;
			}
		}
	}
	//[TargetInd(:,1),TargetInd(:,2)] = find(result==1);
	//Ind_obj = [TargetInd(:,1)+discardCellTop,TargetInd(:,2)+discardCellLeft];
	int(*TargetInd)[2] = (int(*)[2])malloc(sizeof(int) * cout * 2);
	int index = 0;
	for (int j = 0; j < 44; j++){  /*列优先*/
		for (int i = 0; i < 226; i++){
			if (result[i][j] == 1){
				TargetInd[index][0] = i+1;
				TargetInd[index][1] = j+1;
				Ind_obj[index][0] = i+discardCellTop+1;
				Ind_obj[index][1] = j+discardCellLeft+1;
				index++;
			}
		}
	}
	*N_obj = cout;
	//noise_obj = meanCell((TargetInd(:,2)+gaptot-1)*size(meanCell,1)+gaptot+TargetInd(:,1));
	double meanCell_index;
	for (int i = 0; i < cout; i++) {
		meanCell_index = (int(TargetInd[i][1] + gaptot) - 1) * 290 + gaptot + TargetInd[i][0];
		index = 1;
		for (int k = 0; k < 108; k++) { //列优先遍历
			for (int j = 0; j < 290; j++) {
				if (index == meanCell_index) {
					noise_obj[i] = meanCell[j][k];
				}
				index++;
			}
		}
	}
	//CFAR_SNR = sig((Ind_obj(:,2)-1)*size(sig,1)+Ind_obj(:,1))./noise_obj;
	int sig_index;
	for (int i = 0; i < cout; i++) {
		sig_index = ((Ind_obj[i][1] - 1) * 256 + Ind_obj[i][0]);
		index = 1;
		for (int k = 0; k < 64; k++) { //列优先遍历
			for (int j = 0; j < 256; j++) {
				if (index == sig_index) {
					CFAR_SNR[i] = sig_integrate[j][k] / noise_obj[i];
				}
				index++;
			}
		}
	}

	//for (int i = 0; i < 211; i++) {
	//	printf("index:%d CFAR_SNR: %.5lf\n", i+1, CFAR_SNR[i]);
	//}
	//scanf("end");
	free(TargetInd);
	free(result);
	free(meanCell);
	free(sigv);
	free(v);
	free(vecLeft);
	free(vecRight);
	free(sigv_1);
	free(v_1);
	free(vector);
	free(vec);
	free(operator_1);
	free(meanCell_pre);
	free(operator_2);
	TargetInd = NULL;
	result = NULL;
	meanCell = NULL;
	vector = NULL;
	vec = NULL;
	operator_1 = NULL;
	meanCell_pre = NULL;
	operator_2 = NULL;

}


int cmpfunc(const void* a, const void* b)
{
	return (*(int*)a - *(int*)b);
}


void CFAR_CASO_Doppler_overlap(int Ind_obj[][2], int *N_obj, int Ind_obj_Rag[][2], COMPLEX sigCpml[256][64][192], double sig_integ[256][64]) {
	//参数区
	int maxEnable = 0;
	int cellNum0[2] = { 8,4 };
	int gapNum0[2] = { 8,0 };
	int cellNum = 4, gapNum = 0, KO = 3;
	int rangeNumBins = 256; /*size(sig_integ,1)*/
	//工作区
	//detected_Rag_Cell = unique(Ind_obj_Rag(:,1))
	int *detected_Rag_Cell = (int*)malloc(rangeNumBins * sizeof(int));
	int index = 0, flag;
	detected_Rag_Cell[0] = Ind_obj_Rag[0][0];  /*首值*/
	for (int i = 1; i < 211;i++) {
		flag = 1;
		for (int j = 0; j < index + 1; j++) {  /*忽略重复值*/
			if (detected_Rag_Cell[j] == Ind_obj_Rag[i][0]) {
				flag = 0;
				break;
			}
		}
		if (flag == 1) {  /*添加唯一值*/
			detected_Rag_Cell[++index] = Ind_obj_Rag[i][0];
		}
	}
	if (detected_Rag_Cell != 0) {
		qsort(detected_Rag_Cell, index+1, sizeof(detected_Rag_Cell[0]), cmpfunc);  /*快排*/
	}
	else {
		printf("快排失败");
	}
	double(*sig)[64] = (double(*)[64])malloc((index+1)*64 * sizeof(double));//83*64

	//sig = sig_integ(detected_Rag_Cell,:);
	for (int i = 0; i < index+1; i++) {
		for (int j = 0; j < 64; j++) {
			sig[i][j] = sig_integ[detected_Rag_Cell[i]-1][j];
		}
	}
	//for (int k = 0; k < 211; k++) {
	//	printf("index:%d Ind_obj_Rag: %d\n", k + 1, Ind_obj_Rag[k][1]);
	//}
	//scanf("end");
	int M_samp = 83, N_pul = 64;
	int gaptot = gapNum + cellNum;  /* gaptot=4 */
	int detected_Rag_Cell_i, ind1_num, ind_obj_0 = 0,  condition, j0, ind_win;
	double sum = 0;
	double cellave1a, cellave1b, cellave1 = 0, maxInCell;
	int index_a, index_Ind_obj=0;
	double(*sigv) = (double*)malloc(64 * sizeof(double));
	double(*vec) = (double*)malloc((N_pul + gaptot) * 2 * sizeof(double));  //奇怪的断点bug
	int(*cellInd) = (int*)malloc(cellNum * 2 * sizeof(int));
	double(*noiseEst) = (double*)malloc(N_pul * sizeof(double));
	int(*ind_loc_all) = (int*)malloc(N_pul*211 * sizeof(int));
	int(*ind_loc_Dop) = (int*)malloc(N_pul*211 * sizeof(int));
	int(*ind_obj_00)[2] = (int(*)[2])malloc(N_pul*211 *2 * sizeof(int));
	int Ind_obj_FLAG = 0;  /*初值flag*/
	int(*ind1) = (int*)malloc(211 * sizeof(int));
	int(*indR) = (int*)malloc(211 * sizeof(int));
	for (int i = 0; i < M_samp; i++) {
		ind1_num = 0;
		detected_Rag_Cell_i = detected_Rag_Cell[i];
		// ind1 = find(Ind_obj_Rag(:,1) == detected_Rag_Cell_i)
		for (int j = 0; j < 211; j++) {
			if (Ind_obj_Rag[j][0] == detected_Rag_Cell_i) {
				ind1[ind1_num++] = j;
				//printf("index:%d vec: %d\n", i + 1, j);
			}
		}
		// indR = Ind_obj_Rag(ind1, 2)
		for (int k = 0; k < ind1_num; k++) {
			indR[k] = Ind_obj_Rag[ind1[k]][1];
		}

		//sigv=(sig(k,:))
		for (int n = 0; n < 64; n++) {
			sigv[n] = sig[i][n];
		}
		//vec(1:gaptot) = sigv(end - gaptot + 1:end);
		//vec(N_pul + gaptot + 1:end) = sigv(1:gaptot);
		for (int m = 0; m < gaptot; m++) {
			vec[m] = sigv[N_pul - gaptot + m];
			vec[N_pul + gaptot + m] = sigv[m];
		}
		//vec(gaptot + 1: N_pul + gaptot) = sigv;
		for (int p = gaptot; p < N_pul +gaptot; p++) {
			vec[p] = sigv[p-gaptot];
		}

		// ind_loc_all = []; ind_loc_Dop = [];
		index_a = 0;
		for (int j = gaptot; j < N_pul + gaptot; j++) {
			//cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot] -> 0 1 2 3, 5 6 7 8
			for (int k = j-gaptot,i=0; k < j - gapNum; k++,i++) {
				cellInd[i] = k;
				sum += vec[cellInd[i]];
			}
			for (int k = j+gapNum+1, i= gaptot - gapNum; k < j + gaptot + 1; k++, i++) {
				cellInd[i] = k;
				sum += vec[cellInd[i]];
			}
			noiseEst[j - gaptot] = sum;
			sum = 0;
		}

		for (int j = gaptot; j < N_pul + gaptot; j++) {
			j0 = j - gaptot +1; // ！！+1
			//cellInd=[j-gaptot: j-gapNum-1 j+gapNum+1:j+gaptot]
			for (int k = j - gaptot, i=0; k < j  - gapNum; k++, i++) {
				cellInd[i] = k;
			}
			for (int k = j + gapNum + 1, i = gaptot-gapNum; k < j + gaptot + 1; k++, i++) {
				cellInd[i] = k;
			}

			//cellInda = [j - gaptot:j - gapNum - 1]
			//cellIndb = [j + gapNum + 1:j + gaptot]
			//cellave1a =sum(vec(cellInda))/(cellNum)
			//cellave1b =sum(vec(cellIndb))/(cellNum)
			cellave1a = 0;
			cellave1b = 0;
			for (int cellInda = j - gaptot; cellInda < j - gapNum; cellInda++) {
				cellave1a += vec[cellInda];
			}
			cellave1a = cellave1a / cellNum;
			for (int cellIndb = j + gapNum + 1; cellIndb < j + gaptot + 1; cellIndb++) {
				cellave1b += vec[cellIndb];
			}
			cellave1b = cellave1b / cellNum;
			cellave1 = fmin(cellave1a, cellave1b);

			// maxInCell = max(vec(cellInd));
			maxInCell = vec[cellInd[0]];
			for (int i = 0; i < cellNum * 2; i++) {  /*cellInd[]大小: (gaptot-gapNum)*2*/
				if (maxInCell < vec[cellInd[i]]) {
					maxInCell = vec[cellInd[i]];
				}
			}

			if (maxEnable == 1) {
				condition = (vec[j]>KO*cellave1)&&(vec[j]>maxInCell);
			}
			else {
				condition = (vec[j]>KO*cellave1);
			}
			if (condition == 1) {
				//if(find(indR == j0))
				for (int i = 0; i < ind1_num; i++) {
					if (indR[i] == j0) {
						ind_win = detected_Rag_Cell_i;
						//ind_loc_all = [ind_loc_all ind_win]
						ind_loc_all[index_a] = ind_win;
						//ind_loc_Dop = [ind_loc_Dop j0]
						ind_loc_Dop[index_a] = j0;
						index_a++;
					}
				}
			}
		}

		if (index_a > 0) {  /* index_a=length(ind_loc_all) */
			for (int i = 0; i < index_a; i++) {
				ind_obj_00[i][0] = ind_loc_all[i];
				ind_obj_00[i][1] = ind_loc_Dop[i];
			}

			// if size(Ind_obj,1) ==0
			if (Ind_obj_FLAG == 0) {
				//Ind_obj = ind_obj_0
				for (int i = 0; i < index_a; i++) {
					Ind_obj[i][0] = ind_obj_00[i][0];
					Ind_obj[i][1] = ind_obj_00[i][1];
				}
				index_Ind_obj = index_a;
				Ind_obj_FLAG = 1;
			}

			else {
				//ind_obj_0_sum = ind_loc_all + 10000*ind_loc_Dop
				//Ind_obj_sum = Ind_obj(:, 1) + 10000 * Ind_obj(:, 2)
				for (int ii = 0; ii < index_a; ii++) {
					//if (length(find(Ind_obj_sum == ind_obj_0_sum(ii)))==0)
					int find_flag = 0;
					for (int j = 0; j < index_Ind_obj; j++) {
						if (Ind_obj[j][0]+10000*Ind_obj[j][1] == ind_loc_all[ii]+10000*ind_loc_Dop[ii]) {
							find_flag = 1;
							break;
						}
					}
					if (find_flag == 0) {
						//Ind_obj = [Ind_obj ; ind_obj_0(ii,:)];
						Ind_obj[index_Ind_obj][0] = ind_obj_00[ii][0];
						Ind_obj[index_Ind_obj][1] = ind_obj_00[ii][1];
						index_Ind_obj++;
					}
				}

			}
		}
	}
	free(cellInd);
	cellInd = NULL;

	//N_obj = size(Ind_obj,1)
	*N_obj = index_Ind_obj;
	cellNum = cellNum0[0];
	gapNum = gapNum0[0];
	gaptot = gapNum + cellNum;
	int N_obj_valid = 0;
	int ind_range, ind_Dop;
	double min_sig, temp_sig, mean_sig=0, obj_powerThre = 0;
	int obj_numAntenna = 192;
	int(*cellInd_1) = (int*)malloc(cellNum * 2 * sizeof(int));  /*cellNum:4->8*/
	double(*noise_obj_an)[143] = (double(*)[143])malloc(192 * index_Ind_obj * sizeof(double));  /* index_Ind_obj:143 */
	int(*Ind_obj_valid)[2] = (int(*)[2])malloc(143 * 2 * sizeof(int));  /*cellnum:4->8*/

	for (int i_obj = 0; i_obj < *N_obj; i_obj++) {
		ind_range = Ind_obj[i_obj][0]-1;  //索引起点为0
		ind_Dop = Ind_obj[i_obj][1]-1;
		//if (min(abs(sigCpml(ind_range, ind_Dop,:)).^2) < obj.powerThre)
		min_sig = pow(sigCpml[ind_range][ind_Dop][0].re, 2)+ pow(sigCpml[ind_range][ind_Dop][0].im, 2); /*初值*/
		for (int i = 1; i < 192; i++) {
			temp_sig = pow(sigCpml[ind_range][ind_Dop][i].re, 2) + pow(sigCpml[ind_range][ind_Dop][i].im, 2);
			if (min_sig > temp_sig) {
				min_sig = temp_sig;
			}
		}
		if (min_sig < obj_powerThre) {
			continue;
		}

		if (ind_range <= gaptot) {
			// cellInd=[ind_range+gapNum+1:ind_range+gaptot ind_range+gapNum+1:ind_range+gaptot] 两段一样
			for (int k = ind_range + gapNum + 1, i = 0; k < ind_range + gaptot+1; k++, i++) {
				cellInd_1[i] = k;
			}
			for (int k = ind_range + gapNum + 1, i = cellNum; k < ind_range + gaptot + 1; k++, i++) {
				cellInd_1[i] = k;
			}
		}
		else if (ind_range >= rangeNumBins - gaptot + 1) {
			// cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range-gaptot: ind_range-gapNum-1];
			for (int k = ind_range - gaptot, i = 0; k < ind_range - gapNum; k++, i++) {
				cellInd_1[i] = k;
			}
			for (int k = ind_range - gaptot, i = cellNum; k < ind_range - gapNum; k++, i++) {
				cellInd_1[i] = k;
			}
		}
		else {
			// cellInd=[ind_range-gaptot: ind_range-gapNum-1 ind_range+gapNum+1:ind_range+gaptot];
			for (int k = ind_range - gaptot, i = 0; k < ind_range - gapNum; k++, i++) {
				cellInd_1[i] = k;
			}
			for (int k = ind_range + gapNum + 1, i = cellNum; k < ind_range + gaptot+1; k++, i++) {
				cellInd_1[i] = k;
			}
		}


		//noise_obj_an(:, i_obj) = reshape((mean(abs(sigCpml(cellInd, ind_Dop, :)).^2, 1)), obj.numAntenna, 1, 1)
		for (int j = 0; j < 192; j++) { /* step1 mean: 16-1-192 -> 1-1-192; step2 reshape: 1-1-192 -> 192-1-1*/
			mean_sig = 0;
			for (int i = 0; i < cellNum * 2; i++) {
				//printf("i:%d j:%d ind_Dop:%d\n", i,j, ind_Dop);
				//printf("cellInd_1:%d\n", cellInd_1[i]);
				mean_sig += pow(sigCpml[cellInd_1[i]-1][ind_Dop][j].re, 2) + pow(sigCpml[cellInd_1[i]-1][ind_Dop][j].im, 2);
			}
			mean_sig = mean_sig / (cellNum * 2);
			noise_obj_an[j][i_obj] = mean_sig;
		}
		Ind_obj_valid[N_obj_valid][0] = Ind_obj[i_obj][0];
		Ind_obj_valid[N_obj_valid][1] = Ind_obj[i_obj][1];
		N_obj_valid++;
	}
	*N_obj = N_obj_valid;
	//Ind_obj = Ind_obj_valid;
	for (int i = 0; i < *N_obj; i++) {
		Ind_obj[i][0] = Ind_obj_valid[i][0];
		Ind_obj[i][1] = Ind_obj_valid[i][1];
	}
	//for (int k = 0; k < 143; k++) {
	//	printf("Ind_obj:%d  %d  %d\n", k + 1, Ind_obj[k][0], Ind_obj[k][1]);
	//}
	//printf("N_obj: %d", N_obj);
	//scanf("end");
	free(ind1);
	ind1 = NULL;
	free(noiseEst);
	free(indR);
	free(ind_loc_all);
	free(sig);
	free(sigv);
	free(vec);
	free(detected_Rag_Cell);
	free(ind_loc_Dop);
	free(ind_obj_00);
	free(cellInd_1);
	free(Ind_obj_valid);
	free(noise_obj_an);
	noise_obj_an = NULL;
	Ind_obj_valid = NULL;
	cellInd_1 = NULL;
	ind_loc_all = NULL;
	ind_loc_Dop = NULL;
	indR = NULL;
	noiseEst = NULL;
	sig = NULL;
	sigv = NULL;
	vec = NULL;
	detected_Rag_Cell = NULL;
	ind_obj_00 = NULL;
}


void detection_datapath(detectionResult detection_results[], COMPLEX input[256][64][192]) {
	//参数区
	double sig_integrate, abs_input;
	double const angleFFTSize = 128, angleBinSkipLeft = 4, angleBinSkipRight = 4;
	double  obj_rangeBinSize = 0.1465, obj_velocityBinSize = 0.1256;
	int obj_dopplerFFTSize = 64, obj_detectMethod = 1, obj_numAntenna = 192, obj_applyVmaxExtend = 0,
		obj_minDisApplyVmaxExtend = 10, obj_TDM_MIMO_numTX = 12, obj_numRxAnt=16;
	int not_empty_obj_overlapAntenna_ID = 1;
	double(*sig_integrate_1)[64] = (double(*)[64])malloc(sizeof(double) * 256 * 64); /*左右补数据*/
	/*CFAR_CASO_Range_2D的结果*/
	int N_obj_Rag;
	int(*Ind_obj_Rag)[2] = (int(*)[2])malloc(sizeof(int) * 44*226 * 2);
	double(*noise_obj) = (double(*))malloc(sizeof(double) * 211);
	double(*CFAR_SNR) = (double(*))malloc(sizeof(double) * 211);

	//sig_integrate = sum((abs(input)).^2,3) + 1;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 64; j++) {
			sig_integrate_1[i][j] = 0;
			for (int k = 0; k < 192; k++) {
				abs_input = sqrt(pow(input[i][j][k].re, 2) + pow(input[i][j][k].im, 2));
				sig_integrate_1[i][j] += pow(abs_input, 2);
			}
			sig_integrate_1[i][j] += 1;
		}
	}

	if (obj_detectMethod == 1) {
		// [N_obj, Ind_obj] = CFAR_CASO_Doppler_overlap(obj, Ind_obj_Rag, input, sig_integrate);
		CFAR_CASO_Range_2D(sig_integrate_1, &N_obj_Rag, Ind_obj_Rag, noise_obj, CFAR_SNR);
		//for (int i = 0; i < 211; i++) {
		//	printf("index:%d CFAR_SNR: %.5lf\n", i + 1, CFAR_SNR[i]);
		//}
		//scanf("end");
		int N_obj = 0, N_pul = 64, indx1R, indx1D, indR_num, indD_num;
		int(*ind2R) = (int*)malloc(sizeof(int) * 211);  /* cout取211 */
		int ind2D, noiseInd;
		double deltaPhi;
		if (N_obj_Rag > 0) {
			//[N_obj, Ind_obj] = CFAR_CASO_Doppler_overlap(obj, Ind_obj_Rag, input, sig_integrate);
			int(*Ind_obj)[2] = (int(*)[2])malloc(N_pul * 211 * 2 * sizeof(int)); //注意内存分配
			CFAR_CASO_Doppler_overlap(Ind_obj, &N_obj, Ind_obj_Rag, input, sig_integrate_1);
			double(*noise_obj_agg) = (double(*))malloc(sizeof(double) * N_obj);

			for (int i_obj = 0; i_obj < N_obj; i_obj++) {
				indx1R = Ind_obj[i_obj][0];
				indx1D = Ind_obj[i_obj][1];
				indR_num = 0, indD_num = 0;  // reset
				// ind2R = find(Ind_obj_Rag(:,1) == indx1R)
				// ind2D = find(Ind_obj_Rag(ind2R, 2) == indx1D)
				for (int j = 0; j < 211; j++) {
					if (Ind_obj_Rag[j][0] == indx1R) {
						ind2R[indR_num++] = j;
					}
					
				}
				for (int j = 0; j < indR_num; j++) {
					if (Ind_obj_Rag[ind2R[j]][1] == indx1D) {
						ind2D = j;
						break;
					}
				}
				// noiseInd = ind2R(ind2D)
				noiseInd = ind2R[ind2D];
				// noise_obj_agg(i_obj) = noise_obj(noiseInd)
				noise_obj_agg[i_obj] = noise_obj[noiseInd];
			}
			int xind, dopplerInd, i_sig_bin;
			double bin_val_sum=0;
			//detectionResult(*detection_results) = (detectionResult(*))malloc(sizeof(detectionResult) * N_obj);
			int(*RX_ID) = (int(*))malloc(sizeof(int) * obj_numRxAnt);
			COMPLEX(*sig_bin) = (COMPLEX(*))malloc(sizeof(COMPLEX) * obj_TDM_MIMO_numTX * obj_numRxAnt);
			COMPLEX(**bin_val_1) = (COMPLEX**)malloc(sizeof(COMPLEX*) * N_obj);
			for (int i = 0; i < N_obj; i++) {
				bin_val_1[i] = (COMPLEX*)malloc(sizeof(COMPLEX) * obj_TDM_MIMO_numTX * obj_numRxAnt);
			}
			for (int i_obj = 0; i_obj < N_obj;i_obj++) {
				xind = (Ind_obj[i_obj][0] - 1) + 1;
				detection_results[i_obj].rangeInd = Ind_obj[i_obj][0] - 1;
				detection_results[i_obj].range = detection_results[i_obj].rangeInd * obj_rangeBinSize;
				dopplerInd = Ind_obj[i_obj][1] - 1;
				detection_results[i_obj].dopplerInd_org = dopplerInd;
				detection_results[i_obj].dopplerInd = dopplerInd;
				/* velocity estimation */
				detection_results[i_obj].doppler = (double(dopplerInd) - obj_dopplerFFTSize / 2) * obj_velocityBinSize;// something wrong in Matlab
				detection_results[i_obj].doppler_corr = detection_results[i_obj].doppler;
				detection_results[i_obj].noise_var = noise_obj_agg[i_obj];
				detection_results[i_obj].bin_val = bin_val_1[i_obj];  //指向动态数组
				for (int i = 0; i < obj_numAntenna; i++) {
					detection_results[i_obj].bin_val[i]=input[xind-1][Ind_obj[i_obj][1]-1][i];  // index-1
					bin_val_sum += (pow(detection_results[i_obj].bin_val[i].im, 2) + pow(detection_results[i_obj].bin_val[i].re, 2));
				}
				//for (int i = 0; i < obj_numAntenna; i++) {
				//	printf("%d %.5lf + %.5lf i\n",i+1, detection_results[i_obj].bin_val[i].re, detection_results[i_obj].bin_val[i].im);
				//}
				//scanf("end");
				//sum(abs(detection_results (i_obj).bin_val).^2)/sum(detection_results (i_obj).noise_var)
				detection_results[i_obj].estSNR = bin_val_sum / detection_results[i_obj].noise_var; /*bug*/
				bin_val_sum = 0;  //reset
				if (obj_applyVmaxExtend == 1 && (detection_results[i_obj].range > obj_minDisApplyVmaxExtend) && not_empty_obj_overlapAntenna_ID) {
					continue;
					//skip, because obj_applyVmaxExtend is 0
				}
				else {
					/* Doppler phase correction due to TDM MIMO without apply */
					/* Vmax extention algorithm */
					deltaPhi = 2 * PI * (double(dopplerInd) - obj_dopplerFFTSize / 2) / (double(obj_TDM_MIMO_numTX) * obj_dopplerFFTSize);
					i_sig_bin=0;
					for (int i_TX = 0; i_TX < obj_TDM_MIMO_numTX; i_TX++) {
						// RX_ID = (i_TX-1)*obj.numRxAnt+1 : i_TX*obj.numRxAnt
						// sig_bin(RX_ID,: )= sig_bin_org(RX_ID )* exp(-1j*(i_TX-1)*deltaPhi)
						for (int j = 0; j < obj_numRxAnt; j++, i_sig_bin++) {
							RX_ID[j] = i_TX * obj_numRxAnt + j;
							/* (i+r)*(i+r) = (rr-ii)+(ir+ri)j */
							/* exp(-i_TX j * deltaPhi) = cos(i_TX*deltaPhi)-sin(i_TX*deltaPhi)j */
							sig_bin[i_sig_bin].im = detection_results[i_obj].bin_val[RX_ID[j]].im * cos(i_TX * deltaPhi)
								- detection_results[i_obj].bin_val[RX_ID[j]].re * sin(i_TX * deltaPhi);
							sig_bin[i_sig_bin].re = detection_results[i_obj].bin_val[RX_ID[j]].re * cos(i_TX * deltaPhi)
								+ detection_results[i_obj].bin_val[RX_ID[j]].im * sin(i_TX * deltaPhi);

							detection_results[i_obj].bin_val[i_sig_bin] = sig_bin[i_sig_bin];
						}
					}
					//for (int i = 0; i < obj_numAntenna; i++) {
					//	printf("%d %.5lf + %.5lf i\n", i+1, sig_bin[i].re, sig_bin[i].im);
					//}
					//scanf("end");
					//detection_results[i_obj].bin_val = sig_bin;
					detection_results[i_obj].doppler_corr_overlap = detection_results[i_obj].doppler_corr;
					detection_results[i_obj].doppler_corr_FFT = detection_results[i_obj].doppler_corr;
				}
			}
			//for (int i = 0; i < 192; i++) {
			//	//printf("index:%d %d %.5lf %.5lf\n", i + 1, detection_results[i].dopplerInd, detection_results[i].doppler, detection_results[i].estSNR);
			//	printf("i: %d  %.5lf + %.5lf i\n", i+1, detection_results[2].bin_val[i].re, detection_results[2].bin_val[i].im);
			//}
			//scanf("end");
			free(RX_ID);
			free(sig_bin);
			free(Ind_obj);
			free(bin_val_1);
			free(noise_obj_agg);
			noise_obj_agg = NULL;
			bin_val_1 = NULL;
			Ind_obj = NULL;
			sig_bin = NULL;
			RX_ID = NULL;
		}
		free(ind2R);
		ind2R = NULL;
	}
	free(Ind_obj_Rag);
	free(sig_integrate_1);
	free(noise_obj);
	free(CFAR_SNR);
	CFAR_SNR = NULL;
	noise_obj = NULL;
	Ind_obj_Rag = NULL;
	sig_integrate_1 = NULL;
}


void DOA_beamformingFFT_2D(COMPLEX *sig) {
	int angles_DOA_az[2] = { -80, 80 }, angles_DOA_ele[2] = { -20, 20 };
	double d = 0.5046;
	int angleFFTSize = 256, apertureLen_azim, apertureLen_elev;
	int(*D)[2] = (int(*)[2])malloc(sizeof(int) * 192 * 2);  //怎么生成D？
	int D1[192] = { 0,1,2,3,11,12,13,14,46,47,48,49,50,51,52,53,4,5,6,7,15,16,17,18,50,51,52,53,54,55,56,57,8,9,10,11,19,20,21,22,54,55,56,57,58,59,60,61,12,13,14,15,23,24,25,26,58,59,60,61,62,63,64,65,16,17,18,19,27,28,29,30,62,63,64,65,66,67,68,69,20,21,22,23,31,32,33,34,66,67,68,69,70,71,72,73,24,25,26,27,35,36,37,38,70,71,72,73,74,75,76,77,28,29,30,31,39,40,41,42,74,75,76,77,78,79,80,81,32,33,34,35,43,44,45,46,78,79,80,81,82,83,84,85,9,10,11,12,20,21,22,23,55,56,57,58,59,60,61,62,10,11,12,13,21,22,23,24,56,57,58,59,60,61,62,63,11,12,13,14,22,23,24,25,57,58,59,60,61,62,63,64 };
	
	COMPLEX(*sig_sel) = (COMPLEX*)malloc(sizeof(COMPLEX) * 192);
	COMPLEX(*sig_2D)[7] = (COMPLEX(*)[7])malloc(sizeof(COMPLEX) * 86 * 7);
	COMPLEX(*fftOutput)[7] = (COMPLEX(*)[7])malloc(sizeof(COMPLEX) * 256 * 7);
	COMPLEX(*angle_sepc_1D_fft)[7] = (COMPLEX(*)[7])malloc(sizeof(COMPLEX) * 256 * 7);
	COMPLEX(*fftOutput1)[256] = (COMPLEX(*)[256])malloc(sizeof(COMPLEX) * 256 * 256);
	COMPLEX(*angle_sepc_2D_fft)[256] = (COMPLEX(*)[256])malloc(sizeof(COMPLEX) * 256 * 256);
	for (int i = 0; i < 144; i++) {
		D[i][1] = 0 + 1;
	}
	for (int i = 0; i < 16; i++) {
		D[i + 144][1] = 1 + 1;
		D[i + 144 + 16][1] = 4 + 1;
		D[i + 144 + 32][1] = 6 + 1;
	}
	for (int i = 0; i < 192; i++) {
		D[i][0] = D1[i] + 1;
	}
	//apertureLen_azim = max(D(:,1))
	apertureLen_azim = D[0][0];
	for (int i = 1; i < 192; i++) {
		if (D[i][0] > apertureLen_azim) {
			apertureLen_azim = D[i][0];
		}
	}
	//apertureLen_elev = max(D(:,2))
	apertureLen_elev = D[0][1];
	for (int i = 1; i < 192; i++) {
		if (D[i][1] > apertureLen_elev) {
			apertureLen_elev = D[i][1];
		}
	}
	int ind_num = 0, indU_num=1, flag=1, temp, change_flag=0;
	for (int i_line = 0; i_line < apertureLen_elev; i_line++) {
		//初始化
		int(*ind) = (int*)malloc(sizeof(int) * 192);
		int(*indU) = (int*)malloc(sizeof(int) * 192);
		int(*val) = (int*)malloc(sizeof(int) * 192);
		int(*D_sel) = (int*)malloc(sizeof(int) * 192);
		indU_num = 1, ind_num = 0;
		// ind = find(D(:,2) == i_line);
		for (int j = 0; j < 192; j++) {
			if (D[j][1] == i_line+1) {  //i_line+1, 起始值+1
				ind[ind_num++] = j;
			}
		}
		// D_sel = D(ind,1)
		for (int j = 0; j < ind_num; j++) {
			D_sel[j] = D[ind[j]][0];
			sig_sel[j] = sig[ind[j]];
		}
		// [val indU] = unique(D_sel), 取独, then排序
		indU[0] = 1;
		val[0] = D_sel[0];  /*初值*/
		for (int j = 0; j < ind_num; j++) {
			flag = 1;
			for (int t = 0; t < indU_num; t++) {
				if (val[t] == D_sel[j]) {
					flag = 0;  //已有值
					break;
				}
			}
			if (flag == 1) {  //新值
				indU[indU_num] = j;
				val[indU_num++] = D_sel[j];
			}
		}
		//冒泡排序
		while (change_flag == 0) {
			change_flag = 1;
			for (int i = 0; i < indU_num - 1; i++) {
				if (val[i] > val[i + 1]) { //exchange
					temp = val[i];
					val[i] = val[i + 1];
					val[i + 1] = temp;
					temp = indU[i];
					indU[i] = indU[i + 1];
					indU[i + 1] = temp;
					change_flag = 0;
				}
			}
		}
		//  sig_2D(D_sel(indU),i_line) = sig_sel(indU)
		for (int j = 0; j < 86; j++) {  /* 初始化sig_2D */
			sig_2D[j][i_line].im = 0;
			sig_2D[j][i_line].re = 0;
		}
		if (ind_num > 0) {
			for (int j = 0; j < indU_num; j++) {
				sig_2D[D_sel[indU[j]] - 1][i_line] = sig_sel[indU[j]];
			}
		}
		//angle_sepc_1D_fft=fftshift(fft(sig_2D,angleFFTSize,1),1)
		double pr[256], pi[256], fr[256], fi[256];
		for (int j = 0; j < 7; j++) {
			for (int i = 0; i < 86; i++) {
				pr[i] = sig_2D[i][j].re;
				pi[i] = sig_2D[i][j].im;
			}
			for (int i = 86; i < 256; i++) {  // 补零 86->256
				pr[i] = 0;
				pi[i] = 0;
			}
			//每列进行fft
			kfft(pr, pi, 256, 8, fr, fi);  // pr,pi -> input ; fr,fi -> output ; 采样点数 256=2~8
			for (int i = 0; i < 256; i++) { 
				fftOutput[i][j].re = fr[i];
				fftOutput[i][j].im = fi[i];
			}
		}
		fftshift_2_1(fftOutput, angle_sepc_1D_fft);
		//angle_sepc_2D_fft=fftshift(fft(angle_sepc_1D_fft,angleFFTSize,2),2); 
		double pr1[256], pi1[256], fr1[256], fi1[256];
		for (int j = 0; j < 256; j++) {
			for (int i = 0; i < 7; i++) {
				pr1[i] = angle_sepc_1D_fft[j][i].re;
				pi1[i] = angle_sepc_1D_fft[j][i].im;
			}
			for (int i = 7; i < 256; i++) {  // 补零 7->256
				pr1[i] = 0;
				pi1[i] = 0;
			}
			//每行进行fft
			kfft(pr1, pr1, 256, 8, fr1, fi1);
			for (int i = 0; i < 256; i++) {  //
				fftOutput1[j][i].re = fr1[i];
				fftOutput1[j][i].im = fi1[i];
			}
		}

		free(D_sel);
		free(ind);
		free(indU);
		free(val);
		val = NULL;
		indU = NULL;
		ind = NULL;
		D_sel = NULL;
	}
	//for (int i = 0; i < 86; i++) {
	//	printf("%d  re: %.5lf, im: %.5lf\n", i+1, sig_2D[i][2].re, sig_2D[i][2].im);
	//}
	//scanf("end");

	free(D);
	free(sig_sel);
	free(fftOutput);
	free(angle_sepc_1D_fft);
	angle_sepc_1D_fft = NULL;
	fftOutput = NULL;
	sig_sel = NULL;
	D = NULL;
}


void DOA_path(detectionResult detected_obj[], detectionResult out[]) {
	out = detected_obj;
	int numAoAObjCnt = 0, obj_method = 1;
	double estSNR, sum=0;
	detectionResult current_obj;
	COMPLEX(*X)= (COMPLEX*)malloc(sizeof(COMPLEX) * 192);
	COMPLEX(*R)[192] = (COMPLEX(*)[192])malloc(sizeof(COMPLEX) * 192 * 192);
	if (X != NULL && R != NULL) {
		for (int i_obj = 0; i_obj < 143; i_obj++) {
			current_obj = detected_obj[i_obj];
			// estSNR = 10*log10(sum(abs(current_obj.bin_val).^2)/sum(current_obj.noise_var));
			for (int i = 0; i < 192; i++) {
				sum += (pow(current_obj.bin_val[i].im, 2) + pow(current_obj.bin_val[i].re, 2));
			}
			estSNR = 10 * log10(sum / current_obj.noise_var);
			X = current_obj.bin_val;
			//R = X*X'
			for (int i = 0; i < 192; i++) {
				for (int j = 0; j < 192; j++) {
					R[i][j].re = X[i].re * X[j].re + X[i].im * X[j].im;  // X'的虚部取反??
					R[i][j].im = X[i].im * X[j].re - X[i].re * X[j].im;
					//printf("i%d j%d: %.5lf + %.5lfi\n",i+1, j+1, R[i][j].re, R[i][j].im);
				}
			}
			//scanf("end");
			switch (obj_method)
			{
			case 1:
				DOA_beamformingFFT_2D(X);
				break;
			default:
				break;
			}
		}
	}
	free(X);
	free(R);
	R = NULL;
	X = NULL;
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
	double begin_del = 5, end_del = 10;
	//工作区
	while (!feof(fidlist)) {
		fgets(dataFolder_test, 1000, fidlist);
		fgets(dataFolder_calib, 1000, fidlist);
		fgets(module_param_file, 1000, fidlist);
		printf("%s\n%s\n%s\n", dataFolder_test, dataFolder_calib, module_param_file);
		sprintf_s(pathGenParaFile, "%s%s%s%s", input_path, "test", testID, "_param.m");
	}
	fclose(fidlist);
	if (PARAM_FILE_GEN_ON == 1) {
		int parameter_file_gen_json;  /*不明函数*/
	}
	COMPLEX(*adcData)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*输出；256-64-16-12*/
	calibration_datapath(&calibrationObj, adcData);

	COMPLEX(*adcData_slice)[64][16] = (COMPLEX(*)[64][16])malloc(262144 * sizeof(COMPLEX));  /*256-64-16*/
	COMPLEX(*rangeFFTOut)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*256-64-16-12*/
	COMPLEX(*DopplerFFTOut_pre)[64][16][12] = (COMPLEX(*)[64][16][12])malloc(262144 * 12 * sizeof(COMPLEX));  /*256-64-16-12*/
	COMPLEX(*DopplerFFTOut)[64][16 * 12] = (COMPLEX(*)[64][16 * 12])malloc(262144 * 12 * sizeof(COMPLEX));  /*256-64-16x12*/
	COMPLEX(*rangeFFTOut_slice)[64][16] = (COMPLEX(*)[64][16])malloc(sizeof(COMPLEX) * 256 * 64 * 16);
	COMPLEX(*DopplerFFTOut_slice)[64][16] = (COMPLEX(*)[64][16])malloc(sizeof(COMPLEX) * 256 * 64 * 16);
	for (int i_tx = 0; i_tx < 12; i_tx++) {
		//adcData切片
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 64; j++) {
				for (int m = 0; m < 16; m++) {
					adcData_slice[i][j][m] = adcData[i][j][m][i_tx];
				}
			}
		}
		rangeProcCascade_datapath(adcData_slice, rangeFFTOut_slice);
		//rangeFFTOut聚合
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 64; j++) {
				for (int m = 0; m < 16; m++) {
					rangeFFTOut[i][j][m][i_tx] = rangeFFTOut_slice[i][j][m];
				}
			}
		}
		dopplerProcClutterRemove_datapath(rangeFFTOut_slice, DopplerFFTOut_slice);  //地址传参
		//DopplerFFTOut聚合
		for (int i = 0; i < 256; i++) {
			for (int j = 0; j < 64; j++) {
				for (int m = 0; m < 16; m++) {
					DopplerFFTOut_pre[i][j][m][i_tx] = DopplerFFTOut_slice[i][j][m];
				}
			}
		}
		//for (int i = 0; i < 5; i++) {
		//	printf("Doppler re: %.5lf, im: %.5lf\n", DopplerFFTOut[255][63][i][0].re, DopplerFFTOut[255][63][i][0].im);
		//}
		//scanf("end");
	}
	//DopplerFFTOut = reshape(DopplerFFTOut,size(DopplerFFTOut,1), size(DopplerFFTOut,2), size(DopplerFFTOut,3)*size(DopplerFFTOut,4));
	int k;
	for (int i = 0; i < 256; i++) {
		for (int j = 0; j < 64; j++) {
			k = 0;
			for (int p = 0; p < 12; p++) {
				for (int q = 0; q < 16; q++) {
					DopplerFFTOut[i][j][k] = DopplerFFTOut_pre[i][j][q][p];
					k++;
				}
			}
		}
	}
	//sig_integrate = sum((abs(DopplerFFTOut)).^2,3)
	//sig_integrate_LOG = 10*log10(sig_integrate)
	//sig_integrate = sig_integrate(begin_del+1:size(DopplerFFTOut,1)-end_del,:)
	//sig_integrate_LOG = sig_integrate_LOG(begin_del + 1:size(DopplerFFTOut, 1) - end_del,:)
	double(*sig_integrate)[64] = (double(*)[64])malloc(241 * 64 * sizeof(double));  /*241-64*/
	double(*sig_integrate_LOG)[64] = (double(*)[64])malloc(241 * 64 * sizeof(double));  /*241-64*/
	double out_abs;
	for (int i = 5; i < 246; i++) {
		for (int j = 0; j < 64; j++) {
			sig_integrate[i - 5][j] = 0; /*初始化*/
			sig_integrate_LOG[i - 5][j] = 0;
			for (int n = 0; n < 192; n++) {
				out_abs = sqrt(pow(DopplerFFTOut[i][j][n].re, 2) + pow(DopplerFFTOut[i][j][n].im, 2));
				sig_integrate[i - 5][j] += pow(out_abs, 2);
			}
			sig_integrate_LOG[i - 5][j] = 10 * log10(sig_integrate[i - 5][j]);
		}
	}
	//for (int i = 0; i < 10; i++) {
	//	printf("DopplerFFTOut: %.5lf\n", DopplerFFTOut[8][0][0].re);

	//}
	//scanf("end");
	double delta_r = 0.1465;
	double delta_v = 0.1256;
	double range[241], vecolity[64];
	for (int i = 0; i < 241; i++) {
		range[i] = double(i + 5) * delta_r;
	}
	for (int i = 0; i < 64; i++) {
		vecolity[i] = (i + 1 - (64 / 2 + 1)) * delta_v;  /*略有误差*/
	}
	//RD图: Range-Doppler
	writetxt(vecolity, 64, "result/vecolity.txt");
	writetxt(range, 241, "result/range.txt");
	writetxt_2d(sig_integrate, 241, "result/sig_integrate.txt");
	writetxt_2d(sig_integrate_LOG, 241, "result/sig_integrate_LOG.txt");

	//CFAR检测，检测基于192个通道求和数据，数据经过了多普勒维相干积累和多通道非相干积累
	detectionResult(*detection_results) = (detectionResult(*))malloc(sizeof(detectionResult) * 143);
	detection_datapath(detection_results, DopplerFFTOut);
	//for (int i = 0; i < 143; i++) {
	//	printf("index:%d %d %.5lf %.5lf\n", i + 1, detection_results[i].dopplerInd, detection_results[i].doppler, detection_results[i].estSNR);
	//	//printf("i: %d  %.5lf + %.5lf i\n", i+1, detection_results[2].bin_val[i].re, detection_results[2].bin_val[i].im);
	//}
	//scanf("end");
	double(*detect_all_points)[4] = (double(*)[4])malloc(sizeof(double) * 143 * 4);
	for (int iobj = 0; iobj < 143; iobj++) {
		detect_all_points[iobj][0] = detection_results[iobj].rangeInd + 1;
		detect_all_points[iobj][1] = detection_results[iobj].dopplerInd_org + 1;
		detect_all_points[iobj][3] = detection_results[iobj].estSNR;
	}
	//is detection_results empty?
	detectionResult(*angleEst) = (detectionResult(*))malloc(sizeof(detectionResult) * 143);
	DOA_path(detection_results, angleEst);

	free(detection_results);
	free(adcData);
	free(adcData_slice);
	free(rangeFFTOut);
	free(DopplerFFTOut_pre);
	free(sig_integrate);
	free(sig_integrate_LOG);
	free(rangeFFTOut_slice);
	free(DopplerFFTOut_slice);
	detection_results = NULL;
	adcData = NULL;
	adcData_slice = NULL;
	rangeFFTOut = NULL;
	DopplerFFTOut_pre = NULL;
	sig_integrate = NULL;
	sig_integrate_LOG = NULL;
	rangeFFTOut_slice = NULL;
	DopplerFFTOut_slice = NULL;
	return 0;
}

