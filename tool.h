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

auto readBinFile(const char* fileFullPath, int frameIdx, int numSamplePerChirp, int numChirpPerLoop, int numLoops, int numRXPerDevice, int numDevices);
auto read_ADC_bin_TDA2_separateFiles();
auto datapath(CalibrationObj* obj);
auto rangeProcCascade_datapath(COMPLEX input[256][64][16]);
auto DopplerProcClutterRemove_datapath();
#pragma once
