/* LaBr3:Ce SORTING code for TDR format dataframe */
/* Produces ROOT TTrees and histograms for four LaBr3:Ce detectors: L0,L1,L2,L3*/

/* Adapted from P.M Jones code labrsort4 */

/* (c) S. Hart */
/* (10/06/22)  */


/* COMPILE: g++ -std=c++0x labrsortincode.C -o sortingexe `root-config --cflags --libs` -lSpectrum */
/* RUN: ./sortingexe /path/to/raw/experiment/data/ (e.g. /home/shanyn/Data/R50) runXX (PASS THE SPECIFIC NAME OF THE RUN AS ARGV[2] IN ORDER TO CREATE THE CSV FILE WITH THE CORRECT NAME)*/

/* NOTES: First order calibration seems sufficient */

#include<iostream> //INPUT OUTPUT ON SCRREN
#include<string>  //DYNAMIC ARRA#include<iostream> //INPUT OUTPUT ON SCRREN
#include<vector>  //DYNAMIC ARRAYS
#include<string>  //DYNAMIC ARRAYS OF CHAR
#include<dirent.h>
#include<fstream> //INPUT OUTPUT FILES
#include<cmath>
#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <initializer_list>
#include <iomanip>
#include <sstream>
#include <libgen.h>

// ROOT libraries
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TPad.h"
#include "TChain.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH3D.h"
#include "TObject.h"
#include "TString.h"
#include "TRandom.h"
#include "TH2F.h"



#define SIZE  16384

int32_t buffer[SIZE];

uint64_t TStop, TSbot;
uint64_t dataframe;
uint64_t TS, TSlast, TSfirst, counter,  SYNC, SYNClast=0;
double T1, T2;
int64_t TSdiff, SYNCdiff;
int64_t count=0;
uint16_t energy;
unsigned long int pos = 0;
int verbose = 255;

using namespace std;

int main (int argc, char** argv)
{
  FILE *f;
  FILE *csvf;
  int i, n;
  char file_in[50];
  char file_out[50];
  char csv_out[50];
  int  identity, card=0, adcdata;
  double  slowLaBr0, TimeStampLaBr0, slowtimeLaBr0, slowLaBr0calib1stFullRange, slowLaBr0calib1stLowRange, slowLaBr0calib2ndFullRange;
  double  slowLaBr1, TimeStampLaBr1, slowtimeLaBr1, slowLaBr1calib1stFullRange, slowLaBr1calib1stLowRange, slowLaBr1calib2ndFullRange;
  double  slowLaBr2, TimeStampLaBr2, slowtimeLaBr2, slowLaBr2calib1stFullRange, slowLaBr2calib1stLowRange, slowLaBr2calib2ndFullRange;
  double  slowLaBr3, TimeStampLaBr3, slowtimeLaBr3, slowLaBr3calib1stFullRange, slowLaBr3calib1stLowRange, slowLaBr3calib2ndFullRange;
  double  energy0firstFullRange, energy1firstFullRange, energy2firstFullRange, energy3firstFullRange;
  double  energy0firstLowRange, energy1firstLowRange, energy2firstLowRange, energy3firstLowRange;
  double  energy0secondFullRange, energy1secondFullRange, energy2secondFullRange, energy3secondFullRange;
  double RFsignaltime, TimeStampRF, RFenergy;
  double TimeStampPOLARIS, SlowTimePOLARIS, SlowEnergyPOLARIS;
  double rndUnif, rndchannel;
  int tick, tock;
  double sick;
  int twidset=0;
  double Td;
  int RunEnd=999999, Run=0;
  int no_read;
  std::vector<int> time_stamp;
  
  if ((argc < 2) || (argc > 3))
  {
	fprintf(stderr, "Usage: %s dataframe [last_subrun]\n", argv[0]);
	exit(1);
  }
  
  /* ROOT STUFF */
  sprintf(file_out, "%s_rawData.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);

  //first order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a1stFullRange0 = 1.0069, b1stFullRange0 = -52.8029; //L0
  double a1stFullRange1 = 1.0848, b1stFullRange1 = -56.7355; //L1
  double a1stFullRange2 = 0.5532, b1stFullRange2 = -16.1329; //L2
  double a1stFullRange3 = 0.789, b1stFullRange3 = -18.1978; //L3

 //first order polynomial //date: 15/6/22 for 152Eu Source
  double a1stLowRange0 = 0.7517, b1stLowRange0 = -0.1179; //L0
  double a1stLowRange1 = 0.8113, b1stLowRange1 = -0.5323; //L1
  double a1stLowRange2 = 0.5174, b1stLowRange2 = -0.4377; //L2
  double a1stLowRange3 = 0.7139, b1stLowRange3 = 0.588; //L3

  //2nd order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a2ndFullRange0 = 0.0004, b2ndFullRange0 = 0.3737, c2ndFullRange0 = 55.9837; //L0
  double a2ndFullRange1 = 0.0004, b2ndFullRange1 = 0.4026, c2ndFullRange1 = 57.6084; //L1
  double a2ndFullRange2 = 0.0001, b2ndFullRange2 = 0.2726, c2ndFullRange2 = 60.9222; //L2
  double a2ndFullRange3 = 0.0002, b2ndFullRange3 = 0.4864, c2ndFullRange3 = 35.95; //L3

  
  TTree *LaBr0Data = new TTree("LaBr0Data", "LaBr0Data"); 
  LaBr0Data->Branch("EnergySlow_LaBr_Det0_2ndOrder_Calib_FullRange", &slowLaBr0calib2ndFullRange , "slowLaBr0calib2ndFullRange/D");
  LaBr0Data->Branch("EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange", &slowLaBr0calib1stFullRange , "slowLaBr0calib1stFullRange/D");
  LaBr0Data->Branch("EnergySlow_LaBr_Det0_1stOrder_Calib_LowRange", &slowLaBr0calib1stLowRange , "slowLaBr0calib1stLowRange/D");
  LaBr0Data->Branch("EnergySlow_LaBr_Det0", &slowLaBr0 , "slowLaBr0/D");
  LaBr0Data->Branch("TimingSlow_LaBr_Det0", &slowtimeLaBr0 , "slowtimeLaBr0/D");
  LaBr0Data->Branch("TimeStamp_LaBr_Det0", &TimeStampLaBr0 , "TimeStampLaBr0/D");
  LaBr0Data->Branch("TimingSlow_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr0Data->Branch("TimeStamp_RF", &TimeStampRF , "TimeStampRF/D");
  LaBr0Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr0Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");
  LaBr0Data->Branch("TimingPOLARIS", &SlowTimePOLARIS, "SlowTimePOLARIS/D");

  TTree *LaBr1Data = new TTree("LaBr1Data", "LaBr1Data"); 
  LaBr1Data->Branch("EnergySlow_LaBr_Det1_2ndOrder_Calib_FullRange", &slowLaBr1calib2ndFullRange , "slowLaBr1calib2ndFullRange/D");
  LaBr1Data->Branch("EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange", &slowLaBr1calib1stFullRange , "slowLaBr1calib1stFullRange/D");
  LaBr1Data->Branch("EnergySlow_LaBr_Det1_1stOrder_Calib_LowRange", &slowLaBr1calib1stLowRange , "slowLaBr1calib1stLowRange/D");
  LaBr1Data->Branch("EnergySlow_LaBr_Det1", &slowLaBr1 , "slowLaBr1/D");
  LaBr1Data->Branch("TimingSlow_LaBr1", &slowtimeLaBr1 , "slowtimeLaBr1/D");
  LaBr1Data->Branch("TimeStamp_LaBr_Det1", &TimeStampLaBr1 , "TimeStampLaBr1/D");
  LaBr1Data->Branch("TimingSlow_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr1Data->Branch("TimeStamp_RF", &TimeStampRF , "TimeStampRF/D");
  LaBr1Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr1Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");
  LaBr1Data->Branch("TimingPOLARIS", &SlowTimePOLARIS, "SlowTimePOLARIS/D");

  TTree *LaBr2Data = new TTree("LaBr2Data", "LaBr2Data"); 
  LaBr2Data->Branch("EnergySlow_LaBr_Det2_2ndOrder_Calib_FullRange", &slowLaBr2calib2ndFullRange , "slowLaBr2calib2ndFullRange/D");
  LaBr2Data->Branch("EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange", &slowLaBr2calib1stFullRange , "slowLaBr2calib1stFullRange/D");
  LaBr2Data->Branch("EnergySlow_LaBr_Det2_1stOrder_Calib_LowRange", &slowLaBr2calib1stLowRange , "slowLaBr2calib1stLowRange/D");
  LaBr2Data->Branch("EnergySlow_LaBr_Det2", &slowLaBr2 , "slowLaBr2/D");
  LaBr2Data->Branch("TimingSlow_LaBr2", &slowtimeLaBr2 , "slowtimeLaBr2/D");
  LaBr2Data->Branch("TimeStamp_LaBr_Det2", &TimeStampLaBr2 , "TimeStampLaBr2/D");
  LaBr2Data->Branch("TimingSlow_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr2Data->Branch("TimeStamp_RF", &TimeStampRF , "TimeStampRF/D");
  LaBr2Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr2Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");
  LaBr2Data->Branch("TimingPOLARIS", &SlowTimePOLARIS, "SlowTimePOLARIS/D");

  TTree *LaBr3Data = new TTree("LaBr3Data", "LaBr3Data"); 
  LaBr3Data->Branch("EnergySlow_LaBr_Det3_2ndOrder_Calib_FullRange", &slowLaBr3calib2ndFullRange , "slowLaBr3calib2ndFullRange/D");
  LaBr3Data->Branch("EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange", &slowLaBr3calib1stFullRange , "slowLaBr3calib1stFullRange/D");
  LaBr3Data->Branch("EnergySlow_LaBr_Det3_1stOrder_Calib_LowRange", &slowLaBr3calib1stLowRange , "slowLaBr3calib1stLowRange/D");
  LaBr3Data->Branch("EnergySlow_LaBr_Det3", &slowLaBr3 , "slowLaBr3/D");
  LaBr3Data->Branch("TimingSlow_LaBr_Det3", &slowtimeLaBr3 , "slowtimeLaBr3/D");
  LaBr3Data->Branch("TimeStamp_LaBr_Det3", &TimeStampLaBr3 , "TimeStampLaBr3/D");
  LaBr3Data->Branch("TimingSlow_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr3Data->Branch("TimeStamp_RF", &TimeStampRF , "TimeStampRF/D");
  LaBr3Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr3Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");
  LaBr3Data->Branch("TimingPOLARIS", &SlowTimePOLARIS, "SlowTimePOLARIS/D");

  // Uncalibrated Histograms
  TH1D *s0none = new TH1D("L0_Uncalibrated","LaBr3 Det 0 Slow Signal Energy Spectrum - Uncalibrated",8000,0,8000);
  TH1D *s1none = new TH1D("L1_Uncalibrated","LaBr3 Det 1 Slow Signal Energy Spectrum - Uncalibrated",8000,0,8000);
  TH1D *s2none = new TH1D("L2_Uncalibrated","LaBr3 Det 2 Slow Signal Energy Spectrum - Uncalibrated",8000,0,8000);
  TH1D *s3none = new TH1D("L3_Uncalibrated","LaBr3 Det 3 Slow Signal Energy Spectrum - Uncalibrated",8000,0,8000);

  // 1st order calibrations FULL RANGE
  TH1D* s0calib1stFullRange = new TH1D("L0_Calibrated_1st_FullRange", "LaBr3 Det 0 Slow Signal Energy Spectrum - Calibrated 1st Order - Full Range", 8000, 0, 8000);
  TH1D* s1calib1stFullRange = new TH1D("L1_Calibrated_1st_FullRange", "LaBr3 Det 1 Slow Signal Energy Spectrum - Calibrated 1st Order - Full Range", 8000, 0, 8000);
  TH1D* s2calib1stFullRange = new TH1D("L2_Calibrated_1st_FullRange", "LaBr3 Det 2 Slow Signal Energy Spectrum - Calibrated 1st Order - Full Range", 8000, 0, 8000);
  TH1D* s3calib1stFullRange = new TH1D("L3_Calibrated_1st_FullRange", "LaBr3 Det 3 Slow Signal Energy Spectrum - Calibrated 1st Order - Full Range", 8000, 0, 8000);

  // 1st order calibrations LOW RANGE
  TH1D* s0calib1stLowRange = new TH1D("L0_Calibrated_1st_LowRange", "LaBr3 Det 0 Slow Signal Energy Spectrum - Calibrated 1st Order - Low Range", 8000, 0, 8000);
  TH1D* s1calib1stLowRange = new TH1D("L1_Calibrated_1st_LowRange", "LaBr3 Det 1 Slow Signal Energy Spectrum - Calibrated 1st Order - Low Range", 8000, 0, 8000);
  TH1D* s2calib1stLowRange = new TH1D("L2_Calibrated_1st_LowRange", "LaBr3 Det 2 Slow Signal Energy Spectrum - Calibrated 1st Order - Low Range", 8000, 0, 8000);
  TH1D* s3calib1stLowRange = new TH1D("L3_Calibrated_1st_LowRange", "LaBr3 Det 3 Slow Signal Energy Spectrum - Calibrated 1st Order - Low Range", 8000, 0, 8000);

  // 2nd order calibrations FULL RANGE
  TH1D* s0calib2ndFullRange = new TH1D("L0_Calibrated_2nd_FullRange", "LaBr3 Det 0 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s1calib2ndFullRange = new TH1D("L1_Calibrated_2nd_FullRange", "LaBr3 Det 1 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s2calib2ndFullRange = new TH1D("L2_Calibrated_2nd_FullRange", "LaBr3 Det 2 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s3calib2ndFullRange = new TH1D("L3_Calibrated_2nd_FullRange", "LaBr3 Det 3 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);

  // Timing of slow energy signals
  TH1D *timeSlowSig0 = new TH1D("L0_TimeSlow","LaBr3 Det 0 Slow Signal Timing",8000,9e7,9e8);
  TH1D *timeSlowSig1 = new TH1D("L1_TimeSlow","LaBr3 Det 1 Slow Signal Timing",8000,9e7,9e8);
  TH1D *timeSlowSig2 = new TH1D("L2_TimeSlow","LaBr3 Det 2 Slow Signal Timing",8000,9e7,9e8);
  TH1D *timeSlowSig3 = new TH1D("L3_TimeSlow","LaBr3 Det 3 Slow Signal Timing",8000,9e7,9e8);
  
  // Timestamps
  TH1D *Timestamp0 = new TH1D("L0_Timestamp","LaBr3 Det 0 Timestamp",8000,9e7,9e8);
  TH1D *Timestamp1 = new TH1D("L1_Timestamp","LaBr3 Det 1 Timestamp",8000,9e7,9e8);
  TH1D *Timestamp2 = new TH1D("L2_Timestamp","LaBr3 Det 2 Timestamp",8000,9e7,9e8);
  TH1D *Timestamp3 = new TH1D("L3_Timestamp","LaBr3 Det 3 Timestamp",8000,9e7,9e8);

  // POLARIS ENERGY VS TIMESTAMP
  TH2D *energyVtimestampPOLARIS = new TH2D("E_V_Time_POLARIS", "Signal Energy Vs Timestamp - POLARIS", 500,0,250,500,5.5e13,5.7e13);

  // POLARIS ENERGY VS SLOW TIME
  TH2D *energyVslowsigtimePOLARIS = new TH2D("E_V_SlowTime_POLARIS", "Signal Energy Vs Slow Signal Time - POLARIS", 500,0,250,500,5.5e15,5.7e15);

  // POLARIS
  TH1D *POLARISEnergy = new TH1D("POLARIS_Energy","POLARIS Energy",500,0,250);
  TH1D *POLARISTimeStamp = new TH1D("POLARIS_TimeStamp","POLARIS TimeStamp",500,5.5e13,5.7e13);
  TH1D *POLARISSlowTime = new TH1D("POLARIS_SlowTime","POLARIS SlowTime",500,5.5e15,5.7e15);

  // RF Time Walk Correction Information
  TH1D *RFSigTime = new TH1D("RF_Time","RF Signal Time Stamp",8000,9e7,9e8);
  TH1D *RFSigEnergy = new TH1D("RF_Energy","RF Signal Energy",8000,0,6000);
  TH1D *RFTimeStamp = new TH1D("RF_TimeStamp","RF Time Stamp",8000,9e7,9e8);

  // use TRandom to create a random number 
  TRandom* randy = new TRandom();
  double bw = 1;

  // File stuff
  while(Run < RunEnd+1)
  {

      sprintf(file_in, "/%s_%d", argv[1], Run);
      
      f = fopen(file_in, "rb");
      if (!f)
    {
    fprintf(stderr, "Can't open file '%s'\n", file_in);
    goto finish;
    }
      
    printf("File %s opened...\n", file_in);

    
    
    while (!feof(f))
    {

      no_read = fread(buffer, sizeof(buffer[0]), SIZE, f);
    
      if ( (feof(f)) && no_read <= 0 ) goto end;

      pos+=24;

      for (i=6; i< SIZE; i+=2)
      {
        dataframe = buffer[i+1];
        TSbot = buffer[i];

        if ((TSbot == 0) && (dataframe == 0)) goto loop;
        if ((TSbot == 0x5e5e5e5e) && (dataframe == 0x5e5e5e5e)) goto loop;
        if ((TSbot == 0x5e5e5e5e) && (dataframe == 0xffffffff)) goto loop;


        /* dataframe */
        if ((dataframe & 0xc0000000) == 0xc0000000)  
        {
          identity = (dataframe & 0x0fff0000) >> 16;
          adcdata = (dataframe & 0x0000ffff);
          card = (identity / 32);
          TS = (TS & 0x0000fffff0000000ULL);

           
          TS = ((TS | TSbot));
          if (counter == 0)
          {
            TSfirst = TS;
            counter = 1;
          }
          TSdiff = TS-TSlast;

          if ( (TSbot >= 0x0) && (TSbot <= 0x1a) )
          {
            if (twidset == 0) 
            {
              TS=TS+0x10000000;
              twidset=1;
            }

            TSdiff = TS-TSlast;
          }

          //printf(" === i = %d / ident: %d / energy: %d / TS: %ld --- \n", i, ident, adcdata, TS);  // resets at 8161

          /*
            Channel Information:
            ---------------------------------------------------------------------
            Signal  |  Fast Ch   |  Slow Ch  |    HV          |     Slow Time   |
            ---------------------------------------------------------------------
            LaBr 0  |    72      |    64     | -1050          |      80         |
            ---------------------------------------------------------------------
            LaBr 1  |    73      |    65     | -1050          |     81          |
            ---------------------------------------------------------------------
            LaBr 2  |    74      |    66     | -1050          |     82          |
            ---------------------------------------------------------------------
            LaBr 3  |    75      |    67     | -1050          |    83           |
            ---------------------------------------------------------------------
            RF      |    78      |    //   | 5 bunches 100pA  |      94         |
            ---------------------------------------------------------------------
            Polaris |    79      |    //  |     //            |      95         |
            ---------------------------------------------------------------------
          */
          energy = (adcdata);

          // **************************************************  RF Time Walk 
          if (identity == 78)
          {
            RFenergy = energy;
            TimeStampRF = TS;
            LaBr0Data->Fill();
            LaBr1Data->Fill();
            LaBr2Data->Fill();
            LaBr3Data->Fill();
            RFTimeStamp->Fill(TimeStampRF);
            RFSigEnergy->Fill(RFenergy);
            //cout << "RF Energy: " << RFenergy << " RF Time: " << TimeStampRF << endl;
          }

          if (identity == 94) 
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              RFsignaltime = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr0Data->Fill();
              LaBr1Data->Fill();
              LaBr2Data->Fill();
              LaBr3Data->Fill();
              RFSigTime ->Fill(RFsignaltime);
              //cout << "RF Signal Time: " << RFsignaltime << endl;
            }
          } 

          // **************************************************  POLARIS
          if (identity == 79)
          {
            SlowEnergyPOLARIS = energy;
            TimeStampPOLARIS = TS; 
            POLARISTimeStamp->Fill(TimeStampPOLARIS);
            POLARISEnergy->Fill(SlowEnergyPOLARIS);
            LaBr0Data->Fill();
            LaBr1Data->Fill();
            LaBr2Data->Fill();
            LaBr3Data->Fill();
            //cout << "POLARIS Time: " << TimeStampPOLARIS << endl;
          }

          if (identity == 95)
          {
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
                   
            if (tock != 7)
            {
              SlowTimePOLARIS = ((TS * 100.0) + (tock * 20.0) + sick);
              POLARISSlowTime->Fill(SlowTimePOLARIS);
              LaBr0Data->Fill();
              LaBr1Data->Fill();
              LaBr2Data->Fill();
              LaBr3Data->Fill();
              //cout << "POLARIS Slow Time: " << SlowTimePOLARIS << endl;
            }
          }

          
        
          // **************************************************  LaBr0
          if (identity == 64) // energy of slow signal 
          { //
            slowLaBr0 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr0calib1stFullRange = a1stFullRange0*rndchannel + b1stFullRange0;
            slowLaBr0calib1stLowRange = a1stLowRange0*rndchannel + b1stLowRange0;
            slowLaBr0calib2ndFullRange = a2ndFullRange0*(rndchannel*rndchannel) + b2ndFullRange0*rndchannel + c2ndFullRange0;
            TimeStampLaBr0 = TS ; 
            LaBr0Data->Fill();
            s0none->Fill(slowLaBr0);
            s0calib1stFullRange->Fill(slowLaBr0calib1stFullRange);
            s0calib1stLowRange->Fill(slowLaBr0calib1stLowRange);
            s0calib2ndFullRange->Fill(slowLaBr0calib2ndFullRange);
            Timestamp0->Fill(TimeStampLaBr0);
            //cout << "LaBr0 Time: " << TimeStampLaBr0 << endl;
          }

          if (identity == 80) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr0 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr0Data->Fill(); 
              timeSlowSig0->Fill(slowtimeLaBr0);

              //cout << "LaBr0 Slow Signal Time: " << slowtimeLaBr0 << endl;
            }
          }

          // **************************************************  LaBr1
          if (identity == 65) // energy of slow signal
          { 
            slowLaBr1 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr1calib1stFullRange = a1stFullRange1*rndchannel + b1stFullRange1;
            slowLaBr1calib1stLowRange = a1stLowRange1*rndchannel + b1stLowRange1;
            slowLaBr1calib2ndFullRange = a2ndFullRange1*(rndchannel*rndchannel) + b2ndFullRange1*rndchannel + c2ndFullRange1;
            TimeStampLaBr1 = TS ; 
            LaBr1Data->Fill();
            s1none->Fill(slowLaBr1);
            s1calib1stFullRange->Fill(slowLaBr1calib1stFullRange);
            s1calib1stLowRange->Fill(slowLaBr1calib1stLowRange);
            s1calib2ndFullRange->Fill(slowLaBr1calib2ndFullRange);
            Timestamp1->Fill(TimeStampLaBr1);
          }

          if (identity == 81)  // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr1 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr1Data->Fill();
              timeSlowSig1->Fill(slowtimeLaBr1);
            }
          }

          // **************************************************  LaBr2  
          if (identity == 66) // energy of slow signal
          { 
            slowLaBr2 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr2calib1stFullRange = a1stFullRange2*rndchannel + b1stFullRange2;
            slowLaBr2calib1stLowRange = a1stLowRange2*rndchannel + b1stLowRange2;
            slowLaBr2calib2ndFullRange = a2ndFullRange2*(rndchannel*rndchannel) + b2ndFullRange2*rndchannel + c2ndFullRange2;
            TimeStampLaBr2 = TS ; 
            LaBr2Data->Fill();
            s2none->Fill(slowLaBr2);
            s2calib1stFullRange->Fill(slowLaBr2calib1stFullRange);
            s2calib1stLowRange->Fill(slowLaBr2calib1stLowRange);
            s2calib2ndFullRange->Fill(slowLaBr2calib2ndFullRange);
            Timestamp2->Fill(TimeStampLaBr2);
          }

          if (identity == 82) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr2 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr2Data->Fill();
              timeSlowSig2->Fill(slowtimeLaBr2);
            }
          }

          // **************************************************  LaBr3  
          if (identity == 67) // energy of slow signal
          { 
            slowLaBr3 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr3calib1stFullRange = a1stFullRange3*rndchannel + b1stFullRange3;
            slowLaBr3calib1stLowRange = a1stLowRange3*rndchannel + b1stLowRange3;
            slowLaBr3calib2ndFullRange = a2ndFullRange3*(rndchannel*rndchannel) + b2ndFullRange3*rndchannel + c2ndFullRange3;
            TimeStampLaBr3 = TS ; 
            LaBr3Data->Fill();
            s3none->Fill(slowLaBr3);
            s3calib1stFullRange->Fill(slowLaBr3calib1stFullRange);
            s3calib1stLowRange->Fill(slowLaBr3calib1stLowRange);
            s3calib2ndFullRange->Fill(slowLaBr3calib2ndFullRange);
            Timestamp3->Fill(TimeStampLaBr3);
          }

          if (identity == 83) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr3 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr3Data->Fill();
              timeSlowSig3->Fill(slowtimeLaBr3);
            }
          }
          TSlast = TS;
        }

        /* SYNC */
	  
        if ((dataframe & 0xc0f00000) == 0x80400000)
        {
          card = (dataframe & 0x3f000000) >> 24;
          TStop = (dataframe & 0x000fffff);
          TS = (TStop);
          TS = ((TS << 28));
          TS = ((TS | TSbot));

          TSdiff = TS-TSlast;

          SYNC = TS;
          SYNCdiff = SYNC-SYNClast;

          twidset=0;

          SYNClast = SYNC;

        }


        loop:
        pos+=4;
      }

    }

    end:  
    fclose(f);
      
    Run++;
      
  }


  finish:

  int e;
  printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
  printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
  printf("The difference is (min): %" PRIu64 "\n", ((TS-TSfirst)*(e-8))/60);

  LaBr0Data->Write();
  LaBr1Data->Write();
  LaBr2Data->Write();
  LaBr3Data->Write();
  s0none->Write();
  s0calib1stFullRange->Write();
  s0calib1stLowRange->Write();
  s0calib2ndFullRange->Write();
  s1none->Write();
  s1calib1stFullRange->Write();
  s1calib1stLowRange->Write();
  s1calib2ndFullRange->Write();
  s2none->Write();
  s2calib1stFullRange->Write();
  s2calib1stLowRange->Write();
  s2calib2ndFullRange->Write();
  s3none->Write();
  s3calib1stFullRange->Write();
  s3calib1stLowRange->Write();
  s3calib2ndFullRange->Write();
  Timestamp0->Write();
  Timestamp1->Write();
  Timestamp2->Write();
  Timestamp3->Write();
  RFSigTime->Write();
  RFSigEnergy->Write();
  RFTimeStamp->Write();
  timeSlowSig0->Write();
  timeSlowSig1->Write();
  timeSlowSig2->Write();
  timeSlowSig3->Write();
  POLARISEnergy->Write();
  POLARISTimeStamp->Write();
  POLARISSlowTime->Write();

  LaBr0Data->SetBranchAddress("EnergyPOLARIS", &SlowEnergyPOLARIS);
  LaBr0Data->SetBranchAddress("TimeStampPOLARIS", &TimeStampPOLARIS);
  LaBr0Data->SetBranchAddress("TimingPOLARIS", &SlowTimePOLARIS);
  Long64_t nentries = LaBr0Data -> GetEntries();
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr0Data->GetEntry(i);
    energyVtimestampPOLARIS->Fill(SlowEnergyPOLARIS, TimeStampPOLARIS);
    energyVslowsigtimePOLARIS->Fill(SlowEnergyPOLARIS, SlowTimePOLARIS);
  }

  energyVtimestampPOLARIS->Write();
  energyVslowsigtimePOLARIS->Write();

  /*LaBr0Data->SetBranchAddress("TimingSlow_RF", &RFsignaltime);
  LaBr0Data->SetBranchAddress("TimeStamp_RF", &TimeStampRF);
  LaBr0Data->SetBranchAddress("EnergySlow_LaBr_Det0_1stOrder_Calib", &slowLaBr0calib1st);
  LaBr0Data->SetBranchAddress("TimingSlow_LaBr_Det0", &slowtimeLaBr0);
  LaBr0Data->SetBranchAddress("TimeStamp_LaBr_Det0", &TimeStampLaBr0);
  Long64_t nentries = LaBr0Data -> GetEntries();
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr0Data->GetEntry(i);
    energyVslowsigtime0->Fill(slowLaBr0calib1st, slowtimeLaBr0-RFsignaltime);
    energyVtimestamp0->Fill(slowLaBr0calib1st, TimeStampLaBr0-TimeStampRF);
  }

  energyVslowsigtime0->Write();
  energyVtimestamp0->Write();
  */
 
  // **************************************************  CSV Files created from TTrees **************************************************
  /*FILE *fp;

  cout << "Starting to create CSV file for LaBr3 Detector L0" << endl;

  char *dirpath = dirname(strdup(argv[1])); // get the directory path where the argv[1] file is located
  string dirpath_str = dirpath; // create a string containing the directory path where the argv[1] file is located

  // PASS THE SPECIFIC NAME OF THE RUN AS ARGV[2] IN ORDER TO CREATE THE CSV FILE WITH THE CORRECT NAME

  fp = fopen((dirpath_str+"/LaBr3_Det0_"+argv[2]+".csv").c_str(), "w");
  cout << "Created CSV file for LaBr3 Detector L0: " << endl;
  fprintf(fp, "EnergySlow_LaBr_Det0_2ndOrder_Calib, EnergySlow_LaBr_Det0_1stOrder_Calib, EnergySlow_LaBr_Det0, TimingSlow_LaBr_Det0, TimeStamp_LaBr_Det0, TimingSlow_RF, TimeStamp_RF\n");
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr0Data->GetEntry(i);
    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", slowLaBr0calib2nd, slowLaBr0calib1st, slowLaBr0, slowtimeLaBr0, TimeStampLaBr0, RFsignaltime, TimeStampRF);
  }
  fclose(fp);

  fp = fopen((dirpath_str+"/LaBr3_Det1_"+argv[2]+".csv").c_str(), "w");
  cout << "Created CSV file for LaBr3 Detector L1: " << endl;
  fprintf(fp, "EnergySlow_LaBr_Det1_2ndOrder_Calib, EnergySlow_LaBr_Det1_1stOrder_Calib, EnergySlow_LaBr_Det1, TimingSlow_LaBr_Det1, TimeStamp_LaBr_Det1, TimingSlow_RF, TimeStamp_RF\n");
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr1Data->GetEntry(i);
    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", slowLaBr1calib2nd, slowLaBr1calib1st, slowLaBr1, slowtimeLaBr1, TimeStampLaBr1, RFsignaltime, TimeStampRF);
  }
  fclose(fp);

  fp = fopen((dirpath_str+"/LaBr3_Det2_"+argv[2]+".csv").c_str(), "w");
  cout << "Created CSV file for LaBr3 Detector L2: " << endl;
  fprintf(fp, "EnergySlow_LaBr_Det2_2ndOrder_Calib, EnergySlow_LaBr_Det2_1stOrder_Calib, EnergySlow_LaBr_Det2, TimingSlow_LaBr_Det2, TimeStamp_LaBr_Det2, TimingSlow_RF, TimeStamp_RF\n");
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr2Data->GetEntry(i);
    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", slowLaBr2calib2nd, slowLaBr2calib1st, slowLaBr2, slowtimeLaBr2, TimeStampLaBr2, RFsignaltime, TimeStampRF);
  }
  fclose(fp);

  fp = fopen((dirpath_str+"/LaBr3_Det3_"+argv[2]+".csv").c_str(), "w");
  cout << "Created CSV file for LaBr3 Detector L3: " << endl;
  fprintf(fp, "EnergySlow_LaBr_Det3_2ndOrder_Calib, EnergySlow_LaBr_Det3_1stOrder_Calib, EnergySlow_LaBr_Det3, TimingSlow_LaBr_Det3, TimeStamp_LaBr_Det3, TimingSlow_RF, TimeStamp_RF\n");
  for (Long64_t i=0; i < nentries; i++)
  {
    LaBr3Data->GetEntry(i);
    fprintf(fp, "%f, %f, %f, %f, %f, %f, %f\n", slowLaBr3calib2nd, slowLaBr3calib1st, slowLaBr3, slowtimeLaBr3, TimeStampLaBr3, RFsignaltime, TimeStampRF);
  }
  fclose(fp);
  */


  //  Read in Polaris Files 
  /* double energyPol1, xcoord1, ycoord1, zcoord1, energyPol2, xcoord2, ycoord2, zcoord2;
  int cryst1, cryst2;
  uint64_t TSPolaris1, TSPolaris2;
  int num = 0;
  ifstream infile1;
  ifstream infile2;

  // ROOT 
  TTree *PolarisMod1Data = new TTree("PolarisMod1Data", "PolarisModOneData"); 
  PolarisMod1Data->Branch("NoOfCrystInteractions", &cryst1 , "cryst1/I"); // No of interactions recorded in the module within a 10 ns window following the first event trigger (1 = single scatt, 2 = double, etc)
  PolarisMod1Data->Branch("Energy", &energyPol1 , "energyPol1"); // The energy deposited in the crystal in keV. 
  PolarisMod1Data->Branch("Xcoord", &xcoord1 , "xcoord1"); //The event location’s x-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod1Data->Branch("Ycoord", &ycoord1 , "ycoord1"); //The event location’s y-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod1Data->Branch("Zcoord", &zcoord1 , "zcoord1"); //The event location’s z-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod1Data->Branch("TimeStamp", &TSPolaris1 , "TSPolaris1/I"); //Event time of detection in tens of nanoseconds. Two lines with the same timestamp: double scatter event.

  TTree *PolarisMod2Data = new TTree("PolarisMod2Data", "PolarisModTwoData"); 
  PolarisMod2Data->Branch("NoOfCrystInteractions", &cryst2 , "cryst2/I"); // No of interactions recorded in the module within a 10 ns window following the first event trigger (1 = single scatt, 2 = double, etc)
  PolarisMod2Data->Branch("Energy", &energyPol2 , "energyPol2"); // The energy deposited in the crystal in keV. 
  PolarisMod2Data->Branch("Xcoord", &xcoord2 , "xcoord2"); //The event location’s x-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod2Data->Branch("Ycoord", &ycoord2 , "ycoord2"); //The event location’s y-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod2Data->Branch("Zcoord", &zcoord2 , "zcoord2"); //The event location’s z-coordinate in the detector’s coor-dinate system in millimetres.
  PolarisMod2Data->Branch("TimeStamp", &TSPolaris2 , "TSPolaris2/I"); //Event time of detection in tens of nanoseconds. Two lines with the same timestamp: double scatter event.
  
  // Files
  infile1.open("mod51.txt");
  infile2.open("mod52.txt");

  if (infile1.fail()) {cout << "error" << endl;    return 1;}

  if (infile2.fail()) {cout << "error" << endl;    return 1;}

  while (!infile1.eof())
  {
    string line;
    while (getline(infile1, line)) 
    {
      istringstream is(line);

      is >>  cryst1;
      is >>  energyPol1;
      is >>  xcoord1;
      is >>  ycoord1;
      is >>  zcoord1;
      is >>  TSPolaris1;     
      //cout << energyPol1 << endl;

      //printf(" === cryst = %d  / energy: %f / xcoord: %f / ycoord: %f / zcoord: %f / TSPolaris: %ld --- \n", cryst1, energyPol1, xcoord1, ycoord1, zcoord1, TSPolaris1);  // resets at 8161
      PolarisMod1Data->Fill();
      ++num;
    }
  } 

  PolarisMod1Data->Write();
  infile1.close();

  num = 0;
  while (!infile2.eof())
  {
    string line;
    while (getline(infile2, line)) 
    {
      istringstream is(line);

      is >>  cryst2;
      is >>  energyPol2;
      is >>  xcoord2;
      is >>  ycoord2;
      is >>  zcoord2;
      is >>  TSPolaris2;     
      
      //printf(" === cryst = %d  / energy: %f / xcoord: %f / ycoord: %f / zcoord: %f / TSPolaris: %ld --- \n", cryst2, energyPol2, xcoord2, ycoord2, zcoord2, TSPolaris2);  // resets at 8161
      PolarisMod2Data->Fill();
      ++num;
    }
  } 
  PolarisMod2Data->Write();
  infile2.close();
  */
  delete g; // automatically deletes all "owned" histograms, too
  return 0;
}
