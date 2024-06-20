/* LaBr3:Ce SORTING code for TDR format dataframe */
/* Produces ROOT TTrees and histograms for four LaBr3:Ce detectors: L0,L1,L2,L3*/

/* Adapted from P.M Jones code labrsort4 */

/* (c) S. Hart */
/* (27/10/22)  */


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
#include "labrsortincode_withcalib_latest.h"



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

int main (int argc, char* argv[])
{
  FILE *f;
  FILE *csvf;
  int i, n;
  char file_in[50];
  char file_out[50];
  char csv_out[50];
  int  identity, card=0, adcdata;
  double  slowLaBr0, TimeStampLaBr0, slowtimeLaBr0, fasttimeLaBr0, slowLaBr0calib1stFullRange, slowLaBr0calib2ndFullRange;
  double  slowLaBr1, TimeStampLaBr1, slowtimeLaBr1,fasttimeLaBr1, slowLaBr1calib1stFullRange, slowLaBr1calib2ndFullRange;
  double  slowLaBr2, TimeStampLaBr2, slowtimeLaBr2, fasttimeLaBr2, slowLaBr2calib1stFullRange, slowLaBr2calib2ndFullRange;
  double  slowLaBr3, TimeStampLaBr3, slowtimeLaBr3, fasttimeLaBr3, slowLaBr3calib1stFullRange, slowLaBr3calib2ndFullRange;
  double  energy0firstFullRange, energy1firstFullRange, energy2firstFullRange, energy3firstFullRange;
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
  sprintf(file_out, "%s_CalibratedData.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);

  //2nd order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a0 = 0, b0 = 0.74770, c0 = 1.64893;
  double a1 = 0, b1 = 0.80582, c1 = 1.09470;
  double a2 = 0, b2 = 0.51259, c2 = 1.66312;
  double a3 = 0, b3 = 0.70936, c3 = 2.18123;

  // 1st order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a0_1st = 0.75263, b0_1st = 0.71741;
  double a1_1st = 0.81401, b1_1st = -0.54088;
  double a2_1st = 0.51735, b2_1st = 0.41531;
  double a3_1st = 0.71452, b3_1st = 1.14457;


  TTree *LaBr0Data = new TTree("LaBr0Data", "LaBr0Data"); 
  LaBr0Data->Branch("EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange", &slowLaBr0calib1stFullRange , "slowLaBr0calib1stFullRange/D");
  LaBr0Data->Branch("EnergySlow_LaBr_Det0_2ndOrder_Calib_FullRange", &slowLaBr0calib2ndFullRange , "slowLaBr0calib2ndFullRange/D");
  LaBr0Data->Branch("EnergySlow_LaBr_Det0", &slowLaBr0 , "slowLaBr0/D");
  LaBr0Data->Branch("TimingSlow_LaBr0", &slowtimeLaBr0 , "slowtimeLaBr0/D");
  LaBr0Data->Branch("TimingFast_LaBr0", &fasttimeLaBr0 , "fasttimeLaBr0/D");
  LaBr0Data->Branch("TimingFast_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr0Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr0Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");

  TTree *LaBr1Data = new TTree("LaBr1Data", "LaBr1Data"); 
  LaBr1Data->Branch("EnergySlow_LaBr_Det1_1stOrder_Calib_FullRange", &slowLaBr1calib1stFullRange , "slowLaBr1calib1stFullRange/D");
  LaBr1Data->Branch("EnergySlow_LaBr_Det1_2ndOrder_Calib_FullRange", &slowLaBr1calib2ndFullRange , "slowLaBr1calib2ndFullRange/D");
  LaBr1Data->Branch("EnergySlow_LaBr_Det1", &slowLaBr1 , "slowLaBr1/D");
  LaBr1Data->Branch("TimingSlow_LaBr1", &slowtimeLaBr1 , "slowtimeLaBr1/D");
  LaBr1Data->Branch("TimingFast_LaBr1", &fasttimeLaBr1 , "fasttimeLaBr1/D");
  LaBr1Data->Branch("TimingFast_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr1Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr1Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");

  TTree *LaBr2Data = new TTree("LaBr2Data", "LaBr2Data"); 
  LaBr2Data->Branch("EnergySlow_LaBr_Det2_1stOrder_Calib_FullRange", &slowLaBr2calib1stFullRange , "slowLaBr2calib1stFullRange/D");
  LaBr2Data->Branch("EnergySlow_LaBr_Det2_2ndOrder_Calib_FullRange", &slowLaBr2calib2ndFullRange , "slowLaBr2calib2ndFullRange/D");
  LaBr2Data->Branch("EnergySlow_LaBr_Det2", &slowLaBr2 , "slowLaBr2/D");
  LaBr2Data->Branch("TimingSlow_LaBr2", &slowtimeLaBr2 , "slowtimeLaBr2/D");
  LaBr2Data->Branch("TimingFast_LaBr2", &fasttimeLaBr2 , "fasttimeLaBr2/D");
  LaBr2Data->Branch("TimingFast_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr2Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr2Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");

  TTree *LaBr3Data = new TTree("LaBr3Data", "LaBr3Data"); 
  LaBr3Data->Branch("EnergySlow_LaBr_Det3_1stOrder_Calib_FullRange", &slowLaBr3calib1stFullRange , "slowLaBr3calib1stFullRange/D");
  LaBr3Data->Branch("EnergySlow_LaBr_Det3_2ndOrder_Calib_FullRange", &slowLaBr3calib2ndFullRange , "slowLaBr3calib2ndFullRange/D");
  LaBr3Data->Branch("EnergySlow_LaBr_Det3", &slowLaBr3 , "slowLaBr3/D");
  LaBr3Data->Branch("TimingSlow_LaBr3", &slowtimeLaBr3 , "slowtimeLaBr3/D");
  LaBr3Data->Branch("TimingFast_LaBr3", &fasttimeLaBr3 , "fasttimeLaBr3/D");
  LaBr3Data->Branch("TimingFast_RF", &RFsignaltime , "RFsignaltime/D");
  LaBr3Data->Branch("EnergyPOLARIS", &SlowEnergyPOLARIS, "SlowEnergyPOLARIS/D");
  LaBr3Data->Branch("TimeStampPOLARIS", &TimeStampPOLARIS, "TimeStampPOLARIS/D");

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

  // 2nd order calibrations FULL RANGE
  TH1D* s0calib2ndFullRange = new TH1D("L0_Calibrated_2nd_FullRange", "LaBr3 Det 0 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s1calib2ndFullRange = new TH1D("L1_Calibrated_2nd_FullRange", "LaBr3 Det 1 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s2calib2ndFullRange = new TH1D("L2_Calibrated_2nd_FullRange", "LaBr3 Det 2 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);
  TH1D* s3calib2ndFullRange = new TH1D("L3_Calibrated_2nd_FullRange", "LaBr3 Det 3 Slow Signal Energy Spectrum - Calibrated 2nd Order - Full Range", 8000, 0, 8000);


  // use TRandom to create a random number 
  TRandom* randy = new TRandom();
  double bw = 1;

  // File stuff
  while(Run < RunEnd+1)
  {

      sprintf(file_in, "./%s_%d", argv[1], Run);
      
      f = fopen(file_in, "rb"); // open the file as a binary file
      if (!f)
    {
    fprintf(stderr, "Can't open file '%s'\n", file_in);
    goto finish;
    }
      
    printf("File %s opened...\n", file_in);

    
    
    while (!feof(f))
    {

      no_read = fread(buffer, sizeof(buffer[0]), SIZE, f); // read the file into the buffer
    
      if ( (feof(f)) && no_read <= 0 ) goto end; // if the end of the file is reached and no_read is less than or equal to 0, then go to end

      pos+=24; // skip the first 24 bytes of the file

      for (i=6; i< SIZE; i+=2) // loop through the buffer in steps of 2 because each event is 2 bytes
      {
        dataframe = buffer[i+1]; 
        TSbot = buffer[i];

        if ((TSbot == 0) && (dataframe == 0)) goto loop;
        if ((TSbot == 0x5e5e5e5e) && (dataframe == 0x5e5e5e5e)) goto loop; // if the TSbot and dataframe are equal to 0x5e5e5e5e, which is the end of the file, then go to loop
        if ((TSbot == 0x5e5e5e5e) && (dataframe == 0xffffffff)) goto loop; // if the TSbot and dataframe are equal to 0xffffffff, which is the end of the file, then go to loop


        /* dataframe */
        if ((dataframe & 0xc0000000) == 0xc0000000)  
        {
          identity = (dataframe & 0x0fff0000) >> 16; // get the identity of the detector which is the 12 bits in the middle of the dataframe, the 16
          //cout << "Identity: " << identity << endl;
          //if (identity == 76)
          //{
          //  cout << "RF signal" << endl;
          //}
          adcdata = (dataframe & 0x0000ffff); // get the adcdata which is the 16 bits at the end of the dataframe
          card = (identity / 32); // 0-3
          TS = (TS & 0x0000fffff0000000ULL); // get the top 20 bits of the TS

           
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
              TS=TS+0x10000000; // 0x10000000 = 2^28
              twidset=1;
            }

            TSdiff = TS-TSlast;
          }

          //printf(" === i = %d / ident: %d / energy: %d / TS: %ld --- \n", i, ident, adcdata, TS);  // resets at 8161

          /*
            Channel Information:
            ----------------------------------------------------------------------------
            Signal  |  Fast Ch   |  Slow Ch  |    HV                 |     Slow Time   |
            ----------------------------------------------------------------------------
            LaBr 0  |    72      |    64     | -1050                 |      80         |
            ---------------------------------------------------------------------
            LaBr 1  |    73      |    65     | -1050                 |     81          |
            ----------------------------------------------------------------------------
            LaBr 2  |    74      |    66     | -1050                 |     82          |
            ----------------------------------------------------------------------------
            LaBr 3  |    75      |    67     | -1050                 |    83           |
            ----------------------------------------------------------------------------
            RF      |    76      |    //     | 1 in 5 bunches 100pA  |      94         |
            ----------------------------------------------------------------------------
            Polaris |    79      |    //     |     //                |      95         | //sync signal every 2s
            ----------------------------------------------------------------------------
          */
          energy = (adcdata);

          // **************************************************  RF Time Walk 
          if (identity == 76)
          {
            RFenergy = energy;
            TimeStampRF = TS;
            LaBr0Data->Fill();
            LaBr1Data->Fill();
            LaBr2Data->Fill();
            LaBr3Data->Fill();
            //cout << "RF Energy: " << RFenergy << " RF Time: " << TimeStampRF << endl;
          }

          if (identity == 92) 
          { 
            tock = ((adcdata & 0x0000e000) >> 13); // 0-7
            tick = (adcdata & 0x00001fff); // 0-8191
            sick = (tick/8192.0) * 20.0; // 0-20 ns
            if (tock != 7) // 7 is the overflow bit
            {
              RFsignaltime = ((TS * 100.0) + (tock * 20.0) + sick); // in ns
              LaBr0Data->Fill();
              LaBr1Data->Fill();
              LaBr2Data->Fill();
              LaBr3Data->Fill();
              //cout << "RF Signal Time: " << RFsignaltime << endl;
            }
          } 

          // **************************************************  POLARIS
          if (identity == 79)
          {
            SlowEnergyPOLARIS = energy;
            TimeStampPOLARIS = TS;
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
            slowLaBr0calib2ndFullRange = a0*(rndchannel*rndchannel) + b0*(rndchannel) + c0;
            slowLaBr0calib1stFullRange = a0_1st*(rndchannel) + b0_1st;
            TimeStampLaBr0 = TS ;
            LaBr0Data->Fill();
            s0none->Fill(slowLaBr0);
            s0calib1stFullRange->Fill(slowLaBr0calib1stFullRange);
            s0calib2ndFullRange->Fill(slowLaBr0calib2ndFullRange);
            //cout << "LaBr0 Time: " << TimeStampLaBr0 << endl;
          }

          if (identity == 80) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              slowtimeLaBr0 = ((TS * 100.0) + (tock * 20.0) + sick) ;
              LaBr0Data->Fill(); 
              //cout << "LaBr0 Slow Signal Time: " << slowtimeLaBr0 << endl;
            }
          }

          if (identity == 88)
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              fasttimeLaBr0 = ((TS * 100.0) + (tock * 20.0) + sick) ;
              LaBr0Data->Fill(); 
            }
          }

          // **************************************************  LaBr1
          if (identity == 65) // energy of slow signal
          { 
            slowLaBr1 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr1calib2ndFullRange = a1*(rndchannel*rndchannel) + b1*(rndchannel) + c1;
            slowLaBr1calib1stFullRange = a1_1st*(rndchannel) + b1_1st;
            TimeStampLaBr1 = TS ; 
            LaBr1Data->Fill();
            s1none->Fill(slowLaBr1);
            s1calib1stFullRange->Fill(slowLaBr1calib1stFullRange);
            s1calib2ndFullRange->Fill(slowLaBr1calib2ndFullRange);
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
            }
          }

          if (identity == 89) // timing of fast signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              fasttimeLaBr1 = ((TS * 100.0) + (tock * 20.0) + sick) ;
              LaBr1Data->Fill(); 
            }
          }

          // **************************************************  LaBr2  
          if (identity == 66) // energy of slow signal
          { 
            slowLaBr2 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr2calib2ndFullRange = a2*(rndchannel*rndchannel) + b2*(rndchannel) + c2;
            slowLaBr2calib1stFullRange = a2_1st*(rndchannel) + b2_1st;
            TimeStampLaBr2 = TS ; 
            LaBr2Data->Fill();
            s2none->Fill(slowLaBr2);
            s2calib1stFullRange->Fill(slowLaBr2calib1stFullRange);
            s2calib2ndFullRange->Fill(slowLaBr2calib2ndFullRange);
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
            }
          }

          if (identity == 90) // timing of fast signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              fasttimeLaBr2 = ((TS * 100.0) + (tock * 20.0) + sick) ;
              LaBr2Data->Fill(); 
            }
          }

          // **************************************************  LaBr3  
          if (identity == 67) // energy of slow signal
          { 
            slowLaBr3 = energy; 
            rndUnif = randy->Rndm()*bw-(bw/2.);
            rndchannel = energy+rndUnif; 
            slowLaBr3calib2ndFullRange = a3*(rndchannel*rndchannel) + b3*(rndchannel) + c3;
            slowLaBr3calib1stFullRange = a3_1st*(rndchannel) + b3_1st;
            TimeStampLaBr3 = TS ; 
            LaBr3Data->Fill();
            s3none->Fill(slowLaBr3);
            s3calib1stFullRange->Fill(slowLaBr3calib1stFullRange);
            s3calib2ndFullRange->Fill(slowLaBr3calib2ndFullRange);
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
            }
          }
          TSlast = TS;

          if (identity == 91) // timing of fast signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);
            sick = (tick/8192.0) * 20.0;
            if (tock != 7)
            {
              fasttimeLaBr3 = ((TS * 100.0) + (tock * 20.0) + sick) ;
              LaBr3Data->Fill(); 
            }
          }
        }

        /* SYNC */
	  
        if ((dataframe & 0xc0f00000) == 0x80400000)
        {
          card = (dataframe & 0x3f000000) >> 24; // card number
          TStop = (dataframe & 0x000fffff); // 20 bit
          TS = (TStop); // 20 bit
          TS = ((TS << 28)); // shift to make room for the 28 bits of the TStart
          TS = ((TS | TSbot)); // add the 28 bits of the TStart

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

  LaBr0Data->Write();
  LaBr1Data->Write();
  LaBr2Data->Write();
  LaBr3Data->Write();
  s0none->Write();
  s0calib1stFullRange->Write();
  s0calib2ndFullRange->Write();
  s1none->Write();
  s1calib1stFullRange->Write();
  s1calib2ndFullRange->Write();
  s2none->Write();
  s2calib1stFullRange->Write();
  s2calib2ndFullRange->Write();
  s3none->Write();
  s3calib1stFullRange->Write();
  s3calib2ndFullRange->Write();

  printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
  printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
 
  /******************************************************************************/
 /*******************             RF CORRECTION             ********************/
 /******************************************************************************/
/*
  TF1* sbanano = 0;
  int RF_mask = 1, Out_RF_mask = 0;
  double TOF_Sbanano = 0, EnergyNoRF, EnergyRF;
  uint64_t mintime, maxtime, bins;
  
  int N = LaBr0Data -> GetEntries();

 //N = N*0.10;

  LaBr0Data->SetBranchAddress("EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange", &slowLaBr0calib1stFullRange);
  LaBr0Data->SetBranchAddress("TimingSlow_LaBr_Det0", &slowtimeLaBr0);
  LaBr0Data->SetBranchAddress("TimingSlow_RF", &RFsignaltime);
  LaBr0Data->GetEntry(N-1);
  TBranch* Energy_L0_Minus_RF = LaBr0Data->Branch("PromptEnergy", &EnergyNoRF, "EnergyNoRF/I");
  TBranch* RF_Energy = LaBr0Data->Branch("OutOfSyncEnergyRF", &EnergyRF, "EnergyRF/I");

  mintime = TSfirst*10;
  maxtime = TS*10;
  // cut off the first two leading digits from mintime and maxtime and retain the rest
  mintime = mintime;
  maxtime = maxtime;

  bins = (maxtime - mintime)*1e-9;
  cout << "Min time: " << mintime << " s"<< " \nMax time: " << maxtime << " s"<< "\nTime duration: " << ((maxtime-mintime))/60 << " minutes" << "\nbins: " << bins <<endl;

  TH2D *h_tdc_qdc = new TH2D("h_tdc_qdc", "htdcqdc", 8000, 0, 8000, bins, mintime, maxtime); // time in us bins
  TH2D *h_qdc_tdc = new TH2D("h_qdc_tdc", "hqdctdc", bins, mintime, maxtime, 8000, 0, 8000);
  TH1 *h_tdc = new TH1D("h_tdc", "htdc", bins, mintime, maxtime);

  TH1D *h_qdc0_0_Sync = new TH1D("h_qdc0_0_Sync", "In-sync gamma energy calibration", 6000, 0, 6000);
  TH1D *h_qdc0_0_outSync = new TH1D("h_qdc0_0_outSync", "out-sync gamma energy calibration", 6000, 0, 6000);

  for (Long64_t i=0; i < N; i++)
  {
    LaBr0Data->GetEntry(i);
    double slowtimeLaBr0 = slowtimeLaBr0;
    double RFsignaltime = RFsignaltime;
    float TOF = slowtimeLaBr0 - RFsignaltime;
    h_tdc_qdc->Fill(slowLaBr0calib1stFullRange, TOF);
    h_qdc_tdc->Fill(TOF, slowLaBr0calib1stFullRange);
    h_tdc->Fill(TOF);
    
  }
  h_tdc_qdc->Write("TDC_QDC");
  h_qdc_tdc->Write("QDC_TDC"); 
  

  sbanano = Sbanano(h_qdc_tdc, h_tdc_qdc, mintime, maxtime); //fit for sbanano function TF1, see .h file
  TCanvas* plot = new TCanvas();
  h_tdc_qdc->Draw();
  sbanano->Draw("same");
  plot->Write("Sbanano_2D_fitfunction");

  
    for (Long64_t i=0; i < N; i++)
  {
    LaBr0Data->GetEntry(i);
    double slowtimeLaBr0 = slowtimeLaBr0;
    double RFsignaltime = RFsignaltime;
    float TOF = slowtimeLaBr0 - RFsignaltime;
    slowtime0RF -> Fill(TOF);
    TOF_Sbanano = (double)TOF - sbanano->Eval(slowLaBr0calib1stFullRange);

    //I label each event if in-sync or out-of-sync with RF
    RF_mask = 0;
    Out_RF_mask =0;

    //first loop for time > 60 ns is the RF period
    for(double t =mintime; t<(maxtime); t+=60.9) 
    {
      
      if(abs(TOF_Sbanano-t) < 1.5) // if the time of flight is within 1.5 ns of the RF period, then it is in sync with the RF
      {
        RF_mask = 1;
        break;
      }

      if(abs(TOF_Sbanano-(t+5)) < 1.5) // if the time of flight is within 1.5 ns of the RF period + 5 ns, then it is out of sync with the RF
      {
        Out_RF_mask =1;
        break;
      }
    }

    //second loop for time < 0
    for(double t =mintime; t>-(mintime*0.1); t-=610)
    { 

      if(abs(TOF_Sbanano-t) < 15) 
      {
        RF_mask = 1;
        break;
      }

      if(abs(TOF_Sbanano-(t-300)) < 15)
      {
        Out_RF_mask=1;
        break;
      }

    }

    if(RF_mask==1)
    {
      h_qdc0_0_Sync->Fill(slowLaBr0calib1stFullRange);
      EnergyNoRF = slowLaBr0calib1stFullRange;
      LaBr0Data -> Fill();
    }
    if(Out_RF_mask==1)
    {
      h_qdc0_0_outSync->Fill(slowLaBr0calib1stFullRange);
      EnergyRF = slowLaBr0calib1stFullRange;
      LaBr0Data -> Fill();
    }

  }

  h_qdc0_0_Sync -> Write();
  h_qdc0_0_outSync -> Write();  
  LaBr0Data->Write();
  slowtime0RF->Write();

  //SUBTRACTION 
	//ONLY EVENTS PG = SYNC - OutSync
	TH1D* h_qdc0_0_Prompt = (TH1D*)h_qdc0_0_Sync->Clone("h_qdc0_0_Prompt");
  h_qdc0_0_Prompt -> Add(h_qdc0_0_outSync, -1);
  
	//PROMPT GAMMA SPECTRA
	h_qdc0_0_Prompt -> Scale(1., "width");
	h_qdc0_0_Prompt->Write();
*/

  delete g; // automatically deletes all "owned" histograms, too
  return 0;
}
