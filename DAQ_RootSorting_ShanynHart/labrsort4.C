/* Advanced sorting code for TDR format dataframe */
/* Produces ROOT TTrees */

/* Adapted from P.M Jones code labrsort4 */

/* (c) S. Hart */
/* (10/06/22)  */


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
#include "labrsort4.h"

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

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <inttypes.h>
#include <initializer_list>
#include <iomanip>
#include <fstream>
#include <sstream>


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
  int i, n;
  char file_in[50];
  char file_out[50];
  int  identity, card=0, adcdata;
  double  slowLaBr0, TimeStampLaBr0, slowtimeLaBr0;
  double  slowLaBr1, TimeStampLaBr1, slowtimeLaBr1;
  double  slowLaBr2, TimeStampLaBr2, slowtimeLaBr2;
  double  slowLaBr3, TimeStampLaBr3, slowtimeLaBr3;
  double  energy0first, energy1first, energy2first, energy3first;
  double  energy0second, energy1second, energy2second, energy3second;
  double  energy0third, energy1third, energy2third, energy3third;
  int N = 5000;
  double bw = 1;
  double rndchannel, rndUnif, ch;
  int tick, tock;
  double sick;
  int twidset=0;
  double Td;
  int RunEnd=999999, Run=0;
  int no_read;
  std::vector<int> time_stamp;
  bool CalibFlag = false;
  
  	for (int i = 0; i < argc; i++)
    {
		if(strcmp(argv[i],"-calib") == 0) {CalibFlag = true;}
	  }
  
  if ((argc < 2) || (argc > 3))
  {
	fprintf(stderr, "Usage: %s dataframe [last_subrun]\n", argv[0]);
	exit(1);
  }

  /* Calibrations */
  // Energy calibrations
	Calibration LaBrCalib;

  //first order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a10 =  1.5102,  b10 =  -507.8355; 
  double a11 =  1.6303,  b11 =  -509.802;
  double a12 =  1.0404,  b12 =  -509.8517;
  double a13 =  1.4357,  b13 =  -506.9651;

  //2nd order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a20 =  0.00050627,  b20 =  -0.041062,  c20 =  190.9778;
  double a21 =  0.00058954,  b21 =  -0.044864,  c21 =  190.833;
  double a22 =  0.00024151,  b22 =  -0.033294,  c22 =  192.6998;
  double a23 =  0.00045733,  b23 =  -0.03787,  c23 =  190.7451;

  // 3rd order polynomial //date: 15/6/22 for AmBe and 152Eu Sources
  double a30 =  2.3443e-07,  b30 =  -0.00053808,  c30 =  1.0929,  d30 =  -52.7598;
  double a31 =  2.9525e-07,  b31 =  -0.00062997,  c31 =  1.1835,  d31 =  -54.7049;
  double a32 =  7.7626e-08,  b32 =  -0.00026014,  c32 =  0.75698,  d32 =  -54.5558;
  double a33 =  2.0145e-07,  b33 =  0.00048623,  c33 =  1.0389,  d33 =  -52.2604;
  
  /* ROOT STUFF */
  sprintf(file_out, "%s_rawData.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);

  /*
  TTree *LaBr0Data = new TTree("LaBr0Data", "LaBr0RawData"); 
  LaBr0Data->Branch("SlowEnergyLaBr0", &slowLaBr0 , "slowLaBr0");
  LaBr0Data->Branch("SlowTIMINGLaBr0", &slowtimeLaBr0 , "slowtimeLaBr0");
  LaBr0Data->Branch("TSLaBr0", &TimeStampLaBr0 , "TimeStampLaBr0");
  LaBr0Data->Branch("RFTime", &RFsignaltime , "RFsignaltime");
  LaBr0Data->Branch("EnergyCalibSecondOrder", &energy0second , "energy0second");

  TTree *LaBr1Data = new TTree("LaBr1Data", "LaBr1RawData"); 
  LaBr1Data->Branch("SlowEnergyLaBr1", &slowLaBr1 , "slowLaBr1");
  LaBr1Data->Branch("SlowTIMINGLaBr1", &slowtimeLaBr1 , "slowtimeLaBr1");
  LaBr1Data->Branch("TSLaBr1", &TimeStampLaBr1 , "TimeStampLaBr1");
  LaBr1Data->Branch("RFTime", &RFsignaltime , "RFsignaltime");
  LaBr1Data->Branch("EnergyCalibSecondOrder", &energy1second , "energy1second");

  TTree *LaBr2Data = new TTree("LaBr2Data", "LaBr2RawData"); 
  LaBr2Data->Branch("SlowEnergyLaBr2", &slowLaBr2 , "slowLaBr2");
  LaBr2Data->Branch("SlowTIMINGLaBr2", &slowtimeLaBr2 , "slowtimeLaBr2");
  LaBr2Data->Branch("TSLaBr2", &TimeStampLaBr2 , "TimeStampLaBr2");
  LaBr2Data->Branch("RFTime", &RFsignaltime , "RFsignaltime");
  LaBr2Data->Branch("EnergyCalibSecondOrder", &energy2second , "energy2second");

  TTree *LaBr3Data = new TTree("LaBr3Data", "LaBr3RawData"); 
  LaBr3Data->Branch("SlowEnergyLaBr3", &slowLaBr3 , "slowLaBr3");
  LaBr3Data->Branch("SlowTIMINGLaBr3", &slowtimeLaBr3 , "slowtimeLaBr3");
  LaBr3Data->Branch("TSLaBr3", &TimeStampLaBr3 , "TimeStampLaBr3");
  LaBr3Data->Branch("RFTime", &RFsignaltime , "RFsignaltime");
  LaBr3Data->Branch("EnergyCalibSecondOrder", &energy3second , "energy3second");
  */

  // Calibration histograms 152 Eu & AmBe
  TH1D *s0none = new TH1D("L0_Uncalibrated","LaBr3 Det 0 Slow Signal Energy Spectrum - Uncalibrated",1000,0,5000);
  TH1D *s0first = new TH1D("L0_CalibratedFirstOrder","LaBr3 Det 0 Slow Signal Energy Spectrum - Linear Calibrated",1000,0,5000);
  TH1D *s0second = new TH1D("L0_CalibratedSecondOrder","LaBr3 Det 0 Slow Signal Energy Spectrum - Second Order Calibrated",1000,0,5000);
  TH1D *s0third = new TH1D("L0_CalibratedThirdOrder","LaBr3 Det 0 Slow Signal Energy Spectrum - Third Calibrated",1000,0,5000);

  TH1D *s1none = new TH1D("L1_Uncalibrated","LaBr3 Det 1 Slow Signal Energy Spectrum - Uncalibrated",1000,0,5000);
  TH1D *s1first = new TH1D("L1_CalibratedFirstOrder","LaBr3 Det 1 Slow Signal Energy Spectrum - Linear Calibrated",1000,0,5000);
  TH1D *s1second = new TH1D("L1_CalibratedSecondOrder","LaBr3 Det 1 Slow Signal Energy Spectrum - Second Order Calibrated",1000,0,5000);
  TH1D *s1third = new TH1D("L1_CalibratedThirdOrder","LaBr3 Det 1 Slow Signal Energy Spectrum - Third Calibrated",1000,0,5000);

  TH1D *s2none = new TH1D("L2_Uncalibrated","LaBr3 Det 2 Slow Signal Energy Spectrum - Uncalibrated",1000,0,5000);
  TH1D *s2first = new TH1D("L2_CalibratedFirstOrder","LaBr3 Det 2 Slow Signal Energy Spectrum - Linear Calibrated",1000,0,5000);
  TH1D *s2second = new TH1D("L2_CalibratedSecondOrder","LaBr3 Det 2 Slow Signal Energy Spectrum - Second Order Calibrated",1000,0,5000);
  TH1D *s2third = new TH1D("L2_CalibratedThirdOrder","LaBr3 Det 2 Slow Signal Energy Spectrum - Third Calibrated",1000,0,5000);

  TH1D *s3none = new TH1D("L3_Uncalibrated","LaBr3 Det 3 Slow Signal Energy Spectrum - Uncalibrated",1000,0,5000);
  TH1D *s3first = new TH1D("L3_CalibratedFirstOrder","LaBr3 Det 3 Slow Signal Energy Spectrum - Linear Calibrated",1000,0,5000);
  TH1D *s3second = new TH1D("L3_CalibratedSecondOrder","LaBr3 Det 3 Slow Signal Energy Spectrum - Second Order Calibrated",5000,0,5000);
  TH1D *s3third = new TH1D("L3_CalibratedThirdOrder","LaBr3 Det 3 Slow Signal Energy Spectrum - Third Calibrated",1000,0,5000);

  //TRandom* randy = new TRandom(); // random number generator for calibration purposes

  // RF time walk correction histos
  TH1D *RFTiming = new TH1D("RF_Timing","RF Timing", 1000, 1e8, 2.6e10);
  TH1D *LaBr0SlowTime = new TH1D("LaBr0_SlowTime","LaBr0 Slow Time", 1000, 1e8, 2.6e10);
  TH1D *LaBr1SlowTime = new TH1D("LaBr1_SlowTime","LaBr1 Slow Time", 1000, 1e8, 2.6e10);
  TH1D *LaBr2SlowTime = new TH1D("LaBr2_SlowTime","LaBr2 Slow Time", 1000, 1e8, 2.6e10);
  TH1D *LaBr3SlowTime = new TH1D("LaBr3_SlowTime","LaBr3 Slow Time", 1000, 1e8, 2.6e10);
  TH2D *LaBr0LessRF = new TH2D("LaBr0LessRF","LaBr0 2nd Order Calibration Time Walk Correction", 1000, 1e8, 2.6e10, 1000, 0, 5000);
  TH2D *LaBr1LessRF = new TH2D("LaBr1LessRF","LaBr1 2nd Order Calibration Time Walk Correction", 1000, 1e8, 2.6e10, 1000, 0, 5000);
  TH2D *LaBr2LessRF = new TH2D("LaBr2LessRF","LaBr2 2nd Order Calibration Time Walk Correction", 1000, 1e8, 2.6e10, 1000, 0, 5000);
  TH2D *LaBr3LessRF = new TH2D("LaBr3LessRF","LaBr3 2nd Order Calibration Time Walk Correction", 1000, 1e8, 2.6e10, 1000, 0, 5000);

  // File stuff
  while(Run < RunEnd+1)
  {

      sprintf(file_in, "./%s_%d", argv[1], Run);
      
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
          //cout << "Identity: " << identity << endl;
          //printf("  ident: %d /  ", ident);  // resets at 8161
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
            Channel Info:
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
            RF      |    78      |    //   | 5 bunches 100pA  |    94           |
            ---------------------------------------------------------------------
            Polaris |    79      |    //  |     //            |    95           |
            ---------------------------------------------------------------------
          */
          energy = (adcdata);

          // RF Time Walk Corrections
          if (identity == 94) 
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);

            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              double RFsignaltime = ((TS * 100.0) + (tock * 20.0) + sick);
              RFTiming->Fill(RFsignaltime);
            }
          } 

          // LaBr0
          if (identity == 64) // energy of slow signal 
          { 
            slowLaBr0 = energy; TimeStampLaBr0 = (TS & 0x0000fffff0000000ULL); //LaBr0Data->Fill(); 
            s0none->Fill(energy);
          }

          if (identity == 80) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);

            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr0 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr0SlowTime->Fill(slowtimeLaBr0);
            }
              //LaBr0Data->Fill(); 
          }

          // LaBr1
          if (identity == 65) // energy of slow signal
          { 
            slowLaBr1 = energy; TimeStampLaBr1 = (TS & 0x0000fffff0000000ULL); //LaBr1Data->Fill(); 
            s1none->Fill(energy);
          }
          if (identity == 81)  // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);

            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr1 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr1SlowTime->Fill(slowtimeLaBr1);
            }
          }

          // LaBr2  
          if (identity == 66) // energy of slow signal
          { 
            slowLaBr2 = energy; TimeStampLaBr2 = (TS & 0x0000fffff0000000ULL); //LaBr2Data->Fill();
            s2none->Fill(energy);
          }
          if (identity == 82) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);

            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr2 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr2SlowTime->Fill(slowtimeLaBr2);
            }
          }

          // LaBr3  
          if (identity == 67) // energy of slow signal
          { 
            slowLaBr3 = energy; TimeStampLaBr3 = (TS & 0x0000fffff0000000ULL); //LaBr3Data->Fill(); 
            s3none->Fill(energy);
          }
          if (identity == 83) // timing of slow signal
          { 
            tock = ((adcdata & 0x0000e000) >> 13);
            tick = (adcdata & 0x00001fff);

            sick = (tick/8192.0) * 20.0;

            if (tock != 7)
            {
              slowtimeLaBr3 = ((TS * 100.0) + (tock * 20.0) + sick);
              LaBr3SlowTime->Fill(slowtimeLaBr3);
            }
          }
      
        }
        // RF Energy vs Detection Time: Time Walk Corrected
        //LaBr0LessRF -> Fill(energy0second-RFsignaltime, energy0second);  
        //LaBr1LessRF -> Fill(energy1second-RFsignaltime, energy1second);
        //LaBr2LessRF -> Fill(energy2second-RFsignaltime, energy2second);
        //LaBr3LessRF -> Fill(energy3second-RFsignaltime, energy3second);

        loop:
        pos+=4;
      }

      //FIND ENERGY CALIBRATION
      if(!CalibFlag)
      {
        LaBrCalib.enecalib = eneCalib(s0none);
        a10 = LaBrCalib.enecalib[0]; 
        b10 = LaBrCalib.enecalib[1]; 
      }
      else
      {
        LaBrCalib.enecalib[0] = a10;
        LaBrCalib.enecalib[1] = b10;
      }
    }

      //HIST ENE CALIBRATION
      Transform(s0none, LaBrCalib.intercalib01[0], LaBrCalib.intercalib01[1]);
      //cout << "LaBr0: " << LaBrCalib.intercalib01[0] << " " << LaBrCalib.intercalib01[1] << endl;
      //energy calib
			ChannelToEneHisto(s0none, s0first, LaBrCalib);

    end:  

    fclose(f);
  
    Run++;
      
  }


  finish:

  //LaBr0Data->Write();
  //LaBr1Data->Write();
  //LaBr2Data->Write();
  //LaBr3Data->Write();
  s0none->Write();
  s0first->Write();
  s0second->Write();
  s0third->Write();
  s1none->Write();
  s1first->Write();
  s1second->Write();
  s1third->Write();
  s2none->Write();
  s2first->Write();
  s2second->Write();
  s2third->Write();
  s3none->Write();
  s3first->Write();
  s3second->Write();
  s3third->Write();
  LaBr0LessRF->Write();
  LaBr1LessRF->Write();
  LaBr2LessRF->Write();
  LaBr3LessRF->Write();
  RFTiming->Write();
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
