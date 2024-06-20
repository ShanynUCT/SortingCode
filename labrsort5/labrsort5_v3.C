/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (29/03/2022) - RF corrected timing and Trees  */

/* COMPILE: g++ -std=c++0x -O3 labrsort5_v3.C -o exe `root-config --cflags --libs` -lSpectrum */
/* RUN: ./exe /path/to/raw/experiment/data/ (e.g. /home/shanyn/Data/R50) runXX (PASS THE SPECIFIC NAME OF THE RUN AS ARGV[2] IN ORDER TO CREATE THE CSV FILE WITH THE CORRECT NAME)*/
//      ./exe ~/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <inttypes.h> 
#include <time.h>
#include <cstdlib>
#include <pthread.h>

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
uint64_t TStop, TSbot, data;
uint64_t TS, TSinit=0, TSfirst, counter,  SYNC, SYNClast=0;
int64_t item;
uint64_t TimeStamp0[100]={0}, TimeStamp1[100]={0}, TimeStamp2[100]={0} ,TimeStamp3[100]={0}, TimeStamp4[100]={0};
uint16_t Energy[100]={0};
uint16_t QDC[100]={0}, qdc;
uint64_t TS2;
int64_t TSdiff, SYNCdiff;
int64_t count=0;
unsigned long int pos = 0;
int verbose = 255;

int main (int argc, char** argv)
{
  FILE *f;
  int i, j, k;
  int blocks_in=0;
  int twidset=0;
  char file_in[50];
  char file_out[50];
  int ident, card=0, adcdata, detectorID;
  int tick, tock;
  double sick;
  int RunEnd=999999, Run=0;
  
  int m1,m2;
  int sum1, sum2;

  time_t start, end;
  int fastcount[4]={0,0,0,0}, slowcount[4]={0,0,0,0}, fastslowcount[4]={0,0,0,0};
  time(&start);

  /*  W  I  N  D  O  W  */
  // Define coincidence window here in 10's of ns
  int64_t Window = 60; // window = 1/2 * desired window size 22 //60, 65  good //
  /*  W  I  N  D  O  W  */

  // Declare vectors to store the event data
  std::vector<double> energyS[5];
    std::vector<double> energyF[5];
  std::vector<double> timeF[5];
  std::vector<double> timeS[5];
  std::vector<double> slowECalibL[5];
  std::vector<double> slowEPOLARIS(1, 0.0); //Polaris has only one detector
  std::vector<double> slowTPOLARIS(1, 0.0);
  std::vector<int> pop(1, 0);
  double energySlow[5]={0,0,0,0,0}, energyFast[5]={0,0,0,0,0}, timeFast[5]={0,0,0,0,0}, timeSlow[5]={0,0,0,0,0}, slowEnergyPOLARIS[1]={0}, slowTimePOLARIS[1]={0}, insyncEnergySlow[5]={0,0,0,0,0},insyncTimeFast[5]={0,0,0,0,0};
  double tdF, tdS, tdSF;
  int Ed;

  //Energy calibrations 2nd order polynomials (APPLIED) # 15062022 Data
    //   std::vector<double> a2 = {0.00000288786796, 0.00000351711116, 0.00000160595876, 0.00000300000230};
    //   std::vector<double> b2 = {0.747696355, 0.805816609, 0.513776183, 0.709359122};
    //   std::vector<double> c2 = {1.64893477, 1.09470037, 1.37154046, 2.18122723};

  //Energy calibrations 2nd order polynomials (APPLIED) # 28022023 Data
    std::vector<double> a2 = {0.00001942, 0.00010155, 0.00003208, 0.00002200};
    std::vector<double> b2 = {1.81799962, 2.61697133, 2.21345098, 2.05179483};
    std::vector<double> c2 = {-1.65129151, -0.54900373, -4.86476637, -0.97991876};

  double bw = 1, energycalib, rndchannel, rndUnif;
  TRandom3* randy = new TRandom3();

  //Time calibrations
  std::vector<double> tcalibF = {510.00, 621.00, 650.00, 655.00, 652.00};
  std::vector<double> tcalibS = {0, 0, 0, 0, 0};
    
  int no_read;
  std::vector<int> time_stamp;
    
  if ((argc < 2) || (argc > 3))
  {
    fprintf(stderr, "Usage: %s data [last_subrun]\n", argv[0]);
    exit(1);
  }

  /* ROOT STUFF */
  sprintf(file_out, "%s_oldsort.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);
  TFile *g = new TFile(file_out,"recreate");
  printf("ROOT file %s opened...\n", file_out);

    // std::string outputDirectory = argv[1];                        
    // outputDirectory = outputDirectory.substr(0, outputDirectory.length() - 3);
    // std::ofstream outputFile(outputDirectory + "/output.txt");
    // if (!outputFile) 
    // {
    //    std::cout << "Error opening the output file." << std::endl;
    //    return 1;
    // }

  // ################ Create ROOT TTrees ################
  TTree *LaBrData = new TTree("LaBrData","LaBrData");
  LaBrData->Branch("slowECalibL0", &energySlow[0], "slowECalibL0/D");
  LaBrData->Branch("timeSL0", &timeSlow[0], "timeSL0/D");
  LaBrData->Branch("energyFastL0", &energyFast[0], "energyFastL0/D");
  LaBrData->Branch("timeFL0", &timeFast[0], "timeFL0/D");
  LaBrData->Branch("slowECalibL1", &energySlow[1], "slowECalibL1/D");
  LaBrData->Branch("timeSL1", &timeSlow[1], "timeSL1/D");
  LaBrData->Branch("energyFastL1", &energyFast[1], "energyFastL1/D");
  LaBrData->Branch("timeFL1", &timeFast[1], "timeFL1/D");
  LaBrData->Branch("slowECalibL2", &energySlow[2], "slowECalibL2/D");
  LaBrData->Branch("timeSL2", &timeSlow[2], "timeSL2/D");
  LaBrData->Branch("energyFastL2", &energyFast[2], "energyFastL2/D");
  LaBrData->Branch("timeFL2", &timeFast[2], "timeFL2/D");
  LaBrData->Branch("slowECalibL3", &energySlow[3], "slowECalibL3/D");
  LaBrData->Branch("timeSL3", &timeSlow[3], "timeSL3/D");
  LaBrData->Branch("energyFastL3", &energyFast[3], "energyFastL3/D");
  LaBrData->Branch("timeFL3", &timeFast[3], "timeFL3/D");
  LaBrData->Branch("slowEPOLARIS", &slowEnergyPOLARIS[0], "slowEPOLARIS/D");
  LaBrData->Branch("slowTPOLARIS", &slowTimePOLARIS[0], "slowTPOLARIS/D");

  // ################ Create histograms ################
  TH1D** fastE=new TH1D*[5];
  TH1D** slowE=new TH1D*[5];
  TH1D** uncalibE=new TH1D*[5];
  TH1D** coincSlowE1=new TH1D*[5];
  TH1D** coincSlowE2=new TH1D*[5];
  TH1D** fastTD0=new TH1D*[5];
  TH1D** fastTD1=new TH1D*[5];
  TH1D** fastTD2=new TH1D*[5];
  TH1D** fastTD3=new TH1D*[5];
  TH1D** fastTD4=new TH1D*[5];
  TH1D** slowTD0=new TH1D*[5]; 
  TH1D** slowTD1=new TH1D*[5];
  TH1D** slowTD2=new TH1D*[5];
  TH1D** slowTD3=new TH1D*[5];
  TH1D** slowTD4=new TH1D*[5];
  TH1D** slowfastTD=new TH1D*[5];
  TH1D** multSlowE=new TH1D*[5];
  TH2D** fastTDslowE = new TH2D*[5];
  TH1D** insync = new TH1D*[5];
  TH1D** outsync = new TH1D*[5];
  TH1D** oneMeVgate = new TH1D*[5];
  TH1D** twoMeVgate = new TH1D*[5];
  TH1D** threeMeVgate = new TH1D*[5];
  TH1D** fourMeVgate = new TH1D*[5];
  TH1D** fiveMeVgate = new TH1D*[5];
  TH1D** sixMeVgate = new TH1D*[5];
  TH1F** eventCounts = new TH1F*[5];


  for (j = 0; j <=4; j++)
  {
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%02d", j),"Spectrum",8000,0,8000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%02d", j),"Spectrum",8000,0,8000);
    uncalibE[j] = new TH1D(TString::Format("Uncalib_Slow_Energy_L%2d", j),"Spectrum",8000,0,8000);

    coincSlowE1[j] = new TH1D(TString::Format("Coinc1_Energy_L%2d", j),"Spectrum",2000,0,2000);
    coincSlowE2[j] = new TH1D(TString::Format("Coinc2_Energy_L%2d", j),"Spectrum",2000,0,2000);

    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",1000,0,1000);
    fastTD1[j] = new TH1D(TString::Format("Fast Time Difference_L1-L%2d", j),"Spectrum",1000,0,1000);
    fastTD2[j] = new TH1D(TString::Format("Fast Time Difference_L2-L%2d", j),"Spectrum",1000,0,1000);
    fastTD3[j] = new TH1D(TString::Format("Fast Time Difference_L3-L%2d", j),"Spectrum",1000,0,1000);
    fastTD4[j] = new TH1D(TString::Format("Fast Time Difference_RF-L%2d", j),"Spectrum",2000,-1000,1000);

    slowTD0[j] = new TH1D(TString::Format("Slow Time Difference_L0-L%2d", j),"Spectrum",1000,0,1000);
    slowTD1[j] = new TH1D(TString::Format("Slow Time Difference_L1-L%2d", j),"Spectrum",1000,0,1000);
    slowTD2[j] = new TH1D(TString::Format("Slow Time Difference_L2-L%2d", j),"Spectrum",1000,0,1000);
    slowTD3[j] = new TH1D(TString::Format("Slow Time Difference_L3-L%2d", j),"Spectrum",1000,0,1000);
    slowTD4[j] = new TH1D(TString::Format("Slow Time Difference_RF-L%2d", j),"Spectrum",1000,0,1000);

    slowfastTD[j] = new TH1D(TString::Format("Slow-Fast Time Difference_L%2d", j),"Spectrum",2400,-1200,1200);
    multSlowE[j] = new TH1D(TString::Format("Multiplicity_Slow_Energy_%2d", j),"Spectrum",2500,0,2500);
    fastTDslowE[j] = new TH2D(TString::Format("Time_Difference_L%2d_Vs_Slow_Energy", j), TString::Format("Time L%2d - Vs Energy Spectrum", j),2000,-1000,1000,8000,0,8000);
    insync[j] = new TH1D(TString::Format("insync L%2d", j),TString::Format("Time gated energy spectrum of L%2d", j),8000,0,8000);
    outsync[j] = new TH1D(TString::Format("outsync L%2d", j),TString::Format("Time gated energy spectrum of L%2d", j),8000,0,8000);
    oneMeVgate[j] = new TH1D(TString::Format("1MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    twoMeVgate[j] = new TH1D(TString::Format("2MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    threeMeVgate[j] = new TH1D(TString::Format("3MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    fourMeVgate[j] = new TH1D(TString::Format("4MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    fiveMeVgate[j] = new TH1D(TString::Format("5MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    sixMeVgate[j] = new TH1D(TString::Format("6MeVgate L%2d", j),TString::Format("Energy gated time spectrum of L%2d", j),2000,-1000,1000);
    eventCounts[j] = new TH1F(TString::Format("Event Counts L%2d", j), "Event Counts", 3, 0.5, 3.5);
  }
  TH1D *mul1 = new TH1D("Multiplicity 1","spectrum",16,0,15);
  TH1D *mul2 = new TH1D("Multiplicity 2","spectrum",16,0,15);
  TH1D *mul3 = new TH1D("Multiplicity 3","spectrum",16,0,15);
  TH1D *hit1 = new TH1D("Hitpattern 1","spectrum",16,0,15);
  TH1D *hit2 = new TH1D("Hitpattern 2","spectrum",16,0,15);
  TH1D *hit3 = new TH1D("Hitpattern 3","spectrum",16,0,15);
  TH1D *sumspec0 = new TH1D("Sum 0","Spectrum",65536,0,65535);
  TH1D *sumspec00 = new TH1D("Sum 00","Spectrum",65536,0,65535);
  TH1D *sumspec01 = new TH1D("Sum 01","Spectrum",65536,0,65535);
  TH1D *sumspec02 = new TH1D("Sum 02","Spectrum",65536,0,65535);
  TH1D *sumspec03 = new TH1D("Sum 03","Spectrum",65536,0,65535);
  TH1D *sumspec1 = new TH1D("Sum 1","Spectrum",65536,0,65535);

  // ################ Start of data analysis ####################
  //File stuff
  while(Run < RunEnd+1)
   {
        sprintf(file_in, "%s_%d", argv[1], Run);
        
        f = fopen(file_in, "rb");
        if (!f)
        {
        fprintf(stderr, "Can't open file '%s'\n", file_in);
        goto finish;
        }
        
        printf("File %s opened...\n", file_in);

        while (!feof(f))
        {
            //Read a whole block of data in (64kbytes)
        no_read = fread(buffer, sizeof(buffer[0]), SIZE, f);
        blocks_in++;
        
        if ( (feof(f)) && no_read <= 0 ) goto end;

        pos+=24;

            for (i=6; i< SIZE; i+=2)
            {
                //Start of new event data
                //Loop through each 64bit data work and decypher

                data = buffer[i+1];
                TSbot = buffer[i];
            
                if ((TSbot == 0) && (data == 0)) goto loop;
                if ((TSbot == 0x5e5e5e5e) && (data == 0x5e5e5e5e)) goto loop;
                if ((TSbot == 0x5e5e5e5e) && (data == 0xffffffff)) goto loop;
                
                /* DATA */
                if ((data & 0xc0000000) == 0xc0000000) 
                {
                    ident = (data & 0x0fff0000) >> 16;
                    adcdata = (data & 0x0000ffff);
                    card = (ident / 32);
                    if (card == 99) goto skip;
                    TS = (TS & 0x0000fffff0000000ULL);
                    TS = ((TS | TSbot));
                    if (counter == 0)
                    {
                        TSfirst = TS;
                        counter = 1;
                    }
                    if (item == 0) TSinit = TS;
                    // TSdiff is the absolute value of TS-TSinit
                    TSdiff = TS-TSinit;

                    //std::cout << "TSinit " << TSinit << " TS " << TS << " TSdiff " << TSdiff << std::endl;
                    if ( (TSbot >= 0x0) && (TSbot <= 0x1a) )
                    {
                        if (twidset == 0) 
                        {
                        TS=TS+0x10000000;
                        twidset=1;
                        }          
                        TSdiff = TS-TSinit;  
                    }

                    //Time gates
                    if (TSdiff < Window) 
                    {
                        /*
                        Channel Info:
                        ---------------------------------------------------------
                        Signal  |  fastT | fastE  | slowT   |  slowE  |    HV   |    
                        ---------------------------------------------------------
                        LaBr 0  |    88  |   72   |   80    |    64   | -1050   | 
                        ---------------------------------------------------------
                        LaBr 1  |    89  |   73   |    81   |     65  | -1050   |   
                        ---------------------------------------------------------
                        LaBr 2  |    90  |   74   |    82   |    66   | -1050   |  
                        ---------------------------------------------------------
                        LaBr 3  |    91  |   75   |    83   |    67   | -1050   |  
                        --------------------------------------------------------------------
                        RF      |    92  |   76   |   //    |   //    | structure: 5 bunches 100pA  |   
                        --------------------------------------------------------------------
                        Polaris |    95  |   79   |   //    |   //    | structure: ±2s intervals    |   
                        ---------------------------------------------------------------------
                        */
                       // Loop through item values
                        //std::cout << "ident: " << ident << std::endl;
                        // Identify detector and store data accordingly
                        if ((ident > 63) && (ident < 68)) 
                        {
                            detectorID = ident-64;
                            rndUnif = randy->Rndm()*bw-(bw/2.); 
                            rndchannel = adcdata+rndUnif; 
                            energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);
                            if (slowECalibL[detectorID].size() == 0)
                            {
                                slowECalibL[detectorID].push_back(energycalib);
                            }
                            if (energyS[detectorID].size() == 0)
                            {
                                energyS[detectorID].push_back(energycalib);
                            }
                            
                            tock = ((adcdata & 0x0000e000) >> 13);
                            tick = (adcdata & 0x00001fff);
                            sick = (tick/8192.0) * 20.0; 
                            if (tock != 7)
                            {
                                if (timeS[detectorID].size() == 0)
                                {
                                    timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                                }
                            }
                        } 
                        // else if ((ident > 79) && (ident < 84)) 
                        // {
                        //     detectorID = ident-80;
                        //     tock = ((adcdata & 0x0000e000) >> 13);
                        //     tick = (adcdata & 0x00001fff);
                        //     sick = (tick/8192.0) * 20.0; 
                        //     if (tock != 7)
                        //     {
                        //         if (timeS[detectorID].size() == 0)
                        //         {
                        //             timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                        //             //printf("timeS[%d] = %f\n", detectorID, timeS[detectorID][0]);
                        //         }
                        //     }
                        // }
                        else if ((ident > 87) && (ident < 92)) 
                        {
                            detectorID = ident-72;
                            if (energyF[detectorID].size() == 0)
                            {
                                energyF[detectorID].push_back(adcdata);
                            }

                            detectorID = ident-88;
                            tock = ((adcdata & 0x0000e000) >> 13);
                            tick = (adcdata & 0x00001fff);
                            sick = (tick/8192.0) * 20.0; 
                            if (tock != 7)
                            {
                                if (timeF[detectorID].size() == 0)
                                {
                                    timeF[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                                    //printf("timeF[%d] = %f\n", detectorID, timeF[detectorID][0]);
                                }
                            }
                        }
                        else if (ident == 95) 
                        {
                            detectorID = ident-95;
                            tock = ((adcdata & 0x0000e000) >> 13);
                            tick = (adcdata & 0x00001fff);
                            sick = (tick/8192.0) * 20.0; 
                            if (slowEPOLARIS.size() == 0)
                            {
                                slowEPOLARIS[detectorID] = adcdata;
                                slowEPOLARIS.push_back(slowEPOLARIS[0]);
                            }
                            if (slowTPOLARIS.size() == 0)
                            {
                                slowTPOLARIS[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick)/10;
                                slowTPOLARIS.push_back(slowTPOLARIS[0]);
                            }
                            //print the  difference between slowT[0]-fastT[4], slowT[1]-fastT[4], slowT[2]-fastT[4], slowT[3] -fastT[4]                           
                            //std::cout << "slowT[0]-fastT[4] " << timeSlow[0]-timeFast[4] << " slowT[1]-fastT[4] " << timeSlow[1]-timeFast[4] << " slowT[2]-fastT[4] " << timeSlow[2]-timeFast[4] << " slowT[3]-fastT[4] " << timeSlow[3]-timeFast[4] << std::endl;                        
                            //std::cout << "| item " << item << std::endl;
                        }
                        if (item == 0) TSinit=TS;
                        item++;
                    }

                }

                /* SYNC */
                if ((data & 0xc0f00000) == 0x80400000) 
                {
                    card = (data & 0x3f000000) >> 24;
                    TStop = (data & 0x000fffff);
                    TS = (TStop);
                    TS = ((TS << 28));
                    TS = ((TS | TSbot));

                    SYNC = TS;
                    SYNCdiff = SYNC-SYNClast;

                    twidset=0;

                    SYNClast = SYNC;

                }

                /* If we exceed the time different from the first timestamp then do some analysis */
                if (TSdiff > Window)                     
                {   
                    i=i-2; // go back to the start of the event
                    //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;

                    for (j = 0; j <=4; j++)
                    {
                        //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                        //if the size of avector is 0, fill the value with 0
                        if (energyS[j].size() == 0) energyS[j].push_back(0);
                        if (energyF[j].size() == 0) energyF[j].push_back(0);
                        if (timeS[j].size() == 0)  timeS[j].push_back(0);
                        if (timeF[j].size() == 0) timeF[j].push_back(0);
                        if (timeS[j].size() == 0)  timeS[j].push_back(0);

                        timeFast[j] = timeF[j][0];
                        timeSlow[j] = timeS[j][0];
                        timeSlow[j] = timeS[j][0];
                        energySlow[j] = energyS[j][0];
                        energyFast[j] = energyF[j][0];
                    }

                    for (j = 0; j <4; j++)
                    {
                        if (timeFast[j] != 0 && timeSlow[j] == 0) fastcount[j]++;
                        else if (timeSlow[j] != 0 && timeFast[j] == 0) slowcount[j]++;
                        else if (timeSlow[j] != 0 && timeFast[j] != 0) fastslowcount[j]++;

                        eventCounts[j]->Fill(1,fastcount[j]);
                        eventCounts[j]->Fill(2,slowcount[j]);
                        eventCounts[j]->Fill(3,fastslowcount[j]);
                        
                    }

                    slowEnergyPOLARIS[0] = slowEPOLARIS[0];
                    slowTimePOLARIS[0] = slowTPOLARIS[0];

                    LaBrData->Fill();

                    for (j = 0; j <= 4; j++) 
                    { 
                        //Singles spectra 
                        // in order to Fill ROOT histograms we need to convert the vectors to doubles
                        //std::cout << "L" << j << "| energySlow[j] " << energySlow[j] << "| timeFast[j] " << timeFast[j] << "| timeSlow[j] " << timeSlow[j] << "| timeRF" << timeFast[4]<< std::endl;
                        slowE[j]->Fill(energySlow[j]);
                        hit1->Fill(j);
                        
                        // std::cout << "L" << j << "| energySlow[j] " << energySlow[j] << "| timeFast[j] " << timeFast[j] << "| timeSlow[j] " << timeSlow[j] << "| timeRF" << timeFast[4] << std::endl;
                        // ** Uncomment all lines containing outputFile to see what the data looks like **
                        //outputFile << j << " "  << energySlow[j] << " " << timeSlow[j] << " " << timeFast[j] << std::endl; // To store each value on a new line in the text file.
                        
                        for (k = 0; k <= 4; k++) 
                        {
                            if (j != k ) 
                            {
                                tdF=(timeFast[j]-timeFast[k]); // Measured in 1e13s which is 100ns
                                tdF=((tdF+(-tcalibF[j]+tcalibF[k]))/10)+250.0;
                                tdS = (timeSlow[j]-timeSlow[k]);
                                tdS = ((tdS+(-tcalibS[j]+tcalibS[k]))/10)+250.0; // Measured in 1e13s which is 100ns
                                tdSF = (timeSlow[j]-timeFast[j])/10;
                                slowfastTD[j] -> Fill(tdSF);
                                if (j == 0) {fastTD0[k]->Fill(tdF); slowTD0[k]->Fill(tdS);}
                                if (j == 1) {fastTD1[k]->Fill(tdF); slowTD1[k]->Fill(tdS);}
                                if (j == 2) {fastTD2[k]->Fill(tdF); slowTD2[k]->Fill(tdS);}
                                if (j == 3) {fastTD3[k]->Fill(tdF); slowTD3[k]->Fill(tdS);}
                                if (j == 4) 
                                {
                                    fastTD4[k]->Fill(tdF);
                                    slowTD4[k]->Fill(tdS);
                                    fastTDslowE[k]->Fill(tdF,energySlow[k]);

                                    //printf("tdF %f | energySlow[k] %f \n", tdF,energySlow[k]);
                                    
                                    // if ((tdF >= 160) && (tdF <= 240)) 
                                    // if (((tdF >= 137) && (tdF <= 239)) || ((tdF >= 444) && (tdF <= 540)) || ((tdF >= 736) && (tdF <= 804)))
                                    if (((tdF >= 155) && (tdF <= 238)) || ((tdF >= 460) && (tdF <= 802)))
                                    //if (((tdF >= 137) && (tdF <= 239)) || ((tdF >= -154) && (tdF <= -71)) || ((tdF >= 444) && (tdF <= 540)) || ((tdF >= 750) && (tdF <= 804)))
                                    {
                                        insync[k]->Fill(energySlow[k]);
                                        insyncEnergySlow[k] = energySlow[k];
                                        insyncTimeFast[k] = timeFast[k];
                                        //LaBrData->Fill();
                                    }
                                    
                                    // if ((tdF < -154) || (((tdF > -71) && (tdF < 137))) || (((tdF > 239) && (tdF < 444))) || (((tdF > 540) && (tdF < 554))) || ((tdF > 804))) 
                                    // {
                                    //     outsync[k]->Fill(energySlow[k]);
                                    // }

                                    if ((energySlow[k] >= 1000)&&(energySlow[k] < 2000)) 
                                    {
                                        oneMeVgate[k]->Fill(tdF);
                                    }
                                    if ((energySlow[k] >= 2000)&&(energySlow[k] < 3000)) 
                                    {
                                        twoMeVgate[k]->Fill(tdF);
                                    }
                                    if ((energySlow[k] >= 3000)&&(energySlow[k] < 4000)) 
                                    {
                                        threeMeVgate[k]->Fill(tdF);
                                    }
                                    if ((energySlow[k] >= 4000)&&(energySlow[k] < 5000)) 
                                    {
                                        fourMeVgate[k]->Fill(tdF);
                                    }
                                    if ((energySlow[k] >= 5000)&&(energySlow[k] < 6000)) 
                                    {
                                        fiveMeVgate[k]->Fill(tdF);
                                    }
                                    if ((energySlow[k] >= 6000)&&(energySlow[k] < 7000)) 
                                    {
                                        sixMeVgate[k]->Fill(tdF);
                                    }
                                }
                            }
                        }                        
                    }


                    // Reset vectors
                    for (j = 0; j <=4; j++) 
                    {
                        energyS[j].clear();
                        timeF[j].clear();
                        timeS[j].clear();
                        slowECalibL[j].clear();
                    }
                    slowEPOLARIS.clear();
                    slowTPOLARIS.clear();

                    //Reinitialise arrays
                    for (j = 0; j <=4; j++)
                    { 
                        energySlow[j] = 0;
                        timeFast[j] = 0;
                        timeSlow[j] = 0;
                        slowEnergyPOLARIS[j] = 0;
                        slowTimePOLARIS[j] = 0;
                        insyncEnergySlow[j] = 0;
                        insyncTimeFast[j] = 0;
                    }

                    item = 0;
                    //printf("\n");
                    
                }

                
                skip:
                count++;
                loop:
                pos+=8; // 8 bytes per event
            }
        }
        end:  
        //std::cout << "End of file reached" << std::endl;
        fclose(f);
        
        Run++;

    }
    finish:
    
    LaBrData->Write();

    for (j = 0; j <=4; j++)
    {
        slowE[j]->Write();
        slowfastTD[j]->Write();
        fastTD0[j]->Write();
        fastTD1[j]->Write();
        fastTD2[j]->Write();
        fastTD3[j]->Write();
        fastTD4[j]->Write();
        slowTD0[j]->Write();
        slowTD1[j]->Write();
        slowTD2[j]->Write();
        slowTD3[j]->Write();
        slowTD4[j]->Write();
        fastTDslowE[j]->Write();
        insync[j]->Write();
        outsync[j]->Write();

        oneMeVgate[j]->Write();
        twoMeVgate[j]->Write();
        threeMeVgate[j]->Write();
        fourMeVgate[j]->Write();
        fiveMeVgate[j]->Write();
        sixMeVgate[j]->Write();
        eventCounts[j]->Write();
    }
    
    hit1->Write();

    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference is (min): %f\n", TimeDiff);

    time(&end);
    double time_taken = double((end - start)/60.0);
    std::cout << "Time taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " min " << std::endl;
    
    //outputFile.close();
    delete g;
}