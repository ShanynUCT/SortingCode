/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (29/03/2022) - RF corrected timing and Trees  */

/* COMPILE: g++ -std=c++0x labrsort5_v1.C -o exe `root-config --cflags --libs` -lSpectrum */
/* RUN: ./exe /path/to/raw/experiment/data/ (e.g. /home/shanyn/Data/R50) runXX (PASS THE SPECIFIC NAME OF THE RUN AS ARGV[2] IN ORDER TO CREATE THE CSV FILE WITH THE CORRECT NAME)*/

#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <inttypes.h> 
#include <time.h>
#include <cstdlib>

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

#include <pthread.h>
#define SIZE  16384

int32_t buffer[SIZE];

uint64_t TStop, TSbot, data;
uint64_t TS, TSinit=0, TSfirst, counter,  SYNC, SYNClast=0;
int64_t item;
uint64_t TimeStamp0[100]={0}, TimeStamp1[100]={0}, TimeStamp2[100]={0} ,TimeStamp3[100]={0}, TimeStamp4[100]={0};
uint16_t Energy[100]={0};
uint16_t QDC[100]={0}, qdc;
uint64_t TS2;

double T1, T2;

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
  time(&start);

  // Data arrays for events
  double energyF[5] = {-1,-1,-1,-1,-1};
  double energyS[5] = {-1,-1,-1,-1,-1};
  double timeF[5] = {-1,-1,-1,-1,-1};
  double timeS[5] = {-1,-1,-1,-1,-1};
  double tdF, tdS, tdSF;
  int Ed;

  // Define coincidence window here in 10's of ns
  int64_t Window = 23; // window = 1/2 * desired window size 

  //Energy calibrations 1st order polynomials
  double a1[4], b1[4];
  a1[0] = 0.75263259; // 15/06/2022
  b1[0] = 0.71740579; // 15/06/2022

  a1[1] = 0.81401421; // 15/06/2022
  b1[1] = -0.54088457; // 15/06/2022

  a1[2] = 0.51735782; // 15/06/2022
  b1[2] = 0.41341115; // 15/06/2022

  a1[3] = 0.71451586; // 15/06/2022
  b1[3] = 1.14456828; // 15/06/2022

  //Energy calibrations 2nd order polynomials (APPLIED)
  double a2[4], b2[4], c2[4];
  a2[0] = 0.00000288786796; // 15/06/2022 
  b2[0] = 0.747696355; // 15/06/2022
  c2[0] = 1.64893477; // 15/06/2022

  a2[1] = 0.00000351711116; // 15/06/2022
  b2[1] = 0.805816609; // 15/06/2022
  c2[1] = 1.09470037; // 15/06/2022

  a2[2] = 0.00000160595876; // 15/06/2022
  b2[2] = 0.513776183; // 15/06/2022
  c2[2] = 1.37154046; // 15/06/2022

  a2[3] = 0.00000300000230; // 15/06/2022
  b2[3] = 0.709359122; // 15/06/2022
  c2[3] = 2.18122723;  // 15/06/2022

  // a2[0] = 0.00001942 ; // 28/02/2023
  // b2[0] = 1.81799962; // 28/02/2023
  // c2[0] = -1.65129151; // 28/02/2023

  // a2[1] = 0.00010155; // 28/02/2023
  // b2[1] = 2.61697133; // 28/02/2023
  // c2[1] = -0.54900373; // 28/02/2023

  // a2[2] = 0.00003208; // 28/02/2023
  // b2[2] = 2.21345098; // 28/02/2023
  // c2[2] = -4.86476637; // 28/02/2023

  // a2[3] = 0.00002200; // 28/02/2023
  // b2[3] = 2.05179483; // 28/02/2023
  // c2[3] = -0.97991876;  // 28/02/2023

  double bw = 1, energycalib, rndchannel, rndUnif;
  TRandom3* randy = new TRandom3();


  //Time calibrations
  double tcalibF[5];
  tcalibF[0]=   0;
  tcalibF[1]=   120;
  tcalibF[2]=   150;
  tcalibF[3]=   89;
  tcalibF[4]=   0;

  double tcalibS[5];
  tcalibS[0]=   0.00;
  tcalibS[1]=   0;
  tcalibS[2]=   0;
  tcalibS[3]=   0;
  tcalibS[4]=   0.00;
    
  int no_read;
  std::vector<int> time_stamp;


  if ((argc < 2) || (argc > 3))
  {
    fprintf(stderr, "Usage: %s data [last_subrun]\n", argv[0]);
    exit(1);
  }

  /* ROOT STUFF */

  sprintf(file_out, "%s.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);

  TFile *g = new TFile(file_out,"recreate");

  printf("ROOT file %s opened...\n", file_out);


  // ################ Create ROOT TTrees ################

  double insyncSlowECalibL[4] = {-1,-1,-1,-1};
  double slowECalibL[4] = {-1,-1,-1,-1};
  double fastTimeL[4] = {-1,-1,-1,-1};
  double slowTimeL[4] = {-1,-1,-1,-1}; 
  double slowEPOLARIS, slowTimePOLARIS;

  TTree *LaBr0 = new TTree("LaBr0","LaBr0");
  LaBr0->Branch("insyncSlowECalibL0", &insyncSlowECalibL[0]);
  LaBr0->Branch("slowECalibL0", &slowECalibL[0]);
  LaBr0->Branch("fastTimeL0", &fastTimeL[0]);
  LaBr0->Branch("slowTimeL0", &slowTimeL[0]);
  LaBr0->Branch("slowEPOLARIS", &slowEPOLARIS);
  LaBr0->Branch("slowTimePOLARIS", &slowTimePOLARIS);

  TTree *LaBr1 = new TTree("LaBr1","LaBr1");
  LaBr1->Branch("insyncSlowECalibL1", &insyncSlowECalibL[1]);
  LaBr1->Branch("slowECalibL1", &slowECalibL[1]);
  LaBr1->Branch("fastTimeL1", &fastTimeL[1]);
  LaBr1->Branch("slowTimeL1", &slowTimeL[1]);
  LaBr1->Branch("slowEPOLARIS", &slowEPOLARIS);
  LaBr1->Branch("slowTimePOLARIS", &slowTimePOLARIS);

  TTree *LaBr2 = new TTree("LaBr2","LaBr2");
  LaBr2->Branch("insyncSlowECalibL2", &insyncSlowECalibL[2]);
  LaBr2->Branch("slowECalibL2", &slowECalibL[2]);
  LaBr2->Branch("fastTimeL2", &fastTimeL[2]);
  LaBr2->Branch("slowTimeL2", &slowTimeL[2]);
  LaBr2->Branch("slowEPOLARIS", &slowEPOLARIS);
  LaBr2->Branch("slowTimePOLARIS", &slowTimePOLARIS);

  TTree *LaBr3 = new TTree("LaBr3","LaBr3");
  LaBr3->Branch("insyncSlowECalibL3", &insyncSlowECalibL[3]);
  LaBr3->Branch("slowECalibL3", &slowECalibL[3]);
  LaBr3->Branch("fastTimeL3", &fastTimeL[3]);
  LaBr3->Branch("slowTimeL3", &slowTimeL[3]);
  LaBr3->Branch("slowEPOLARIS", &slowEPOLARIS);
  LaBr3->Branch("slowTimePOLARIS", &slowTimePOLARIS);

  // ################ Create histograms ################
  TH1D** fastE=new TH1D*[5];
  TH1D** slowE=new TH1D*[5];

  TH1D** coincSlowE1=new TH1D*[5];
  TH1D** coincSlowE2=new TH1D*[5];

  TH1D** fastTD0=new TH1D*[5];
  TH1D** fastTD1=new TH1D*[5];
  TH1D** fastTD2=new TH1D*[5];
  TH1D** fastTD3=new TH1D*[5];
  TH1D** fastTD4=new TH1D*[5];

  TH1D** slowTD0=new TH1D*[5]; //NB: 4 is the number of detectors, RF has no slow signal but is included in the loop
  TH1D** slowTD1=new TH1D*[5];
  TH1D** slowTD2=new TH1D*[5];
  TH1D** slowTD3=new TH1D*[5];
  TH1D** slowTD4=new TH1D*[5];
  TH1D** insync=new TH1D*[5];

  TH1D** multSlowE=new TH1D*[5];
  TH2D** fastTDslowE = new TH2D*[5];

  for (j = 0; j <=4; j++)
  {
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%02d", j),"Spectrum",8000,0,8000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%02d", j),"Spectrum",8000,0,8000);

    coincSlowE1[j] = new TH1D(TString::Format("Coinc1_Energy_L%2d", j),"Spectrum",2000,0,2000);
    coincSlowE2[j] = new TH1D(TString::Format("Coinc2_Energy_L%2d", j),"Spectrum",2000,0,2000);

    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",600,0,600);
    fastTD1[j] = new TH1D(TString::Format("Fast Time Difference_L1-L%2d", j),"Spectrum",600,0,600);
    fastTD2[j] = new TH1D(TString::Format("Fast Time Difference_L2-L%2d", j),"Spectrum",600,0,600);
    fastTD3[j] = new TH1D(TString::Format("Fast Time Difference_L3-L%2d", j),"Spectrum",600,0,600);
    fastTD4[j] = new TH1D(TString::Format("Fast Time Difference_RF-L%2d", j),"Spectrum",600,0,600);

    slowTD0[j] = new TH1D(TString::Format("Slow Time Difference_L0-L%2d", j),"Spectrum",600,0,600);
    slowTD1[j] = new TH1D(TString::Format("Slow Time Difference_L1-L%2d", j),"Spectrum",600,0,600);
    slowTD2[j] = new TH1D(TString::Format("Slow Time Difference_L2-L%2d", j),"Spectrum",600,0,600);
    slowTD3[j] = new TH1D(TString::Format("Slow Time Difference_L3-L%2d", j),"Spectrum",600,0,600);
    slowTD4[j] = new TH1D(TString::Format("Slow Time Difference_RF-L%2d", j),"Spectrum",600,0,600);

    multSlowE[j] = new TH1D(TString::Format("Multiplicity_Slow_Energy_%2d", j),"Spectrum",2500,0,2500);
    insync[j] = new TH1D(TString::Format("Insync_Slow_Energy_%2d_range = 166-248ns", j),"Spectrum",2500,0,2500);
    fastTDslowE[j] = new TH2D(TString::Format("Time_Difference_L%2d_Vs_Slow_Energy", j), TString::Format("Time L%2d - Vs Energy Spectrum", j),600,0,600,8000,0,8000);
  }
    
  //  TH1D *BumpL3E = new TH1D("Bump L3","Bump: 380-513 ns - Detector L3",12000,0,12000);   
  //  TH1D *TailGauss = new TH1D("Tail Gauss","Tail: 227-245 ns - Detector L0",8000,0,8000);
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
        
        // print what the file contains
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

          TSdiff = TS-TSinit;

          //std::cout << "ident: " << ident << " | adcdata: " << adcdata << " | TS: " << TS << " | TSdiff: " << TSdiff << " | TSinit: " << TSinit << std::endl;

          if ( (TSbot >= 0x0) && (TSbot <= 0x1a) )
          {
            if (twidset == 0) 
            {
              TS=TS+0x10000000;
              twidset=1;
            }
            
            TSdiff = TS-TSinit;
          }

          //if ((TS - TSfirst)*(1e-8)/60.0 > 7)	  goto finish; // Here we cut data after x minute
      
          //Time gates
          if (TSdiff < Window) // Commented out to obtain full energy and time spectra before time gates
          {
            //if (item == 0) printf("\n");
            //printf("Event %ld -> %d : %" PRIx64 " %3" PRIu64 " %" PRIx64 "", item, ident, TS, TSdiff, TSinit);
            //if (ident == 76) printf("  ******");
            //printf("\n");

            //Ident 64,65,66,67,68,69,70,71 are slow signals 
            if ((ident > 63) && (ident < 68))
            {
              detectorID = ident-64;
              rndUnif = randy->Rndm()*bw-(bw/2.); 
              rndchannel = adcdata+rndUnif; 
              //std::cout << "energy " << adcdata << " rndchannel " << rndchannel << std::endl;
              energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);
              if (energyS[detectorID] == -1) 
              {
                energyS[detectorID] = energycalib;
                slowECalibL[detectorID] = energycalib;
                if (detectorID==0) LaBr0->Fill();
                if (detectorID==1) LaBr1->Fill();
                if (detectorID==2) LaBr2->Fill();
                if (detectorID==3) LaBr3->Fill();
              }
            }

            //Ident 72,73,74,75,76,77,78,79 are fast signals
            if ((ident > 71) && (ident < 77))
            {
              detectorID = ident-72;
              rndUnif = randy->Rndm()*bw-(bw/2.);
              rndchannel = adcdata+rndUnif; 
              energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);		  
              if (energyF[detectorID] == -1) 
              {
                energyF[detectorID] = energycalib;
              }
            }

            if (ident == 79)
            {
              slowTimePOLARIS = TS;
              slowEPOLARIS = adcdata;
              LaBr0->Fill();
              LaBr1->Fill();
              LaBr2->Fill();
              LaBr3->Fill();
            }

            //Ident 80,81,82,83,84,85,86,87 are timing signals of slow signals (RF has no slow signal)
            if ((ident > 79) && (ident < 84))
            {
              detectorID = ident-80;
              
              tock = ((adcdata & 0x0000e000) >> 13);
              tick = (adcdata & 0x00001fff);

              sick = (tick/8192.0) * 20.0; // here the 20ns is the time resolution of the ADC (8192 channels)

              if (tock != 7)
              {
                if (timeS[detectorID] == -1) 
                {
                  timeS[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick);
                  slowTimeL[detectorID] = timeS[detectorID];
                }
                if (detectorID == 0) LaBr0->Fill();
                if (detectorID == 1) LaBr1->Fill();
                if (detectorID == 2) LaBr2->Fill();
                if (detectorID == 3) LaBr3->Fill();
              }
            }

            //Ident 88,89,90,91,92,93,94,95 are timing signals of fast signals
            if ((ident > 87) && (ident < 93))
            {
              detectorID = ident-88;
              
              tock = ((adcdata & 0x0000e000) >> 13);
              tick = (adcdata & 0x00001fff);

              sick = (tick/8192.0) * 20.0; 

              if (tock != 7)
              {
                if (timeF[detectorID] == -1) 
                {
                  timeF[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick);
                  fastTimeL[detectorID] = timeF[detectorID];
                }
                if (detectorID == 0) LaBr0->Fill();
                if (detectorID == 1) LaBr1->Fill();
                if (detectorID == 2) LaBr2->Fill();
                if (detectorID == 3) LaBr3->Fill();
              }
            }
          }
          if (item == 0) TSinit=TS;
          item++;
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

          //Be kind, rewind one
          i=i-2;

          T1=0;
          T2=0;

          //printf("\n");
          for (j = 0; j <= 4; j++)
          { 
          //Singles spectra
          {            
            //std::cout << "j = " << j << "| energyF = " << energyF[j] <<  " | timeF = " << timeF[j] << " | energyS = " << energyS[j] << " | timeS = " << timeS[j] << " | timeF - timeS = " << timeF[j]-timeS[j] << std::endl;
            fastE[j]->Fill(energyF[j]);
            slowE[j]->Fill(energyS[j]);
            hit1->Fill(j);
            

            //Time difference spectra
            for (k = 0; k <= 4; k++)
            {
              {
                if (j != k )
                {
                  tdF=(timeF[j]-timeF[k]); //here
                  tdF= (((tdF+(-tcalibF[j]+tcalibF[k]))/10)+250.0);
                  //std::cout << "k = " << k << "| energyF = " << energyF[k] <<  " | timeF = " << timeF[k] << " | energyS = " << energyS[k] << " | timeS = " << timeS[k] << " | timeF - timeS = " << timeF[k]-timeS[k] << std::endl;
                  // std::cout << "j = " << j << "| k = " << k << "| timeF[j] = " << timeF[j] << " | timeF[k] = " << timeF[k] << "| tdF = " << tdF << std::endl;
                  
                  if (j == 0) fastTD0[k]->Fill(tdF);
                  if (j == 1) fastTD1[k]->Fill(tdF);
                  if (j == 2) fastTD2[k]->Fill(tdF);
                  if (j == 3) fastTD3[k]->Fill(tdF);
                  if (j == 4) 
                  { 
                    fastTD4[k]->Fill(tdF); 
                    fastTDslowE[k]->Fill(tdF,energyS[k]); 
                    if ((tdF >= 166)&&(tdF<=248))
                    {
                      insync[k]->Fill(energyS[k]);
                    }
                  }

                  tdS = (timeS[j]-timeS[k]);
                  tdS = (((tdS+(-tcalibS[k]+tcalibS[j]))/10)+250.0); // Measured in 1e13s which is 100ns

                  if (j == 0) slowTD0[k]->Fill(tdS);
                  if (j == 1) slowTD1[k]->Fill(tdS);
                  if (j == 2) slowTD2[k]->Fill(tdS);
                  if (j == 3) slowTD3[k]->Fill(tdS);
                  if (j == 4)
                  {
                    tdS = (timeF[j]-timeS[k]);
                    tdS = ((tdS+(-tcalibF[j]+tcalibS[k]))/10)+250.0; // Measured in 1e13s which is 100ns
                    slowTD4[k]->Fill(tdS);
                  }
                }
              }
            }
          }
        }

          // More Data analysis here

          m1=0;
          m2=0;
          sum1=0;
          sum2=0;
        
          for (j = 0; j <=4; j++)
          { 
            for (k = j+1; k <=4; k++)
            { 
              tdF=(timeF[j]-timeF[k]); // here
              tdF=((tdF+(-tcalibF[j]+tcalibF[k]))/10)+250.0;
              //std::cout << "tdF = " << tdF << std::endl;
              
              if ( (tdF >= 249) && (tdF <= 251) )
              {
                coincSlowE1[j]->Fill(energyS[j]);
                coincSlowE2[j]->Fill(energyS[k]);

                hit2->Fill(j);
                hit2->Fill(k);

                m1=m1+1;

                sum1=sum1+energyS[j]+energyS[k];
            
                if ((energyS[j] > 471) && (energyS[j] < 551))
                {
                  if ((energyS[k] > 471) && (energyS[k] < 551))
                  {
                    hit3->Fill(j);
                    hit3->Fill(k);

                    m2=m2+1;
                    sum2=sum2+energyS[j]+energyS[k];
                  }
                }
              }
            }
          }

          mul1->Fill(m1);
          mul2->Fill(m2);

          if ((sum1 > 942) && (sum1 < 1102))
          {
            mul3->Fill(m1);
          }

        
          if (sum1 > 0) sumspec0->Fill(sum1);
          if (sum2 > 0) sumspec1->Fill(sum2);

          if (m1 == 0)  sumspec00->Fill(sum1);
          if (m1 == 1)  sumspec01->Fill(sum1);
          if (m1 == 2)  sumspec02->Fill(sum1);
          if (m1 == 3)  sumspec03->Fill(sum1);

        
          for (j = 0; j <=4; j++)
          { 
            if ((m1>0) && (m1<5)) multSlowE[m1]->Fill(energyS[j]);
          }

          //Reinitialise arrays
          for (j = 0; j <=4; j++)
          { 
            energyS[j] = -1;
            energyF[j] = -1;
            timeF[j] = -1;
            timeS[j] = -1;
            slowECalibL[j] = -1;
            insyncSlowECalibL[j] = -1;
            fastTimeL[j] = -1;
            slowTimeL[j] = -1;
          }
        
          //Not so clever
          item=0;
        }
    
        skip:

        count++;

        loop:
        pos+=8; // 8 bytes per event
      }
    }

    end:  

    fclose(f);
    
    Run++;
      
  }

  finish:

  double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;

  printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
  printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
  printf("The time difference is (min): %f\n", TimeDiff);


  for (j = 0; j <=4; j++)
  {
    fastE[j]->Write();
    slowE[j]->Write();

    coincSlowE1[j]->Write();
    coincSlowE2[j]->Write();

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
    insync[j]->Write();

    multSlowE[j]->Write();
    fastTDslowE[j]->Write();
  }

  LaBr0->Write();
  LaBr1->Write();
  LaBr2->Write();
  LaBr3->Write();

  mul1->Write();
  mul2->Write();
  mul3->Write();
  hit1->Write();
  hit2->Write();

  sumspec0->Write();
  sumspec00->Write();
  sumspec01->Write();
  sumspec02->Write();
  sumspec03->Write();
  sumspec1->Write();

  // take the ratio of the number of events in the slowE[j] histogram and the number of events in the insync[j] histogram
  //print this value to the screen
  for (j = 0; j <=4; j++)
  {
    double ratio = (insync[j]->GetEntries()/slowE[j]->GetEntries())*100;
    std::cout << "Ratio of insync[" << j << "] to total[" << j << "] is " << ratio << std::endl;
  }


  time(&end);

  // Calculating total time taken by the program.
  double time_taken = double((end - start)/60.0);
  std::cout << "Time taken by program is : " << time_taken << std::setprecision(5);
  std::cout << " min " << std::endl;

  delete g;
}
