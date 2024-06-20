/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (29/03/2022) - RF corrected timing and Trees  */

/* COMPILE: g++ -std=c++0x labrsort5_v3.2.C -o exe `root-config --cflags --libs` -lSpectrum */
/* RUN: ./exe /path/to/raw/experiment/data/ (e.g. /home/shanyn/Data/R50) runXX (PASS THE SPECIFIC NAME OF THE RUN AS ARGV[2] IN ORDER TO CREATE THE CSV FILE WITH THE CORRECT NAME)*/
//      ./exe /Users/shanyn/Documents/PhD/ExperimentResults/2022/220615/analysis/ExperimentAngleAnalysis/WaterTarget/run27-66MeV_Proton_beam-combo_run-D0_0deg_18cm-Water_Target_2cm_depth-20min/R63
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
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

  /*  W  I  N  D  O  W  */
  // Define coincidence window here in 10's of ns
  int64_t Window = 65; // window = 1/2 * desired window size 22 //60, 65  good //
  /*  W  I  N  D  O  W  */

  // Declare vectors to store the event data
  std::vector<double> energyS[5];
  std::vector<double> timeF[5];
  std::vector<double> timeS[5];
  std::vector<double> slowECalibL[5];
  std::vector<double> slowEPOLARIS;
  std::vector<double> slowTimePOLARIS;

  double tdF, tdS, tdSF;
  int Ed;
  
  //Energy calibrations 2nd order polynomials (APPLIED)
  std::vector<double> a2 = {0.00000288786796, 0.00000351711116, 0.00000160595876, 0.00000300000230};
  std::vector<double> b2 = {0.747696355, 0.805816609, 0.513776183, 0.709359122};
  std::vector<double> c2 = {1.64893477, 1.09470037, 1.37154046, 2.18122723};
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

  sprintf(file_out, "%s.root", argv[1]);
  if (argc > 2) sscanf(argv[2], "%d", &RunEnd);
  TFile *g = new TFile(file_out,"recreate");
  printf("ROOT file %s opened...\n", file_out);

  // ################ Create ROOT TTrees ################
  TTree *LaBrData = new TTree("LaBrData","LaBrData");
  LaBrData->Branch("slowECalibL0", &slowECalibL[0], "slowECalibL0/D");
  LaBrData->Branch("fastTimeL0", &timeF[0], "fastTimeL0/D");
  LaBrData->Branch("slowTimeL0", &timeS[0], "slowTimeL0/D");
  LaBrData->Branch("slowECalibL1", &slowECalibL[1], "slowECalibL1/D");
  LaBrData->Branch("fastTimeL1", &timeF[1], "fastTimeL1/D");
  LaBrData->Branch("slowTimeL1", &timeS[1], "slowTimeL1/D");
  LaBrData->Branch("slowECalibL2", &slowECalibL[2], "slowECalibL2/D");
  LaBrData->Branch("fastTimeL2", &timeF[2], "fastTimeL2/D");
  LaBrData->Branch("slowTimeL2", &timeS[2], "slowTimeL2/D");
  LaBrData->Branch("slowECalibL3", &slowECalibL[3], "slowECalibL3/D");
  LaBrData->Branch("fastTimeL3", &timeF[3], "fastTimeL3/D");
  LaBrData->Branch("slowTimeL3", &timeS[3], "slowTimeL3/D");
  LaBrData->Branch("slowEPOLARIS", &slowEPOLARIS[0], "slowEPOLARIS/D");
  LaBrData->Branch("slowTimePOLARIS", &slowTimePOLARIS[0], "slowTimePOLARIS/D");
  LaBrData->Branch("timeRF", &timeF[4], "timeRF/D");

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

  TH1D** slowTD0=new TH1D*[5]; //NB: 4 is the number of detectors, RF has no slow signal but is included in the loop
  TH1D** slowTD1=new TH1D*[5];
  TH1D** slowTD2=new TH1D*[5];
  TH1D** slowTD3=new TH1D*[5];
  TH1D** slowTD4=new TH1D*[5];
  TH1D** slowfastTD=new TH1D*[5];
  TH1D** multSlowE=new TH1D*[5];

  TH2D** fastTDslowE = new TH2D*[5];
  TH1D** insync = new TH1D*[5];

  for (j = 0; j <=4; j++)
  {
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%02d", j),"Spectrum",8000,0,8000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%02d", j),"Spectrum",8000,0,8000);
    uncalibE[j] = new TH1D(TString::Format("Uncalib_Slow_Energy_L%2d", j),"Spectrum",8000,0,8000);

    coincSlowE1[j] = new TH1D(TString::Format("Coinc1_Energy_L%2d", j),"Spectrum",2000,0,2000);
    coincSlowE2[j] = new TH1D(TString::Format("Coinc2_Energy_L%2d", j),"Spectrum",2000,0,2000);

    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",600,0,600);
    fastTD1[j] = new TH1D(TString::Format("Fast Time Difference_L1-L%2d", j),"Spectrum",600,0,600);
    fastTD2[j] = new TH1D(TString::Format("Fast Time Difference_L2-L%2d", j),"Spectrum",600,0,600);
    fastTD3[j] = new TH1D(TString::Format("Fast Time Difference_L3-L%2d", j),"Spectrum",600,0,600);
    fastTD4[j] = new TH1D(TString::Format("Fast Time Difference_RF-L%2d", j),"Spectrum",1000,0,1000);

    slowTD0[j] = new TH1D(TString::Format("Slow Time Difference_L0-L%2d", j),"Spectrum",600,0,600);
    slowTD1[j] = new TH1D(TString::Format("Slow Time Difference_L1-L%2d", j),"Spectrum",600,0,600);
    slowTD2[j] = new TH1D(TString::Format("Slow Time Difference_L2-L%2d", j),"Spectrum",600,0,600);
    slowTD3[j] = new TH1D(TString::Format("Slow Time Difference_L3-L%2d", j),"Spectrum",600,0,600);
    slowTD4[j] = new TH1D(TString::Format("Slow Time Difference_RF-L%2d", j),"Spectrum",600,0,600);

    slowfastTD[j] = new TH1D(TString::Format("Slow-Fast Time Difference_L%2d", j),"Spectrum",2400,-1200,1200);
    multSlowE[j] = new TH1D(TString::Format("Multiplicity_Slow_Energy_%2d", j),"Spectrum",2500,0,2500);
    fastTDslowE[j] = new TH2D(TString::Format("Time_Difference_L%2d_Vs_Slow_Energy", j), TString::Format("Time L%2d - Vs Energy Spectrum", j),600,0,600,8000,0,8000);
    insync[j] = new TH1D(TString::Format("insync L%2d", j),TString::Format("In sync energy spectrum for TD gate 200-300 ns - Detector L%2d", j),8000,0,8000);
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

          //if ((TS - TSfirst)*(1e-8)/60.0 > 7)	  goto finish; // Here we cut data after x minute
      
          //Time gates
          if (TSdiff < Window) // Commented out to obtain full energy and time spectra before time gates
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
            for (item = 0; item < 10; ++item)  // 6 is the number of items per event
            {
              //std::cout << "ident: " << ident << std::endl;
              // Identify detector and store data accordingly
              if ((ident > 63) && (ident < 68)) 
              {
                detectorID = ident-64;
                rndUnif = randy->Rndm()*bw-(bw/2.); 
                rndchannel = adcdata+rndUnif; 
                energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);
                //if (energyS[detectorID] == std::vector<double>(-1)) 
                {
                  slowECalibL[detectorID].push_back(energycalib); 
                  energyS[detectorID].push_back(energycalib);
                }
              } 
              else if ((ident > 79) && (ident < 84)) 
              {
                detectorID = ident-80;
                tock = ((adcdata & 0x0000e000) >> 13);
                tick = (adcdata & 0x00001fff);
                sick = (tick/8192.0) * 20.0; 
                if (tock != 7)
                {
                  //if (timeS[detectorID] == std::vector<double>(-1)) 
                  {
                    timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                  }
                }
              } 
              else if ((ident > 87) && (ident < 93)) 
              {
                detectorID = ident-88;
                tock = ((adcdata & 0x0000e000) >> 13);
                tick = (adcdata & 0x00001fff);
                sick = (tick/8192.0) * 20.0; 
                if (tock != 7)
                {
                  //if (timeF[detectorID] == std::vector<double>(-1)) 
                  {
                    timeF[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                  }
                }
              }
              else if (ident == 95) 
              {
                detectorID = ident-95;
                tock = ((adcdata & 0x0000e000) >> 13);
                tick = (adcdata & 0x00001fff);
                sick = (tick/8192.0) * 20.0; 
                slowTimePOLARIS[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick)/10;
                slowEPOLARIS[detectorID] = adcdata;
                slowEPOLARIS.push_back(slowEPOLARIS[0]);
                slowTimePOLARIS.push_back(slowTimePOLARIS[0]);
              }

            }
            // Store event data in the tree
            // for (int l = 0; l <= 4; ++l) 
            // {
            //   energyS[l] = std::vector<double>(1, (energyS[l].empty()) ? -1 : energyS[l][0]);// Store only the first value
            //   slowECalibL[l] = std::vector<double>(1, (slowECalibL[l].empty()) ? -1 : slowECalibL[l][0]); // Store only the first value
            //   timeF[l] = std::vector<double>(1, (timeF[l].empty()) ? -1 : timeF[l][0]);
            //   timeS[l] = std::vector<double>(1, (timeS[l].empty()) ? -1 : timeS[l][0]);
            // }
            // slowEPOLARIS = std::vector<double>(1, (slowEPOLARIS.empty()) ? -1 : slowEPOLARIS[0]);
            // slowTimePOLARIS = std::vector<double>(1, (slowTimePOLARIS.empty()) ? -1 : slowTimePOLARIS[0]);
            LaBrData->Fill();
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
            //if ((timeS[j] != -1) || (timeF[j] != -1))
            //Singles spectra 
            // in order to Fill ROOT histograms we need to convert the vectors to doubles
            double energyS[j], timeF[j], timeS[j];
            //std::cout << "L" << j << "| energyS[j] " << energyS[j] << "| timeF[j] " << timeF[j] << "| timeS[j] " << timeS[j] << "| timeRF" << timeF[4]<< std::endl;
            slowE[j]->Fill(energyS[j]);
            hit1->Fill(j);

            //Time difference spectra
            for (k = 0; k <= 4; k++)
            { 
              if (j != k )
              {
                tdF=(timeF[j]-timeF[k]); //here
                tdF=((tdF+(-tcalibF[j]+tcalibF[k]))/10);

                tdS = (timeS[j]-timeS[k]);
                tdS = ((tdS+(-tcalibS[j]+tcalibS[k]))/10)+250.0; // Measured in 1e13s which is 100ns

                tdSF = (timeS[j]-timeF[j])/10;
                slowfastTD[j] -> Fill(tdSF);

                if (j == 0) {fastTD0[k]->Fill(tdF);}
                if (j == 1) fastTD1[k]->Fill(tdF);
                if (j == 2) fastTD2[k]->Fill(tdF);
                if (j == 3) fastTD3[k]->Fill(tdF);
                if (j == 4) 
                {
                  fastTD4[k]->Fill(tdF);
                  fastTDslowE[k]->Fill(tdF,energyS[k]);

                  if ((tdF >= 200) && (tdF <= 300)) 
                  {
                    insync[k]->Fill(energyS[k]);
                  }
                }
                if (j == 0) slowTD0[k]->Fill(tdS);
                if (j == 1) slowTD1[k]->Fill(tdS);
                if (j == 2) slowTD2[k]->Fill(tdS);
                if (j == 3) slowTD3[k]->Fill(tdS);
                if (j == 4) slowTD4[k]->Fill(tdS);
              }
            }
            //std::cout << "L" << j << "| energyS[j] " << energyS[j] << "| timeF[j] " << timeF[j] << "| timeS[j] " << timeS[j] << "| timeRF" << timeF[4]<< std::endl;
          
            // ** Uncomment to see what the data looks like **
            // std::string outputDirectory = argv[1];
            // // Remove the last 3 characters from the output directory string
            // outputDirectory = outputDirectory.substr(0, outputDirectory.length() - 3);

            // std::ofstream outputFile(outputDirectory + "/output.txt");
            // if (!outputFile) 
            // {
            //   std::cout << "Error opening the output file." << std::endl;
            //   return 1;
            // }

            // Store the first 10000 values of the variables in the text file
            // int numValues = 1000;
            // for (int i = 0; i < numValues; ++i) 
            // {
            //   for (int j = 0; j <= 4; ++j) 
            //   {
            //     outputFile << j << " "  << energyS[j] << " " << timeF[j] << " " << timeS[j] << std::endl; // To store each value on a new line in the text file.
            //   }
            // }

            // outputFile.close();
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
              double energyS[j], timeF[j], timeS[j];
              tdF=(timeF[j]-timeF[k]); // here
              tdF=((tdF+(-tcalibF[j]+tcalibF[k]))/10)+250.0;
              
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
            double energyS[j], timeF[j], timeS[j];
            if ((m1>0) && (m1<5)) multSlowE[m1]->Fill(energyS[j]);
          }
          // Initialize vectors
          for (int j = 0; j <= 4; ++j) 
          {
            energyS[j].clear();
            timeF[j].clear();
            timeS[j].clear();
            slowECalibL[j].clear();
          }
          slowEPOLARIS.clear();
          slowTimePOLARIS.clear();
            
          //Not so clever
          item=0;
          //printf("\n");
          //std::cout << "TSdiff " << TSdiff << std::endl;
        
        }
  
        skip:

        count++;

        loop:
        pos+=8; // 8 bytes per event
      }
      //std::cout << "End of block " << blocks_in << std::endl;
    }
    end:  
    //std::cout << "End of file reached" << std::endl;
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
    uncalibE[j]->Write();

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

    multSlowE[j]->Write();
    fastTDslowE[j]->Write();

    insync[j]->Write();
    slowfastTD[j] ->Write();
  }
  
  LaBrData->Write();
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

  TCanvas* plot51 = new TCanvas();
  plot51->Divide(1,2);
  plot51->cd(1);
  fastTD4[0]->Draw();
  fastTD4[0]->GetYaxis()->SetTitle( "Counts" );
  fastTD4[0]->SetLineColor(kBlue);
  fastTD4[0]->SetLineWidth(2);
  fastTD4[0]->SetStats(kFALSE);
  fastTD4[0]->SetTitle( "" );
  plot51->cd(2);
  fastTDslowE[0]->Draw("colz");
  fastTDslowE[0]->GetZaxis()->SetRangeUser(0,20);
  fastTDslowE[0]->GetYaxis()->SetTitle( "Energy (keV)" );
  fastTDslowE[0]->GetXaxis()->SetTitle( "Time Difference L0 - RF (ns)" );
  fastTDslowE[0]->SetStats(kFALSE);
  fastTDslowE[0]->SetTitle( "" );
  plot51->Write( "TimeDiffEnergyL0" );
  
  TCanvas* plot52 = new TCanvas();
  plot52->Divide(1,2);
  plot52->cd(1);
  fastTD4[1]->Draw();
  fastTD4[1]->GetYaxis()->SetTitle( "Counts" );
  fastTD4[1]->SetLineColor(kBlue);
  fastTD4[1]->SetLineWidth(2);
  fastTD4[1]->SetStats(kFALSE);
  fastTD4[1]->SetTitle( "" );
  plot52->cd(2);
  fastTDslowE[1]->Draw("colz");
  fastTDslowE[1]->GetZaxis()->SetRangeUser(0,20);
  fastTDslowE[1]->GetYaxis()->SetTitle( "Energy (keV)" );
  fastTDslowE[1]->GetXaxis()->SetTitle( "Time Difference L1 - RF (ns)" );
  fastTDslowE[1]->SetStats(kFALSE);
  fastTDslowE[1]->SetTitle( "" );
  plot52->Write( "TimeDiffEnergyL1" );

  TCanvas* plot53 = new TCanvas();
  plot53->Divide(1,2);
  plot53->cd(1);
  fastTD4[2]->Draw();
  fastTD4[2]->GetYaxis()->SetTitle( "Counts" );
  fastTD4[2]->SetLineColor(kBlue);
  fastTD4[2]->SetLineWidth(2);
  fastTD4[2]->SetTitle( "" );
  fastTD4[2]->SetStats(kFALSE);
  plot53->cd(2);
  fastTDslowE[2]->Draw("colz");
  fastTDslowE[2]->GetZaxis()->SetRangeUser(0,20);
  fastTDslowE[2]->GetYaxis()->SetTitle( "Energy (keV)" );
  fastTDslowE[2]->GetXaxis()->SetTitle( "Time Difference L2 - RF (ns)" );
  fastTDslowE[2]->SetStats(kFALSE);
  fastTDslowE[2]->SetTitle( "" );
  plot53->Write( "TimeDiffEnergyL2" );

  TCanvas* plot54 = new TCanvas();
  plot54->Divide(1,2);
  plot54->cd(1);
  fastTD4[3]->Draw();
  fastTD4[3]->GetYaxis()->SetTitle( "Counts" );
  fastTD4[3]->SetLineColor(kBlue);
  fastTD4[3]->SetLineWidth(2);
  fastTD4[3]->SetStats(kFALSE);
  fastTD4[3]->SetTitle( "" );
  plot54->cd(2);
  fastTDslowE[3]->Draw("colz");
  fastTDslowE[3]->GetZaxis()->SetRangeUser(0,20);
  fastTDslowE[3]->GetYaxis()->SetTitle( "Energy (keV)" );
  fastTDslowE[3]->GetXaxis()->SetTitle( "Time Difference L3 - RF (ns)" );
  fastTDslowE[3]->SetTitle( "" );
  fastTDslowE[3]->SetStats(kFALSE);
  plot54->Write( "TimeDiffEnergyL3" );

  TCanvas* plot71 = new TCanvas();
  plot71->Divide(2,1);
  plot71->cd(1);
  fastTDslowE[0]->Draw("colz");
  fastTDslowE[0]->GetZaxis()->SetRangeUser(0,40);
  fastTDslowE[0]->GetXaxis()->SetTitle( "Time Difference L0 - RF (ns)" );
  fastTDslowE[0]->GetYaxis()->SetTitle( "Energy (keV)" );
  plot71->cd(2);
  slowE[0]->Draw("hbar");
  slowE[0]->GetXaxis()->SetTitle( "Counts" );
  slowE[0]->GetXaxis()->SetLabelSize(0.07);
  slowE[0]->GetXaxis()->SetTitleSize(0.07);
  slowE[0]->SetStats(kFALSE);
  slowE[0]->SetTitle( "");
  plot71->Write( "L0_EnergyTimeDiff_Energy" );

  TCanvas* plot72 = new TCanvas();
  plot72->Divide(2,1);
  plot72->cd(1);
  fastTDslowE[1]->Draw("colz");
  fastTDslowE[1]->GetZaxis()->SetRangeUser(0,40);
  fastTDslowE[1]->GetXaxis()->SetTitle( "Time Difference L1 - RF (ns)" );
  fastTDslowE[1]->GetYaxis()->SetTitle( "Energy (keV)" );
  plot72->cd(2);
  slowE[1]->Draw("hbar");
  slowE[1]->GetXaxis()->SetTitle( "Counts" );
  slowE[1]->GetXaxis()->SetLabelSize(0.07);
  slowE[1]->GetXaxis()->SetTitleSize(0.07);
  slowE[1]->SetStats(kFALSE);
  slowE[1]->SetTitle( "");
  plot72->Write( "L1_EnergyTimeDiff_Energy" );

  TCanvas* plot73 = new TCanvas();
  plot73->Divide(2,1);
  plot73->cd(1);
  fastTDslowE[2]->Draw("colz");
  fastTDslowE[2]->GetZaxis()->SetRangeUser(0,40);
  fastTDslowE[2]->GetXaxis()->SetTitle( "Time Difference L2 - RF (ns)" );
  fastTDslowE[2]->GetYaxis()->SetTitle( "Energy (keV)" );
  plot73->cd(2);
  slowE[2]->Draw("hbar");
  slowE[2]->GetXaxis()->SetTitle( "Counts" );
  slowE[2]->GetXaxis()->SetLabelSize(0.07);
  slowE[2]->GetXaxis()->SetTitleSize(0.07);
  slowE[2]->SetStats(kFALSE);
  slowE[2]->SetTitle( "");
  plot73->Write( "L2_EnergyTimeDiff_Energy" );

  TCanvas* plot74 = new TCanvas();
  plot74->Divide(2,1);
  plot74->cd(1);
  fastTDslowE[3]->Draw("colz");
  fastTDslowE[3]->GetZaxis()->SetRangeUser(0,40);
  fastTDslowE[3]->GetXaxis()->SetTitle( "Time Difference L3 - RF (ns)" );
  fastTDslowE[3]->GetYaxis()->SetTitle( "Energy (keV)" );
  fastTDslowE[3]->GetXaxis()->SetRangeUser(240,260);
  plot74->cd(2);
  slowE[3]->Draw("hbar");
  slowE[3]->GetXaxis()->SetTitle( "Counts" );
  slowE[3]->GetXaxis()->SetLabelSize(0.07);
  slowE[3]->GetXaxis()->SetTitleSize(0.07);
  slowE[3]->SetStats(kFALSE);
  slowE[3]->SetTitle( "");
  plot74->Write( "L3_EnergyTimeDiff_Energy" );
  time(&end);
  
  // Calculating total time taken by the program.
  double time_taken = double((end - start)/60.0);
  std::cout << "Time taken by program is : " << time_taken << std::setprecision(5);
  std::cout << " min " << std::endl;

  delete g;

}
