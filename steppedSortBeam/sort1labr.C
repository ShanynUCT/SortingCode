/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (01/06/2023) */

/* This code reads in the list mode data from the buffer. It sorts the time and (calibrated) energy for each detector into a tree. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort1labr.C -o exe1 `root-config --cflags --libs` -lSpectrum 
     ./exe1 ~/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63
*/
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

  // Declare vectors to store the event data
  std::vector<double> timeF[5];
  std::vector<double> timeS[5];
  std::vector<double> slowECalibL[5];
  std::vector<double> fastEPOLARIS(1, 0.0); //Polaris has only one detector
  std::vector<double> fastTPOLARIS(1, 0.0);
  std::vector<int> pop(1, 0);
  double energySlow[5]={0,0,0,0,0}, timeFast[5]={0,0,0,0,0}, timeSlow[5]={0,0,0,0,0}, slowEnergyCalib[5]={0,0,0,0,0}, fastEnergyPOLARIS[1]={0}, fastTimePOLARIS[1]={0}, insyncEnergySlow[5]={0,0,0,0,0};
  double tdF, tdS, tdSF;
  int Ed;

  //Energy calibrations 2nd order polynomials (APPLIED) # 15062022 Data
  std::vector<double> a2 = {0.00000288786796, 0.00000351711116, 0.00000160595876, 0.00000300000230};
  std::vector<double> b2 = {0.747696355, 0.805816609, 0.513776183, 0.709359122};
  std::vector<double> c2 = {1.64893477, 1.09470037, 1.37154046, 2.18122723};

  //Energy calibrations 2nd order polynomials (APPLIED) # 28022023 Data
    // std::vector<double> a2 = {0.00001942, 0.00010155, 0.00003208, 0.00002200};
    // std::vector<double> b2 = {1.81799962, 2.61697133, 2.21345098, 2.05179483};
    // std::vector<double> c2 = {-1.65129151, -0.54900373, -4.86476637, -0.97991876};

  double bw = 1, energycalib, rndchannel, rndUnif;
  TRandom3* randy = new TRandom3();

    
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
  LaBrData->Branch("slowECalibL0", &energySlow[0], "slowECalibL0/D");
  LaBrData->Branch("timeSL0", &timeSlow[0], "timeSL0/D");
  LaBrData->Branch("timeFL0", &timeFast[0], "timeFL0/D");
  LaBrData->Branch("slowECalibL1", &energySlow[1], "slowECalibL1/D");
  LaBrData->Branch("timeSL1", &timeSlow[1], "timeSL1/D");
  LaBrData->Branch("timeFL1", &timeFast[1], "timeFL1/D");
  LaBrData->Branch("slowECalibL2", &energySlow[2], "slowECalibL2/D");
  LaBrData->Branch("timeSL2", &timeSlow[2], "timeSL2/D");
  LaBrData->Branch("timeFL2", &timeFast[2], "timeFL2/D");
  LaBrData->Branch("slowECalibL3", &energySlow[3], "slowECalibL3/D");
  LaBrData->Branch("timeSL3", &timeSlow[3], "timeSL3/D");
  LaBrData->Branch("timeFL3", &timeFast[3], "timeFL3/D");
  LaBrData->Branch("timeRF", &timeFast[4], "timeRF/D");
  LaBrData->Branch("fastEPOLARIS", &fastEnergyPOLARIS[0], "fastEPOLARIS/D");
  LaBrData->Branch("fastTPOLARIS", &fastTimePOLARIS[0], "fastTPOLARIS/D");

  // ################ Create histograms ################
  TH1D** slowE=new TH1D*[5];
  TH1D** fastTD4=new TH1D*[5];
  TH1D** fastTD0=new TH1D*[5];
  TH1F** eventCounts = new TH1F*[5];

  for (j = 0; j <=4; j++)
  {
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%02d", j),"Spectrum",8000,0,8000);
    fastTD4[j] = new TH1D(TString::Format("Fast Time Difference_RF-L%2d", j),"Spectrum",2000,-1000,1000);
    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",2000,-1000,1000);
    eventCounts[j] = new TH1F(TString::Format("Event Counts L%2d", j), "Event Counts", 3, 0.5, 3.5);
  }

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

                    // Identify detector and store data accordingly
                    if ((ident > 63) && (ident < 68)) 
                    {
                        detectorID = ident-64;
                        rndUnif = randy->Uniform(-2, 2); 
                        rndchannel = adcdata; 
                        energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);

                        tock = ((adcdata & 0x0000e000) >> 13);
                        tick = (adcdata & 0x00001fff);
                        sick = (tick/8192.0) * 20.0; 
                        if (tock != 7)
                        {
                            if ((timeS[detectorID].size() == 0) && (slowECalibL[detectorID].size() == 0))
                            {
                                timeS[detectorID].push_back(((TS * 100.0) + (tock * 20.0) + sick)); // time in 1e12s which is 10ns
                                slowECalibL[detectorID].push_back(energycalib+rndUnif);
                                //printf("timeS[%d] = %f\n", detectorID, timeS[detectorID][0]);
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
                            if (timeF[detectorID].size() == 0)
                            {
                                timeF[detectorID].push_back(((TS * 100.0) + (tock * 20.0) + sick)); // time in 1e12s which is 10ns
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
                        if ((fastEPOLARIS.size() == 0) && (fastTPOLARIS.size() == 0))
                        {
                            fastEPOLARIS.push_back(adcdata);
                            fastTPOLARIS.push_back((TS * 100.0) + (tock * 20.0) + sick);
                        }
                    }


                    if (TSdiff >= 1 )                     
                    {   
                        i=i-2; 
                        //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;
                        for (j = 0; j <=4; j++)
                        {
                            //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                            //if the size of avector is 0, fill the value with 0
                            if (slowECalibL[j].size() == 0) slowECalibL[j].push_back(0);
                            if (timeF[j].size() == 0) timeF[j].push_back(0);
                            if (timeS[j].size() == 0)  timeS[j].push_back(0);

                            timeFast[j] = timeF[j][0]/10;
                            timeSlow[j] = timeS[j][0]/10;
                            energySlow[j] = slowECalibL[j][0];
                        }
                        fastEnergyPOLARIS[0] = fastEPOLARIS[0];
                        fastTimePOLARIS[0] = fastTPOLARIS[0]/10;


                        LaBrData->Fill();

                        for (j = 0; j <4; j++)
                        {
                            if (energySlow[j] >= 1000)
                            {
                                if (timeFast[j] != 0 && timeSlow[j] == 0) fastcount[j]++;
                                else if (timeSlow[j] != 0 && timeFast[j] == 0) slowcount[j]++;
                                else if (timeSlow[j] != 0 && timeFast[j] != 0) fastslowcount[j]++;

                                eventCounts[j]->Fill(1,fastcount[j]);
                                eventCounts[j]->Fill(2,slowcount[j]);
                                eventCounts[j]->Fill(3,fastslowcount[j]);

                                eventCounts[j]->GetXaxis()->SetBinLabel(1,"Only Fast");
                                eventCounts[j]->GetXaxis()->SetBinLabel(2,"Only Slow");
                                eventCounts[j]->GetXaxis()->SetBinLabel(3,"Coinc Fast+Slow");
                            }                        
                        }

                        for (j = 0; j <= 3; j++) 
                        { 
                            //Singles spectra 
                            //std::cout << "L" << j << "| energySlow[j] " << energySlow[j] << "| timeFast[j] " << timeFast[j] << "| timeSlow[j] " << timeSlow[j] << "| timeRF" << timeFast[4]<< std::endl;
                            slowE[j]->Fill(energySlow[j]);
                        }

                        for (j = 0; j <= 4; j++) 
                        { 
                            for (k = 0; k <= 4; k++) 
                            {
                                if (j != k ) 
                                {
                                    tdF=(timeFast[j]-timeFast[k]); // Measured in 1e13s which is 100ns
                                    if (j == 4) fastTD4[k]->Fill(tdF);
                                    if (j == 0) fastTD0[k]->Fill(tdF);

                                }
                            }                        
                        }


                        // Reset vectors
                        for (j = 0; j <=4; j++) 
                        {
                            timeF[j].clear();
                            timeS[j].clear();
                            slowECalibL[j].clear();
                        }
                        fastEPOLARIS.clear();
                        fastTPOLARIS.clear();

                        //Reinitialise arrays
                        for (j = 0; j <4; j++)
                        { 
                            energySlow[j] = 0;
                            timeFast[j] = 0;
                            timeSlow[j] = 0;
                        }                   
                        fastEnergyPOLARIS[0] = 0;
                        fastTimePOLARIS[0] = 0; 
                        
                    }

                    TSinit=TS;

                }

                /* SYNC */
                if ((data & 0xc0f00000) == 0x80400000) 
                {
                    card = (data & 0x3f000000) >> 24;
                    TStop = (data & 0x000fffff);
                    TS = (TStop);
                    TS = ((TS << 28));
                    TS = ((TS | TSbot));

                    TSdiff = TS-TSinit;

                    SYNC = TS;
                    SYNCdiff = SYNC-SYNClast;

                    twidset=0;

                    SYNClast = SYNC;

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
        fastTD4[j]->Write();
        fastTD0[j]->Write();
        eventCounts[j]->Write();
    }
    
    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference is (min): %f\n", TimeDiff);

    time(&end);
    double time_taken = double((end - start)/60.0);
    std::cout << "Time taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " min " << std::endl;
    
    delete g;
}