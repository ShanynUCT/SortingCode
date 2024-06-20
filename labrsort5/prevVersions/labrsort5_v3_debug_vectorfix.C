/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (29/03/2022) - RF corrected timing and Trees  */

/* COMPILE: g++ -std=c++0x labrsort5_v3_debug.C -o exe `root-config --cflags --libs` -lSpectrum */
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
                       if (item < 6)
                        {
                            //std::cout << "ident: " << ident << std::endl;
                            // Identify detector and store data accordingly
                            if ((ident > 63) && (ident < 68)) 
                            {
                                detectorID = ident-64;
                                rndUnif = randy->Rndm()*bw-(bw/2.); 
                                rndchannel = adcdata+rndUnif; 
                                energycalib = (a2[detectorID]*(rndchannel*rndchannel))+(b2[detectorID]*rndchannel)+(c2[detectorID]);
                                slowECalibL[detectorID].push_back(energycalib); 
                                energyS[detectorID].push_back(energycalib);
                            } 
                            else if ((ident > 79) && (ident < 84)) 
                            {
                                detectorID = ident-80;
                                tock = ((adcdata & 0x0000e000) >> 13);
                                tick = (adcdata & 0x00001fff);
                                sick = (tick/8192.0) * 20.0; 
                                if (tock != 7)
                                {
                                    timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
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
                                    timeF[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
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
                        LaBrData->Fill();
                        //printf("Filling Tree \n");
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
                    i=i-2; // go back to the start of the event
                    //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;

                    for (j = 0; j <= 4; j++)
                    { 
                        //Singles spectra 
                        // in order to Fill ROOT histograms we need to convert the vectors to doubles
                        //std::cout << "L" << j << "| energyS[j] " << energyS[j] << "| timeF[j] " << timeF[j] << "| timeS[j] " << timeS[j] << "| timeRF" << timeF[4]<< std::endl;
                        for (const auto& energySj : energyS[j]) 
                        {
                            for (const auto& timeFj : timeF[j]) 
                            {
                                for (const auto& timeSj : timeS[j]) 
                                { 
                                    slowE[j]->Fill(energySj);
                                
                                    hit1->Fill(j);

                                    for (k = 0; k <= 4; k++)
                                    {
                                        for (const auto& energySk : energyS[k]) 
                                        {
                                            for (const auto& timeFk : timeF[k]) 
                                            { 
                                                for (const auto& timeSk : timeS[k]) 
                                                { 
                                                    if (j != k )
                                                    {
                                                        tdF=(timeFj-timeFk); //here
                                                        tdF=((tdF+(-tcalibF[j]+tcalibF[k]))/10);
                                                        tdS = (timeSj-timeSk);
                                                        tdS = ((tdS+(-tcalibS[j]+tcalibS[k]))/10)+250.0; // Measured in 1e13s which is 100ns
                                                        tdSF = (timeSj-timeFj)/10;
                                                        slowfastTD[j] -> Fill(tdSF);
                                                        if (j == 0) {fastTD0[k]->Fill(tdF); slowTD0[k]->Fill(tdS);}
                                                        if (j == 1) {fastTD1[k]->Fill(tdF); slowTD1[k]->Fill(tdS);}
                                                        if (j == 2) {fastTD2[k]->Fill(tdF); slowTD2[k]->Fill(tdS);}
                                                        if (j == 3) {fastTD3[k]->Fill(tdF); slowTD3[k]->Fill(tdS);}
                                                        if (j == 4) 
                                                        {
                                                            fastTD4[k]->Fill(tdF);
                                                            slowTD4[k]->Fill(tdS);
                                                            fastTDslowE[k]->Fill(tdF,energySk);

                                                            if ((tdF >= 200) && (tdF <= 300)) 
                                                            {
                                                                insync[k]->Fill(energySk);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                   item = 0;
                   // Reset vectors
                    for (j = 0; j <=4; j++)
                    {
                    energyS[j].clear();
                    timeF[j].clear();
                    timeS[j].clear();
                    slowECalibL[j].clear();
                    }
                    slowEPOLARIS.clear();
                    slowTimePOLARIS.clear();
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