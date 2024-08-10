/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (01/06/2023) */

/* This code reads in the list mode data from the buffer. It sorts the time and (calibrated) energy for each detector into a tree. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort1labr_timingshapingjune2024.C -o exe1timingshaping `root-config --cflags --libs` -lSpectrum 
     ./exe1timingshaping ~/Users/shanyn/Documents/PhD/Exp/2024/timing_shaping_characterisation_june2024/CC_polaris_2inch/run8/R8
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
  int RunEnd=999999, Run=0;
  
  int m1,m2;
  int sum1, sum2;

  time_t start, end;
  int fastcount[4]={0,0,0,0}, slowcount[4]={0,0,0,0}, fastslowcount[4]={0,0,0,0};
  time(&start);

  // Declare vectors to store the event data
  std::vector<double> energyS[5];
  std::vector<double> energycalibS[5];
  std::vector<double> energyF[5];
  std::vector<double> timeF[5];
  std::vector<double> timeS[5];
  std::vector<double> fastEPOLARIS(1, 0.0); //Polaris has only one detector
  std::vector<double> fastTPOLARIS(1, 0.0);
  std::vector<int> pop(1, 0);
  double energySlow[5]={0,0,0,0,0}, timeFast[5]={0,0,0,0,0}, timeSlow[5]={0,0,0,0,0},  fastEnergyPOLARIS[1]={0}, fastTimePOLARIS[1]={0}, energyFast[5]={0,0,0,0,0}, energy[5]={0,0,0,0,0};
  double tdF, tdSF0, tdSF1;
  int Ed;

    //Energy calibrations 2nd order polynomials (APPLIED) # timing_shaping_characterisation_june_2024 Data (L0 has 152Eu and 137Cs and 60Co peaks 2nd order, L2 has 152Eu 2nd order)
    std::vector<double> a2slow = {3.82160000E-07, 5.0248710E-07};
    std::vector<double> b2slow = {2.80362000E-01, 2.7453690E-01};
    std::vector<double> c2slow = {3.73025000E-01, 9.07E-01};

    std::vector<double> a2fast = {1.56438000E-06, 2.19003000E-07};
    std::vector<double> b2fast = {2.69853000E-01, 1.90136900E-01};
    std::vector<double> c2fast = {1.09687000E+01, 0.00000000E+00};

    //Energy calibrations 2nd order polynomials (APPLIED) # timing_shaping_characterisation_june_2024 Data (L0 has 22Na peaks linear fit, L2 has 152Eu 2nd order)
    /* std::vector<double> a2slow = {0, 4.1046890E-07};
    std::vector<double> b2slow = {2.76500000E-01, 2.7537890E-01};
    std::vector<double> c2slow = {0.00000000E+00, 0.00000000E+00};

    std::vector<double> a2fast = {0, 2.19003000E-07};
    std::vector<double> b2fast = {1.85123000E-01, 1.90136900E-01};
    std::vector<double> c2fast = {0.00000000E+00, 0.00000000E+00}; */

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
  LaBrData->Branch("fastECalibL0", &energyFast[0], "fastECalibL0/D");
  LaBrData->Branch("timeSL0", &timeSlow[0], "timeSL0/D");
  LaBrData->Branch("timeFL0", &timeFast[0], "timeFL0/D");
  LaBrData->Branch("slowECalibL1", &energySlow[1], "slowECalibL1/D");
  LaBrData->Branch("fastECalibL1", &energyFast[1], "fastECalibL1/D");
  LaBrData->Branch("timeSL1", &timeSlow[1], "timeSL1/D");
  LaBrData->Branch("timeFL1", &timeFast[1], "timeFL1/D");
  LaBrData->Branch("slowECalibL2", &energySlow[2], "slowECalibL2/D");
  LaBrData->Branch("fastECalibL2", &energyFast[2], "fastECalibL2/D");
  LaBrData->Branch("timeSL2", &timeSlow[2], "timeSL2/D");
  LaBrData->Branch("timeFL2", &timeFast[2], "timeFL2/D");
  LaBrData->Branch("slowECalibL3", &energySlow[3], "slowECalibL3/D");
  LaBrData->Branch("fastECalibL3", &energyFast[3], "fastECalibL3/D");
  LaBrData->Branch("timeSL3", &timeSlow[3], "timeSL3/D");
  LaBrData->Branch("timeFL3", &timeFast[3], "timeFL3/D");
  LaBrData->Branch("fastEPOLARIS", &fastEnergyPOLARIS[0], "fastEPOLARIS/D");
  LaBrData->Branch("fastTPOLARIS", &fastTimePOLARIS[0], "fastTPOLARIS/D");

  // ################ Create histograms ################
  TH1D** slowcalibE=new TH1D*[5];
  TH1D** slowE=new TH1D*[5];
  TH1D** fastE=new TH1D*[5];
  TH1D** fastTD0=new TH1D*[5];
  TH1D** slowfastTD=new TH1D*[5];
  TH2D** slowfastTDSE=new TH2D*[5];
  TH2D** slowfastTDFE=new TH2D*[5];

  for (j = 0; j <=4; j++)
  {
    slowcalibE[j] = new TH1D(TString::Format("Slow_Energy_Calib_L%2d", j),"Spectrum",8000,0,8000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%2d", j),"Spectrum",16380,0,16380);
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%2d", j),"Spectrum",8000,0,8000);
    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",2000,-1000,1000);
    slowfastTD[j] = new TH1D(TString::Format("Slow-Fast Time Difference_L0-L%2d", j),"Spectrum",2000,-1000,1000);
    slowfastTDSE[j] =new TH2D(TString::Format("Slow-Fast Time Difference vs Slow Energy_-L%2d", j),"Spectrum",2000,-1000,1000,8000,0,8000);
    slowfastTDFE[j] =new TH2D(TString::Format("Slow-Fast Time Difference vs Fast Energy_-L%2d", j),"Spectrum",2000,-1000,1000,8000,0,8000);
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
                    TS = (TS & 0x0000fffff0000000ULL); // 100 MHz sampling speed
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
                    Polaris |    95  |   79   |   //    |   //    | structure: 2s intervals    |   
                    ---------------------------------------------------------------------
                    */

                    // Identify detector and store data accordingly
                    if ((ident > 63) && (ident < 66)) 
                    {
                        detectorID = ident-64;
                        energycalib = 0;
                        rndUnif = randy->Uniform(-3.5, 3.5);
                        energycalib = (a2slow[detectorID]*(adcdata*adcdata))+(b2slow[detectorID]*adcdata)+(c2slow[detectorID]);
                        if (energycalibS[detectorID].size() == 0) energycalibS[detectorID].push_back(energycalib + rndUnif);
                        if (energyS[detectorID].size() == 0) energyS[detectorID].push_back(adcdata);
                    }
                    else if ((ident >79) && (ident < 82)) 
                    {
                        detectorID = ident-80;
                        int tick, tock;
                        double sick;
                        tock = ((adcdata & 0x0000e000) >> 13);
                        tick = (adcdata & 0x00001fff);
                        sick = (tick/8192.0) * 20.0; 
                        // if (tock != 7) // The time slow doesnt have the tock information -> there is no fractonal range in the byte information between 0 and 7 for the 2ns sampling speed interval. Uncommenting this gives very low statistics.
                        {
                            if (timeS[detectorID].size() == 0)  timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick); // time in 1e12s which is 10ns
                        }
                    }
                    else if ((ident > 71) && (ident < 74)) 
                    {
                        detectorID = ident-72;
                        energycalib = 0;
                        rndUnif = randy->Uniform(-5.2, 5.2);
                        energycalib = (a2fast[detectorID]*(adcdata*adcdata))+(b2fast[detectorID]*adcdata)+(c2fast[detectorID]);
                        if (energyF[detectorID].size() == 0) energyF[detectorID].push_back(energycalib+rndUnif);
                        //if (energyF[detectorID].size() == 0) energyF[detectorID].push_back(adcdata);
                    } 
                    else if ((ident > 87) && (ident < 90)) 
                    {
                        detectorID = ident-88;
                        int tick, tock;
                        double sick;
                        tock = ((adcdata & 0x0000e000) >> 13);// 500 MHz sampling speed (2 ns) bits 13-15
                        tick = (adcdata & 0x00001fff);// 500 MHz sampling speed. Time difference between 2ns samples bits 0-12. 0 is start of tick, 8191 is start of next tick (2,4,6,8)
                        sick = (tick/8192.0) * 20.0; // between the 2 ns samples
                        if (tock != 7)
                        {
                            if (timeF[detectorID].size() == 0)
                            {
                                timeF[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick);
                            }
                        }
                    }
                    else if (ident == 95) 
                    {
                        detectorID = ident-95;
                        int tick, tock;
                        double sick;
                        tock = ((adcdata & 0x0000e000) >> 13);
                        tick = (adcdata & 0x00001fff);
                        sick = (tick/8192.0) * 20.0; 
                        if (fastEPOLARIS.size() == 0)
                        {
                            fastEPOLARIS[detectorID] = adcdata;
                            fastEPOLARIS.push_back(fastEPOLARIS[0]);
                        }
                        if (fastTPOLARIS.size() == 0)
                        {
                            fastTPOLARIS[detectorID] = ((TS * 100.0) + (tock * 20.0) + sick);
                            fastTPOLARIS.push_back(fastTPOLARIS[0]);
                        }
                    }
                    

                    // printf("TSdiff %lld \n", TSdiff);
                    if (TSdiff >= 90 )                     
                    {   
                        i=i-2; 
                        //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;
                        for (j = 0; j <=4; j++)
                        {
                            //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                            if (energyS[j].size() == 0) energyS[j].push_back(0);
                            if (energycalibS[j].size() == 0) energycalibS[j].push_back(0);
                            if (timeF[j].size() == 0) timeF[j].push_back(0);
                            if (timeS[j].size() == 0)  timeS[j].push_back(0);
                            if (energyF[j].size() == 0) energyF[j].push_back(0);

                            timeFast[j] = timeF[j][0]/10;
                            timeSlow[j] = timeS[j][0]/10;
                            energySlow[j] = energycalibS[j][0];
                            energy[j] = energyS[j][0];
                            energyFast[j] = energyF[j][0];

                            // printf("j %d | energySlow %f | timeFast %f | timeSlow %f | energyFast %f \n", j, energySlow[j], timeFast[j], timeSlow[j], energyFast[j]);
                        }
                        fastEnergyPOLARIS[0] = fastEPOLARIS[0];
                        fastTimePOLARIS[0] = fastTPOLARIS[0]/10;
                       // printf("timeF %f | timeS %f | energyS %f\n", timeFast[0]  , timeSlow[0]  , energySlow[0]);
                        LaBrData->Fill();

                        for (j = 0; j <= 3; j++) 
                        { 
                            //Singles spectra 
                            //std::cout << "L" << j << "| energySlow[j] " << energySlow[j] << "| timeFast[j] " << timeFast[j] << "| timeSlow[j] " << timeSlow[j] << "| timeRF" << timeFast[4]<< std::endl;
                            slowcalibE[j]->Fill(energySlow[j]);
                            fastE[j]->Fill(energyFast[j]);
                            slowE[j]->Fill(energy[j]);
                        }

                        for (j = 0; j <= 0; j++) 
                        { 
                            tdSF0 = (timeSlow[j]-timeFast[j]); // Measured in nanoseconds
                            if (tdSF0 == 0) continue;
                            slowfastTD[j]->Fill(tdSF0);
                            slowfastTDFE[j]->Fill(tdSF0, energyFast[j]);
                            slowfastTDSE[j]->Fill(tdSF0, energySlow[j]);

                            for (k = 1; k <= 1; k++) 
                            {
                                tdSF1 = (timeSlow[k]-timeFast[k]); // Measured in nanoseconds
                                if (tdSF1 == 0) continue;
                                slowfastTD[k]->Fill(tdSF1);
                                slowfastTDFE[k]->Fill(tdSF1, energyFast[k]);
                                slowfastTDSE[k]->Fill(tdSF1, energySlow[k]);
                                
                                if (j != k ) 
                                {
                                    if ((timeFast[j] == 0) &&(timeFast[k]== 0)) break;
                                    tdF=(timeFast[j]-timeFast[k]); // Measured in nanoseconds
                                    if (j == 0) 
                                    {
                                        fastTD0[k]->Fill(tdF); 
                                        //printf("tdF %f\n", tdF);
                                    }
                                }
                            }                        
                        }

                        // Reset vectors
                        for (j = 0; j <=4; j++) 
                        {
                            energyS[j].clear();
                            energycalibS[j].clear();
                            energyF[j].clear();
                            timeF[j].clear();
                            timeS[j].clear();
                            energySlow[j] = 0;
                            timeFast[j] = 0;
                            timeSlow[j] = 0;
                            energyFast[j] = 0;
                        }
                        fastEPOLARIS.clear();
                        fastTPOLARIS.clear();                 
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

    // calculate the time difference between the first and last event in LaBrData
    LaBrData->SetBranchAddress("timeFL0", &timeFast[0]);
    int nEnt = LaBrData->GetEntries();
    double firstTime = 0;
    // find first non-zero entry
    for (int entryIndex = 0; entryIndex < nEnt; ++entryIndex) 
    {
        LaBrData->GetEntry(entryIndex);
        if ((timeFast[0] >0)) 
        {
            firstTime = timeFast[0];
            break;
        }
    }

    // find last non-zero entry
    double lastTime = 0;
    for (int entryIndex = nEnt-1; entryIndex >= 0; --entryIndex) 
    {
        LaBrData->GetEntry(entryIndex);
        if ((timeFast[0] >0)) 
        {
            lastTime = timeFast[0];
            break;
        }
    }

    printf("\nFirst time: %.2f \n", firstTime);
    printf("\nLast time: %.2f \n", lastTime);

    printf("\nLaBrData TTree run duration: %.2f min (JUST A CHECK)\n", (lastTime - firstTime)*1e-9/ 60.0);
    
    LaBrData->Write();

    for (j = 0; j <=4; j++)
    {
        slowE[j]->Write();
        slowcalibE[j]->Write();
        fastE[j]->Write();
        fastTD0[j]->Write();
        slowfastTD[j]->Write();
        slowfastTDSE[j]->Write();
        slowfastTDFE[j]->Write();
    }
    
    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("\nThe first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference of the RXX_ file is (min): %f\n", TimeDiff);

    time(&end);
    double time_taken = double((end - start)/60.0);
    std::cout << "\nTime taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " min " << std::endl;

    
    delete g;
}