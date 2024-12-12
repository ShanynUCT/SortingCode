/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (01/06/2023) */

/* This code reads in the list mode data from the buffer. It sorts the time and (calibrated) energy for each detector into a tree. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort1labr_13082024.C -o exe13082024 `root-config --cflags --libs` -lSpectrum 
     ./exe13082024 /Users/shanyn/Documents/PhD/PANGoLINS/Measurements/DetectorAssemblies/LaBr3/sn230824-05/13082024/R
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
  time(&start);

  // Declare vectors to store the event data
  std::vector<double> energyS[1];
  std::vector<double> energycalibS[1];
  std::vector<double> timeS[1];
  double energySlow[1]={0}, timeSlow[1]={0}, energy[1]={0};
  int Ed;

    //LaBr3Ce sn230824-05 pol2
    /* double a2slow = 0;	
    double b2slow = 4.06270E-06;
    double c2slow = 6.51000E-02; */

    //SrI2Eu sn230912-04  pol3
    double a2slow = 0;
    double b2slow = 1.68760E-06;
    double c2slow = 1.53180E-01;

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

  // ################ Create histograms ################
  TH1D** slowcalibE=new TH1D*[5];
  TH1D** slowE=new TH1D*[5];

  for (j = 0; j <1; j++)
  {
    slowcalibE[j] = new TH1D(TString::Format("Slow_Energy_Calib_L%2d", j),"Spectrum",8000,0,8000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%2d", j),"Spectrum",16380,0,16380);
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
                    LaBr 0  |   -   |   -   |   80    |    64   | -1050   | 
                    ---------------------------------------------------------
                    */

                    // Identify detector and store data accordingly
                    if ((ident > 63) && (ident < 65)) 
                    {
                        detectorID = ident-64;
                        energycalib = 0;
                        rndUnif = randy->Uniform(-4,4);
                        energycalib = a2slow*(adcdata*adcdata*adcdata)+b2slow*(adcdata*adcdata)+c2slow*adcdata;
                        if (energycalibS[detectorID].size() == 0) energycalibS[detectorID].push_back(energycalib + rndUnif);
                        if (energyS[detectorID].size() == 0) energyS[detectorID].push_back(adcdata);
                    }
                    else if ((ident >79) && (ident < 81)) 
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
                    

                    // printf("TSdiff %lld \n", TSdiff);
                    if (TSdiff >= 90 )                     
                    {   
                        i=i-2; 
                        //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;
                        for (j = 0; j <1; j++)
                        {
                            //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                            if (energyS[j].size() == 0) energyS[j].push_back(0);
                            if (energycalibS[j].size() == 0) energycalibS[j].push_back(0);
                            if (timeS[j].size() == 0)  timeS[j].push_back(0);

                            timeSlow[j] = timeS[j][0]/10;
                            energySlow[j] = energycalibS[j][0];
                            energy[j] = energyS[j][0];
                            // printf("j %d | energySlow %f | timeFast %f | timeSlow %f | energyFast %f \n", j, energySlow[j], timeFast[j], timeSlow[j], energyFast[j]);
                        }
                       // printf("timeF %f | timeS %f | energyS %f\n", timeFast[0]  , timeSlow[0]  , energySlow[0]);
                        LaBrData->Fill();

                        for (j = 0; j <1; j++) 
                        { 
                            //Singles spectra 
                            slowcalibE[j]->Fill(energySlow[j]);
                            slowE[j]->Fill(energy[j]);
                        }

                        
                        // Reset vectors
                        for (j = 0; j <1; j++) 
                        {
                            energyS[j].clear();
                            energycalibS[j].clear();
                            timeS[j].clear();
                            energySlow[j] = 0;
                            timeSlow[j] = 0;
                        }                      
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

    for (j = 0; j <1; j++)
    {
        slowE[j]->Write();
        slowcalibE[j]->Write();
    } 

    std::string inputfile = argv[1];
    std::string delimiter = "/";
    std::string directory = inputfile.substr(0, inputfile.find_last_of(delimiter));
    std::string filename = inputfile.substr(inputfile.find_last_of(delimiter) + 1);

    std::ofstream uncalibrated;
    uncalibrated.open(directory + "/"+filename+"_hist.asc");
    for (int i = 0; i < 17000; i++)
    {
        uncalibrated << i << "," << slowE[0]->GetBinContent(i) << std::endl;
    }
    uncalibrated.close();

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