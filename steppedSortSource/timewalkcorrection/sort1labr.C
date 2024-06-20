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
     ./exe1 ../../../../Exp/2024/timing_shaping_characterisation_june2024/R2  
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
#include "TH2D.h"


#define SIZE  16384
int32_t buffer[SIZE];
uint64_t TStop, TSbot, data;
uint64_t TS, TSinit=0, TSfirst, counter,  SYNC, SYNClast=0, TSfirstpolaris = -1, TSlastpolaris;
int64_t TSdiff, SYNCdiff;
uint64_t TSdiffpolaris;
int64_t count=0, countpolaris = 0;
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
  std::vector<double> energycalibS[2];
  std::vector<double> energycalibF[2];
  std::vector<double> timeF[2];
  std::vector<double> timeS[2];
  std::vector<double> timeP[1];
  std::vector<double> energyP[1];
  std::vector<double> energyuncalibS[2];
  std::vector<double> energyuncalibF[2];

  double energySlow[2]={0,0}, timeFast[2]={0,0}, timeSlow[2]={0,0}, energyFast[2]={0,0}, energySlowUncalib[2]={0,0}, energyFastUncalib[2]={0,0}, energyPolaris[1]={0}, timePolaris[1]={0};

  double tdF, tdSF0, tdSF1;
  int Ed;
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
  LaBrData->Branch("timeP", &timePolaris[0], "timeP/D");
  LaBrData->Branch("energyP", &energyPolaris[0], "energyP/D");

  // ################ Create histograms ################
  TH1D** slowcalibE=new TH1D*[2];
  TH1D** fastcalibE=new TH1D*[2];
  TH1D** slowE=new TH1D*[2];
  TH1D** fastE=new TH1D*[2];
  TH1D** fastTD0=new TH1D*[2];
  TH1D** slowfastTD=new TH1D*[2];
  TH2D** slowfastTDSE=new TH2D*[2];
  TH2D** slowfastTDFE=new TH2D*[2];
  TH2D** slowfastE=new TH2D*[2];

  for (j = 0; j <2; j++)
  {
    slowcalibE[j] = new TH1D(TString::Format("Slow_Energy_Calib_L%2d", j),"Spectrum", 2000,0,2000);
    slowE[j] = new TH1D(TString::Format("Slow_Energy_Uncalib_L%2d", j),"Spectrum", 16384,0,16384);
    fastcalibE[j] = new TH1D(TString::Format("Fast_Energy_Calib_L%2d", j),"Spectrum", 2000,0,2000);
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%2d", j),"Spectrum", 16384,0,16384);
    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum", 100,-5,5);
    slowfastTD[j] = new TH1D(TString::Format("Slow-Fast Time Difference_L%2d", j),"Spectrum", 1000,0,1000);
    slowfastTDSE[j] =new TH2D(TString::Format("Slow-Fast Time Difference vs Slow Energy_-L%2d", j),"Spectrum", 1000,0,1000, 2000,0,2000);
    slowfastTDFE[j] =new TH2D(TString::Format("Slow-Fast Time Difference vs Fast Energy_-L%2d", j),"Spectrum", 1000,0,1000, 2000,0,2000);
    slowfastE[j] =new TH2D(TString::Format("Slow-Fast Energy_-L%2d", j),"Spectrum", 2000,0,2000, 2000,0,2000);
  }

    TH2D** polaris_energytime = new TH2D*[1];
    polaris_energytime[0] = new TH2D("Polaris_Energy_Time","Spectrum", 2000,0,2000, 2000,0,2000);

  //Time calibrations
  double Ta[2];
  //measurements in SHANY9
  //R1
  /* Ta[0]=   0.00;
  Ta[1]=   -11.39+1.8; */
  //R12
  Ta[0]=   0.00;
  Ta[1]=   -11.39+1.8;

  //Energy cablibrations
  double as[2], bs[2], cs[2], ds[2], af[2], bf[2], cf[2], df[2];
    /* // LaBr slow 0 (Labelled 8)
    as[0] =  0.000002190000;
    bs[0] =  0.433381500000;
    cs[0] =  0.000000000000;

        // LaBr slow 1  (Labelled 2)
    as[1] =  0.000001726752;
    bs[1] =  0.430148300000;
    cs[1] =  0.000000000000;

    // LaBr fast 0 (Labelled 9)
    af[0] =  0.000000105753;
    bf[0] =  0.180526000000;
    cf[0] =  0.000000000000;

    // LaBr fast 1 (Labelled 3)
    af[1] =  0.000000224606;
    bf[1] =  0.186934400000;
    cf[1] =  0.000000000000; */


    // LaBr slow 0 (Labelled 8)  152Eu calibration SHANY9/R12 factor ch/e = 3.5
    as[0] =  7.34444200E-07;
    bs[0] =  2.80468500E-01;
    cs[0] =  0.;

    // LaBr fast 0 (Labelled 8) 152Eu calibration SHANY9/R12 factor ch/e = 5.4
    af[0] =  3.47720100E-07;
    bf[0] =  1.81586700E-01;
    cf[0] =  0.;

    // LaBr slow 1  (Labelled 2) 152Eu calibration SHANY9/R7 factor ch/e = 3.6
    as[1] =  4.1046890E-07; 
    bs[1] =  2.7537890E-01;
    cs[1] =  0.000000000000;

    // LaBr fast 1 (Labelled 3) 152Eu calibration SHANY9/R7 factor ch/e = 5.2
    af[1] =  2.19003000E-07;
    bf[1] =  1.90136900E-01;
    cf[1] =  0.;
        

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
                    ----------------------------------------------------------------------
                    Signal               |  fastT | fastE  | slowT   |  slowE  |    HV   |    
                    ----------------------------------------------------------------------
                    LaBr 0 (Labelled 8)  |    88  |   72   |   80    |    64   | -1200   | 
                    ----------------------------------------------------------------------
                    LaBr 1  (Labelled 2) |    89  |   73   |    81   |     65  | -1200   |   
                    --------------------------------------------------------------------------------------
                    Polaris              |    95  |   79   |   //    |   //    | structure: 2s intervals |   
                    --------------------------------------------------------------------------------------

                    */
                    //if (ident == 64) printf("adcdata %d\n", adcdata);
                    // DETECTOR 0 slow signal is split between channels 0 and 1 but I'm just extracting 0 (64) 
                    if (ident > 63 && ident < 66)
                    {
                        detectorID = ident-64;
                        energycalib = 0;
                        rndUnif = 0;
                        rndUnif = randy->Uniform(-3, 3);
                        rndchannel = adcdata;
                        energycalib = (as[detectorID]*(rndchannel*rndchannel))+(bs[detectorID]*rndchannel)+(cs[detectorID]);
                        if (energycalibS[detectorID].size() == 0) energycalibS[detectorID].push_back(energycalib + rndUnif);
                        if (energyuncalibS[detectorID].size() == 0) energyuncalibS[detectorID].push_back(adcdata);


                        //if (ident == 64) printf("adcdata %d\n", adcdata);
                    }
                    else if (ident > 79 && ident < 82)
                    {
                        detectorID = ident-80;
                        if (timeS[detectorID].size() == 0) 
                        {
                            timeS[detectorID].push_back(TS*100);
                            //printf("timeS[%d] = %f\n", detectorID, timeS[detectorID][0]);
                        }
                    }
                    else if (ident > 71 && ident < 74)
                    {
                        detectorID = ident-72;
                        rndUnif = 0;
                        energycalib = 0;
                        rndUnif = randy->Uniform(-5, 5);
                        rndchannel = adcdata;
                        energycalib = (af[detectorID]*(rndchannel*rndchannel))+(bf[detectorID]*rndchannel)+(cf[detectorID]);
                        if (energycalibF[detectorID].size() == 0) energycalibF[detectorID].push_back(energycalib + rndUnif);
                        if (energyuncalibF[detectorID].size() == 0) energyuncalibF[detectorID].push_back(adcdata);
                    } 
                    else if (ident > 87 && ident < 90)
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
                    else if (ident == 79)
                    {
                        detectorID = ident-79;
                        if (timeP[detectorID].size() == 0) timeP[detectorID].push_back(TS*100);
                        if (energyP[detectorID].size() == 0) energyP[detectorID].push_back(adcdata);
                    }
                    else
                    {
                        goto skip;
                    }

                    // printf("TSdiff %lld \n", TSdiff);
                    if (TSdiff > 70 )                     
                    {   
                        i=i-2; 
                        //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;
                        for (j = 0; j <2; j++)
                        {
                            //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                            if (energycalibS[j].size() == 0) energycalibS[j].push_back(0);
                            if (timeF[j].size() == 0) timeF[j].push_back(0);
                            if (timeS[j].size() == 0)  timeS[j].push_back(0);
                            if (energycalibF[j].size() == 0) energycalibF[j].push_back(0);
                            if (energyuncalibS[j].size() == 0) energyuncalibS[j].push_back(0);
                            if (energyuncalibF[j].size() == 0) energyuncalibF[j].push_back(0);

                            timeFast[j] = timeF[j][0]/10;
                            timeSlow[j] = timeS[j][0]/10;
                            energySlow[j] = energycalibS[j][0];
                            energyFast[j] = energycalibF[j][0];
                            energySlowUncalib[j] = energyuncalibS[j][0];
                            energyFastUncalib[j] = energyuncalibF[j][0];

                            if (j == 0) 
                            {
                                if (timeP[0].size() == 0) timeP[0].push_back(0);
                                if (energyP[0].size() == 0) energyP[0].push_back(0);
                                timePolaris[0] = timeP[0][0]/10;
                                energyPolaris[0] = energyP[0][0];
                            }
                            //if (energySlow[0] !=0)  std::cout << "L" << 0 << "| timeSlow " << timeSlow[0] << "| energySlow " << energySlow[0] << "| timeFast " << timeFast[0] << "| energyFast " << energyFast[0] << std::endl;
                            //if (energySlow[1] !=0) std::cout << "L" << 1 << "| timeSlow " << timeSlow[1] << "| energySlow " << energySlow[1] << "| timeFast " << timeFast[1] << "| energyFast " << energyFast[1] << std::endl;
                        }

                       //printf("timeF %f | timeS %f | energyS %f\n", timeFast[0]  , timeSlow[0]  , energySlow[0]);
                        LaBrData->Fill();
                                            

                        for (j = 0; j < 2; j++) 
                        { 
                            if ((timeSlow[j] > 0) && (timeFast[j] > 0)) 
                            //if (timeSlow[j] > 0)
                            {
                                // std::cout << "L" << j << "| timeSlow " << timeSlow[j] << "| energySlow " << energySlow[j] << "| timeFast " << timeFast[j] << "| energyFast " << energyFast[j] << std::endl;
                                slowcalibE[j]->Fill(energySlow[j]);
                                slowE[j]->Fill(energySlowUncalib[j]);
                                fastcalibE[j]->Fill(energyFast[j]);
                                fastE[j]->Fill(energyFastUncalib[j]);
                                tdSF0 = (timeSlow[j]-timeFast[j]); // Measured in nanoseconds
                                //printf("tdSF0 %f\n", tdSF0);
                                slowfastTD[j]->Fill(tdSF0);
                                slowfastTDFE[j]->Fill(tdSF0, energyFast[j]);
                                slowfastTDSE[j]->Fill(tdSF0, energySlow[j]);
                                slowfastE[j]->Fill(energySlow[j], energyFast[j]);
                            }

                            for (k = 1; k < 2; k++) 
                            {
                                if (j != k ) 
                                {
                                    if ((timeFast[j] > 0) && (timeFast[k] > 0)) 
                                    {

                                        tdF=(timeFast[j]-timeFast[k]); // Measured in nanoseconds
                                        tdF = tdF+(Ta[j]-Ta[k]);
                                        if (j == 0) 
                                        {
                                            fastTD0[k]->Fill(tdF); 
                                            //printf("tdF %f\n", tdF);
                                        }
                                    }
                                }
                            }                        
                        }

                        // Reset vectors
                        for (j = 0; j < 2; j++) 
                        {
                            energycalibS[j].clear();
                            energycalibF[j].clear();
                            energyuncalibS[j].clear();
                            energyuncalibF[j].clear();
                            timeF[j].clear();
                            timeS[j].clear();
                            energySlow[j] = 0;
                            timeFast[j] = 0;
                            timeSlow[j] = 0;
                            energyFast[j] = 0;
                            energySlowUncalib[j] = 0;
                            energyFastUncalib[j] = 0;
                        }
                        timePolaris[0] = 0;
                        energyPolaris[0] = 0;
                        energyP[0].clear();
                        timeP[0].clear();
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

    for (j = 0; j < 2; j++)
    {
        slowcalibE[j]->Write();
        slowE[j]->Write();
        fastcalibE[j]->Write();
        fastE[j]->Write();
        fastTD0[j]->Write();
        slowfastTD[j]->Write();
        slowfastTDSE[j]->Write();
        slowfastTDFE[j]->Write();
    }
    polaris_energytime[0]->Write();

    std::string directory = argv[1];
    size_t found = directory.find_last_of("/");
    directory = directory.substr(0, found);

    // the run name is everything after the last '/' in the input argument
    std::string RunString = argv[1];
    found = RunString.find_last_of("/");
    RunString = RunString.substr(found+1);

    
    // create a text file of slowcalibE where the first column is the energy and the second column is the counts
    std::ofstream slowcalibEfile;
    slowcalibEfile.open(directory + "/" + RunString + "_slowcalibE_L0.asc");
   // get bin 1,2,3..... to 2000 and hte counts
    for (int i = 0; i < 2001; i++)
    {
        slowcalibEfile << i << "," << slowcalibE[0]->GetBinContent(i) << std::endl;
    }
    slowcalibEfile.close();

    // create a text file of slowcalibE where the first column is the energy and the second column is the counts
    std::ofstream slowcalibEfile1;
    slowcalibEfile1.open(directory + "/" + RunString + "_slowcalibE_L2.asc");
    // get bin 1,2,3..... to 2000 and hte counts
    for (int i = 0; i < 2001; i++)
    {
        slowcalibEfile1 << i << "," << slowcalibE[1]->GetBinContent(i) << std::endl;
    }
    slowcalibEfile1.close();
    
    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("\nThe first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference of the RXX_ file is (min): %f\n", TimeDiff);

    printf("Total number of sync pulses: %lld\n", countpolaris);

    time(&end);
    double time_taken = double((end - start));
    std::cout << "\nTime taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " s " << std::endl;

    
    delete g;
}