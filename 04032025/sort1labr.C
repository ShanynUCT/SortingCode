/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms */

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1 (base): (c) P.M. Jones (13/07/18)*/
/* Adapted:  (c) S. Hart (01/06/2023) */

/* This code reads in the list mode data from the buffer. It sorts the time and (rated) energy for each detector into a tree. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort1labr.C -o exe `root-config --cflags --libs` -lSpectrum 
     ./exe ~/Documents/PhD/exp/2025/04032025/R17
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
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TRandom.h"
#include "TRandom3.h"

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
  int fastcount[2]={0,0}, slowcount[2]={0,0}, fastslowcount[2]={0,0};
  time(&start);

  // Declare vectors to store the event data
  std::vector<double> energyS[2];
  std::vector<double> energyF[2];
  std::vector<double> timeF[2];
  std::vector<double> timeS[2];
  std::vector<int> pop(1, 0);
  double energySlow[2]={0,0}, timeFast[2]={0,0}, timeSlow[2]={0,0}, energyFast[2]={0,0};
  double tdF, tdS;
  int Ed;

    // 2nd order polynomials (APPLIED) only using 60Co peaks but works well

    // LaBr3 DA-DA coincidence
    std::vector<double> p3fast = {1.2297879694477388e-10, 2.082460415289269e-09};
    std::vector<double> p2fast = {2.6215774005579137e-06, -1.9763963523846775e-05};
    std::vector<double> p1fast = {0.10936703529359502, 0.21787276320993212};
    std::vector<double> p0fast = {0.4453158269436617, 0.39106855912559646};

    std::vector<double> p3slow = {0,0};
    std::vector<double> p2slow = {4.564683211159151e-06, 4.281828161925052e-06};
    std::vector<double> p1slow = {0.06555801446690429, 0.06837424859671552};
    std::vector<double> p0slow = {0.4672198515954785, 0.465811805244691};

  double bw = 1, energy, rndchannel, rndUnif;
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
  LaBrData->Branch("slowEL0", &energySlow[0], "slowEL0/D");
  LaBrData->Branch("fastEL0", &energyFast[0], "fastEL0/D");
  LaBrData->Branch("timeSL0", &timeSlow[0], "timeSL0/D");
  LaBrData->Branch("timeFL0", &timeFast[0], "timeFL0/D");
  LaBrData->Branch("slowEL1", &energySlow[1], "slowEL1/D");
  LaBrData->Branch("fastEL1", &energyFast[1], "fastEL1/D");
  LaBrData->Branch("timeSL1", &timeSlow[1], "timeSL1/D");
  LaBrData->Branch("timeFL1", &timeFast[1], "timeFL1/D");

  // ################ Create histograms ################
  TH1D** slowE=new TH1D*[2];
  TH1D** fastE=new TH1D*[2];
  TH1D** fastTD0=new TH1D*[2];
  TH1D** fastTD1=new TH1D*[2];
  TH1D** slowTD0=new TH1D*[2];
  TH1D** slowTD1=new TH1D*[2];
  for (j = 0; j <=1; j++)
  {
    slowE[j] = new TH1D(TString::Format("Slow_Energy_L%2d", j),"Spectrum",16380,0,16380);
    fastE[j] = new TH1D(TString::Format("Fast_Energy_L%2d", j),"Spectrum",16350,0,16350);
    fastTD0[j] = new TH1D(TString::Format("Fast Time Difference_L0-L%2d", j),"Spectrum",10000,0,100); // 0 to 100 ns with 10 ps binning
    fastTD1[j] = new TH1D(TString::Format("Fast Time Difference_L1-L%2d", j),"Spectrum",10000,0,100); 
    slowTD0[j] = new TH1D(TString::Format("Slow Time Difference_L0-L%2d", j),"Spectrum",10000,0,100); 
    slowTD1[j] = new TH1D(TString::Format("Slow Time Difference_L1-L%2d", j),"Spectrum",10000,0,100);
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
                    LaBr 0  |    88  |   72   |   80    |    64   | -1100   | 
                    ---------------------------------------------------------
                    LaBr 1  |    89  |   73   |    81   |     65  | -1100   |   
                    ---------------------------------------------------------
                    */

                    // Identify detector and store data accordingly
                    if ((ident > 63) && (ident < 66)) 
                    {
                        detectorID = ident-64;
                        if (energyS[detectorID].size() == 0)
                        {
                            double calibslowE = p3slow[detectorID]*pow(adcdata,3) + p2slow[detectorID]*pow(adcdata,2) + p1slow[detectorID]*pow(adcdata,1) + p0slow[detectorID];
                            rndUnif = randy->Uniform(-2,2);
                            energyS[detectorID].push_back(calibslowE+rndUnif);
                            //printf("energyS[%d] = %f\n", detectorID, energyS[detectorID][0]);
                            //printf("adcdata = %d\n\n", adcdata);
                        }

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
                            if (timeS[detectorID].size() == 0) 
                            {
                                timeS[detectorID].push_back((TS * 100.0) + (tock * 20.0) + sick); // time in 1e12s which is 10ns
                                // printf("timeS[%d] = %f\n", detectorID, timeS[detectorID][0]);
                            }
                        }
                    }
                    else if ((ident > 71) && (ident < 76)) 
                    {
                        detectorID = ident-72;
                        if (energyF[detectorID].size() == 0)
                        {
                            double calibfastE = p3fast[detectorID]*pow(adcdata,3) + p2fast[detectorID]*pow(adcdata,2) + p1fast[detectorID]*pow(adcdata,1) + p0fast[detectorID];
                            rndUnif = randy->Uniform(-2,2);
                            energyF[detectorID].push_back(calibfastE+rndUnif);
                        }
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
                                // printf("timeF[%d] = %f\n", detectorID, timeF[detectorID][0]);
                            }
                        }
                    }
                    

                    // printf("TSdiff %lld \n", TSdiff);
                    if (TSdiff >= 10 ) //90                     
                    {   
                        i=i-2; 
                        //std::cout <<"LOOP2 " << "i" << i << "| TSdiff " << TSdiff << "| TS " << TS << "| TSinit " << TSinit << std::endl;
                        for (j = 0; j <=1; j++)
                        {
                            //printf("j %d | energyS[j].size() %lu | timeF[j].size() %lu | timeS[j].size() %lu \n", j, energyS[j].size(), timeF[j].size(), timeS[j].size());
                            if (energyS[j].size() == 0) energyS[j].push_back(0);
                            if (timeF[j].size() == 0) timeF[j].push_back(0);
                            if (timeS[j].size() == 0)  timeS[j].push_back(0);
                            if (energyF[j].size() == 0) energyF[j].push_back(0);

                            timeFast[j] = timeF[j][0]/10;
                            timeSlow[j] = timeS[j][0]/10;
                            energySlow[j] = energyS[j][0];
                            energyFast[j] = energyF[j][0];

                            // printf("j %d | energySlow %f | timeFast %f | timeSlow %f | energyFast %f \n", j, energySlow[j], timeFast[j], timeSlow[j], energyFast[j]);
                        }
                       // printf("timeF %f | timeS %f | energyS %f\n", timeFast[0]  , timeSlow[0]  , energySlow[0]);
                        LaBrData->Fill();

                        for (j = 0; j <= 1; j++) 
                        { 
                            //Singles spectra 
                            //std::cout << "L" << j << "| energySlow[j] " << energySlow[j] << "| timeFast[j] " << timeFast[j] << "| timeSlow[j] " << timeSlow[j] << "| timeRF" << timeFast[4]<< std::endl;
                            slowE[j]->Fill(energySlow[j]);
                            fastE[j]->Fill(energyFast[j]);

                            for (k = 0; k <= 1; k++) 
                            {
                                if (j != k ) 
                                {
                                    tdF=(timeFast[j]-timeFast[k]); // Measured in nanoseconds
                                    tdS = (timeSlow[j]-timeSlow[k]);
                                    if (j == 1) 
                                    {
                                        if (energyFast[0] > 1122 && energyFast[0] < 1252) // 1173 keV peak
                                        {
                                            if (energyFast[1] > 1252 && energyFast[1] < 1402) // 1332 keV peak
                                            {
                                                fastTD1[k]->Fill(tdF);
                                                fastTD0[j]->Fill(tdF);
                                                //if (tdF>0 && tdF<1000) printf("tdF %f\n", tdF);
                                            }
                                        }
                                        // print the slowEnergy[0] and slowEnergy[1] values
                                        //printf("slowE[0] %f | slowE[1] %f\n", energySlow[0], energySlow[1]);
                                        if (energySlow[0] > 1122 && energySlow[0] < 1252) // 1173 keV peak
                                        {
                                            if (energySlow[1] > 1252 && energySlow[1] < 1402) // 1332 keV peak
                                            {
                                                slowTD1[k]->Fill(tdS);
                                                slowTD0[j]->Fill(tdS);
                                            }
                                        }
                                    }
                                }
                            }                        
                        }

                        // Reset vectors
                        for (j = 0; j <=1; j++) 
                        {
                            energyS[j].clear();
                            energyF[j].clear();
                            timeF[j].clear();
                            timeS[j].clear();
                            energySlow[j] = 0;
                            timeFast[j] = 0;
                            timeSlow[j] = 0;
                            energyFast[j] = 0;
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

    // plot the fastTD0[] histogram, fit a gaussian to the peak and draw the gaussian in red with a line thickness of 4. then write the fit parameters and the time resolution as (2.35*sigma)/sqrt(2) to a text block on the canvas
    TF1 *fit0 = new TF1("fit0","gaus",0,100);
    fastTD0[0]->Fit("fit0","Q");
    fastTD0[0]->Draw();
    fit0->SetLineColor(kRed);
    fit0->SetLineWidth(4);
    fit0->Draw("same");
    double mean0 = fit0->GetParameter(1);
    double sigma0 = fit0->GetParameter(2);
    double timeRes0 = (2.35*sigma0)/sqrt(2);
    TLatex *tex0 = new TLatex(0.5,0.5,Form("Time Resolution = %f", timeRes0));
    tex0->SetNDC();
    tex0->Draw();
    // save the canvas as a .root file and as a .png file to the directory of the .root file
    g->Write();
    
    
    

    LaBrData->Write();

    for (j = 0; j <=1; j++)
    {
        slowE[j]->Write();
        fastE[j]->Write();
        fastTD0[j]->Write();
        fastTD1[j]->Write();
        slowTD0[j]->Write();
        slowTD1[j]->Write();
    }


    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("\nThe first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference of the RXX_ file is (min): %f\n", TimeDiff);

    time(&end);

    double time_taken = double((end - start));
    std::cout << "\nTime taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " seconds " << std::endl;

    
    delete g;
}