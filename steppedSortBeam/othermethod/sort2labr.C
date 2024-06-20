/* Advanced sorting code for TDR format data for source data (ONLY ONE EXPONENTIAL TIME WALK FIT)*/
/* Produces ROOT histograms and TTrees*/

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1:  (c) S. Hart (01/06/2023) */

/* This code reads in the R63.root file produced by sort1labr.c */
/* The tree has the following stucture:   
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
  LaBrData->Branch("fastTPOLARIS", &fastTimePOLARIS[0], "fastTPOLARIS/D");*/

/* The tree of the ROOT file contains information about the fast time, slow time & slow energy for the detectors, */
/* as well as the fast time of the RF (detector4) & the fast energy and time of the POLARIS sync pulse. */
/* This code calculates the cabling/PMT/... time offset between the slow signal and the fast signal of each detector.*/
/* The code aligns the slow and fast signals to more acturately track slow energy events (good energy resolution) */
/* and fast time events (good timing resolution). */
/* It then produces a root file called R63_aligned.root (that is not ordered in time) for event sorting. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort2labr.C -o exe2 `root-config --cflags --libs` -lSpectrum 
     ./exe2 ~/Documents/PhD/exp/2023/230228/run4_22Na_138uCi_10min/R6.root source/beam
*/


#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TTreeIndex.h>
#include <TLegend.h>
#include <TMath.h>
#include <TThread.h>
#include <Math/WrappedTF1.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <vector>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <inttypes.h> 
#include <stdexcept>
#include <fstream>
#include <map>
#include <thread>

const int NumDetectors = 4;
const int RFIndex = 4;
const int POLARISIndex = 5;

time_t start, end;

// variables for input tree data
std::vector<Double_t> energySlow(NumDetectors);
std::vector<Double_t> timeSlow(NumDetectors);
std::vector<Double_t> timeFast(NumDetectors+1); // Additional branch for RF
std::vector<Double_t> fastEnergyPOLARIS(1);
std::vector<Double_t> fastTimePOLARIS(1);
std::vector<Double_t> alignedSlowECalibL(NumDetectors);

// variables for copytree data
std::vector<Double_t> energySlowCp(NumDetectors);
std::vector<Double_t> timeSlowCp(NumDetectors);
std::vector<Double_t> timeSlowCorrected(NumDetectors);
std::vector<Double_t> timeSlowCorrectedCp(NumDetectors);
std::vector<Double_t> timeFastCp(NumDetectors+1); 
std::vector<Double_t> fastEnergyPOLARISCp(1);
std::vector<Double_t> fastTimePOLARISCp(1);
std::vector<Double_t> timeOffsetBeforeCorrection(NumDetectors);
std::vector<Double_t> timeOffsetBeforeCorrectionCp(NumDetectors);

// Variables for aligned tree data
Int_t detectorID;
std::vector<Double_t> alignedTimeFast(1);
std::vector<Double_t> alignedEnergyFast(1);
std::vector<Double_t> alignedTimeSlow(1);
std::vector<Double_t> alignedEnergySlow(1);
std::vector<Double_t> timeGlobal(1);
 std::vector<TH1D*> slowfasttimediff_hists(NumDetectors);
std::vector<TH1D*> timeOffsetCorrectedHistsNext(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists(NumDetectors); 
std::vector<TH2D*> slowfasttimediff_energy_hists1(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists1(NumDetectors); 
std::vector<TH2D*> slowfasttimediff_energy_hists2(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists2(NumDetectors); 
std::vector<TH2D*> slowfasttimediff_energy_hists3(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists3(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists4(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists4(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists5(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists5(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists6(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists6(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists7(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists7(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists8(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists8(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists9(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists9(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists10(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists10(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists11(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists11(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists12(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists12(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists13(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists13(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists14(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists14(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists15(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists15(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists16(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists16(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists17(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists17(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists18(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists18(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists19(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists19(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists20(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists20(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists21(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists21(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists22(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists22(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists23(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists23(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists24(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists24(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists25(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists25(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists26(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists26(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists27(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists27(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists28(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists28(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists29(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists29(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists30(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists30(NumDetectors);
std::vector<TH2D*> timeEnergyOffsetCorrectedHistsNext(NumDetectors);
std::vector<TH1D*> energyhist(NumDetectors);

std::vector<Double_t> mean_fasttimediff(NumDetectors-1);
std::vector<Double_t> sigma_fasttimediff(NumDetectors-1);
bool tensOfNanoSeconds = false;
bool hundredsOfNanoSeconds = false;

// initialise the preComputedValues map
std::map<std::pair<int, double>, double> precomputedValues1;
std::map<std::pair<int, double>, double> precomputedValues2;
std::map<std::pair<int, double>, double> precomputedValues3;
std::map<std::pair<int, double>, double> precomputedValues4;
std::map<std::pair<int, double>, double> precomputedValues5;
std::map<std::pair<int, double>, double> precomputedValues6;
std::map<std::pair<int, double>, double> precomputedValues7;
std::map<std::pair<int, double>, double> precomputedValues8;
std::map<std::pair<int, double>, double> precomputedValues9;
std::map<std::pair<int, double>, double> precomputedValues10;
std::map<std::pair<int, double>, double> precomputedValues11;
std::map<std::pair<int, double>, double> precomputedValues12;
std::map<std::pair<int, double>, double> precomputedValues13;
std::map<std::pair<int, double>, double> precomputedValues14;
std::map<std::pair<int, double>, double> precomputedValues15;
std::map<std::pair<int, double>, double> precomputedValues16;
std::map<std::pair<int, double>, double> precomputedValues17;
std::map<std::pair<int, double>, double> precomputedValues18;
std::map<std::pair<int, double>, double> precomputedValues19;
std::map<std::pair<int, double>, double> precomputedValues20;
std::map<std::pair<int, double>, double> precomputedValues21;
std::map<std::pair<int, double>, double> precomputedValues22;
std::map<std::pair<int, double>, double> precomputedValues23;
std::map<std::pair<int, double>, double> precomputedValues24;
std::map<std::pair<int, double>, double> precomputedValues25;
std::map<std::pair<int, double>, double> precomputedValues26;
std::map<std::pair<int, double>, double> precomputedValues27;
std::map<std::pair<int, double>, double> precomputedValues28;
std::map<std::pair<int, double>, double> precomputedValues29;
std::map<std::pair<int, double>, double> precomputedValues30;


// BANANA FIT FUNCTIONS //
TF1* banana1(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots);
TF1* banana1fit[NumDetectors];
TF1* banana2(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists2, TH2D* energy_slowfasttimediff_hists2, std::string dirPathplots);
TF1* banana2fit[NumDetectors];
TF1* banana3(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists3, TH2D* energy_slowfasttimediff_hists3, std::string dirPathplots);
TF1* banana3fit[NumDetectors];
TF1* banana4(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists4, TH2D* energy_slowfasttimediff_hists4, std::string dirPathplots);
TF1* banana4fit[NumDetectors];
TF1* banana5(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots);
TF1* banana5fit[NumDetectors];
TF1* banana6(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6, TH2D* energy_slowfasttimediff_hists6, std::string dirPathplots);
TF1* banana6fit[NumDetectors];
TF1* banana7(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots);
TF1* banana7fit[NumDetectors];
TF1* banana8(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists8, TH2D* energy_slowfasttimediff_hists8, std::string dirPathplots);
TF1* banana8fit[NumDetectors];
TF1* banana9(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists9, TH2D* energy_slowfasttimediff_hists9, std::string dirPathplots);
TF1* banana9fit[NumDetectors];
TF1* banana10(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists10, TH2D* energy_slowfasttimediff_hists10, std::string dirPathplots);
TF1* banana10fit[NumDetectors];
TF1* banana11(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists11, TH2D* energy_slowfasttimediff_hists11, std::string dirPathplots);
TF1* banana11fit[NumDetectors];
TF1* banana12(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists12, TH2D* energy_slowfasttimediff_hists12, std::string dirPathplots);
TF1* banana12fit[NumDetectors];
TF1* banana13(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists13, TH2D* energy_slowfasttimediff_hists13, std::string dirPathplots);
TF1* banana13fit[NumDetectors];
TF1* banana14(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists14, TH2D* energy_slowfasttimediff_hists14, std::string dirPathplots);
TF1* banana14fit[NumDetectors];
TF1* banana15(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists15, TH2D* energy_slowfasttimediff_hists15, std::string dirPathplots);
TF1* banana15fit[NumDetectors];
TF1* banana16(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists16, TH2D* energy_slowfasttimediff_hists16, std::string dirPathplots);
TF1* banana16fit[NumDetectors];
TF1* banana17(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists17, TH2D* energy_slowfasttimediff_hists17, std::string dirPathplots);
TF1* banana17fit[NumDetectors];
TF1* banana18(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists18, TH2D* energy_slowfasttimediff_hists18, std::string dirPathplots);
TF1* banana18fit[NumDetectors];
TF1* banana19(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists19, TH2D* energy_slowfasttimediff_hists19, std::string dirPathplots);
TF1* banana19fit[NumDetectors];
TF1* banana20(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists20, TH2D* energy_slowfasttimediff_hists20, std::string dirPathplots);
TF1* banana20fit[NumDetectors];
TF1* banana21(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists21, TH2D* energy_slowfasttimediff_hists21, std::string dirPathplots);
TF1* banana21fit[NumDetectors];
TF1* banana22(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists22, TH2D* energy_slowfasttimediff_hists22, std::string dirPathplots);
TF1* banana22fit[NumDetectors];
TF1* banana23(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists23, TH2D* energy_slowfasttimediff_hists23, std::string dirPathplots);
TF1* banana23fit[NumDetectors];
TF1* banana24(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists24, TH2D* energy_slowfasttimediff_hists24, std::string dirPathplots);
TF1* banana24fit[NumDetectors];
TF1* banana25(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists25, TH2D* energy_slowfasttimediff_hists25, std::string dirPathplots);
TF1* banana25fit[NumDetectors];
TF1* banana26(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists26, TH2D* energy_slowfasttimediff_hists26, std::string dirPathplots);
TF1* banana26fit[NumDetectors];
TF1* banana27(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists27, TH2D* energy_slowfasttimediff_hists27, std::string dirPathplots);
TF1* banana27fit[NumDetectors];
TF1* banana28(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists28, TH2D* energy_slowfasttimediff_hists28, std::string dirPathplots);
TF1* banana28fit[NumDetectors];
TF1* banana29(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists29, TH2D* energy_slowfasttimediff_hists29, std::string dirPathplots);
TF1* banana29fit[NumDetectors];
TF1* banana30(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists30, TH2D* energy_slowfasttimediff_hists30, std::string dirPathplots);
TF1* banana30fit[NumDetectors];

// END BANANA FIT FUNCTIONS //

// banana1 good fit
TF1* banana1(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots)
{
    std::cout <<"\nFitting banana1 (0-8000 keV) for detector L \n" << i << std::endl;
    std::cout <<"\nFitting banana1 time range (235-245 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists1->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists1->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists1->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists1->Draw();

    yBins = slowfasttimediff_energy_hists1->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists1->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana1Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana1s for detector %d...\n", i);


    TCanvas* cb = new TCanvas();
	TGraphErrors* gr1 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr1 -> SetTitle(Form("Banana1 Fit L%d", i));
    //gr1 -> SetStats(0);
	gr1 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr1 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr1 -> SetMarkerStyle(20);
    gr1 -> SetMarkerSize(0.8);
    gr1 -> SetMarkerColor(kBlue);
    gr1->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit1 = new TF1("fit1", "pol0(0)", 0, 8000);  
    
    gr1->Draw("AL*");
    gr1->Fit("fit1", "");
    fit1->Draw("same");
    //gr1->GetXaxis()->SetRangeUser(0, 8000);
    cb -> SaveAs(Form("%s/banana1Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr1;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana1 fit for detector %d...\n", i);
    return fit1;
}

TF1* banana2(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists2, TH2D* energy_slowfasttimediff_hists2, std::string dirPathplots)
{
    std::cout <<"Fitting banana2 (0-8000 keV) - detector L \n" << i << std::endl;
    std::cout <<"Fitting banana2 time range (245-255 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists2->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists2->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists2->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists2->Draw();

    yBins = slowfasttimediff_energy_hists2->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists2->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana2Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana2s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr2 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr2 -> SetTitle(Form("Banana2 Fit L%d", i));
    //gr2 -> SetStats(0);
	gr2 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr2 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr2 -> SetMarkerStyle(20);
    gr2 -> SetMarkerSize(0.8);
    gr2 -> SetMarkerColor(kBlue);
    gr2->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit2 = new TF1("fit2", "pol0(0)", 0, 8000); 

    gr2->Draw("AL*");
    gr2->Fit("fit2", "");
    fit2->Draw("same");
    cb -> SaveAs(Form("%s/banana2Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr2;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana2 fit for detector %d...\n", i);
    return fit2;
}

TF1* banana3(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists3, TH2D* energy_slowfasttimediff_hists3, std::string dirPathplots)
{
    std::cout <<"Fitting banana3 (0-8000 keV) for detector L \n" << i << std::endl;
    std::cout <<"Fitting banana3 time range (255-265 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists3->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists3->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists3->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists3->Draw();

    yBins = slowfasttimediff_energy_hists3->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists3->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana3Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana3s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr3 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr3 -> SetTitle(Form("Banana3 Fit L%d", i));
    //gr3 -> SetStats(0);
	gr3 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr3 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr3 -> SetMarkerStyle(20);
    gr3 -> SetMarkerSize(0.8);
    gr3 -> SetMarkerColor(kBlue);
    gr3->GetXaxis()->SetRangeUser(0, 8000);

    TF1* fit3 = new TF1("fit3", " pol0(0)", 0, 8000); 
    gr3->Draw("AL*");
    gr3->Fit("fit3", "");
    fit3->Draw("same");
    cb -> SaveAs(Form("%s/banana3Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr3;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana3 fit for detector %d...\n", i);
    return fit3;
}

TF1* banana4(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists4, TH2D* energy_slowfasttimediff_hists4, std::string dirPathplots)
{
    std::cout <<"Fitting banana4 (0-8000 keV) for detector L \n" << i << std::endl;
    std::cout <<"Fitting banana4 time range (265-275 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists4->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists4->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists4->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists4->Draw();

    yBins = slowfasttimediff_energy_hists4->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        // if (bin > 600) slicewidth = 150;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists4->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana4Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana4s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr4 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr4 -> SetTitle(Form("Banana4 Fit L%d", i));
    //gr4 -> SetStats(0);
	gr4 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr4 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr4 -> SetMarkerStyle(20);
    gr4 -> SetMarkerSize(0.8);
    gr4 -> SetMarkerColor(kBlue);
    gr4->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit4 = new TF1("fit4", "pol0(0)", 0, 8000);   
    gr4->Draw("AL*");
    gr4->Fit("fit4", "");
    fit4->Draw("same");
    cb -> SaveAs(Form("%s/banana4Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr4;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana4 fit for detector %d...\n", i);
    return fit4;
}

TF1* banana5(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots)
{
    std::cout <<"Fitting banana5 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana5 time range (275-285 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists5->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists5->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists5->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists5->Draw();

    yBins = slowfasttimediff_energy_hists5->GetYaxis()->GetNbins();

    while (bin < yBins)
    {

        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists5->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana5Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana5s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr5 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr5 -> SetTitle(Form("Banana5 Fit L%d", i));
    //gr5 -> SetStats(0);
	gr5 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr5 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr5 -> SetMarkerStyle(20);
    gr5 -> SetMarkerSize(0.8);
    gr5 -> SetMarkerColor(kBlue);
    gr5->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit5 = new TF1("fit5", "pol0(0)", 0, 8000);
    gr5->Draw("AL*");
    gr5->Fit("fit5", "");
    fit5->Draw("same");
    cb -> SaveAs(Form("%s/banana5Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr5;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana5 fit for detector %d...\n", i);
    return fit5;
}

TF1* banana6(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6, TH2D* energy_slowfasttimediff_hists6, std::string dirPathplots)
{
    // This is a Leading Edge region 
    std::cout <<"Fitting banana6 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana6 time range (285-295 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists6->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists6->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists6->Draw();

    yBins = slowfasttimediff_energy_hists6->GetYaxis()->GetNbins();

    while (bin < yBins)
    {if (bin <= 400) slicewidth = 5;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists6->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana6Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana6s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr6 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr6 -> SetTitle(Form("Banana6 Fit L%d", i));
    //gr6 -> SetStats(0);
	gr6 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr6 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr6 -> SetMarkerStyle(20);
    gr6 -> SetMarkerSize(0.8);
    gr6 -> SetMarkerColor(kBlue);
    TF1* fit6 = new TF1("fit6", "pol0(0)", 0, 8000);
    gr6->Draw("AL*");
    gr6->Fit("fit6", "");
    fit6->Draw("same");

    cb -> SaveAs(Form("%s/banana6Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr6;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana6 fit for detector %d...\n", i);
    return fit6;
}

TF1* banana7(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots)
{
    std::cout <<"Fitting banana7 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana7 time range (295-305 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists7->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists7->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists7->Draw();

    yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();

    while (bin < yBins)
    {

        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists7->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana7Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana7s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr7 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr7 -> SetTitle(Form("Banana7 Fit L%d", i));
    //gr7 -> SetStats(0);
	gr7 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr7 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr7 -> SetMarkerStyle(20);
    gr7 -> SetMarkerSize(0.8);
    gr7 -> SetMarkerColor(kBlue);
    gr7->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit7 = new TF1("fit7", "pol0(0)", 0, 8000);
    gr7->Draw("AL*");
    gr7->Fit("fit7", "");
    fit7->Draw("same");
    cb -> SaveAs(Form("%s/banana7Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr7;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana7 fit for detector %d...\n", i);
    return fit7;
}

TF1* banana8(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists8, TH2D* energy_slowfasttimediff_hists8, std::string dirPathplots)
{
    // This is a Leading Edge region 
    std::cout <<"Fitting banana8 (0-8000 keV) for detector  \n" << i << std::endl;
    std::cout <<"Fitting banana8 time range (305-315 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists8->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists8->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists8->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists8->Draw();

    yBins = slowfasttimediff_energy_hists8->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 400) slicewidth = 5;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists8->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana8Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana8s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr8 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr8 -> SetTitle(Form("Banana8 Fit L%d", i));
    //gr8 -> SetStats(0);
	gr8 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr8 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr8 -> SetMarkerStyle(20);
    gr8 -> SetMarkerSize(0.8);
    gr8 -> SetMarkerColor(kBlue);
    gr8->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit8 = new TF1("fit8", "pol0(0)", 0, 8000);
    gr8->Draw("AL*");
    gr8->Fit("fit8", "");
    fit8->Draw("same");
    cb -> SaveAs(Form("%s/banana8Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr8;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana8 fit for detector %d...\n", i);
    return fit8;
}

TF1* banana9(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists9, TH2D* energy_slowfasttimediff_hists9, std::string dirPathplots)
{
    std::cout <<"Fitting banana9 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana9 time range (315-325 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists9->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists9->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists9->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists9->Draw();

    yBins = slowfasttimediff_energy_hists9->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists9->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana9Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana9s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr9 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr9 -> SetTitle(Form("Banana9 Fit L%d", i));
    //gr9 -> SetStats(0);
	gr9 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr9 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr9 -> SetMarkerStyle(20);
    gr9 -> SetMarkerSize(0.8);
    gr9 -> SetMarkerColor(kBlue);
    gr9->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit9 = new TF1("fit9", "pol0(0)", 0, 8000);
    gr9->Draw("AL*");
    gr9->Fit("fit9", "");
    fit9->Draw("same");
    cb -> SaveAs(Form("%s/banana9Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr9;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana9 fit for detector %d...\n", i);
    return fit9;
}

TF1* banana10(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists10, TH2D* energy_slowfasttimediff_hists10, std::string dirPathplots)
{
    std::cout <<"Fitting banana10 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana10 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists10->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists10->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists10->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists10->Draw();

    yBins = slowfasttimediff_energy_hists10->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists10->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana10Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana10s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr10 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr10 -> SetTitle(Form("Banana10 Fit L%d", i));
    //gr10 -> SetStats(0);
	gr10 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr10 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr10 -> SetMarkerStyle(20);
    gr10 -> SetMarkerSize(0.8);
    gr10 -> SetMarkerColor(kBlue);
    gr10->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit10 = new TF1("fit10", "pol0(0)", 0, 8000);
    gr10->Draw("AL*");
    gr10->Fit("fit10", "");
    fit10->Draw("same");
    cb -> SaveAs(Form("%s/banana10Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr10;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana10 fit for detector %d...\n", i);
    return fit10;
}

TF1* banana11(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists11, TH2D* energy_slowfasttimediff_hists11, std::string dirPathplots)
{
    std::cout <<"Fitting banana11 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana11 time range (335-345 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists11->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists11->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists11->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists11->Draw();

    yBins = slowfasttimediff_energy_hists11->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists11->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana11Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana11s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr11 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr11 -> SetTitle(Form("Banana11 Fit L%d", i));
    //gr11 -> SetStats(0);
	gr11 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr11 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr11 -> SetMarkerStyle(20);
    gr11 -> SetMarkerSize(0.8);
    gr11 -> SetMarkerColor(kBlue);
    gr11->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit11 = new TF1("fit11", "pol0(0)", 0, 8000);
    gr11->Draw("AL*");
    gr11->Fit("fit11", "");
    fit11->Draw("same");
    cb -> SaveAs(Form("%s/banana11Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr11;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana11 fit for detector %d...\n", i);
    return fit11;
}

TF1* banana12(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists12, TH2D* energy_slowfasttimediff_hists12, std::string dirPathplots)
{
    std::cout <<"Fitting banana12 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana12 time range (345-355 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists12->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists12->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists12->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists12->Draw();

    yBins = slowfasttimediff_energy_hists12->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists12->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana12Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana12s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr12 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr12 -> SetTitle(Form("Banana12 Fit L%d", i));
    //gr12 -> SetStats(0);
	gr12 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr12 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr12 -> SetMarkerStyle(20);
    gr12 -> SetMarkerSize(0.8);
    gr12 -> SetMarkerColor(kBlue);
    gr12->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit12 = new TF1("fit12", "pol0(0)", 0, 8000);
    gr12->Draw("AL*");
    gr12->Fit("fit12", "");
    fit12->Draw("same");
    cb -> SaveAs(Form("%s/banana12Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr12;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana12 fit for detector %d...\n", i);
    return fit12;
}

TF1* banana13(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists13, TH2D* energy_slowfasttimediff_hists13, std::string dirPathplots)
{
    std::cout <<"Fitting banana13 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana13 time range (355-365 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists13->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists13->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists13->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists13->Draw();

    yBins = slowfasttimediff_energy_hists13->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists13->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana13Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana13s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr13 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr13 -> SetTitle(Form("Banana13 Fit L%d", i));
    //gr13 -> SetStats(0);
	gr13 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr13 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr13 -> SetMarkerStyle(20);
    gr13 -> SetMarkerSize(0.8);
    gr13 -> SetMarkerColor(kBlue);
    gr13->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit13 = new TF1("fit13", "pol0(0)", 0, 8000);
    gr13->Draw("AL*");
    gr13->Fit("fit13", "");
    fit13->Draw("same");
    cb -> SaveAs(Form("%s/banana13Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr13;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana13 fit for detector %d...\n", i);
    return fit13;
}

TF1* banana14(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists14, TH2D* energy_slowfasttimediff_hists14, std::string dirPathplots)
{
    std::cout <<"Fitting banana14 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana14 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists14->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists14->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists14->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists14->Draw();

    yBins = slowfasttimediff_energy_hists14->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists14->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana14Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana14s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr14 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr14 -> SetTitle(Form("Banana14 Fit L%d", i));
    //gr14 -> SetStats(0);
	gr14 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr14 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr14 -> SetMarkerStyle(20);
    gr14 -> SetMarkerSize(0.8);
    gr14 -> SetMarkerColor(kBlue);
    gr14->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit14 = new TF1("fit14", "pol0(0)", 0, 8000);
    gr14->Draw("AL*");
    gr14->Fit("fit14", "");
    fit14->Draw("same");
    cb -> SaveAs(Form("%s/banana14Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr14;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana14 fit for detector %d...\n", i);
    return fit14;
}

TF1* banana15(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists15, TH2D* energy_slowfasttimediff_hists15, std::string dirPathplots)
{
    std::cout <<"Fitting banana15 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana15 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists15->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists15->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists15->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists15->Draw();

    yBins = slowfasttimediff_energy_hists15->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists15->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana15Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana15s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr15 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr15 -> SetTitle(Form("Banana15 Fit L%d", i));
    //gr15 -> SetStats(0);
	gr15 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr15 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr15 -> SetMarkerStyle(20);
    gr15 -> SetMarkerSize(0.8);
    gr15 -> SetMarkerColor(kBlue);
    gr15->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit15 = new TF1("fit15", "pol0(0)", 0, 8000);
    gr15->Draw("AL*");
    gr15->Fit("fit15", "");
    fit15->Draw("same");
    cb -> SaveAs(Form("%s/banana15Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr15;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana15 fit for detector %d...\n", i);
    return fit15;
}

TF1* banana16(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists16, TH2D* energy_slowfasttimediff_hists16, std::string dirPathplots)
{
    std::cout <<"Fitting banana16 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana16 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists16->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists16->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists16->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists16->Draw();

    yBins = slowfasttimediff_energy_hists16->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists16->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana16Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana16s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr16 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr16 -> SetTitle(Form("Banana16 Fit L%d", i));
    //gr16 -> SetStats(0);
	gr16 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr16 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr16 -> SetMarkerStyle(20);
    gr16 -> SetMarkerSize(0.8);
    gr16 -> SetMarkerColor(kBlue);
    gr16->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit16 = new TF1("fit16", "pol0(0)", 0, 8000);
    gr16->Draw("AL*");
    gr16->Fit("fit16", "");
    fit16->Draw("same");
    cb -> SaveAs(Form("%s/banana16Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr16;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana16 fit for detector %d...\n", i);
    return fit16;
}

TF1* banana17(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists17, TH2D* energy_slowfasttimediff_hists17, std::string dirPathplots)
{
    std::cout <<"Fitting banana17 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana17 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists17->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists17->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists17->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists17->Draw();

    yBins = slowfasttimediff_energy_hists17->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists17->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana17Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana17s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr17 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr17 -> SetTitle(Form("Banana17 Fit L%d", i));
    //gr17 -> SetStats(0);
	gr17 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr17 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr17 -> SetMarkerStyle(20);
    gr17 -> SetMarkerSize(0.8);
    gr17 -> SetMarkerColor(kBlue);
    gr17->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit17 = new TF1("fit17", "pol0(0)", 0, 8000);
    gr17->Draw("AL*");
    gr17->Fit("fit17", "");
    fit17->Draw("same");
    cb -> SaveAs(Form("%s/banana17Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr17;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana17 fit for detector %d...\n", i);
    return fit17;
}

TF1* banana18(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists18, TH2D* energy_slowfasttimediff_hists18, std::string dirPathplots)
{
    std::cout <<"Fitting banana18 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana18 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists18->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists18->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists18->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists18->Draw();

    yBins = slowfasttimediff_energy_hists18->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists18->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana18Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana18s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr18 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr18 -> SetTitle(Form("Banana18 Fit L%d", i));
    //gr18 -> SetStats(0);
	gr18 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr18 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr18 -> SetMarkerStyle(20);
    gr18 -> SetMarkerSize(0.8);
    gr18 -> SetMarkerColor(kBlue);
    gr18->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit18 = new TF1("fit18", "pol0(0)", 0, 8000);
    gr18->Draw("AL*");
    gr18->Fit("fit18", "");
    fit18->Draw("same");
    cb -> SaveAs(Form("%s/banana18Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr18;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana18 fit for detector %d...\n", i);
    return fit18;
}

TF1* banana19(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists19, TH2D* energy_slowfasttimediff_hists19, std::string dirPathplots)
{
    std::cout <<"Fitting banana19 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana19 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists19->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists19->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists19->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists19->Draw();

    yBins = slowfasttimediff_energy_hists19->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists19->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana19Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana19s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr19 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr19 -> SetTitle(Form("Banana19 Fit L%d", i));
    //gr19 -> SetStats(0);
	gr19 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr19 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr19 -> SetMarkerStyle(20);
    gr19 -> SetMarkerSize(0.8);
    gr19 -> SetMarkerColor(kBlue);
    gr19->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit19 = new TF1("fit19", "pol0(0)", 0, 8000);
    gr19->Draw("AL*");
    gr19->Fit("fit19", "");
    fit19->Draw("same");
    cb -> SaveAs(Form("%s/banana19Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr19;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana19 fit for detector %d...\n", i);
    return fit19;
}

TF1* banana20(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists20, TH2D* energy_slowfasttimediff_hists20, std::string dirPathplots)
{
    std::cout <<"Fitting banana20 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana20 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists20->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists20->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists20->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists20->Draw();

    yBins = slowfasttimediff_energy_hists20->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists20->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana20Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana20s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr20 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr20 -> SetTitle(Form("Banana20 Fit L%d", i));
    //gr20 -> SetStats(0);
	gr20 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr20 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr20 -> SetMarkerStyle(20);
    gr20 -> SetMarkerSize(0.8);
    gr20 -> SetMarkerColor(kBlue);
    gr20->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit20 = new TF1("fit20", "pol0(0)", 0, 8000);
    gr20->Draw("AL*");
    gr20->Fit("fit20", "");
    fit20->Draw("same");
    cb -> SaveAs(Form("%s/banana20Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr20;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana20 fit for detector %d...\n", i);
    return fit20;
}

TF1* banana21(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists21, TH2D* energy_slowfasttimediff_hists21, std::string dirPathplots)
{
    std::cout <<"Fitting banana21 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana21 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists21->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists21->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists21->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists21->Draw();

    yBins = slowfasttimediff_energy_hists21->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists21->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana21Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana21s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr21 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr21 -> SetTitle(Form("Banana21 Fit L%d", i));
    //gr21 -> SetStats(0);
	gr21 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr21 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr21 -> SetMarkerStyle(20);
    gr21 -> SetMarkerSize(0.8);
    gr21 -> SetMarkerColor(kBlue);
    gr21->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit21 = new TF1("fit21", "pol0(0)", 0, 8000);
    gr21->Draw("AL*");
    gr21->Fit("fit21", "");
    fit21->Draw("same");
    cb -> SaveAs(Form("%s/banana21Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr21;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana21 fit for detector %d...\n", i);
    return fit21;
}

TF1* banana22(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists22, TH2D* energy_slowfasttimediff_hists22, std::string dirPathplots)
{
    std::cout <<"Fitting banana22 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana22 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists22->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists22->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists22->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists22->Draw();

    yBins = slowfasttimediff_energy_hists22->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists22->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana22Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana22s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr22 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr22 -> SetTitle(Form("Banana22 Fit L%d", i));
    //gr22 -> SetStats(0);
	gr22 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr22 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr22 -> SetMarkerStyle(20);
    gr22 -> SetMarkerSize(0.8);
    gr22 -> SetMarkerColor(kBlue);
    gr22->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit22 = new TF1("fit22", "pol0(0)", 0, 8000);
    gr22->Draw("AL*");
    gr22->Fit("fit22", "");
    fit22->Draw("same");
    cb -> SaveAs(Form("%s/banana22Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr22;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana22 fit for detector %d...\n", i);
    return fit22;
}

TF1* banana23(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists23, TH2D* energy_slowfasttimediff_hists23, std::string dirPathplots)
{
    std::cout <<"Fitting banana23 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana23 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists23->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists23->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists23->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists23->Draw();

    yBins = slowfasttimediff_energy_hists23->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists23->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana23Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana23s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr23 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr23 -> SetTitle(Form("Banana23 Fit L%d", i));
    //gr23 -> SetStats(0);
	gr23 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr23 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr23 -> SetMarkerStyle(20);
    gr23 -> SetMarkerSize(0.8);
    gr23 -> SetMarkerColor(kBlue);
    gr23->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit23 = new TF1("fit23", "pol0(0)", 0, 8000);
    gr23->Draw("AL*");
    gr23->Fit("fit23", "");
    fit23->Draw("same");
    cb -> SaveAs(Form("%s/banana23Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr23;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana23 fit for detector %d...\n", i);
    return fit23;
}

TF1* banana24(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists24, TH2D* energy_slowfasttimediff_hists24, std::string dirPathplots)
{
    std::cout <<"Fitting banana24 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana24 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists24->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists24->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists24->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists24->Draw();

    yBins = slowfasttimediff_energy_hists24->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists24->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana24Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana24s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr24 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr24 -> SetTitle(Form("Banana24 Fit L%d", i));
    //gr24 -> SetStats(0);
	gr24 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr24 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr24 -> SetMarkerStyle(20);
    gr24 -> SetMarkerSize(0.8);
    gr24 -> SetMarkerColor(kBlue);
    gr24->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit24 = new TF1("fit24", "pol0(0)", 0, 8000);
    gr24->Draw("AL*");
    gr24->Fit("fit24", "");
    fit24->Draw("same");
    cb -> SaveAs(Form("%s/banana24Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr24;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana24 fit for detector %d...\n", i);
    return fit24;
}

TF1* banana25(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists25, TH2D* energy_slowfasttimediff_hists25, std::string dirPathplots)
{
    std::cout <<"Fitting banana25 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana25 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists25->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists25->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists25->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists25->Draw();

    yBins = slowfasttimediff_energy_hists25->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists25->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana25Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana25s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr25 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr25 -> SetTitle(Form("Banana25 Fit L%d", i));
    //gr25 -> SetStats(0);
	gr25 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr25 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr25 -> SetMarkerStyle(20);
    gr25 -> SetMarkerSize(0.8);
    gr25 -> SetMarkerColor(kBlue);
    gr25->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit25 = new TF1("fit25", "pol0(0)", 0, 8000);
    gr25->Draw("AL*");
    gr25->Fit("fit25", "");
    fit25->Draw("same");
    cb -> SaveAs(Form("%s/banana25Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr25;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana25 fit for detector %d...\n", i);
    return fit25;
}

TF1* banana26(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists26, TH2D* energy_slowfasttimediff_hists26, std::string dirPathplots)
{
    std::cout <<"Fitting banana26 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana26 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists26->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists26->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists26->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists26->Draw();

    yBins = slowfasttimediff_energy_hists26->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists26->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana26Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana26s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr26 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr26 -> SetTitle(Form("Banana26 Fit L%d", i));
    //gr26 -> SetStats(0);
	gr26 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr26 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr26 -> SetMarkerStyle(20);
    gr26 -> SetMarkerSize(0.8);
    gr26 -> SetMarkerColor(kBlue);
    gr26->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit26 = new TF1("fit26", "pol0(0)", 0, 8000);
    gr26->Draw("AL*");
    gr26->Fit("fit26", "");
    fit26->Draw("same");
    cb -> SaveAs(Form("%s/banana26Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr26;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana26 fit for detector %d...\n", i);
    return fit26;
}

TF1* banana27(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists27, TH2D* energy_slowfasttimediff_hists27, std::string dirPathplots)
{
    std::cout <<"Fitting banana27 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana27 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists27->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists27->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists27->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists27->Draw();

    yBins = slowfasttimediff_energy_hists27->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists27->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana27Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana27s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr27 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr27 -> SetTitle(Form("Banana27 Fit L%d", i));
    //gr27 -> SetStats(0);
	gr27 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr27 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr27 -> SetMarkerStyle(20);
    gr27 -> SetMarkerSize(0.8);
    gr27 -> SetMarkerColor(kBlue);
    gr27->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit27 = new TF1("fit27", "pol0(0)", 0, 8000);
    gr27->Draw("AL*");
    gr27->Fit("fit27", "");
    fit27->Draw("same");
    cb -> SaveAs(Form("%s/banana27Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr27;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana27 fit for detector %d...\n", i);
    return fit27;
}

TF1* banana28(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists28, TH2D* energy_slowfasttimediff_hists28, std::string dirPathplots)
{
    std::cout <<"Fitting banana28 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana28 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists28->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists28->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists28->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists28->Draw();

    yBins = slowfasttimediff_energy_hists28->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists28->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana28Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana28s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr28 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr28 -> SetTitle(Form("Banana28 Fit L%d", i));
    //gr28 -> SetStats(0);
	gr28 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr28 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr28 -> SetMarkerStyle(20);
    gr28 -> SetMarkerSize(0.8);
    gr28 -> SetMarkerColor(kBlue);
    gr28->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit28 = new TF1("fit28", "pol0(0)", 0, 8000);
    gr28->Draw("AL*");
    gr28->Fit("fit28", "");
    fit28->Draw("same");
    cb -> SaveAs(Form("%s/banana28Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr28;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana28 fit for detector %d...\n", i);
    return fit28;
}

TF1* banana29(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists29, TH2D* energy_slowfasttimediff_hists29, std::string dirPathplots)
{
    std::cout <<"Fitting banana29 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana29 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists29->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists29->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists29->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists29->Draw();

    yBins = slowfasttimediff_energy_hists29->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists29->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana29Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana29s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr29 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr29 -> SetTitle(Form("Banana29 Fit L%d", i));
    //gr29 -> SetStats(0);
	gr29 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr29 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr29 -> SetMarkerStyle(20);
    gr29 -> SetMarkerSize(0.8);
    gr29 -> SetMarkerColor(kBlue);
    gr29->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit29 = new TF1("fit29", "pol0(0)", 0, 8000);
    gr29->Draw("AL*");
    gr29->Fit("fit29", "");
    fit29->Draw("same");
    cb -> SaveAs(Form("%s/banana29Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr29;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana29 fit for detector %d...\n", i);
    return fit29;
}

TF1* banana30(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists30, TH2D* energy_slowfasttimediff_hists30, std::string dirPathplots)
{
    std::cout <<"Fitting banana30 (0-8000 keV) for detector \n" << i << std::endl;
    std::cout <<"Fitting banana30 time range (325-335 keV) "<< std::endl;

    int yBins = slowfasttimediff_energy_hists30->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists30->GetXaxis()->GetNbins();

    int slicewidth = 80, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists30->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists30->Draw();

    yBins = slowfasttimediff_energy_hists30->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin <= 500) slicewidth = 50;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists30->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);
        
        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");
        
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back((slicewidth)/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/banana30Slices_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("Fitting banana30s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr30 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr30 -> SetTitle(Form("Banana30 Fit L%d", i));
    //gr30 -> SetStats(0);
	gr30 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr30 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr30 -> SetMarkerStyle(20);
    gr30 -> SetMarkerSize(0.8);
    gr30 -> SetMarkerColor(kBlue);
    gr30->GetXaxis()->SetRangeUser(0, 8000);
    TF1* fit30 = new TF1("fit30", "pol0(0)", 0, 8000);
    gr30->Draw("AL*");
    gr30->Fit("fit30", "");
    fit30->Draw("same");
    cb -> SaveAs(Form("%s/banana30Fit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr30;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();
    printf("Returning banana30 fit for detector %d...\n", i);
    return fit30;
}
void runduration(TTree* tree, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds,unsigned long long numEntries)
{
    for (int i = 0; i < 1; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
    tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
    tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

    Double_t firstTime = 0;
    Double_t lastTime = 0;
    for (unsigned long long entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);
        if (timeFast[0] > 0)
        {
            firstTime = timeFast[0];
            break;
        }
    }

    // Find the last non-zero entry in the timeFL0 branch of the tree
    for (unsigned long long entry = numEntries - 1; entry >= 0; entry--) 
    {
        tree->GetEntry(entry);
        if (timeFast[0] > 0)
        {
            lastTime = timeFast[0];
            break;
        }
    }

    Double_t timeDuration = (lastTime - firstTime)*1e-9/60;

    if (timeDuration/10 > 3) // if its more than 90 minutes (our runs are 10 and 20 min so 10/10 is 1, etc.)
    { 
        printf("\033[1;31mWARNING: The time is not in nanoseconds. Must divide time in tree if not already done. \033[0m\n");
        if (timeDuration/100 > 3 && timeDuration/100 < 33) 
        {
            hundredsOfNanoSeconds = true;
            printf("\033[1;31mWARNING: The time is in the order of hundreds of nanoseconds. \033[0m\n");
        }
        else if (timeDuration/10 > 3 && timeDuration/10 < 33) 
        {
            tensOfNanoSeconds = true;
            printf("\033[1;31mWARNING: The time is in the order of tens of nanoseconds. \033[0m\n");
        }
        
    }
    printf("\n\nTime duration of the run = %f minutes\n", timeDuration);
}

std::vector<Double_t> fasttime_calibration_constants(bool beam, TFile *outputFile ,std::string dirPathplots, TDirectory *dir_fasttimecalib, TTree* tree, bool tensOfNanoSeconds, bool hundredsOfNanoSeconds, unsigned long long numEntries)
{
    if (beam)
    {
        printf("\n\nBeam run: Calculating inter-detector fast time calibration constants (used in next function)...\n");
        // Set branch addresses for each detector
        for (int i = 0; i < 1; i++) 
        {
            tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
            tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
            tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
        }
        tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
        tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
        tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

        // Initialise the histograms
        std::vector<TH1D*> fasttimecalib_hists(NumDetectors-1);
        for (int i = 1; i < 1; i++) 
        {
            fasttimecalib_hists[i] = new TH1D(Form("fasttimecalib_hists%d", i), Form("Fast time calibration L0 - L%d", i), 1000 , 0, 1000);
        }

        for (unsigned long long entry = 0; entry < numEntries; entry++) 
        {
            tree->GetEntry(entry);

            Double_t timeFastDet0 = 0;  // Initialize
            if (timeFast[0] == 0) continue;

            if (tensOfNanoSeconds) timeFastDet0 = timeFast[0] / 10;
            else if (hundredsOfNanoSeconds) timeFastDet0 = timeFast[0] / 100;
            else timeFastDet0 = timeFast[0];

            for (int i = 1; i < 1; i++) 
            {
                tree->GetEntry(entry + 1); // Get the next entry for timeFastDetX

                if (timeFast[i] != 0) 
                {
                    Double_t timeFastDetX = 0;  // Initialize
                    if (tensOfNanoSeconds) timeFastDetX = timeFast[i] / 10;
                    else if (hundredsOfNanoSeconds) timeFastDetX = timeFast[i] / 100;
                    else timeFastDetX = timeFast[i];

                    Double_t deltaTime = (timeFastDetX - timeFastDet0);

                    //printf("timeFastL0: %f, timeFastL%d: %f\n", timeFastDet0, i, timeFastDetX);
                    //printf("Detector %d - 0: Fast time diff = %f\n", i, deltaTime);

                    fasttimecalib_hists[i]->Fill(deltaTime);
                }
            }
        }

        // Create a separate fit function for each histogram
        std::vector<TF1*> fasttimecalib_fit(NumDetectors-1);
        for (int i = 1; i < 1; i++) 
        {
            fasttimecalib_fit[i] = new TF1(Form("fasttimecalib_fit%d", i), "gaus(0)", 260, 360);
        }

        // Fit the histograms and extract parameters
        for (int i = 1; i < 1; i++) 
        {
            fasttimecalib_hists[i]->GetXaxis()->SetRangeUser(260, 360);
            fasttimecalib_hists[i]->Fit(Form("fasttimecalib_fit%d", i), "Q");
            fasttimecalib_hists[i]->GetXaxis()->SetRangeUser(0, 1000);

            mean_fasttimediff[i] = fasttimecalib_fit[i]->GetParameter(1);
            sigma_fasttimediff[i] = fasttimecalib_fit[i]->GetParameter(2);
            printf("Detector %d: Fast time diff mean = %f\n", i, mean_fasttimediff[i]);
            printf("Detector %d: Fast time diff sigma = %f\n", i, sigma_fasttimediff[i]);
        }

        // Create a canvas and plot the histograms with fits
        TCanvas* c1 = new TCanvas("c1", "Fast time calibration", 800, 600);
        c1->Divide(2, 2);
        for (int i = 1; i < 1; i++) 
        {
            c1->cd(i);
            fasttimecalib_hists[i]->SetTitle(Form("Fast time calibration L0 - L%d", i));
            fasttimecalib_hists[i]->SetStats(0);
            fasttimecalib_hists[i]->GetXaxis()->SetTitle("Fast time difference (ns)");
            fasttimecalib_hists[i]->GetYaxis()->SetTitle("Counts/ns");
            fasttimecalib_hists[i]->Draw();
            fasttimecalib_fit[i]->Draw("same");
        }

        // Save canvas and clean up
        c1->Write();
        c1->SaveAs(Form("%s/fasttimecalib.root", dirPathplots.c_str()));
        delete c1; 
        
        mean_fasttimediff[0] = 0;
    }
    else 
    {
        printf("Source run: not calculating inter-detector fast time calibration constants ...\n");
        for (int i = 1; i < 1; i++) 
        {
            mean_fasttimediff[i] = 0;
            sigma_fasttimediff[i] = 0;
        }
    }
    
    return mean_fasttimediff;
}

void align_offset(TFile *outputFile, std::string dirPathplots, TTree* tree, TDirectory *dir_banana, TDirectory *dir_aligned, std::vector<Double_t> mean_fasttimediff, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds, unsigned long long numEntries) {
    
    dir_aligned->cd();
    printf("\nCloning the TTree...\n");
    // Set branch addresses for each detector
    for (int i = 0; i < 1; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
    tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
    tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

    TTree* copyTree = new TTree("LaBrDataCopy", "LaBrDataCopy");
    for (int i = 0; i < 1; i++) 
    {
        copyTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        copyTree->Branch(Form("timeSL%d", i), &timeSlowCp[i], Form("timeSL%d/D", i));
        copyTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
    }
    copyTree->Branch("timeRF", &timeFastCp[RFIndex], "timeRF/D");
    copyTree->Branch("fastEPOLARIS", &fastEnergyPOLARISCp[0], "fastEPOLARIS/D");
    copyTree->Branch("fastTPOLARIS", &fastTimePOLARISCp[0], "fastTPOLARIS/D");


    for (unsigned long long entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);
        for (int i = 0; i < 1; i++) 
        {
            if (energySlow[i]< 14000) 
            {
                if (timeSlow[i] > 0) 
                {
                    energySlowCp[i] = energySlow[i];
                    timeOffsetBeforeCorrection[i] = 0;
                    if (tensOfNanoSeconds) timeSlowCp[i] = timeSlow[i]/10;
                    else if (hundredsOfNanoSeconds) timeSlowCp[i] = timeSlow[i]/100;
                    else timeSlowCp[i] = timeSlow[i];
                }
            }
            else 
            {
                energySlowCp[i] = 0;
                timeSlowCp[i] = 0;
            }
            // printf("energySlowCp[%d] = %f, timeSlowCp[%d] = %f, timeFastCp[%d] = %f\n", i, energySlowCp[i], i, timeSlowCp[i], i, timeFastCp[i]);
        }
        if (tensOfNanoSeconds) 
        {
            if (timeFast[0] > 0) timeFastCp[0] = timeFast[0]/10;
            if (timeFast[1] > 0) timeFastCp[1] = timeFast[1]/10 - mean_fasttimediff[1];
            if (timeFast[2] > 0) timeFastCp[2] = timeFast[2]/10 - mean_fasttimediff[2];
            if (timeFast[3] > 0) timeFastCp[3] = timeFast[3]/10 - mean_fasttimediff[3];
            timeFastCp[RFIndex] = timeFast[RFIndex]/10;
            fastTimePOLARISCp[0] = fastTimePOLARIS[0]/10;
            fastEnergyPOLARISCp[0] = fastEnergyPOLARIS[0];
        }
        else if (hundredsOfNanoSeconds) 
        {
            if (timeFast[0] > 0) timeFastCp[0] = timeFast[0]/100;
            if (timeFast[1] > 0) timeFastCp[1] = timeFast[1]/100 - mean_fasttimediff[1];
            if (timeFast[2] > 0) timeFastCp[2] = timeFast[2]/100 - mean_fasttimediff[2];
            if (timeFast[3] > 0) timeFastCp[3] = timeFast[3]/100 - mean_fasttimediff[3];
            timeFastCp[RFIndex] = timeFast[RFIndex]/100;
            fastTimePOLARISCp[0] = fastTimePOLARIS[0]/100;
            fastEnergyPOLARISCp[0] = fastEnergyPOLARIS[0];
        }
        else 
        {
            if (timeFast[0] > 0) timeFastCp[0] = timeFast[0];
            if (timeFast[1] > 0) timeFastCp[1] = timeFast[1] - mean_fasttimediff[1];
            if (timeFast[2] > 0) timeFastCp[2] = timeFast[2] - mean_fasttimediff[2];
            if (timeFast[3] > 0) timeFastCp[3] = timeFast[3] - mean_fasttimediff[3];
            timeFastCp[RFIndex] = timeFast[RFIndex];
            fastTimePOLARISCp[0] = fastTimePOLARIS[0];
            fastEnergyPOLARISCp[0] = fastEnergyPOLARIS[0];
        }
        copyTree->Fill();
    }


    // Initialise the histograms
    for (int i = 0; i < 1; i++) 
    {
        slowfasttimediff_hists[i] = new TH1D(Form("slowfasttimediff_hists_%d", i), Form("Uncorrected time offset: L%d", i),  600, 0, 600);

        // Full
        slowfasttimediff_energy_hists[i] = new TH2D(Form("slowfasttimediff_energy_hists_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  600, 0, 600, 8000, 0, 8000);
        energy_slowfasttimediff_hists[i] = new TH2D(Form("energy_slowfasttimediff_hists_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 600, 0, 600);

        // banana 1
        slowfasttimediff_energy_hists1[i] = new TH2D(Form("slowfasttimediff_energy_hists1_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 235, 245, 8000, 0, 8000);
        energy_slowfasttimediff_hists1[i] = new TH2D(Form("energy_slowfasttimediff_hists1_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 235, 245);

        // banana 2
        slowfasttimediff_energy_hists2[i] = new TH2D(Form("slowfasttimediff_energy_hists2_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 245, 255, 8000, 0, 8000);
        energy_slowfasttimediff_hists2[i] = new TH2D(Form("energy_slowfasttimediff_hists2_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 245, 255);

        // banana 3
        slowfasttimediff_energy_hists3[i] = new TH2D(Form("slowfasttimediff_energy_hists3_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 255, 265, 8000, 0, 8000);
        energy_slowfasttimediff_hists3[i] = new TH2D(Form("energy_slowfasttimediff_hists3_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 255, 265);

        // banana 4
        slowfasttimediff_energy_hists4[i] = new TH2D(Form("slowfasttimediff_energy_hists4_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 265, 275, 8000, 0, 8000);
        energy_slowfasttimediff_hists4[i] = new TH2D(Form("energy_slowfasttimediff_hists4_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 265, 275);

        // banana 5
        slowfasttimediff_energy_hists5[i] = new TH2D(Form("slowfasttimediff_energy_hists5_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 275, 285, 8000, 0, 8000);
        energy_slowfasttimediff_hists5[i] = new TH2D(Form("energy_slowfasttimediff_hists5_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 275, 285);

        // banana 6
        slowfasttimediff_energy_hists6[i] = new TH2D(Form("slowfasttimediff_energy_hists6_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 285, 295, 8000, 0, 8000);
        energy_slowfasttimediff_hists6[i] = new TH2D(Form("energy_slowfasttimediff_hists6_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 285, 295);
        
        // banana 7
        slowfasttimediff_energy_hists7[i] = new TH2D(Form("slowfasttimediff_energy_hists7_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 295, 305, 8000, 0, 8000);
        energy_slowfasttimediff_hists7[i] = new TH2D(Form("energy_slowfasttimediff_hists7_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 295, 305);

        // banana 8
        slowfasttimediff_energy_hists8[i] = new TH2D(Form("slowfasttimediff_energy_hists8_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 305, 315, 8000, 0, 8000);
        energy_slowfasttimediff_hists8[i] = new TH2D(Form("energy_slowfasttimediff_hists8_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 305, 315);

        // banana 9
        slowfasttimediff_energy_hists9[i] = new TH2D(Form("slowfasttimediff_energy_hists9_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 315, 325, 8000, 0, 8000);
        energy_slowfasttimediff_hists9[i] = new TH2D(Form("energy_slowfasttimediff_hists9_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 315, 325);

        // banana 10
        slowfasttimediff_energy_hists10[i] = new TH2D(Form("slowfasttimediff_energy_hists10_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 325, 335, 8000, 0, 8000);
        energy_slowfasttimediff_hists10[i] = new TH2D(Form("energy_slowfasttimediff_hists10_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 325, 335);

        // banana 11
        slowfasttimediff_energy_hists11[i] = new TH2D(Form("slowfasttimediff_energy_hists11_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 335, 345, 8000, 0, 8000);
        energy_slowfasttimediff_hists11[i] = new TH2D(Form("energy_slowfasttimediff_hists11_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 335, 345);

        // banana 12
        slowfasttimediff_energy_hists12[i] = new TH2D(Form("slowfasttimediff_energy_hists12_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 345, 355, 8000, 0, 8000);
        energy_slowfasttimediff_hists12[i] = new TH2D(Form("energy_slowfasttimediff_hists12_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 345, 355);

        // banana 13
        slowfasttimediff_energy_hists13[i] = new TH2D(Form("slowfasttimediff_energy_hists13_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 355, 365, 8000, 0, 8000);
        energy_slowfasttimediff_hists13[i] = new TH2D(Form("energy_slowfasttimediff_hists13_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 355, 365);

        // banana 14
        slowfasttimediff_energy_hists14[i] = new TH2D(Form("slowfasttimediff_energy_hists14_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 365, 375, 8000, 0, 8000);
        energy_slowfasttimediff_hists14[i] = new TH2D(Form("energy_slowfasttimediff_hists14_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 365, 375);

        // banana 15
        slowfasttimediff_energy_hists15[i] = new TH2D(Form("slowfasttimediff_energy_hists15_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 375, 385, 8000, 0, 8000);
        energy_slowfasttimediff_hists15[i] = new TH2D(Form("energy_slowfasttimediff_hists15_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 375, 385);

        // banana 16
        slowfasttimediff_energy_hists16[i] = new TH2D(Form("slowfasttimediff_energy_hists16_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 385, 395, 8000, 0, 8000);
        energy_slowfasttimediff_hists16[i] = new TH2D(Form("energy_slowfasttimediff_hists16_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 385, 395);

        // banana 17
        slowfasttimediff_energy_hists17[i] = new TH2D(Form("slowfasttimediff_energy_hists17_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 395, 405, 8000, 0, 8000);
        energy_slowfasttimediff_hists17[i] = new TH2D(Form("energy_slowfasttimediff_hists17_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 395, 405);

        // banana 18
        slowfasttimediff_energy_hists18[i] = new TH2D(Form("slowfasttimediff_energy_hists18_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 405, 415, 8000, 0, 8000);
        energy_slowfasttimediff_hists18[i] = new TH2D(Form("energy_slowfasttimediff_hists18_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 405, 415);

        // banana 19
        slowfasttimediff_energy_hists19[i] = new TH2D(Form("slowfasttimediff_energy_hists19_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 415, 425, 8000, 0, 8000);
        energy_slowfasttimediff_hists19[i] = new TH2D(Form("energy_slowfasttimediff_hists19_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 415, 425);

        // banana 20
        slowfasttimediff_energy_hists20[i] = new TH2D(Form("slowfasttimediff_energy_hists20_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 425, 435, 8000, 0, 8000);
        energy_slowfasttimediff_hists20[i] = new TH2D(Form("energy_slowfasttimediff_hists20_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 425, 435);

        // banana 21
        slowfasttimediff_energy_hists21[i] = new TH2D(Form("slowfasttimediff_energy_hists21_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 435, 445, 8000, 0, 8000);
        energy_slowfasttimediff_hists21[i] = new TH2D(Form("energy_slowfasttimediff_hists21_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 435, 445);

        // banana 22
        slowfasttimediff_energy_hists22[i] = new TH2D(Form("slowfasttimediff_energy_hists22_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 445, 455, 8000, 0, 8000);
        energy_slowfasttimediff_hists22[i] = new TH2D(Form("energy_slowfasttimediff_hists22_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 445, 455);

        // banana 23
        slowfasttimediff_energy_hists23[i] = new TH2D(Form("slowfasttimediff_energy_hists23_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 455, 465, 8000, 0, 8000);
        energy_slowfasttimediff_hists23[i] = new TH2D(Form("energy_slowfasttimediff_hists23_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 455, 465);

        // banana 24
        slowfasttimediff_energy_hists24[i] = new TH2D(Form("slowfasttimediff_energy_hists24_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 465, 475, 8000, 0, 8000);
        energy_slowfasttimediff_hists24[i] = new TH2D(Form("energy_slowfasttimediff_hists24_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 465, 475);

        // banana 25
        slowfasttimediff_energy_hists25[i] = new TH2D(Form("slowfasttimediff_energy_hists25_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 475, 485, 8000, 0, 8000);
        energy_slowfasttimediff_hists25[i] = new TH2D(Form("energy_slowfasttimediff_hists25_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 475, 485);

        // banana 26
        slowfasttimediff_energy_hists26[i] = new TH2D(Form("slowfasttimediff_energy_hists26_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 485, 495, 8000, 0, 8000);
        energy_slowfasttimediff_hists26[i] = new TH2D(Form("energy_slowfasttimediff_hists26_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 485, 495);

        // banana 27
        slowfasttimediff_energy_hists27[i] = new TH2D(Form("slowfasttimediff_energy_hists27_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 495, 505, 8000, 0, 8000);
        energy_slowfasttimediff_hists27[i] = new TH2D(Form("energy_slowfasttimediff_hists27_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 495, 505);

        // banana 28
        slowfasttimediff_energy_hists28[i] = new TH2D(Form("slowfasttimediff_energy_hists28_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 505, 515, 8000, 0, 8000);
        energy_slowfasttimediff_hists28[i] = new TH2D(Form("energy_slowfasttimediff_hists28_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 505, 515);

        // banana 29
        slowfasttimediff_energy_hists29[i] = new TH2D(Form("slowfasttimediff_energy_hists29_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 515, 525, 8000, 0, 8000);
        energy_slowfasttimediff_hists29[i] = new TH2D(Form("energy_slowfasttimediff_hists29_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 515, 525);

        // banana 30
        slowfasttimediff_energy_hists30[i] = new TH2D(Form("slowfasttimediff_energy_hists30_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  10, 525, 535, 8000, 0, 8000);
        energy_slowfasttimediff_hists30[i] = new TH2D(Form("energy_slowfasttimediff_hists30_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 10, 525, 535);

        timeOffsetCorrectedHistsNext[i] = new TH1D(Form("timeOffsetCorrectedHistsNext%d", i), Form("Corrected time offset - Next Event: L%d", i),  800, -400, 400);

        timeEnergyOffsetCorrectedHistsNext[i] = new TH2D(Form("timeEnergyOffsetCorrectedHistsNext%d", i), Form("Corrected time offset vs energy - Next Event: L%d", i),  800, -400, 400, 8000, 0, 8000);
        energyhist[i] = new TH1D(Form("energyhist%d", i), Form(" energy for corrected time: L%d", i),  8000, 0, 8000);

    }
    
    //________________________________________________________________________________________________________________________
    // Initialize new TTree
    TTree* offsetTree = new TTree("LaBrDataOffset", "LaBrDataOffset"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < 1; i++) 
    {
        offsetTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        offsetTree->Branch(Form("timeSL%d", i), &timeSlowCp[i], Form("timeSL%d/D", i));
        offsetTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        offsetTree->Branch(Form("timeOffsetL%d", i), &timeOffsetBeforeCorrection[i], Form("timeOffsetL%d/D", i)); //new Branch
    }
    offsetTree->Branch("timeRF", &timeFastCp[RFIndex], "timeRF/D");
    offsetTree->Branch("fastEPOLARIS", &fastEnergyPOLARISCp[0], "fastEPOLARIS/D");
    offsetTree->Branch("fastTPOLARIS", &fastTimePOLARISCp[0], "fastTPOLARIS/D");

    // Initialize new TTree
    TTree* correctedTree = new TTree("LaBrDataCorrected", "LaBrDataCorrected"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < 1; i++) 
    {
        correctedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrected[i], Form("timeSCorrectedL%d/D", i)); //new Branch
        correctedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
    }
    correctedTree->Branch("timeRF", &timeFastCp[RFIndex], "timeRF/D");
    correctedTree->Branch("fastEPOLARIS", &fastEnergyPOLARISCp[0], "fastEPOLARIS/D");
    correctedTree->Branch("fastTPOLARIS", &fastTimePOLARISCp[0], "fastTPOLARIS/D");

    // _________________________________________________________________________________________________________________________________________________________
    for (int i = 0; i < 1; i++) 
    {
        printf("\n\nCalculating time offset for detector %d...\n", i);
        for (unsigned long long entry = 0; entry < numEntries; entry++) 
        {
            copyTree->GetEntry(entry);
            //printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCp[i], i, energySlowCp[i]);              

            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;
            Double_t deltaTime = 0;
            Double_t energy = 0;

            for (unsigned long long nextEntry = entry + 1; nextEntry < numEntries; nextEntry++) 
            {
                copyTree->GetEntry(nextEntry);  
                // printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCp[i], i, energySlowCp[i]);       
                       
                if (timeSlowCp[i] == 0) continue;
                nextTimeSlow = timeSlowCp[i];

                deltaTime = (nextTimeSlow - refTimeFast); 
                // printf("Detector %d: Time offset = %f\n", i, deltaTime);
                energy = energySlowCp[i];
                if ((deltaTime > 600.0) || (deltaTime < 0.0)) break;
                else
                {
                    timeOffsetBeforeCorrection[i] = static_cast<Double_t>(static_cast<int>(deltaTime));
                    // printf("Detector %d: Time offset = %f\n", i, timeOffsetBeforeCorrection[i]);
                    // printf("Detector %d: deltaTime = %f\n", i, deltaTime);
                    
                    slowfasttimediff_hists[i]->Fill(deltaTime);
                    slowfasttimediff_energy_hists[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists1[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists1[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists2[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists2[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists3[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists3[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists4[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists4[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists5[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists5[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists6[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists6[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists7[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists7[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists8[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists8[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists9[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists9[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists10[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists10[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists11[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists11[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists12[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists12[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists13[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists13[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists14[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists14[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists15[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists15[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists16[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists16[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists17[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists17[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists18[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists18[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists19[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists19[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists20[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists20[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists21[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists21[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists22[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists22[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists23[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists23[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists24[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists24[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists25[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists25[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists26[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists26[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists27[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists27[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists28[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists28[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists29[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists29[i]->Fill(energy, deltaTime);
                    slowfasttimediff_energy_hists30[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists30[i]->Fill(energy, deltaTime);
                    
                    entry = nextEntry;
                    break;
                }

            }
            offsetTree->Fill();
        }
        
        TCanvas* c1 = new TCanvas(Form("c1_%d", i), Form("Uncorrected time offset: L%d", i), 800, 600);
        slowfasttimediff_hists[i]->GetXaxis()->SetRangeUser(0, 600);
        slowfasttimediff_hists[i]->SetFillColor(kBlue);
        slowfasttimediff_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_hists[i]->GetYaxis()->SetTitle("Counts/ns");
        slowfasttimediff_hists[i]->Draw();
        c1->SaveAs(Form("%s/uncorrected_time_offset_hist_L%d.root", dirPathplots.c_str(), i));
        c1->Write(Form("uncorrected_time_offset_hist_L%d", i));
        delete c1;

        TCanvas* c12 = new TCanvas(Form("c12_%d", i), Form("Uncorrected time offset vs energy: L%d", i), 800, 600);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetRangeUser(0, 600);
        slowfasttimediff_energy_hists[i]->GetYaxis()->SetRangeUser(0, 8000);
        slowfasttimediff_energy_hists[i]->SetStats(0);
        slowfasttimediff_energy_hists[i]->SetMaximum(500);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_energy_hists[i]->GetYaxis()->SetTitle("Energy (keV)");
        slowfasttimediff_energy_hists[i]->Draw("colz");
        c12->SaveAs(Form("%s/uncorrected_time_offset_energy_L%d.root", dirPathplots.c_str(), i));
        c12->Write(Form("uncorrected_time_offset_energy_L%d", i));
        delete c12;

        energy_slowfasttimediff_hists[i]->GetYaxis()->SetRangeUser(0, 600);

        dir_banana->cd();
        banana1fit[i] = banana1(outputFile, i, slowfasttimediff_energy_hists1[i], energy_slowfasttimediff_hists1[i], dirPathplots);
        banana2fit[i] = banana2(outputFile, i, slowfasttimediff_energy_hists2[i], energy_slowfasttimediff_hists2[i], dirPathplots);
        banana3fit[i] = banana3(outputFile, i, slowfasttimediff_energy_hists3[i], energy_slowfasttimediff_hists3[i], dirPathplots);
        banana4fit[i] = banana4(outputFile, i, slowfasttimediff_energy_hists4[i], energy_slowfasttimediff_hists4[i], dirPathplots);
        banana5fit[i] = banana5(outputFile, i, slowfasttimediff_energy_hists5[i], energy_slowfasttimediff_hists5[i], dirPathplots);
        banana6fit[i] = banana6(outputFile, i, slowfasttimediff_energy_hists6[i], energy_slowfasttimediff_hists6[i], dirPathplots);
        banana7fit[i] = banana7(outputFile, i, slowfasttimediff_energy_hists7[i], energy_slowfasttimediff_hists7[i], dirPathplots);
        banana8fit[i] = banana8(outputFile, i, slowfasttimediff_energy_hists8[i], energy_slowfasttimediff_hists8[i], dirPathplots);
        banana9fit[i] = banana9(outputFile, i, slowfasttimediff_energy_hists9[i], energy_slowfasttimediff_hists9[i], dirPathplots);
        banana10fit[i] = banana10(outputFile, i, slowfasttimediff_energy_hists10[i], energy_slowfasttimediff_hists10[i], dirPathplots);
        banana11fit[i] = banana11(outputFile, i, slowfasttimediff_energy_hists11[i], energy_slowfasttimediff_hists11[i], dirPathplots);
        banana12fit[i] = banana12(outputFile, i, slowfasttimediff_energy_hists12[i], energy_slowfasttimediff_hists12[i], dirPathplots);
        banana13fit[i] = banana13(outputFile, i, slowfasttimediff_energy_hists13[i], energy_slowfasttimediff_hists13[i], dirPathplots);
        banana14fit[i] = banana14(outputFile, i, slowfasttimediff_energy_hists14[i], energy_slowfasttimediff_hists14[i], dirPathplots);
        banana15fit[i] = banana15(outputFile, i, slowfasttimediff_energy_hists15[i], energy_slowfasttimediff_hists15[i], dirPathplots);
        banana16fit[i] = banana16(outputFile, i, slowfasttimediff_energy_hists16[i], energy_slowfasttimediff_hists16[i], dirPathplots);
        banana17fit[i] = banana17(outputFile, i, slowfasttimediff_energy_hists17[i], energy_slowfasttimediff_hists17[i], dirPathplots);
        banana18fit[i] = banana18(outputFile, i, slowfasttimediff_energy_hists18[i], energy_slowfasttimediff_hists18[i], dirPathplots);
        banana19fit[i] = banana19(outputFile, i, slowfasttimediff_energy_hists19[i], energy_slowfasttimediff_hists19[i], dirPathplots);
        banana20fit[i] = banana20(outputFile, i, slowfasttimediff_energy_hists20[i], energy_slowfasttimediff_hists20[i], dirPathplots);
        banana21fit[i] = banana21(outputFile, i, slowfasttimediff_energy_hists21[i], energy_slowfasttimediff_hists21[i], dirPathplots);
        banana22fit[i] = banana22(outputFile, i, slowfasttimediff_energy_hists22[i], energy_slowfasttimediff_hists22[i], dirPathplots);
        banana23fit[i] = banana23(outputFile, i, slowfasttimediff_energy_hists23[i], energy_slowfasttimediff_hists23[i], dirPathplots);
        banana24fit[i] = banana24(outputFile, i, slowfasttimediff_energy_hists24[i], energy_slowfasttimediff_hists24[i], dirPathplots);
        banana25fit[i] = banana25(outputFile, i, slowfasttimediff_energy_hists25[i], energy_slowfasttimediff_hists25[i], dirPathplots);
        banana26fit[i] = banana26(outputFile, i, slowfasttimediff_energy_hists26[i], energy_slowfasttimediff_hists26[i], dirPathplots);
        banana27fit[i] = banana27(outputFile, i, slowfasttimediff_energy_hists27[i], energy_slowfasttimediff_hists27[i], dirPathplots);
        banana28fit[i] = banana28(outputFile, i, slowfasttimediff_energy_hists28[i], energy_slowfasttimediff_hists28[i], dirPathplots);
        banana29fit[i] = banana29(outputFile, i, slowfasttimediff_energy_hists29[i], energy_slowfasttimediff_hists29[i], dirPathplots);
        banana30fit[i] = banana30(outputFile, i, slowfasttimediff_energy_hists30[i], energy_slowfasttimediff_hists30[i], dirPathplots);
        
        TCanvas * c2 = new TCanvas("c2", "banana", 800, 600);
        energy_slowfasttimediff_hists[i]->Draw();
        banana1fit[i]->SetLineColor(kRed);
        banana2fit[i]->SetLineColor(kGreen);
        banana3fit[i]->SetLineColor(kBlue);
        banana4fit[i]->SetLineColor(kYellow);
        banana5fit[i]->SetLineColor(kMagenta);
        banana6fit[i]->SetLineColor(kCyan);
        banana7fit[i]->SetLineColor(kOrange);
        banana8fit[i]->SetLineColor(kWhite);
        banana9fit[i]->SetLineColor(kRed);
        banana10fit[i]->SetLineColor(kGreen);
        banana11fit[i]->SetLineColor(kBlue);
        banana12fit[i]->SetLineColor(kYellow);
        banana13fit[i]->SetLineColor(kMagenta);
        banana14fit[i]->SetLineColor(kCyan);
        banana15fit[i]->SetLineColor(kOrange);
        banana16fit[i]->SetLineColor(kWhite);
        banana17fit[i]->SetLineColor(kRed);
        banana18fit[i]->SetLineColor(kGreen);
        banana19fit[i]->SetLineColor(kBlue);
        banana20fit[i]->SetLineColor(kYellow);
        banana21fit[i]->SetLineColor(kMagenta);
        banana22fit[i]->SetLineColor(kCyan);
        banana23fit[i]->SetLineColor(kOrange);
        banana24fit[i]->SetLineColor(kWhite);
        banana25fit[i]->SetLineColor(kRed);
        banana26fit[i]->SetLineColor(kGreen);
        banana27fit[i]->SetLineColor(kBlue);
        banana28fit[i]->SetLineColor(kYellow);
        banana29fit[i]->SetLineColor(kMagenta);
        banana30fit[i]->SetLineColor(kCyan);
        banana1fit[i]->Draw("same");
        banana2fit[i]->Draw("same");
        banana3fit[i]->Draw("same");
        banana4fit[i]->Draw("same");
        banana5fit[i]->Draw("same");
        banana6fit[i]->Draw("same");
        banana7fit[i]->Draw("same");
        banana8fit[i]->Draw("same");
        banana9fit[i]->Draw("same");
        banana10fit[i]->Draw("same");
        banana11fit[i]->Draw("same");
        banana12fit[i]->Draw("same");
        banana13fit[i]->Draw("same");
        banana14fit[i]->Draw("same");
        banana15fit[i]->Draw("same");
        banana16fit[i]->Draw("same");
        banana17fit[i]->Draw("same");
        banana18fit[i]->Draw("same");
        banana19fit[i]->Draw("same");
        banana20fit[i]->Draw("same");
        banana21fit[i]->Draw("same");
        banana22fit[i]->Draw("same");
        banana23fit[i]->Draw("same");
        banana24fit[i]->Draw("same");
        banana25fit[i]->Draw("same");
        banana26fit[i]->Draw("same");
        banana27fit[i]->Draw("same");
        banana28fit[i]->Draw("same");
        banana29fit[i]->Draw("same");
        banana30fit[i]->Draw("same");
        c2->SaveAs(Form("%s/bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), i));
        c2->Write(Form("bananas_2Dfitfunction_L%d", i));
        delete c2;
    
        dir_aligned->cd();

        // ________________________________________________________________
        Double_t minEnergy = 1;
        Double_t maxEnergy = 8000.0;
        Double_t energyStep = 1.0;

        // Iterate over the energy values and compute the precomputed values

        // printf("\n\nCalculating precomputed values for detector %d...\n", i);
        for (Double_t e = minEnergy; e <= maxEnergy; e += energyStep)
        {
            precomputedValues1[std::make_pair(i, e)] = banana1fit[i]->Eval(e);
            precomputedValues2[std::make_pair(i, e)] = banana2fit[i]->Eval(e);
            precomputedValues3[std::make_pair(i, e)] = banana3fit[i]->Eval(e);
            precomputedValues4[std::make_pair(i, e)] = banana4fit[i]->Eval(e);
            precomputedValues5[std::make_pair(i, e)] = banana5fit[i]->Eval(e);
            precomputedValues6[std::make_pair(i, e)] = banana6fit[i]->Eval(e);
            precomputedValues7[std::make_pair(i, e)] = banana7fit[i]->Eval(e);
            precomputedValues8[std::make_pair(i, e)] = banana8fit[i]->Eval(e);
            precomputedValues9[std::make_pair(i, e)] = banana9fit[i]->Eval(e);
            precomputedValues10[std::make_pair(i, e)] = banana10fit[i]->Eval(e);
            precomputedValues11[std::make_pair(i, e)] = banana11fit[i]->Eval(e);
            precomputedValues12[std::make_pair(i, e)] = banana12fit[i]->Eval(e);
            precomputedValues13[std::make_pair(i, e)] = banana13fit[i]->Eval(e);
            precomputedValues14[std::make_pair(i, e)] = banana14fit[i]->Eval(e);
            precomputedValues15[std::make_pair(i, e)] = banana15fit[i]->Eval(e);
            precomputedValues16[std::make_pair(i, e)] = banana16fit[i]->Eval(e);
            precomputedValues17[std::make_pair(i, e)] = banana17fit[i]->Eval(e);
            precomputedValues18[std::make_pair(i, e)] = banana18fit[i]->Eval(e);
            precomputedValues19[std::make_pair(i, e)] = banana19fit[i]->Eval(e);
            precomputedValues20[std::make_pair(i, e)] = banana20fit[i]->Eval(e);
            precomputedValues21[std::make_pair(i, e)] = banana21fit[i]->Eval(e);
            precomputedValues22[std::make_pair(i, e)] = banana22fit[i]->Eval(e);
            precomputedValues23[std::make_pair(i, e)] = banana23fit[i]->Eval(e);
            precomputedValues24[std::make_pair(i, e)] = banana24fit[i]->Eval(e);
            precomputedValues25[std::make_pair(i, e)] = banana25fit[i]->Eval(e);
            precomputedValues26[std::make_pair(i, e)] = banana26fit[i]->Eval(e);
            precomputedValues27[std::make_pair(i, e)] = banana27fit[i]->Eval(e);
            precomputedValues28[std::make_pair(i, e)] = banana28fit[i]->Eval(e);
            precomputedValues29[std::make_pair(i, e)] = banana29fit[i]->Eval(e);
            precomputedValues30[std::make_pair(i, e)] = banana30fit[i]->Eval(e);
        }
    
        // ________________________________________________________________
        printf("\n\nCalculating the corrected time offset from the fit precomputed values for detector %d...\n", i);     
        for (unsigned long long entry = 0; entry < numEntries; entry++) 
        {
            offsetTree->GetEntry(entry);     
            timeSlowCorrected[i] = 0;
            Double_t t = timeOffsetBeforeCorrection[i];

            if (t == 0) continue;

            Double_t e = static_cast<Double_t>(static_cast<int>(energySlowCp[i]));
            Double_t timeOffset = 0.0;

            if ((t > 235) && (t<=245))
            {
                timeOffset = precomputedValues1[std::make_pair(i, e)];
            }
            else if ((t > 245) && (t<=255))
            {
                timeOffset = precomputedValues2[std::make_pair(i, e)];
            }
            else if ((t > 255) && (t<=265))
            {
                timeOffset = precomputedValues3[std::make_pair(i, e)];
            }
            else if ((t > 265) && (t<=275))
            {
                timeOffset = precomputedValues4[std::make_pair(i, e)];
            }
            else if ((t > 275) && (t<=285))
            {
                timeOffset = precomputedValues5[std::make_pair(i, e)];
            }
            else if ((t > 285) && (t<=295))
            {
                timeOffset = precomputedValues6[std::make_pair(i, e)];
            }
            else if ((t > 295) && (t<=305))
            {
                timeOffset = precomputedValues7[std::make_pair(i, e)];
            }
            else if ((t > 305) && (t<=315))
            {
                timeOffset = precomputedValues8[std::make_pair(i, e)];
            }
            else if ((t > 315) && (t<=325))
            {
                timeOffset = precomputedValues9[std::make_pair(i, e)];
            }
            else if ((t > 325) && (t<=335))
            {
                timeOffset = precomputedValues10[std::make_pair(i, e)];
            }
            else if ((t > 335) && (t<=345))
            {
                timeOffset = precomputedValues11[std::make_pair(i, e)];
            }
            else if ((t > 345) && (t<=355))
            {
                timeOffset = precomputedValues12[std::make_pair(i, e)];
            }
            else if ((t > 355) && (t<=365))
            {
                timeOffset = precomputedValues13[std::make_pair(i, e)];
            }
            else if ((t > 365) && (t<=375))
            {
                timeOffset = precomputedValues14[std::make_pair(i, e)];
            }
            else if ((t > 375) && (t<=385))
            {
                timeOffset = precomputedValues15[std::make_pair(i, e)];
            }
            else if ((t > 385) && (t<=395))
            {
                timeOffset = precomputedValues16[std::make_pair(i, e)];
            }
            else if ((t > 395) && (t<=405))
            {
                timeOffset = precomputedValues17[std::make_pair(i, e)];
            }
            else if ((t > 405) && (t<=415))
            {
                timeOffset = precomputedValues18[std::make_pair(i, e)];
            }
            else if ((t > 415) && (t<=425))
            {
                timeOffset = precomputedValues19[std::make_pair(i, e)];
            }
            else if ((t > 425) && (t<=435))
            {
                timeOffset = precomputedValues20[std::make_pair(i, e)];
            }
            else if ((t > 435) && (t<=445))
            {
                timeOffset = precomputedValues21[std::make_pair(i, e)];
            }
            else if ((t > 445) && (t<=455))
            {
                timeOffset = precomputedValues22[std::make_pair(i, e)];
            }
            else if ((t > 455) && (t<=465))
            {
                timeOffset = precomputedValues23[std::make_pair(i, e)];
            }
            else if ((t > 465) && (t<=475))
            {
                timeOffset = precomputedValues24[std::make_pair(i, e)];
            }
            else if ((t > 475) && (t<=485))
            {
                timeOffset = precomputedValues25[std::make_pair(i, e)];
            }
            else if ((t > 485) && (t<=495))
            {
                timeOffset = precomputedValues26[std::make_pair(i, e)];
            }
            else if ((t > 495) && (t<=505))
            {
                timeOffset = precomputedValues27[std::make_pair(i, e)];
            }
            else if ((t > 505) && (t<=515))
            {
                timeOffset = precomputedValues28[std::make_pair(i, e)];
            }
            else if ((t > 515) && (t<=525))
            {
                timeOffset = precomputedValues29[std::make_pair(i, e)];
            }
            else if ((t > 525) && (t<=535))
            {
                timeOffset = precomputedValues30[std::make_pair(i, e)];
            }
            else
            {
                timeOffset = 0;
            }

            // Calculate timeSlowCorrected and fill correctedTree using timeOffset
            if (timeOffset != 0.0) 
            {
                timeSlowCorrected[i] = timeSlowCp[i] - timeOffset;
                correctedTree->Fill();
            }
        }

        delete offsetTree;

        // ________________________________________________________________
        printf("\n\n Making corrected histograms for detector %d...\n", i);

        for (unsigned long long entry = 0; entry < numEntries; entry++) 
        {
            correctedTree->GetEntry(entry);

            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t energy = 0;
            Double_t nextTimeSlow = 0;
            Double_t deltaTime = 0;

            for (unsigned long long nextEntry = entry + 1; nextEntry < numEntries; nextEntry++) 
            {
                correctedTree->GetEntry(nextEntry);     
                if (timeSlowCorrected[i] == 0) continue;
                energy = energySlowCp[i];
                nextTimeSlow = timeSlowCorrected[i];
                
                deltaTime = (nextTimeSlow - refTimeFast);
               //  printf("Detector %d: Time Offset Forward = %f\n", i, deltaTime);

                // printf("Detector %d: Time Offset Forward = %f\n", i, deltaTime2);
                if ((deltaTime > 600.0) || (deltaTime < 0.0)) break;
                else
                {
                    timeOffsetCorrectedHistsNext[i]->Fill(deltaTime);
                    timeEnergyOffsetCorrectedHistsNext[i]->Fill(deltaTime, energy);
                    energyhist[i]->Fill(energy);
                    entry = nextEntry;
                    break;
                }
            }
        }

        // ________________________________________________________________
    
        timeOffsetCorrectedHistsNext[i]->Write();
        timeEnergyOffsetCorrectedHistsNext[i]->Write();
        energyhist[i]->Write();

        TCanvas * c5 = new TCanvas("c5", "CorrectedTimeOffsetNext", 800, 600);
        timeOffsetCorrectedHistsNext[i]->GetXaxis()->SetRangeUser(-400, 400);
        timeOffsetCorrectedHistsNext[i]->SetFillColor(kBlue);
        timeOffsetCorrectedHistsNext[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        timeOffsetCorrectedHistsNext[i]->GetYaxis()->SetTitle("Counts/ns");
        timeOffsetCorrectedHistsNext[i]->Draw();
        c5->SaveAs(Form("%s/corrected_time_offset_hist_L%d.root", dirPathplots.c_str(), i));
        delete c5;

        TCanvas * c6 = new TCanvas("c6", "CorrectedTimeOffsetEnergyNext", 800, 600);
        timeEnergyOffsetCorrectedHistsNext[i]->GetXaxis()->SetRangeUser(-400, 400);
        timeEnergyOffsetCorrectedHistsNext[i]->GetYaxis()->SetRangeUser(0, 8000);
        timeEnergyOffsetCorrectedHistsNext[i]->SetStats(0);
        timeEnergyOffsetCorrectedHistsNext[i]->SetMaximum(500);
        timeEnergyOffsetCorrectedHistsNext[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        timeEnergyOffsetCorrectedHistsNext[i]->GetYaxis()->SetTitle("Energy (keV)");
        timeEnergyOffsetCorrectedHistsNext[i]->Draw("colz");
        c6->SaveAs(Form("%s/corrected_time_energy_offset_hist_L%d.root", dirPathplots.c_str(), i));
        delete c6;

        TCanvas *c15 = new TCanvas("c15", "Energy", 800, 600);
        energyhist[i]->GetXaxis()->SetRangeUser(0, 8000);
        energyhist[i]->GetXaxis()->SetTitle("Energy (keV)");
        energyhist[i]->GetYaxis()->SetTitle("Counts");
        energyhist[i]->Draw();
        c15->SaveAs(Form("%s/energyhist_L%d.root", dirPathplots.c_str(), i));
        delete c15;

        // ________________________________________________________________
        // Write the histograms to the output file
        slowfasttimediff_hists[i]->Write();
        slowfasttimediff_energy_hists[i]->Write();
        energy_slowfasttimediff_hists[i]->Write();
    }

    // ________________________________________________________________
    // Set the branches for the aligned data in the new copyTree
    TTree *alignedTree = new TTree("AlignedData", "AlignedData");
    alignedTree->Branch("detectorID", &detectorID, "detectorID/I");
    alignedTree->Branch("globalTime", &timeGlobal[0], "globalTime/D");
    alignedTree->Branch("timeF", &alignedTimeFast[0], "timeF/D");
    alignedTree->Branch("energyF", &alignedEnergyFast[0], "energyF/D");
    alignedTree->Branch("timeS", &alignedTimeSlow[0], "timeS/D");
    alignedTree->Branch("energyS", &alignedEnergySlow[0], "energyS/D");

    // Align the slow time and energy data for each detector based on the calculated offsets
    for (int i = 0; i < 4 +2; i++) 
    {
        printf("\n\nAligning data for detector %d...\n", i);
        for (unsigned long long entry = 0; entry < numEntries; entry++) 
        {
            correctedTree->GetEntry(entry);

            if (i == RFIndex) 
            {
                detectorID = RFIndex;
                if (timeFastCp[i] == 0) continue;
                timeGlobal[0] = timeFastCp[i];
                alignedTimeFast[0] = timeFastCp[i];
                alignedEnergyFast[0] = 0;
                alignedTimeSlow[0] = 0;
                alignedEnergySlow[0] = 0;
                // printf("RF Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i == POLARISIndex) 
            {
                detectorID = POLARISIndex;
                if (fastTimePOLARISCp[0] == 0) continue;
                timeGlobal[0] = fastTimePOLARISCp[0];
                alignedTimeFast[0] = fastTimePOLARISCp[0];
                alignedEnergyFast[0] = fastEnergyPOLARISCp[0];
                alignedTimeSlow[0] = 0;
                alignedEnergySlow[0] = 0;
                // printf("POLARIS Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i < RFIndex)
            {
                detectorID = i;
                // printf("Detector %d: timeFastCp = %f | , timeSlowCorrectedCp = %f\n", i, timeFastCp[i], timeSlowCorrectedCp[i]);
                if ((timeFastCp[i] == 0) && (timeSlowCorrectedCp[i] != 0)) 
                {
                    timeGlobal[0] = timeSlowCorrectedCp[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrectedCp[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                    // printf("Loop1 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrectedCp = %f, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrectedCp[i] == 0)) 
                {
                    timeGlobal[0] = timeFastCp[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrectedCp[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                    // printf("Loop2 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrectedCp = %f\n, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrectedCp[i] != 0)) 
                {
                    // The ratio in the print below should be 1
                    //printf("Error: Both fast and slow times are non-zero for detector %d\n, timeFastCp = %f | , timeSlowCorrectedCp + timeOffsets = %f\n, ratio of timeFastCp to timeSlowCp + timeOffsets = %f\n", i, timeFastCp[i], timeSlowCp[i] + timeOffsets[i], (timeFastCp[i]/(timeSlowCp[i] + timeOffsets[i])));
                    timeGlobal[0] = timeFastCp[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrectedCp[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                }
                else 
                {
                    // The slow energy in the print below should be zero
                    //printf("Error: Both fast and slow times are zero for detector %d\n, Slow Energy = %f\n", i, energySlowCp[i]);
                    continue;
                }
                alignedTree->Fill();

            }
        }
    }

    alignedTree->Write();
}

// ________________________________________ Main function ________________________________________
int main(int argc, char* argv[])
{
    // Start timer
    time(&start);

    // Extract the run number and directory path from the terminal input
    // example of input: clear && g++ -std=c++0x sort2labr.C -o exe `root-config --cflags --libs` -lSpectrum && ./exe ~/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63.root
    if (argc < 2) {
        std::cout << "Usage: ./program <dir_file>" << std::endl;
        return 1;
    }
    std::string inputFileName = argv[1];
    //std::string inputFileName = "/Users/shanyn/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63.root";
    std::string dataType = argv[2];
    //std::string dataType = "source";

    // the below bool is used to determine whether the data is beam or angle data
    // this is used to determine if the fast time detector intercalibration should be applied 
    // it is not applied for source data
    bool beam = false;
    if (dataType == "beam") beam = true;

    std::string runStr = inputFileName.substr(inputFileName.find_last_of("/") + 2);
    std::string dirPath = inputFileName.substr(0, inputFileName.find_last_of("/") + 1);
    int run = std::atoi(runStr.c_str());
    std::string outputFileName = dirPath + "R" + std::to_string(run) + "_aligned.root";

    std::cout << "Run: " << run << std::endl;
    std::cout << "Directory Path: " << dirPath << std::endl;
    std::cout << "Output File Name: " << outputFileName << std::endl;
    std::cout << "_______________________________________________" << std::endl;


    std::string dirPathplots = dirPath + "plots/sort2labr/";
    // check to see if dirPath contains a folder called plots. if not, create it
    if (dirPath.find("plots/sort2labr/") == std::string::npos)
    {
        std::string dirPath1 = dirPath + "plots/";
        std::string dirPathplots = dirPath1 + "sort2labr/";
        std::string command = "mkdir " + dirPath1;
        std::system(command.c_str());
        command = "mkdir " + dirPathplots;
        system(command.c_str());
    }

    // ________________________________________________________________________________________________________________________
    // Reserve memory for heavy arrays
    energySlow.reserve(NumDetectors);
    timeSlow.reserve(NumDetectors);
    timeFast.reserve(NumDetectors + 1);
    fastEnergyPOLARIS.reserve(1);
    fastTimePOLARIS.reserve(1);
    
    // Open the ROOT file
    TFile*file = new TFile(inputFileName.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file->Get("LaBrData"));
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    unsigned long long numEntries = (tree->GetEntries());
    printf("Number of entries: %lld\n", numEntries);

    if (!file->IsOpen()) 
    {
    std::cerr << "Error: File not open" << std::endl;
    return 1;
    }

    if (!outputFile->IsOpen()) 
    {
        std::cerr << "Error: Output file not open" << std::endl;
        return 1;
    }

    TDirectory *dir_fasttimecalib = outputFile->mkdir("dir_fasttimecalib");
    TDirectory *dir_banana = outputFile->mkdir("dir_banana");
    TDirectory *dir_aligned = outputFile->mkdir("dir_aligned");

    // ________________________________________________________________________________________________________________________
    
    // FUNCTIONS
    runduration(tree, tensOfNanoSeconds, hundredsOfNanoSeconds,numEntries);
    fasttime_calibration_constants(beam, outputFile, dirPathplots, dir_fasttimecalib, tree, tensOfNanoSeconds, hundredsOfNanoSeconds,numEntries); // this function is only applied to beam data
    align_offset(outputFile, dirPathplots, tree, dir_banana, dir_aligned, mean_fasttimediff, tensOfNanoSeconds, hundredsOfNanoSeconds, numEntries);

    // ________________________________________________________________________________________________________________________

    // Cleanup
    file->Close();
    outputFile->Close();

    time(&end);
    double time_taken = double((end - start)/60.0);
    std::cout << "Time taken by program is : " << time_taken << std::setprecision(5);
    std::cout << " min " << std::endl;

    return 0;
}
