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
  LaBrData->Branch("energyP", &energyPolaris[0], "energyP/D");
  LaBrData->Branch("timeP", &timePolaris[0], "timeP/D");*/

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
#include <TMath.h>
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

const int NumDetectors = 2;
const int RFIndex = 4;
const int POLARISIndex = 5;

time_t start, end;

// variables for input tree data
std::vector<Double_t> energySlow(NumDetectors);
std::vector<Double_t> timeSlow(NumDetectors);
std::vector<Double_t> timeFast(NumDetectors); // Additional branch needs to be added for RF for beamq
std::vector<Double_t> energyPolaris(1);
std::vector<Double_t> timePolaris(1);
std::vector<Double_t> alignedSlowECalibL(NumDetectors);

// variables for copytree data
std::vector<Double_t> energySlowCp(NumDetectors);
std::vector<Double_t> timeSlowCp(NumDetectors);
std::vector<Double_t> timeSlowCorrected(NumDetectors);
std::vector<Double_t> timeSlowCorrectedCp(NumDetectors);
std::vector<Double_t> timeFastCp(NumDetectors); 
std::vector<Double_t> energyPolarisCp(1);
std::vector<Double_t> timePolarisCp(1);
std::vector<Double_t> timeOffsetBeforeCorrection(NumDetectors);
std::vector<Double_t> timeOffsetBeforeCorrectionCp(NumDetectors);

// Variables for aligned tree data
Int_t detectorID;
std::vector<Double_t> alignedTimeFast(1);
std::vector<Double_t> alignedEnergyFast(1);
std::vector<Double_t> alignedTimeSlow(1);
std::vector<Double_t> alignedEnergySlow(1);

 std::vector<TH1D*> slowfasttimediff_hists(NumDetectors);
std::vector<TH1D*> timeOffsetCorrectedHists(NumDetectors);
std::vector <TH1D*> timefastslowratio_hists(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists(NumDetectors); 
std::vector<TH2D*> slowfasttimediff_energy_hists1(1);
std::vector<TH2D*> energy_slowfasttimediff_hists1(1); 
std::vector<TH2D*> slowfasttimediff_energy_hists2(1);
std::vector<TH2D*> energy_slowfasttimediff_hists2(1); 
std::vector<TH2D*> slowfasttimediff_energy_hists3(1);
std::vector<TH2D*> energy_slowfasttimediff_hists3(1);
std::vector<TH2D*> slowfasttimediff_energy_hists4(1);
std::vector<TH2D*> energy_slowfasttimediff_hists4(1);
std::vector<TH2D*> slowfasttimediff_energy_hists5(1);
std::vector<TH2D*> energy_slowfasttimediff_hists5(1);
std::vector<TH2D*> slowfasttimediff_energy_hists6(1);
std::vector<TH2D*> energy_slowfasttimediff_hists6(1);
std::vector<TH2D*> slowfasttimediff_energy_hists7(1);
std::vector<TH2D*> energy_slowfasttimediff_hists7(1);
std::vector<TH2D*> slowfasttimediff_energy_hists8(1);
std::vector<TH2D*> energy_slowfasttimediff_hists8(1);
std::vector<TH2D*> slowfasttimediff_energy_hists9(1);
std::vector<TH2D*> energy_slowfasttimediff_hists9(1);
std::vector<TH2D*> slowfasttimediff_energy_hists10(1);
std::vector<TH2D*> energy_slowfasttimediff_hists10(1);
std::vector<TH2D*> slowfasttimediff_energy_hists11(1);
std::vector<TH2D*> slowfasttimediff_energy_hists6_2(1);
std::vector<TH2D*> energy_slowfasttimediff_hists6_2(1);
std::vector<TH2D*> slowfasttimediff_energy_hists6_3(1);
std::vector<TH2D*> energy_slowfasttimediff_hists6_3(1);
std::vector<TH2D*> energy_slowfasttimediff_hists11(1);
std::vector<TH2D*> timeEnergyOffsetCorrectedHists(NumDetectors);

std::vector<TH2D*> slowfasttimediff_energy_hists_again(1);
std::vector<TH2D*> energy_slowfasttimediff_hists_again(1);


std::vector<Double_t> mean_fasttimediff(NumDetectors-1);
std::vector<Double_t> sigma_fasttimediff(NumDetectors-1);
bool tensOfNanoSeconds = false;
bool hundredsOfNanoSeconds = false;

// initialise the preComputedValues map
std::map<std::pair<int, double>, double> precomputedValues1;
std::map<std::pair<int, double>, double> precomputedValues2;
std::map<std::pair<int, double>, double> precomputedValues3;
std::map<std::pair<int, double>, double> precomputedValues4;
std::map<std::pair<int, double>, double> precomputedValues5_1;
std::map<std::pair<int, double>, double> precomputedValues5_2;
std::map<std::pair<int, double>, double> precomputedValues6;
std::map<std::pair<int, double>, double> precomputedValues6_2;
std::map<std::pair<int, double>, double> precomputedValues6_3;
std::map<std::pair<int, double>, double> precomputedValues7_1;
std::map<std::pair<int, double>, double> precomputedValues7_2;
std::map<std::pair<int, double>, double> precomputedValues8;
std::map<std::pair<int, double>, double> precomputedValues9;
std::map<std::pair<int, double>, double> precomputedValues10;
std::map<std::pair<int, double>, double> precomputedValues11;
std::map<std::pair<int, double>, double> precomputedValues_again;

// BANANA FIT FUNCTIONS //
TF1* banana1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots);
TF1* banana1fit[1];

TF1* banana2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists2, TH2D* energy_slowfasttimediff_hists2, std::string dirPathplots);
TF1* banana2fit[1];

TF1* banana3(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists3, TH2D* energy_slowfasttimediff_hists3, std::string dirPathplots);
TF1* banana3fit[1];

TF1* banana4(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists4, TH2D* energy_slowfasttimediff_hists4, std::string dirPathplots);
TF1* banana4fit[1];

TF1* banana5_1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots);
TF1* banana5_1fit[1];

TF1* banana5_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots);
TF1* banana5_2fit[1];

TF1* banana6(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6, TH2D* energy_slowfasttimediff_hists6, std::string dirPathplots);
TF1* banana6fit[1];

TF1* banana7_1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots);
TF1* banana7_1fit[1];

TF1* banana7_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots);
TF1* banana7_2fit[1];

TF1* banana8(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists8, TH2D* energy_slowfasttimediff_hists8, std::string dirPathplots);
TF1* banana8fit[1];

TF1* banana9(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists9, TH2D* energy_slowfasttimediff_hists9, std::string dirPathplots);
TF1* banana9fit[1];

TF1* banana10(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists10, TH2D* energy_slowfasttimediff_hists10, std::string dirPathplots);
TF1* banana10fit[1];

TF1* banana11(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists11, TH2D* energy_slowfasttimediff_hists11, std::string dirPathplots);
TF1* banana11fit[1];

TF1* banana6_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6_2, TH2D* energy_slowfasttimediff_hists6_2, std::string dirPathplots);
TF1* banana6_2fit[1];

TF1* banana6_3(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6_3, TH2D* energy_slowfasttimediff_hists6_3, std::string dirPathplots);
TF1* banana6_3fit[1];

TF1* bananafinal(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists_again, TH2D* energy_slowfasttimediff_hists_again, std::string dirPathplots);
TF1* bananafinalfit[1];
// END BANANA FIT FUNCTIONS //

// banana1 good fit
TF1* banana1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots)
{
    std::cout <<"\nFitting banana1 (< 40 ns  offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists1->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists1->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana1Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana1s for detector %d...\n", i);


    TCanvas* cb = new TCanvas();
	TGraphErrors* gr1 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr1 -> SetTitle(Form("Banana1 Fit L%d", i));
    //gr1 -> SetStats(0);
	gr1 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr1 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr1 -> SetMarkerStyle(20);
    gr1 -> SetMarkerSize(0.8);
    gr1 -> SetMarkerColor(kBlue);
    gr1->GetYaxis()->SetRangeUser(26, 46);
    gr1->GetXaxis()->SetRangeUser(400, 2000);
    TF1* fit1 = new TF1("fit1", "expo(0)+pol2(2)", 450, 2000);
    gr1->Draw("AL*");
    gr1->Fit("fit1", ""); // Perform the fit using the initial parameters
    fit1 -> SetParameters( 5.17847, -0.00871973, 30.7512, 0.0054537, -2.12777e-06);
    fit1->Draw("same");
    //gr1->GetXaxis()->SetRangeUser(0, 2000);
    cb -> SaveAs(Form("%s/R%s_banana1Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists2, TH2D* energy_slowfasttimediff_hists2, std::string dirPathplots)
{
    std::cout <<"Fitting banana2 (>40<56 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists2->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists2->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana2Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana2s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr2 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr2 -> SetTitle(Form("Banana2 Fit L%d", i));
    //gr2 -> SetStats(0);
	gr2 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr2 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr2 -> SetMarkerStyle(20);
    gr2 -> SetMarkerSize(0.8);
    gr2 -> SetMarkerColor(kBlue);
    gr2->GetYaxis()->SetRangeUser(42, 54);
    gr2->GetXaxis()->SetRangeUser(212, 2000);
    TF1* fit2 = new TF1("fit2", "expo(0)+pol2(2)", 200, 1600); 
    gr2->Draw("AL*");
    gr2->Fit("fit2", "");
    fit2->SetParameters(6.48699,-0.0213749,47.7895,-0.00057454,4.54381e-07);
    fit2->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana2Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana3(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists3, TH2D* energy_slowfasttimediff_hists3, std::string dirPathplots)
{
    std::cout <<"Fitting banana3 (>56 < 73 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists3->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists3->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana3Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana3s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr3 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr3 -> SetTitle(Form("Banana3 Fit L%d", i));
    //gr3 -> SetStats(0);
	gr3 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr3 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr3 -> SetMarkerStyle(20);
    gr3 -> SetMarkerSize(0.8);
    gr3 -> SetMarkerColor(kBlue);
    gr3->GetYaxis()->SetRangeUser(0, 120);
    gr3->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit3 = new TF1("fit3", "pol3(0)", 47, 1700);
    gr3->Draw("AL*");
    gr3->Fit("fit3", "");
    fit3->SetParameters(70.9007, -0.0225665,-1.03825e-05,1.41102e-08 );
    fit3->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana3Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana4(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists4, TH2D* energy_slowfasttimediff_hists4, std::string dirPathplots)
{
    std::cout <<"Fitting banana4 (>73 < 102 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists4->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists4->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana4Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana4s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr4 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr4 -> SetTitle(Form("Banana4 Fit L%d", i));
    //gr4 -> SetStats(0);
	gr4 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr4 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr4 -> SetMarkerStyle(20);
    gr4 -> SetMarkerSize(0.8);
    gr4 -> SetMarkerColor(kBlue);
    gr4->GetYaxis()->SetRangeUser(0, 102);
    gr4->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit4 = new TF1("fit4", "pol2(0)+expo(3)", 50, 400);
    gr4->Draw("AL*");
    gr4->Fit("fit4", "");
    fit4->SetParameters(-9.38799e+02,2.74694e+00,-2.25030e-03,6.97939e+00,-3.27849e-03);
    fit4->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana4Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana5_1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots)
{
    std::cout <<"Fitting banana5 (>102 < 300 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists5->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists5->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana5_1Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana5s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr5 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr5 -> SetTitle(Form("Banana5 Fit L%d", i));
    //gr5 -> SetStats(0);
	gr5 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr5 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr5 -> SetMarkerStyle(20);
    gr5 -> SetMarkerSize(0.8);
    gr5 -> SetMarkerColor(kBlue);
    gr5->GetYaxis()->SetRangeUser(0, 300);
    gr5->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit5 = new TF1("fit5", "expo(0)", 0, 150);
    gr5->Draw("AL*");
    gr5->Fit("fit5", "");
    fit5->SetParameters(5.52032e+00,-1.73912e-02);
    fit5->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana5_1Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana5_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists5, TH2D* energy_slowfasttimediff_hists5, std::string dirPathplots)
{
    std::cout <<"Fitting banana5 (>102 < 300 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists5->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists5->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana5_2Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana5s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
	TGraphErrors* gr5 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr5 -> SetTitle(Form("Banana5 Fit L%d", i));
    //gr5 -> SetStats(0);
	gr5 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
	gr5 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr5 -> SetMarkerStyle(20);
    gr5 -> SetMarkerSize(0.8);
    gr5 -> SetMarkerColor(kBlue);
    gr5->GetYaxis()->SetRangeUser(0, 300);
    gr5->GetXaxis()->SetRangeUser(0, 2000);
    gr5->Draw("AL*");
    TF1 * fitt = new TF1("fitt", "expo(0)", 180, 450);
    gr5->Fit("fitt", "");
    fitt->SetParameters(4.85465e+00,-2.86881e-17);
    fitt->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana5_2Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    cb->Write();

    delete gr5;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana5 fit for detector %d...\n", i);
    return fitt;
}

TF1* banana6(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6, TH2D* energy_slowfasttimediff_hists6, std::string dirPathplots)
{
    std::cout <<"Fitting banana6 (< 18 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists6->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
    {
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana6Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana6s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr6 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr6 -> SetTitle(Form("Banana6 Fit L%d", i));
    //gr6 -> SetStats(0);
    gr6 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr6 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr6 -> SetMarkerStyle(20);
    gr6 -> SetMarkerSize(0.8);
    gr6 -> SetMarkerColor(kBlue);
    gr6->GetYaxis()->SetRangeUser(0, 30);
    gr6->GetXaxis()->SetRangeUser(550, 2000);
    TF1* fit6 = new TF1("fit6", "pol2(0)", 330, 2000);
    gr6->Draw("AL*");
    gr6->Fit("fit6", "");
    fit6->SetParameters(2.63595e+01, -1.24823e-02,3.80862e-06);
    fit6->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana6Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    cb->Write();

    delete gr6;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana3 fit for detector %d...\n", i);
    return fit6;
}

TF1* banana7_1(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots)
{
    std::cout <<"Fitting banana7 (>18 < 32 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists7->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana7_2Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana7s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr7 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr7 -> SetTitle(Form("Banana7 Fit L%d", i));
    //gr7 -> SetStats(0);
    gr7 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr7 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr7 -> SetMarkerStyle(20);
    gr7 -> SetMarkerSize(0.8);
    gr7 -> SetMarkerColor(kBlue);
    gr7->GetYaxis()->SetRangeUser(18, 32);
    gr7->GetXaxis()->SetRangeUser(0, 900);
    TF1* fit7 = new TF1("fit7", "pol2(0)", 222, 900);
    gr7->Draw("AL*");
    gr7->Fit("fit7", "");
    fit7->SetParameters(37.8859,-0.0182283,9.5443e-06);
    fit7->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana7_1Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana7_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots)
{
    std::cout <<"Fitting banana7 (>18 < 32 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists7->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana7_2Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana7s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr7 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr7 -> SetTitle(Form("Banana7 Fit L%d", i));
    //gr7 -> SetStats(0);
    gr7 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr7 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr7 -> SetMarkerStyle(20);
    gr7 -> SetMarkerSize(0.8);
    gr7 -> SetMarkerColor(kBlue);
    gr7->GetYaxis()->SetRangeUser(18, 32);
    gr7->GetXaxis()->SetRangeUser(9000, 2000);
    TF1* fitt7 = new TF1("fitt7", "pol2(0)", 900, 1970);
    gr7->Draw("AL*");
    gr7->Fit("fitt7", "");
    fitt7->SetParameters(23.2273, -0.00423945, 1.47897e-06);
    fitt7->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana7_2Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    cb->Write();

    delete gr7;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana7 fit for detector %d...\n", i);
    return fitt7;
}

TF1* banana8(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists8, TH2D* energy_slowfasttimediff_hists8, std::string dirPathplots)
{
    std::cout <<"Fitting banana8 (>32 < 56ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists8->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists8->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));
        
        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana8Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana8s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr8 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr8 -> SetTitle(Form("Banana8 Fit L%d", i));
    //gr8 -> SetStats(0);
    gr8 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr8 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr8 -> SetMarkerStyle(20);
    gr8 -> SetMarkerSize(0.8);
    gr8 -> SetMarkerColor(kBlue);
    gr8->GetYaxis()->SetRangeUser(32, 56);
    gr8->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit8 = new TF1("fit8", "pol3(0)", 111,1500);
    gr8->Draw("AL*");
    gr8->Fit("fit8", "");
    fit8->SetParameters(71.4321,-0.0943406,7.96066e-05,-2.21986e-08);
    fit8->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana8Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana9(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists9, TH2D* energy_slowfasttimediff_hists9, std::string dirPathplots)
{
    std::cout <<"Fitting banana9(>56 < 74 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists9->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists9->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth/2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;

    }

    p -> SaveAs(Form("%s/R%s_banana9Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana9s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr9 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr9 -> SetTitle(Form("Banana9 Fit L%d", i));
    //gr9 -> SetStats(0);
    gr9 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr9 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr9 -> SetMarkerStyle(20);
    gr9 -> SetMarkerSize(0.8);
    gr9 -> SetMarkerColor(kBlue);
    gr9->GetYaxis()->SetRangeUser(55, 60);
    gr9->GetXaxis()->SetRangeUser(400, 2000);
    TF1* fit9 = new TF1("fit9", "expo(0)", 340,1580);
    gr9->Draw("AL*");
    gr9->Fit("fit9", "");
    fit9->SetParameters(4.03405e+00,8.59436e-06);
    fit9->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana9Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana10(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists10, TH2D* energy_slowfasttimediff_hists10, std::string dirPathplots)
{
    std::cout <<"Fitting banana10 (>74 < 102 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists10->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists10->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth / 2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana10Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;
    
    printf("Fitting banana10s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr10 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr10 -> SetTitle(Form("Banana10 Fit L%d", i));
    //gr10 -> SetStats(0);
    gr10 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr10 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr10 -> SetMarkerStyle(20);
    gr10 -> SetMarkerSize(0.8);
    gr10 -> SetMarkerColor(kBlue);
    gr10->GetYaxis()->SetRangeUser(74, 102);
    gr10->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit10 = new TF1("fit10", "pol0(0)", 0, 2000);
    gr10->Draw("AL*");
    gr10->Fit("fit10", "");
    fit10->SetParameters(87.6);
    fit10->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana10Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana11(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists11, TH2D* energy_slowfasttimediff_hists11, std::string dirPathplots)
{
    std::cout <<"Fitting banana11 (>102 < 300 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists11->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists11->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

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
        //if (bin > 4200) slicewidth = 150;

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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth / 2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana11Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana11s for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr11 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr11 -> SetTitle(Form("Banana11 Fit L%d", i));
    //gr11 -> SetStats(0);
    gr11 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr11 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr11 -> SetMarkerStyle(20);
    gr11 -> SetMarkerSize(0.8);
    gr11 -> SetMarkerColor(kBlue);
    gr11->GetYaxis()->SetRangeUser(0, 300);
    gr11->GetXaxis()->SetRangeUser(0, 140);
    TF1* fit11 = new TF1("fit11", "expo(0)", 0, 2000);
    gr11->Draw("AL*");
    gr11->Fit("fit11", "");
    fit11->SetParameters(5.33894e+00,-1.19345e-02);
    fit11->Draw("same");
    cb -> SaveAs(Form("%s/R%s_banana11Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
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

TF1* banana6_2(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6_2, TH2D* energy_slowfasttimediff_hists6_2, std::string dirPathplots)
{
    std::cout <<"Fitting banana11 (>102 < 300 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists6_2->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6_2->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists6_2->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists6_2->Draw();

    yBins = slowfasttimediff_energy_hists6_2->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;

        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists6_2->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth / 2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana6_2Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana6_2s for detector %d...\n", i);

    TCanvas* c6_2 = new TCanvas();
    TGraphErrors* gr6_2 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr6_2 -> SetTitle(Form("Banana6_2 Fit L%d", i));
    //gr6_2 -> SetStats(0);
    gr6_2 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr6_2 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr6_2 -> SetMarkerStyle(20);
    gr6_2 -> SetMarkerSize(0.8);
    gr6_2 -> SetMarkerColor(kBlue);
    gr6_2->GetYaxis()->SetRangeUser(0, 200);
    gr6_2->GetXaxis()->SetRangeUser(0, 300);
    TF1* fit6_2 = new TF1("fit6_2", "pol2(0)+expo(3)", 9, 300);
    gr6_2->Draw("AL*");
    gr6_2->Fit("fit6_2", "");
    fit6_2->SetParameters(-294.674,1.86449,-0.00292699,6.22,-0.00791423); 
    fit6_2->Draw("same");
    c6_2 -> SaveAs(Form("%s/R%s_banana6_2Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    c6_2->Write();

    delete gr6_2;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana6_2 fit for detector %d...\n", i);
    return fit6_2;
}

TF1* banana6_3(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists6_3, TH2D* energy_slowfasttimediff_hists6_3, std::string dirPathplots)
{
    std::cout <<"Fitting banana11 (>102 < 300 ns offset) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists6_3->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6_3->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists6_3->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists6_3->Draw();

    yBins = slowfasttimediff_energy_hists6_3->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;

        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists6_3->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth / 2.);
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_banana6_3Slices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting banana6_3s for detector %d...\n", i);

    TCanvas* c6_3 = new TCanvas();
    TGraphErrors* gr6_3 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr6_3 -> SetTitle(Form("Banana6_3 Fit L%d", i));
    //gr6_3 -> SetStats(0);
    gr6_3 -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr6_3 -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr6_3 -> SetMarkerStyle(20);
    gr6_3 -> SetMarkerSize(0.8);
    gr6_3 -> SetMarkerColor(kBlue);
    gr6_3->GetYaxis()->SetRangeUser(0, 200);
    gr6_3->GetXaxis()->SetRangeUser(0, 300);
    TF1* fit6_3 = new TF1("fit6_3", "pol2(0)+expo(3)", 100, 200);
    gr6_3->Draw("AL*");
    gr6_3->Fit("fit6_3", "");
    fit6_3->SetParameters(-294.674,1.86449,-0.00292699,6.22,-0.00791423); 
    fit6_3->Draw("same");
    c6_3 -> SaveAs(Form("%s/R%s_banana6_3Fit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    c6_3->Write();

    delete gr6_3;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning banana6_3 fit for detector %d...\n", i);
    return fit6_3;
}

TF1* bananafinal(std::string run,TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists_again, TH2D* energy_slowfasttimediff_hists_again, std::string dirPathplots)
{
    std::cout <<"Fitting final banana for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists_again->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists_again->GetXaxis()->GetNbins();

    int slicewidth = 20, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge(yBins);
    std::vector<Double_t> slice_charge_err(yBins);
    std::vector<Double_t> slice_time(yBins);
    std::vector<Double_t> slice_time_err(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists_again->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists_again->Draw();

    yBins = slowfasttimediff_energy_hists_again->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 4200) slicewidth = 150;

        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists_again->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
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

        if (fitSlice->GetParameter(1) <= 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        slice_time.push_back(fitSlice->GetParameter(1));
        slice_charge.push_back((bin+binmax)/2.);
        slice_charge_err.push_back(slicewidth / 2.);
        slice_time_err.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/R%s_bananafinalSlices_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    delete p;

    printf("Fitting final banana for detector %d...\n", i);

    TCanvas* cb = new TCanvas();
    TGraphErrors* gr = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
    gr -> SetTitle(Form("Final Banana Fit L%d", i));
    //gr -> SetStats(0);
    gr -> GetYaxis()->SetTitle("Centroid Time Difference (ns)");
    gr -> GetXaxis()->SetTitle("Centroid Energy (keV)");
    gr -> SetMarkerStyle(20);
    gr -> SetMarkerSize(0.8);
    gr -> SetMarkerColor(kBlue);
    gr->GetYaxis()->SetRangeUser(0,30);
    gr->GetXaxis()->SetRangeUser(0, 2000);
    TF1* fit = new TF1("fit", "pol0(0)", 0, 1500);
    gr->Draw("AL*");
    gr->Fit("fit", "");
    fit->SetParameters(12.6376);
    fit->Draw("same");
    cb -> SaveAs(Form("%s/R%s_bananafinalFit_L%d.root", dirPathplots.c_str(), run.c_str(), i));
    cb->Write();

    delete gr;
    printf("Clearing variables for detector...\n");
    slice_charge.clear();
    slice_charge_err.clear();
    slice_time.clear();
    slice_time_err.clear();

    printf("Returning final banana fit for detector %d...\n", i);
    return fit;
}

void runduration(TTree* tree, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds,uint64_t numEntries)
{
    for (int i = 0; i < NumDetectors; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("energyP", &energyPolaris[0]);
    tree->SetBranchAddress("timeP", &timePolaris[0]);

    Double_t firstTime = 0;
    Double_t lastTime = 0;

    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);
        if (timeFast[0] > 0)
        {
            firstTime = timeFast[0];
            break;
        }
    }

    // Find the last non-zero entry in the timeFL0 branch of the tree
    for (uint64_t entry = numEntries - 1; entry >= 0; entry--) 
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

std::vector<Double_t> fasttime_calibration_constants(std::string run,bool beam, TFile *outputFile ,std::string dirPathplots, TDirectory *dir_fasttimecalib, TTree* tree, bool tensOfNanoSeconds, bool hundredsOfNanoSeconds, uint64_t numEntries)
{
    if (beam)
    {
        printf("\n\nBeam run: Calculating inter-detector fast time calibration constants (used in next function)...\n");
        // Set branch addresses for each detector
        for (int i = 0; i < NumDetectors; i++) 
        {
            tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
            tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
            tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
        }
        tree->SetBranchAddress("energyP", &energyPolaris[0]);
        tree->SetBranchAddress("timeP", &timePolaris[0]);

        // Initialise the histograms
        std::vector<TH1D*> fasttimecalib_hists(NumDetectors-1);
        for (int i = 1; i < NumDetectors; i++) 
        {
            fasttimecalib_hists[i] = new TH1D(Form("fasttimecalib_hists%d", i), Form("Fast time calibration L0 - L%d", i), 1000 , 0, 1000);
        }

        for (uint64_t entry = 0; entry < numEntries; entry++) 
        {
            tree->GetEntry(entry);

            Double_t timeFastDet0 = 0;  // Initialize
            if (timeFast[0] == 0) continue;

            if (tensOfNanoSeconds) timeFastDet0 = timeFast[0] / 10;
            else if (hundredsOfNanoSeconds) timeFastDet0 = timeFast[0] / 100;
            else timeFastDet0 = timeFast[0];

            for (int i = 1; i < NumDetectors; i++) 
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
        for (int i = 1; i < NumDetectors; i++) 
        {
            fasttimecalib_fit[i] = new TF1(Form("fasttimecalib_fit%d", i), "gaus(0)", 260, 360);
        }

        // Fit the histograms and extract parameters
        for (int i = 1; i < NumDetectors; i++) 
        {
            fasttimecalib_hists[i]->GetXaxis()->SetRangeUser(260, 360);
            fasttimecalib_hists[i]->Fit(Form("fasttimecalib_fit%d", i), "Q");
            fasttimecalib_hists[i]->GetXaxis()->SetRangeUser(0, 1000);

            mean_fasttimediff[i] = fasttimecalib_fit[i]->GetParameter(1);
            sigma_fasttimediff[i] = fasttimecalib_fit[i]->GetParameter(2);
            printf("Detector %d: Fast time diff Centroid = %f\n", i, mean_fasttimediff[i]);
            printf("Detector %d: Fast time diff sigma = %f\n", i, sigma_fasttimediff[i]);
        }

        // Create a canvas and plot the histograms with fits
        TCanvas* c1 = new TCanvas("c1", "Fast time calibration", 800, 600);
        c1->Divide(2, 2);
        for (int i = 1; i < NumDetectors; i++) 
        {
            c1->cd(i);
            fasttimecalib_hists[i]->SetTitle(Form("Fast time calibration L0 - L%d", i));
            fasttimecalib_hists[i]->SetStats(0);
            fasttimecalib_hists[i]->GetXaxis()->SetTitle("Fast time difference (ns)");
            fasttimecalib_hists[i]->GetYaxis()->SetTitle("Counts (1/ns)");
            fasttimecalib_hists[i]->Draw();
            fasttimecalib_fit[i]->Draw("same");
        }

        // Save canvas and clean up
        c1->Write();
        c1->SaveAs(Form("%s/R%s_fasttimecalib.root", dirPathplots.c_str(), run.c_str()));
        delete c1; 
        
        mean_fasttimediff[0] = 0;
    }
    else 
    {
        printf("Source run: not calculating inter-detector fast time calibration constants ...\n");
        for (int i = 1; i < NumDetectors; i++) 
        {
            mean_fasttimediff[i] = 0;
            sigma_fasttimediff[i] = 0;
        }
    }
    
    return mean_fasttimediff;
}

void align_offset(std::string run,TFile *outputFile, std::string dirPathplots, TTree* tree, TDirectory *dir_banana, TDirectory *dir_aligned, std::vector<Double_t> mean_fasttimediff, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds, uint64_t numEntries) {
    
    printf("\n\nAligning offsets...\n");
    
    dir_aligned->cd();
    // Set branch addresses for each detector
    for (int i = 0; i < NumDetectors; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("energyP", &energyPolaris[0]);
    tree->SetBranchAddress("timeP", &timePolaris[0]);

    TTree* copyTree = new TTree("LaBrDataCopy", "LaBrDataCopy");
    for (int i = 0; i < NumDetectors; i++) 
    {
        copyTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        copyTree->Branch(Form("timeSL%d", i), &timeSlowCp[i], Form("timeSL%d/D", i));
        copyTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
    }
    copyTree->Branch("energyP", &energyPolarisCp[0], "energyP/D");
    copyTree->Branch("timeP", &timePolarisCp[0], "timeP/D");

    double slowTimeScaleFactor = 1.0;
    if (tensOfNanoSeconds) slowTimeScaleFactor = 0.1;
    else if (hundredsOfNanoSeconds) slowTimeScaleFactor = 0.01;

    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);
        for (int i = 0; i < NumDetectors; i++) 
        {
            energySlowCp[i] = energySlow[i];
            timeOffsetBeforeCorrection[i] = 0;
            timeSlowCp[i] = timeSlow[i] * slowTimeScaleFactor;
        }
        timeFastCp[0] = timeFast[0] ;
        for (int i = 1; i < NumDetectors; i++) 
        {
            timeFastCp[i] = (timeFast[i] - mean_fasttimediff[i]) * slowTimeScaleFactor;
        }
        timePolarisCp[0] = timePolaris[0] * slowTimeScaleFactor;
        energyPolarisCp[0] = energyPolaris[0];
        // printf("energySlowCp[%d] = %f, timeSlowCp[%d] = %f, timeFastCp[%d] = %f\n", i, energySlowCp[i], i, timeSlowCp[i], i, timeFastCp[i]);
        
        copyTree->Fill();
    }

    // Initialise the histograms
    for (int i = 0; i < NumDetectors; i++) 
    {
        slowfasttimediff_hists[i] = new TH1D(Form("slowfasttimediff_hists_%d", i), Form("Uncorrected time offset: L%d", i),  600, 0, 600);
        
        // Full
        slowfasttimediff_energy_hists[i] = new TH2D(Form("slowfasttimediff_energy_hists_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  600, 0, 600, 2000, 0, 2000);
        energy_slowfasttimediff_hists[i] = new TH2D(Form("energy_slowfasttimediff_hists_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 600, 0, 600);
        /* timeOffsetCorrectedHistsfirst[i] = new TH1D(Form("timeOffsetCorrectedHists%d", i), Form("Corrected time offset : L%d", i),  600, -300, 300);
        timeEnergyOffsetCorrectedHistsfirst[i] = new TH2D(Form("timeEnergyOffsetCorrectedHists%d", i), Form("Corrected time offset vs energy: L%d", i),  600, -300, 300, 2000, 0, 2000);
         */
        timeOffsetCorrectedHists[i] = new TH1D(Form("timeOffsetCorrectedHists_%d", i), Form("Corrected time offset : L%d", i),  600, -300, 300);
        timeEnergyOffsetCorrectedHists[i] = new TH2D(Form("timeEnergyOffsetCorrectedHists_%d", i), Form("Corrected time offset vs energy: L%d", i),  600, -300, 300, 2000, 0, 2000);
        timefastslowratio_hists[i] = new TH1D(Form("timefastslowratio_hists_%d", i), Form("Time fast/slow ratio: L%d", i),  20, -10, 10);
    }

    for (int i = 0; i < 1; i++) 
    {
        // banana 1 for L0 (Det8)
        slowfasttimediff_energy_hists1[i] = new TH2D(Form("slowfasttimediff_energy_hists1_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  40, 0, 40, 2000, 0, 2000);
        energy_slowfasttimediff_hists1[i] = new TH2D(Form("energy_slowfasttimediff_hists1_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 40, 0, 40);

        // banana 2 for L0 (Det8)
        slowfasttimediff_energy_hists2[i] = new TH2D(Form("slowfasttimediff_energy_hists2_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  16, 40, 56,  2000, 0, 2000);
        energy_slowfasttimediff_hists2[i] = new TH2D(Form("energy_slowfasttimediff_hists2_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 16, 40, 56);

        // banana 3 for L0 (Det8)
        slowfasttimediff_energy_hists3[i] = new TH2D(Form("slowfasttimediff_energy_hists3_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  17, 56, 73, 2000, 0, 2000);
        energy_slowfasttimediff_hists3[i] = new TH2D(Form("energy_slowfasttimediff_hists3_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 17, 56, 73);

        // banana 4 for L0 (Det8)
        slowfasttimediff_energy_hists4[i] = new TH2D(Form("slowfasttimediff_energy_hists4_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  29, 73, 102,  2000, 0, 2000);
        energy_slowfasttimediff_hists4[i] = new TH2D(Form("energy_slowfasttimediff_hists4_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 29, 73, 102);

        // banana 5 for L0 (Det8)
        slowfasttimediff_energy_hists5[i] = new TH2D(Form("slowfasttimediff_energy_hists5_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  198, 102, 300, 2000, 0, 2000);
        energy_slowfasttimediff_hists5[i] = new TH2D(Form("energy_slowfasttimediff_hists5_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 198, 102, 300);

        slowfasttimediff_energy_hists6_2[i] = new TH2D(Form("slowfasttimediff_energy_hists6_2_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  160, 40, 200, 2000, 0, 2000);
        energy_slowfasttimediff_hists6_2[i] = new TH2D(Form("energy_slowfasttimediff_hists6_2_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 160, 40, 200);

        slowfasttimediff_energy_hists6_3[i] = new TH2D(Form("slowfasttimediff_energy_hists6_3_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  160, 40, 200, 2000, 0, 2000);
        energy_slowfasttimediff_hists6_3[i] = new TH2D(Form("energy_slowfasttimediff_hists6_3_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 160, 40, 200);

        /*  slowfasttimediff_energy_hists_again[i] = new TH2D(Form("slowfasttimediff_energy_hists_again_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  30, 10, 40, 2000, 0, 2000);
        energy_slowfasttimediff_hists_again[i] = new TH2D(Form("energy_slowfasttimediff_hists_again_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 30, 10, 40); */
    }
    for (int i = 1; i < NumDetectors; i++)
    {
     //________________________________________________________________________________________________________
        // banana 6 for L2 (Det3)
        slowfasttimediff_energy_hists6[i] = new TH2D(Form("slowfasttimediff_energy_hists6_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  18, 0, 18, 2000, 0, 2000);
        energy_slowfasttimediff_hists6[i] = new TH2D(Form("energy_slowfasttimediff_hists6_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 18, 0, 18);

        // banana 7 for L2 (Det3)
        slowfasttimediff_energy_hists7[i] = new TH2D(Form("slowfasttimediff_energy_hists7_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  14, 18, 32, 2000, 0, 2000);
        energy_slowfasttimediff_hists7[i] = new TH2D(Form("energy_slowfasttimediff_hists7_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 14, 18, 32);

        // banana 8 for L2 (Det3)
        slowfasttimediff_energy_hists8[i] = new TH2D(Form("slowfasttimediff_energy_hists8_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  24, 32, 56, 2000, 0, 2000);
        energy_slowfasttimediff_hists8[i] = new TH2D(Form("energy_slowfasttimediff_hists8_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 24, 32, 56);

        // banana 9 for L2 (Det3)
        slowfasttimediff_energy_hists9[i] = new TH2D(Form("slowfasttimediff_energy_hists9_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  18, 56, 74, 2000, 0, 2000);
        energy_slowfasttimediff_hists9[i] = new TH2D(Form("energy_slowfasttimediff_hists9_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 18, 56, 74);

        // banana 10 for L2 (Det3)
        slowfasttimediff_energy_hists10[i] = new TH2D(Form("slowfasttimediff_energy_hists10_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  28, 74, 102, 2000, 0, 2000);
        energy_slowfasttimediff_hists10[i] = new TH2D(Form("energy_slowfasttimediff_hists10_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 28, 74, 102);

        // banana 11 for L2 (Det3)
        slowfasttimediff_energy_hists11[i] = new TH2D(Form("slowfasttimediff_energy_hists11_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  198, 102, 300, 2000, 0, 2000);
        energy_slowfasttimediff_hists11[i] = new TH2D(Form("energy_slowfasttimediff_hists11_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 198, 102, 300);
    }
        
    // Initialize new TTree
    TTree* offsetTree = new TTree("LaBrDataOffset", "LaBrDataOffset"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
    {
        offsetTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        offsetTree->Branch(Form("timeSL%d", i), &timeSlowCp[i], Form("timeSL%d/D", i));
        offsetTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        offsetTree->Branch(Form("timeOffsetL%d", i), &timeOffsetBeforeCorrection[i], Form("timeOffsetL%d/D", i)); //new Branch
    }
    offsetTree->Branch("energyP", &energyPolarisCp[0], "energyP/D");
    offsetTree->Branch("timeP", &timePolarisCp[0], "timeP/D");

    // Initialize new TTree
    TTree* correctedTree = new TTree("LaBrDataCorrected", "LaBrDataCorrected"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
    {
        correctedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrected[i], Form("timeSCorrectedL%d/D", i)); //new Branch
        correctedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
    }
    correctedTree->Branch("energyP", &energyPolarisCp[0], "energyP/D");
    correctedTree->Branch("timeP", &timePolarisCp[0], "timeP/D");

    // Initialize new TTree
   /*  TTree* correctedcorrectedTree = new TTree("LaBrDataCorrected2", "LaBrDataCorrected2"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
    {
        correctedcorrectedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedcorrectedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        correctedcorrectedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrected[i], Form("timeSCorrectedL%d/D", i)); //new Branch
        correctedcorrectedTree->Branch(Form("timeOffsetL%d", i), &timeOffsetBeforeCorrectionCp[i], Form("timeOffsetL%d/D", i)); //new Branch

    }
    correctedcorrectedTree->Branch("fastEPOLARIS", &energyPolarisCp[0], "fastEPOLARIS/D");
    correctedcorrectedTree->Branch("fastTPOLARIS", &timePolarisCp[0], "fastTPOLARIS/D");

    // Initialize new TTree
    TTree* correctedcorrectedcorrectedTree = new TTree("LaBrDataCorrected3", "LaBrDataCorrected3"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
    {
        correctedcorrectedcorrectedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedcorrectedcorrectedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        correctedcorrectedcorrectedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrectedCp[i], Form("timeSCorrectedL%d/D", i)); //new Branch
    }
    correctedcorrectedcorrectedTree->Branch("fastEPOLARIS", &energyPolarisCp[0], "fastEPOLARIS/D");
    correctedcorrectedcorrectedTree->Branch("fastTPOLARIS", &timePolarisCp[0], "fastTPOLARIS/D"); */


    // Variables for aligned tree data
    Int_t detectorID;
    Double_t alignedTimeFast;
    Double_t alignedEnergyFast;
    Double_t alignedTimeSlow;
    Double_t alignedEnergySlow;
    Double_t timeGlobal;
    
    struct EntryData 
    {
        Int_t detectorID;
        Double_t timeGlobal;
        Double_t alignedTimeFast;
        Double_t alignedEnergyFast;
        Double_t alignedTimeSlow;
        Double_t alignedEnergySlow;
    };
    std::vector<EntryData> entryDataVec;

    // _________________________________________________________________________________________________________________________________________________________
    
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {
            if (entry==0) printf("\n\nCalculating time offset for detector %d...\n", i);
            copyTree->GetEntry(entry);

            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;
            Double_t deltaTime = 0;
            Double_t energy = 0;
            uint64_t nextEntry = entry;

            while (nextEntry < numEntries) 
            {
                copyTree->GetEntry(nextEntry);
                // printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCp[i], i, energySlowCp[i]);       
                       
                if (timeSlowCp[i] == 0) 
                {
                    nextEntry++; // Move to the next entry if timeSlowCp[i] is zero
                    continue;
                }
                nextTimeSlow = timeSlowCp[i];

                deltaTime = (nextTimeSlow - refTimeFast); 
                energy = energySlowCp[i];
                if (deltaTime > 600) break;
                else
                {
                    timeOffsetBeforeCorrection[i] = static_cast<Double_t>(static_cast<int>(deltaTime));
                    // printf("Detector %d: deltaTime = %f\n", i, deltaTime);
                    
                    if (nextEnergySlow == 0) continue;
                    slowfasttimediff_hists[i]->Fill(deltaTime);
                    slowfasttimediff_energy_hists[i]->Fill(deltaTime, energy);
                    energy_slowfasttimediff_hists[i]->Fill(energy, deltaTime);

                    if (i ==0) 
                    {
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
                        slowfasttimediff_energy_hists6_2[i]->Fill(deltaTime, energy);
                        energy_slowfasttimediff_hists6_2[i]->Fill(energy, deltaTime);
                        slowfasttimediff_energy_hists6_3[i]->Fill(deltaTime, energy);
                        energy_slowfasttimediff_hists6_3[i]->Fill(energy, deltaTime);
                    }
                    else if (i == 1)
                    {
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
                    }
                    break;
                }
                
            }
            offsetTree->Fill();
        }
    }

    printf("\n\nFinished calculating time offsets...\n");

    for (int i = 0; i < NumDetectors; i++) 
    {    
        TCanvas* c1 = new TCanvas(Form("c1_%d", i), Form("Uncorrected time offset: L%d", i), 800, 600);
        slowfasttimediff_hists[i]->GetXaxis()->SetRangeUser(0, 600);
        slowfasttimediff_hists[i]->SetStats(0);
        slowfasttimediff_hists[i]->SetFillColor(kBlue);
        slowfasttimediff_hists[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        slowfasttimediff_hists[i]->GetYaxis()->SetTitle("Counts (1/ns)");
        slowfasttimediff_hists[i]->Draw();
        c1->SaveAs(Form("%s/R%s_time_offset_hist_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        c1->Write(Form("time_offset_hist_L%d", i));
        delete c1;
              
        TCanvas* c18 = new TCanvas(Form("c18_%d", i), Form("Uncorrected time offset energy matrix: L%d", i), 800, 600);
        slowfasttimediff_energy_hists[i]->SetStats(0);
        slowfasttimediff_energy_hists[i]->SetFillColor(kBlue);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        slowfasttimediff_energy_hists[i]->GetYaxis()->SetTitle("Energy (keV)");
        slowfasttimediff_energy_hists[i]->Draw();
        c18->SaveAs(Form("%s/R%s_time_offset_energy_hist_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        c18->Write(Form("time_offset_energy_hist_L%d", i));
        delete c18;

        dir_banana->cd();

        if (i==0)
        {   
            printf("\n\nPerforming the banana fit for detector %d...\n", i);
            banana1fit[i] = banana1(run,outputFile, i, slowfasttimediff_energy_hists1[i], energy_slowfasttimediff_hists1[i], dirPathplots);
            banana2fit[i] = banana2(run,outputFile, i, slowfasttimediff_energy_hists2[i], energy_slowfasttimediff_hists2[i], dirPathplots);
            banana3fit[i] = banana3(run,outputFile, i, slowfasttimediff_energy_hists3[i], energy_slowfasttimediff_hists3[i], dirPathplots);
            banana4fit[i] = banana4(run,outputFile, i, slowfasttimediff_energy_hists4[i], energy_slowfasttimediff_hists4[i], dirPathplots);
            banana5_1fit[i] = banana5_1(run,outputFile, i, slowfasttimediff_energy_hists5[i], energy_slowfasttimediff_hists5[i], dirPathplots);
            banana5_2fit[i] = banana5_2(run,outputFile, i, slowfasttimediff_energy_hists5[i], energy_slowfasttimediff_hists5[i], dirPathplots);
            banana6_2fit[i] = banana6_2(run,outputFile, i, slowfasttimediff_energy_hists6_2[i], energy_slowfasttimediff_hists6_2[i], dirPathplots);
            banana6_3fit[i] = banana6_3(run,outputFile, i, slowfasttimediff_energy_hists6_3[i], energy_slowfasttimediff_hists6_3[i], dirPathplots);

            TCanvas * c2 = new TCanvas("c2", "banana", 800, 600);
            energy_slowfasttimediff_hists[i]->Draw();
            banana1fit[i]->SetLineColor(kRed);
            banana2fit[i]->SetLineColor(kBlue);
            banana3fit[i]->SetLineColor(kGreen);
            banana4fit[i]->SetLineColor(kMagenta);
            banana5_1fit[i]->SetLineColor(kYellow);
            banana5_2fit[i]->SetLineColor(kYellow);
            banana6_2fit[i]->SetLineColor(kWhite);
            banana6_3fit[i]->SetLineColor(kWhite);
            banana1fit[i]->Draw("same");
            banana2fit[i]->Draw("same");
            banana3fit[i]->Draw("same");
            banana4fit[i]->Draw("same");
            banana5_1fit[i]->Draw("same");
            banana5_2fit[i]->Draw("same");
            banana6_2fit[i]->Draw("same");
            banana6_3fit[i]->Draw("same");
            TLegend *leg1 = new TLegend(0.1,0.7,0.48,0.9);
            leg1->AddEntry(banana1fit[i],"banana1fit","l");
            leg1->AddEntry(banana2fit[i],"banana2fit","l");
            leg1->AddEntry(banana3fit[i],"banana3fit","l");
            leg1->AddEntry(banana4fit[i],"banana4fit","l");
            leg1->AddEntry(banana5_1fit[i],"banana5_1fit","l");
            leg1->AddEntry(banana5_2fit[i],"banana5_2fit","l");
            leg1->AddEntry(banana6_2fit[i],"banana6_2fit","l");
            leg1->AddEntry(banana6_3fit[i],"banana6_3fit","l");
            leg1->Draw();
            c2->SaveAs(Form("%s/R%s_bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), run.c_str(), i));
            c2->Write(Form("bananas_2Dfitfunction_L%d", i));
            delete c2;
        }
        else 
        {   
            printf("\n\nPerforming the banana fit for detector %d...\n", i);
            banana6fit[i] = banana6(run,outputFile, i, slowfasttimediff_energy_hists6[i], energy_slowfasttimediff_hists6[i], dirPathplots);
            banana7_1fit[i] = banana7_1(run,outputFile, i, slowfasttimediff_energy_hists7[i], energy_slowfasttimediff_hists7[i], dirPathplots);
            banana7_2fit[i] = banana7_2(run,outputFile, i, slowfasttimediff_energy_hists7[i], energy_slowfasttimediff_hists7[i], dirPathplots);
            banana8fit[i] = banana8(run,outputFile, i, slowfasttimediff_energy_hists8[i], energy_slowfasttimediff_hists8[i], dirPathplots);
            banana9fit[i] = banana9(run,outputFile, i, slowfasttimediff_energy_hists9[i], energy_slowfasttimediff_hists9[i], dirPathplots);
            banana10fit[i] = banana10(run,outputFile, i, slowfasttimediff_energy_hists10[i], energy_slowfasttimediff_hists10[i], dirPathplots);
            banana11fit[i] = banana11(run,outputFile, i, slowfasttimediff_energy_hists11[i], energy_slowfasttimediff_hists11[i], dirPathplots);
              
            TCanvas * c3 = new TCanvas("c3", "banana", 800, 600);
            energy_slowfasttimediff_hists[i]->Draw();
            banana6fit[i]->SetLineColor(kRed);
            banana7_1fit[i]->SetLineColor(kBlue);
            banana7_2fit[i]->SetLineColor(kBlue);
            banana8fit[i]->SetLineColor(kGreen);
            banana9fit[i]->SetLineColor(kMagenta);
            banana10fit[i]->SetLineColor(kYellow);
            banana11fit[i]->SetLineColor(kOrange);
            banana8fit[i]->Draw("same");
            banana6fit[i]->Draw("same");
            banana7_1fit[i]->Draw("same");
            banana7_2fit[i]->Draw("same");
            banana9fit[i]->Draw("same");
            banana10fit[i]->Draw("same");
            banana11fit[i]->Draw("same");
            TLegend *leg2 = new TLegend(0.1,0.7,0.48,0.9);
            leg2->AddEntry(banana6fit[i],"banana6fit","l");
            leg2->AddEntry(banana7_1fit[i],"banana7_1fit","l");
            leg2->AddEntry(banana7_2fit[i],"banana7_2fit","l");
            leg2->AddEntry(banana8fit[i],"banana8fit","l");
            leg2->AddEntry(banana9fit[i],"banana9fit","l");
            leg2->AddEntry(banana10fit[i],"banana10fit","l");
            leg2->AddEntry(banana11fit[i],"banana11fit","l");
            leg2->Draw();
            c3->SaveAs(Form("%s/R%s_bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), run.c_str(), i));
            c3->Write(Form("bananas_2Dfitfunction_L%d", i));
            delete c3;
        }
    }
    dir_aligned->cd();

    // ________________________________________________________________
    Double_t minEnergy = 1, maxEnergy = 2000.0, energyStep = 1.0;

    printf("\nPrecomputing the fit values for the time offset correction...\n");
    
    for (Double_t e = minEnergy; e <= maxEnergy; e += energyStep)
    {
        for (int i = 0; i < NumDetectors; i++) 
        {
            if (i==0)
            {
                precomputedValues1[std::make_pair(i, e)] = banana1fit[i]->Eval(e);
                precomputedValues2[std::make_pair(i, e)] = banana2fit[i]->Eval(e);
                precomputedValues3[std::make_pair(i, e)] = banana3fit[i]->Eval(e);
                precomputedValues4[std::make_pair(i, e)] = banana4fit[i]->Eval(e);
                precomputedValues5_1[std::make_pair(i, e)] = banana5_1fit[i]->Eval(e);
                precomputedValues5_2[std::make_pair(i, e)] = banana5_2fit[i]->Eval(e);
                precomputedValues6_2[std::make_pair(i, e)] = banana6_2fit[i]->Eval(e);
                precomputedValues6_3[std::make_pair(i, e)] = banana6_3fit[i]->Eval(e);
            }
            else
            {
                precomputedValues6[std::make_pair(i, e)] = banana6fit[i]->Eval(e);
                precomputedValues7_1[std::make_pair(i, e)] = banana7_1fit[i]->Eval(e);
                precomputedValues7_2[std::make_pair(i, e)] = banana7_2fit[i]->Eval(e);
                precomputedValues8[std::make_pair(i, e)] = banana8fit[i]->Eval(e);
                precomputedValues9[std::make_pair(i, e)] = banana9fit[i]->Eval(e);
                precomputedValues10[std::make_pair(i, e)] = banana10fit[i]->Eval(e);
                precomputedValues11[std::make_pair(i, e)] = banana11fit[i]->Eval(e);
            }
        }
    }

    // ________________________________________________________________
    
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {  if (entry ==0) printf("\n\nCalculating the time offset from the fit precomputed values for detector %d...\n", i);       
            offsetTree->GetEntry(entry);     
            timeSlowCorrected[i] = 0;
            Double_t t = timeOffsetBeforeCorrection[i];

            if (t == 0) continue;

            Double_t e = static_cast<Double_t>(static_cast<int>(energySlowCp[i]));
            Double_t timeOffset = 0.0;

            if (i==0)
            {
                if ((t <= 40))
                {
                    timeOffset = precomputedValues1[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 1 = %f\n", i, timeOffset);
                } 
                else if (t>40 && t<=56) 
                {
                    timeOffset = precomputedValues2[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 2 = %f\n", i, timeOffset);
                }
                else if (t>56 && t<=73) 
                {
                    timeOffset = precomputedValues3[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 3 = %f\n", i, timeOffset);
                } 
                else if (t>73 && t<=102) 
                {
                    timeOffset = precomputedValues4[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 4 = %f\n", i, timeOffset);
                }
                /* else if (t>102 && e< 150) 
                {
                    timeOffset = precomputedValues5_1[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
                }  */
                else if (t>102 && e> 180) 
                {
                    timeOffset = precomputedValues5_2[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
                } 
                else if (t > 60 && e< 150) 
                {
                    timeOffset = precomputedValues6_2[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
                    if (t > 100)timeOffset = precomputedValues6_3[std::make_pair(i, e)];
                } 
                
                else continue;
            }
            else if (i == 1)
            {
                if (t <= 18)
                {
                    timeOffset = precomputedValues6[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 6 = %f\n", i, timeOffset);
                }
                else if (t<=32 && t>18 && e < 900)
                {
                    timeOffset = precomputedValues7_1[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 7 = %f\n", i, timeOffset);
                }
                else if (t<=32 && t>18 && e > 900)
                {
                    timeOffset = precomputedValues7_2[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 7 = %f\n", i, timeOffset);
                }
                else if (t<=56 && t>32)
                {
                    timeOffset = precomputedValues8[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 8 = %f\n", i, timeOffset);
                }
                else if (t<=74 && t>56)
                {
                    timeOffset = precomputedValues9[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 9 = %f\n", i, timeOffset);
                }
                else if (t<=102 && t>74)
                {
                    timeOffset = precomputedValues10[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 10 = %f\n", i, timeOffset);
                }
                else if (t>102)
                {
                    timeOffset = precomputedValues11[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 11 = %f\n", i, timeOffset);
                }
                else continue;
            }
            
            // Calculate timeSlowCorrected and fill correctedTree using timeOffset
            if (timeOffset > 0) 
            {
                timeSlowCorrected[i] = timeSlowCp[i] - timeOffset;
                //printf("Detector %d: Time slow corrected = %f\n", i, timeSlowCorrected[i]);
                //printf("Detector %d: Time slow = %f\n", i, timeSlowCp[i]/1e14);
                if (timeSlowCorrected[i]/1e14 < 1) 
                {
                    continue;
                }
                correctedTree->Fill();
            }
        }
    }
/*     for (int i = 0; i < NumDetectors; i++) 
    { 
        for (unsigned long long entry = 0; entry < numEntries; entry++)
        {
            correctedTree->GetEntry(entry);
            if (i == 1) 
            {
                timeFastCp[i] = timeFastCp[0];
                timeSlowCorrected[i] = timeSlowCorrected[0];
                energySlowCp[i] = energySlowCp[0];
                timeOffsetBeforeCorrection[i] = 0;
            }
            else 
            {
                if (timeFastCp[i] == 0) continue;
                Double_t refTimeFast = timeFastCp[i];
                Double_t nextTimeSlow = 0;
                Double_t deltaTime = 0;
                Double_t energy = 0;

                uint64_t nextEntry = entry;

                while (nextEntry < numEntries) 
                {
                    correctedTree->GetEntry(nextEntry);
                    // printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCorrected[i], i, energySlowCp[i]);       
                        
                    if (timeSlowCorrected[i] == 0) 
                    {
                        nextEntry++; // Move to the next entry if timeSlowCp[i] is zero
                        continue;
                    }
                    nextTimeSlow = timeSlowCorrected[i];

                    deltaTime = (nextTimeSlow - refTimeFast); 
                    energy = energySlowCp[i];
                    if (deltaTime > 600) break;
                    else
                    {
                        timeOffsetBeforeCorrectionCp[i] = static_cast<Double_t>(static_cast<int>(deltaTime));
                        slowfasttimediff_energy_hists_again[i]->Fill(deltaTime, energy);
                        energy_slowfasttimediff_hists_again[i]->Fill(energy, deltaTime);
                        timeOffsetCorrectedHistsfirst[i]->Fill(deltaTime);
                        timeEnergyOffsetCorrectedHistsfirst[i]->Fill(deltaTime, energy);
                        break;
                    }
                }    
            }
            correctedcorrectedTree->Fill();
        }
    }
    delete correctedTree;
    delete offsetTree;

    // Draw the histograms
    for (int i = 0; i < NumDetectors; i++) 
    {
        TCanvas * c44 = new TCanvas("c44", "CorrectedTimeOffsetFirst", 800, 600);
        timeOffsetCorrectedHistsfirst[i]->GetXaxis()->SetRangeUser(-100, 100);
        timeOffsetCorrectedHistsfirst[i]->SetFillColor(kBlue);
        timeOffsetCorrectedHistsfirst[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        timeOffsetCorrectedHistsfirst[i]->GetYaxis()->SetTitle("Counts (1/ns)");
        timeOffsetCorrectedHistsfirst[i]->Draw();
        c44->SaveAs(Form("%s/R%s_time_offset_corrected_first_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        delete c44;

        TCanvas * c47 = new TCanvas("c47", "CorrectedTimeOffsetEnergyFirst", 800, 600);
        timeEnergyOffsetCorrectedHistsfirst[i]->GetXaxis()->SetRangeUser(-100, 100);
        timeEnergyOffsetCorrectedHistsfirst[i]->GetYaxis()->SetRangeUser(0, 600);
        timeEnergyOffsetCorrectedHistsfirst[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        timeEnergyOffsetCorrectedHistsfirst[i]->GetYaxis()->SetTitle("Energy (keV)");
        timeEnergyOffsetCorrectedHistsfirst[i]->Draw("colz");
        c47->SaveAs(Form("%s/R%s_time_offset_energy_corrected_first_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        delete c47;
    }

    dir_banana->cd();

    for (int i = 0; i < 1; i++) 
    {
        bananafinalfit[i] = bananafinal(run,outputFile, i, slowfasttimediff_energy_hists_again[i], energy_slowfasttimediff_hists_again[i], dirPathplots);
        TCanvas * c8 = new TCanvas("c8", "banana", 800, 600);
        energy_slowfasttimediff_hists_again[i]->Draw();
        bananafinalfit[i]->SetLineColor(kRed);
        bananafinalfit[i]->Draw("same");
        c8->SaveAs(Form("%s/R%s_bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        c8->Write(Form("bananas_2Dfitfunction_L%d", i));
        delete c8;
    }
    dir_aligned->cd();

    for (Double_t e = minEnergy; e <= maxEnergy; e += energyStep)
    {
        for (int i = 0; i < 1; i++) 
        {
            if (i==0)
            {
                precomputedValues_again[std::make_pair(i, e)] = bananafinalfit[i]->Eval(e);
            }
            else continue;
        }
    }

    // ________________________________________________________________
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {  if (entry ==0) printf("\n\nCalculating the final corrected time offset from the fit precomputed values for detector %d...\n", i);       
            
            correctedcorrectedTree->GetEntry(entry);    
            timeSlowCorrectedCp[i] = 0; 
            Double_t t = timeOffsetBeforeCorrectionCp[i];

            Double_t e = static_cast<Double_t>(static_cast<int>(energySlowCp[i]));
            Double_t timeOffset = 0.0;

            if (i==0)
            {
                if ((t >=10))
                {
                    timeOffset = precomputedValues_again[std::make_pair(i, e)];
                    //printf("Detector %d: Time Offset function 1 = %f\n", i, timeOffset);
                } 
                
            }
            else 
            {
                timeOffset = 0;
                continue;
            }

            timeSlowCorrectedCp[i] = timeSlowCorrected[i] - timeOffset;
            if (timeSlowCorrectedCp[i] < 1e14) continue;
            correctedcorrectedcorrectedTree->Fill();
        }
    }
    delete correctedcorrectedTree; */

    /* // ________________________________________________________________
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < 1; i++) 
        {
            if (entry == 0) printf("\n\nCalculating the corrected time offset for detector %d...\n", i);
            
            correctedcorrectedcorrectedTree->GetEntry(entry);
            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;  // Initialize nextTimeSlow
            Double_t deltaTime = 0; // Initialize deltaTime
            Double_t nextEnergySlow = 0; // Initialize nextEnergySlow

            uint64_t nextEntry = entry;

            while (nextEntry < numEntries) 
            {
                correctedcorrectedcorrectedTree->GetEntry(nextEntry);

                if (timeSlowCorrectedCp[i] == 0) 
                {
                    nextEntry++; // Move to the next entry if timeSlowCp[i] is zero
                    continue;
                }
                nextEnergySlow = energySlowCp[i];
                nextTimeSlow = timeSlowCorrectedCp[i];
                //printf("\nDetector %d: Fast Time = %f, Slow Time = %f\n", i, refTimeFast, nextTimeSlow);
                deltaTime = (nextTimeSlow - refTimeFast);
                if (deltaTime > 1000.0) break; // Stop searching if time difference exceeds 1 microsecond
                //printf("\nDetector %d: Corrected Time Offset = %f\n", i, deltaTime);
                timeOffsetCorrectedHists[i]->Fill(deltaTime);
                timeEnergyOffsetCorrectedHists[i]->Fill(deltaTime, nextEnergySlow);
                break;
            }
        }
    } */

   // ________________________________________________________________
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {
            if (entry == 0) printf("\n\nCalculating the corrected time offset for detector %d...\n", i);
            
            correctedTree->GetEntry(entry);
            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;  // Initialize nextTimeSlow
            Double_t deltaTime = 0; // Initialize deltaTime
            Double_t nextEnergySlow = 0; // Initialize nextEnergySlow

            uint64_t nextEntry = entry;

            while (nextEntry < numEntries) 
            {
                correctedTree->GetEntry(nextEntry);

                if (timeSlowCorrected[i] == 0) 
                {
                    nextEntry++; // Move to the next entry if timeSlowCp[i] is zero
                    continue;
                }
                nextEnergySlow = energySlowCp[i];
                nextTimeSlow = timeSlowCorrected[i];
                //printf("\nDetector %d: Fast Time = %f, Slow Time = %f\n", i, refTimeFast, nextTimeSlow);
                deltaTime = (nextTimeSlow - refTimeFast);
                if (deltaTime > 1000.0) break; // Stop searching if time difference exceeds 1 microsecond
                //printf("\nDetector %d: Corrected Time Offset = %f\n", i, deltaTime);
                if (nextEnergySlow == 0) continue;
                timeOffsetCorrectedHists[i]->Fill(deltaTime);
                timeEnergyOffsetCorrectedHists[i]->Fill(deltaTime, nextEnergySlow);
                break;
            }
        }
    }

    for (int i = 0; i < NumDetectors; i++) 
    {

        timeOffsetCorrectedHists[i] ->Write();
        timeEnergyOffsetCorrectedHists[i] ->Write();

        TCanvas * c5 = new TCanvas("c5", "CorrectedTimeOffsetFirst", 800, 600);
        timeOffsetCorrectedHists[i]->GetXaxis()->SetRangeUser(-100, 100);
        timeOffsetCorrectedHists[i]->SetFillColor(kBlue);
        timeOffsetCorrectedHists[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        timeOffsetCorrectedHists[i]->GetYaxis()->SetTitle("Counts (1/ns)");
        timeOffsetCorrectedHists[i]->Draw();
        c5->SaveAs(Form("%s/R%s_time_offset_corrected_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        delete c5;

        TCanvas * c6 = new TCanvas("c6", "CorrectedTimeOffsetEnergyFirst", 800, 600);
        timeEnergyOffsetCorrectedHists[i]->GetXaxis()->SetRangeUser(-100, 100);
        timeEnergyOffsetCorrectedHists[i]->GetYaxis()->SetRangeUser(0, 2000);
        timeEnergyOffsetCorrectedHists[i]->SetFillColor(kBlue);
        timeEnergyOffsetCorrectedHists[i]->GetXaxis()->SetTitle("Time slow - fast (ns)");
        timeEnergyOffsetCorrectedHists[i]->GetYaxis()->SetTitle("Energy (keV)");
        timeEnergyOffsetCorrectedHists[i]->Draw("colz");
        c6->SaveAs(Form("%s/R%s_time_energy_offset_corrected_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        delete c6;

        // ________________________________________________________________
        // Write the histograms to the output file
        slowfasttimediff_hists[i]->Write();
        slowfasttimediff_energy_hists[i]->Write();
        energy_slowfasttimediff_hists[i]->Write();

    }

    // ________________________________________________________________

    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {
            if (entry == 0) printf("\n\nAligning the data for detector %d...\n", i);
            correctedTree->GetEntry(entry);
            if ((timePolarisCp[0] !=0) && (energyPolarisCp[0]!=0))
            {
                detectorID = POLARISIndex;
                timeGlobal = timePolarisCp[0];
                alignedTimeFast = timePolarisCp[0];
                alignedEnergyFast = energyPolarisCp[0];
                alignedTimeSlow = 0;
                alignedEnergySlow = 0;
                // printf("\nPOLARIS Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i < RFIndex)
            {
                detectorID = i;
                if ((timeFastCp[i] != 0) && (timeSlowCorrected[i] != 0)) 
                {
                    // The ratio in the print below should be 1
                    //printf("\nError: Both fast and slow times are non-zero for detector %d\n, timeFastCp = %f | , timeSlowCorrected + timeOffsets = %f\n, ratio of timeFastCp to timeSlowCp + timeOffsets = %f\n", i, timeFastCp[i], timeSlowCp[i] + timeOffsets[i], (timeFastCp[i]/(timeSlowCp[i] + timeOffsets[i])));
                    timeGlobal = timeFastCp[i];
                    alignedTimeFast = timeFastCp[i];
                    alignedEnergyFast = 0;
                    alignedTimeSlow = timeSlowCorrected[i];
                    alignedEnergySlow = energySlowCp[i];
                    //alignedTree1->Fill();
                }
                else 
                {
                    // The slow energy in the print below should be zero
                    //printf("\nError: Both fast and slow times are zero for detector %d\n, Slow Energy = %f | , timeFastCp = %f | , timeSlowCorrected = %f\n", i, energySlowCp[i], timeFastCp[i], timeSlowCorrected[i]);
                    continue;
                }
                double ratio = alignedTimeFast/alignedTimeSlow;
                if (ratio != 0) timefastslowratio_hists[i]->Fill(alignedTimeFast/alignedTimeSlow);
            }

            // Prepare the EntryData struct and fill it
            EntryData entryData;
            entryData.detectorID = detectorID;
            entryData.timeGlobal = timeGlobal;
            entryData.alignedTimeFast = alignedTimeFast;
            entryData.alignedEnergyFast = alignedEnergyFast;
            entryData.alignedTimeSlow = alignedTimeSlow;
            entryData.alignedEnergySlow = alignedEnergySlow;

            entryDataVec.push_back(entryData);
        }
    }    

    for (int i = 0; i < NumDetectors; i++) 
    {
        timefastslowratio_hists[i]->Write();

        TCanvas * c7 = new TCanvas("c7", "TimeFastSlowRatio", 800, 600);
        timefastslowratio_hists[i]->GetXaxis()->SetRangeUser(0, 2);
        timefastslowratio_hists[i]->SetFillColor(kBlue);
        timefastslowratio_hists[i]->GetXaxis()->SetTitle("Time fast /slow");
        timefastslowratio_hists[i]->GetYaxis()->SetTitle("Counts (1/ns)");
        timefastslowratio_hists[i]->Draw();
        c7->SaveAs(Form("%s/R%s_time_fast_slow_ratio_L%d.root", dirPathplots.c_str(), run.c_str(), i));
        delete c7;
    }

    // Sort the entryDataVec based on timeGlobal in ascending order
    std::sort(entryDataVec.begin(), entryDataVec.end(), [](const EntryData& a, const EntryData& b) {return a.alignedTimeFast < b.alignedTimeFast;});

    // Set the branches for the aligned data in the new copyTree
    TTree *alignedTree1 = new TTree("AlignedData", "AlignedData");
    alignedTree1->Branch("detectorID", &detectorID, "detectorID/I");
    //alignedTree1->Branch("globalTime", &timeGlobal, "globalTime/D");
    alignedTree1->Branch("timeF", &alignedTimeFast, "timeF/D");
    alignedTree1->Branch("energyF", &alignedEnergyFast, "energyF/D");
    alignedTree1->Branch("timeS", &alignedTimeSlow, "timeS/D");
    alignedTree1->Branch("energyS", &alignedEnergySlow, "energyS/D");
    alignedTree1->Reset();

    // Iterate through the sorted vector and fill the alignedTree1
    printf("\nFilling the aligned tree in order of timeGlobal...\n");
    for (const EntryData& entryData : entryDataVec) 
    {
        detectorID = entryData.detectorID;
        //timeGlobal = entryData.timeGlobal;
        alignedTimeFast = entryData.alignedTimeFast;
        alignedEnergyFast = entryData.alignedEnergyFast;
        alignedTimeSlow = entryData.alignedTimeSlow;
        alignedEnergySlow = entryData.alignedEnergySlow;
        alignedTree1->Fill();
    }

    uint64_t nEntries = (alignedTree1->GetEntries());

    uint64_t lastEntry;
    uint64_t firstEntry;
    Double_t lastTime;
    Double_t firstTime;

    for (uint64_t entry = 0; entry < nEntries; entry++) 
    {
        alignedTree1->GetEntry(entry);
        if (alignedTimeFast <= 0 )continue;
        firstTime = alignedTimeFast;
        break;
    }
   
    lastEntry = nEntries - 1;
    for (uint64_t entry = lastEntry; entry >= 0; entry--) 
    {
        alignedTree1->GetEntry(entry);
        lastTime = alignedTimeFast;
        break;
    }
    
    alignedTree1->Write();


    Double_t runLength = (lastTime - firstTime)*1e-9/60;
    printf("\n First entry alignedTimeFast = %f\n", firstTime);
    printf("\n Last entry alignedTimeFast = %f\n", lastTime);
    printf("\nRun Length (According to alignedTimeFast variable. Used to check if it's in sequential time) = %f\n", runLength);

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
    int run1 = std::atoi(runStr.c_str());
    std::string outputFileName = dirPath + "R" + std::to_string(run1) + "_aligned.root";

    std::cout << "Run: " << run1 << std::endl;
    std::cout << "Directory Path: " << dirPath << std::endl;
    std::cout << "Output File Name: " << outputFileName << std::endl;
    std::cout << "_______________________________________________" << std::endl;

    std::string run = std::to_string(run1);

    std::string dirPathplots = dirPath + "plots/sort2labr/";
    // check to see if dirPath contains a folder called plots. if not, create it
    if (dirPath.find("sort2labr_plots/") == std::string::npos)
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
    timeFast.reserve(NumDetectors);
    energyPolaris.reserve(1);
    timePolaris.reserve(1);
    
    // Open the ROOT file
    TFile*file = new TFile(inputFileName.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file->Get("LaBrData"));
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    uint64_t numEntries = tree->GetEntries();

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
    fasttime_calibration_constants(run, beam, outputFile, dirPathplots, dir_fasttimecalib, tree, tensOfNanoSeconds, hundredsOfNanoSeconds,numEntries); // this function is only applied to beam data
    align_offset(run,outputFile, dirPathplots, tree, dir_banana, dir_aligned, mean_fasttimediff, tensOfNanoSeconds, hundredsOfNanoSeconds, numEntries);

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
