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

const int NumDetectors = 2;
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
std::vector<Double_t> timeFastCp(NumDetectors+1); 
std::vector<Double_t> fastEnergyPOLARISCp(1);
std::vector<Double_t> fastTimePOLARISCp(1);
std::vector<Double_t> timeOffsetBeforeCorrection(NumDetectors);

// Variables for aligned tree data
Int_t detectorID;
std::vector<Double_t> alignedTimeFast(1);
std::vector<Double_t> alignedEnergyFast(1);
std::vector<Double_t> alignedTimeSlow(1);
std::vector<Double_t> alignedEnergySlow(1);
std::vector<Double_t> timeGlobal(1);
std::vector<TH1D*> slowfasttimediff_hists(NumDetectors);
std::vector<TH1D*> timeOffsetCorrectedHists(NumDetectors);
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
// std::vector<TH2D*> slowfasttimediff_energy_hists7(NumDetectors);
// std::vector<TH2D*> energy_slowfasttimediff_hists7(NumDetectors);
std::vector<TH2D*> timeEnergyOffsetCorrectedHists(NumDetectors);

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

// TF1* banana7(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots);
// TF1* banana7fit[NumDetectors];


// END BANANA FIT FUNCTIONS //

// banana1 good fit
TF1* banana1(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots)
{
    std::cout <<"Fitting banana1 (2000-8000 keV) for CFD - detector L \n" << i << std::endl;

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

        slice_charge.push_back((bin+binmax)/2.+2000);
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
    gr1->GetYaxis()->SetRangeUser(0, 500);
    gr1->GetXaxis()->SetRangeUser(2000, 8000);
    TF1* fit1 = new TF1("fit1", "pol0(0)", 2000, 8000);  
    
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
    std::cout <<"Fitting banana2 (1000-8000 keV) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists2->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists2->GetXaxis()->GetNbins();

    int slicewidth = 10, bin = 1, binmax = 1;

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

        slice_charge.push_back((bin+binmax)/2.+1000);
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
    gr2->GetYaxis()->SetRangeUser(0, 500);
    gr2->GetXaxis()->SetRangeUser(1000, 8000);
    TF1* fit2 = new TF1("fit2", "pol0(0)", 1000, 8000); 

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
    std::cout <<"Fitting banana3 (0-8000 keV) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists3->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists3->GetXaxis()->GetNbins();

    int slicewidth = 10, bin = 1, binmax = 1;

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
    gr3->GetYaxis()->SetRangeUser(0, 500);
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
    std::cout <<"Fitting banana4 (2000-8000 keV) for CFD - detector L \n" << i << std::endl;

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
        if (bin > 600) slicewidth = 150;

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

        slice_charge.push_back((bin+binmax)/2.+2000);
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
    gr4->GetYaxis()->SetRangeUser(0, 500);
    gr4->GetXaxis()->SetRangeUser(2000, 8000);
    TF1* fit4 = new TF1("fit4", "pol0(0)", 2000, 8000);   
    //fit4->SetParameter(0, 3.41957e+02);
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
    std::cout <<"Fitting banana5 (0-1000 keV) for detector Leading Edge \n" << i << std::endl;

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
    gr5->GetYaxis()->SetRangeUser(0, 600);
    gr5->GetXaxis()->SetRangeUser(0, 1000);
    //TF1* fit5 = new TF1("fit5", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)+exp([p4]+[p5]*x)+exp([p6]+[p7]*x)", 0, 1000);
    TF1* fit5 = new TF1("fit5", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)", 0, 1000);
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
    std::cout <<"Fitting banana6 (0-2000 keV) for detector Leading Edge \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists6->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6->GetXaxis()->GetNbins();

    int slicewidth = 10, bin = 1, binmax = 1;

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
    gr6->GetYaxis()->SetRangeUser(0, 500);
    gr6->GetXaxis()->SetRangeUser(0, 2000);

    TF1* fit6 = new TF1("fit6", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)", 0, 2000 ); 
    fit6->SetParameters(5.86287e+00, -7.12939e-06, 4.48340e+00, -1.96397e-03);
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

// TF1* banana7(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists7, TH2D* energy_slowfasttimediff_hists7, std::string dirPathplots)
// {
//     std::cout <<"Fitting banana7 (700-1000 keV) for detector Leading Edge " << i << std::endl;

//     int yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();
//     int xBins = slowfasttimediff_energy_hists7->GetXaxis()->GetNbins();

//     int slicewidth = 10, bin = 1, binmax = 1;

//     std::vector<Double_t> slice_charge(yBins);
//     std::vector<Double_t> slice_charge_err(yBins);
//     std::vector<Double_t> slice_time(yBins);
//     std::vector<Double_t> slice_time_err(yBins);

//     TCanvas * p = new TCanvas("p", "slices", 800, 600);
//     p->Divide(5,5);
//     int l = 3;
//     p->cd(1);
//     slowfasttimediff_energy_hists7->Draw();
//     p->cd(2);
//     energy_slowfasttimediff_hists7->Draw();

//     yBins = slowfasttimediff_energy_hists7->GetYaxis()->GetNbins();

//     while (bin < yBins)
//     {

//         binmax = bin + slicewidth;
//         if (binmax > yBins) binmax = yBins;

//         TH1D* ProjX = slowfasttimediff_energy_hists7->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
//         if (ProjX->GetEntries() == 0) 
//         {
//             delete ProjX;
//             bin = binmax;
//             continue;
//         }
//         TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  35, 60);

//         p->cd(l);
//         ProjX->Fit("fitSlice", "Q");
//         ProjX->Draw("same");
//         fitSlice->Draw("same");

//         slice_charge.push_back((bin+binmax)/2.+700);
//         slice_charge_err.push_back((slicewidth)/2.);
//         slice_time.push_back(fitSlice->GetParameter(1));
//         slice_time_err.push_back(fitSlice->GetParameter(2));
//         l++;
//         bin = binmax;
//     }

//     p -> SaveAs(Form("%s/banana7Slices_L%d.root", dirPathplots.c_str(), i));
//     delete p;

//     printf("Fitting banana7s for detector %d...\n", i);

//     TCanvas* cb = new TCanvas();
// 	TGraphErrors* gr7 = new TGraphErrors(slice_charge.size(), &slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]);
//     gr7 -> SetTitle(Form("Banana7 Fit L%d", i));
//     gr7 -> SetStats(0);
// 	gr7 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
// 	gr7 -> GetXaxis()->SetTitle("Mean Energy (keV)");
//     gr7 -> SetMarkerStyle(20);
//     gr7 -> SetMarkerSize(0.8);
//     gr7 -> SetMarkerColor(kBlue);
//     gr7->GetYaxis()->SetRangeUser(0, 600);
//     gr7->GetXaxis()->SetRangeUser(700, 1000);
//     TF1* fit7 = new TF1("fit7", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)", 700, 1000);
//     gr7->Draw("AL*");
//     gr7->Fit("fit7", "");
//     fit7->Draw("same");
//     cb -> SaveAs(Form("%s/banana7Fit_L%d.root", dirPathplots.c_str(), i));
//     cb->Write();

//     delete gr7;
//     printf("Clearing variables for detector...\n");
//     slice_charge.clear();
//     slice_charge_err.clear();
//     slice_time.clear();
//     slice_time_err.clear();

//     printf("Returning banana7 fit for detector %d...\n", i);
//     return fit7;
// }

//       //        //           //           //          //          // 

void runduration(TTree* tree, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds,unsigned long long numEntries)
{
    for (int i = 0; i < NumDetectors; i++) 
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
    printf("Time duration of the run = %f minutes\n", timeDuration);
}

std::vector<Double_t> fasttime_calibration_constants(bool beam, TFile *outputFile ,std::string dirPathplots, TDirectory *dir_fasttimecalib, TTree* tree, bool tensOfNanoSeconds, bool hundredsOfNanoSeconds, unsigned long long numEntries)
{
    if (beam)
    {
        printf("Beam run: Calculating inter-detector fast time calibration constants (used in next function)...\n");
        // Set branch addresses for each detector
        for (int i = 0; i < NumDetectors; i++) 
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
        for (int i = 1; i < NumDetectors; i++) 
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
            printf("Detector %d: Fast time diff mean = %f\n", i, mean_fasttimediff[i]);
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
        for (int i = 1; i < NumDetectors; i++) 
        {
            mean_fasttimediff[i] = 0;
            sigma_fasttimediff[i] = 0;
        }
    }
    
    return mean_fasttimediff;
}

// Function to process a range of entries using multiple threads
void process_entries_range(TFile *outputFile, std::string dirPathplots, TTree* copyTree, TDirectory *dir_banana, TDirectory *dir_aligned, unsigned long long threadStart, unsigned long long threadEnd) 
{

    // Initialize new TTree
    TTree* offsetTree = new TTree("LaBrDataOffset", "LaBrDataOffset"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
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
    for (int i = 0; i < NumDetectors; i++) 
    {
        correctedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        correctedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrected[i], Form("timeSCorrectedL%d/D", i)); //new Branch
    }
    correctedTree->Branch("timeRF", &timeFastCp[RFIndex], "timeRF/D");
    correctedTree->Branch("fastEPOLARIS", &fastEnergyPOLARISCp[0], "fastEPOLARIS/D");
    correctedTree->Branch("fastTPOLARIS", &fastTimePOLARISCp[0], "fastTPOLARIS/D");
    // ________________________________________________________________
    for (int i = 0; i < NumDetectors; i++) 
    {
        printf("Calculating time offset for detector %d...\n", i);
        for (unsigned long long entry = threadStart; entry < threadEnd; entry++) 
        {
            copyTree->GetEntry(entry);
            //printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCp[i], i, energySlowCp[i]);              

            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;
            Double_t deltaTime = 0;
            Double_t energy = 0;

            for (unsigned long long nextEntry = entry + 1; nextEntry < threadEnd; nextEntry++) 
            {
                copyTree->GetEntry(nextEntry);  
                // printf("timeFastL%d: %f, timeSlowL%d: %f, energySlowL%d: %f\n", i, timeFastCp[i], i, timeSlowCp[i], i, energySlowCp[i]);       
                       
                if (timeSlowCp[i] == 0) continue;
                nextTimeSlow = timeSlowCp[i];

                deltaTime = (nextTimeSlow - refTimeFast)/10; 
                energy = energySlowCp[i];

                if (deltaTime > 1000.0) break; // Stop searching if time difference exceeds 1 microsecond
                timeOffsetBeforeCorrection[i] = static_cast<Double_t>(static_cast<int>(deltaTime));
                //printf("Detector %d: Time Offset Forward = %f\n", i, timeOffsetBeforeCorrection[i]);

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
                //slowfasttimediff_energy_hists7[i]->Fill(deltaTime, energy);
                //energy_slowfasttimediff_hists7[i]->Fill(energy, deltaTime);
            }
            offsetTree->Fill();
        }
        
        TCanvas* c1 = new TCanvas(Form("c1_%d", i), Form("Uncorrected time offset: L%d", i), 800, 600);
        slowfasttimediff_hists[i]->GetXaxis()->SetRangeUser(0, 600);
        slowfasttimediff_hists[i]->SetStats(0);
        slowfasttimediff_hists[i]->SetFillColor(kBlue);
        slowfasttimediff_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_hists[i]->GetYaxis()->SetTitle("Counts/ns");
        slowfasttimediff_hists[i]->Draw();
        c1->SaveAs(Form("%s/time_offset_hist_L%d.root", dirPathplots.c_str(), i));
        c1->Write(Form("time_offset_hist_L%d", i));
        delete c1;
        energy_slowfasttimediff_hists[i]->GetYaxis()->SetRangeUser(0, 600);

        dir_banana->cd();
        banana1fit[i] = banana1(outputFile, i, slowfasttimediff_energy_hists1[i], energy_slowfasttimediff_hists1[i], dirPathplots);
        banana2fit[i] = banana2(outputFile, i, slowfasttimediff_energy_hists2[i], energy_slowfasttimediff_hists2[i], dirPathplots);
        banana3fit[i] = banana3(outputFile, i, slowfasttimediff_energy_hists3[i], energy_slowfasttimediff_hists3[i], dirPathplots);
        banana4fit[i] = banana4(outputFile, i, slowfasttimediff_energy_hists4[i], energy_slowfasttimediff_hists4[i], dirPathplots);
        banana5fit[i] = banana5(outputFile, i, slowfasttimediff_energy_hists5[i], energy_slowfasttimediff_hists5[i], dirPathplots);
        banana6fit[i] = banana6(outputFile, i, slowfasttimediff_energy_hists6[i], energy_slowfasttimediff_hists6[i], dirPathplots);
        //banana7fit[i] = banana7(outputFile, i, slowfasttimediff_energy_hists7[i], energy_slowfasttimediff_hists7[i], dirPathplots);
        
        TCanvas * c = new TCanvas("c", "banana", 800, 600);
        energy_slowfasttimediff_hists[i]->Draw();
        banana1fit[i]->SetLineColor(kRed);
        banana2fit[i]->SetLineColor(kGreen);
        banana3fit[i]->SetLineColor(kBlue);
        banana4fit[i]->SetLineColor(kYellow);
        banana5fit[i]->SetLineColor(kMagenta);
        banana6fit[i]->SetLineColor(kCyan);
        //banana7fit[i]->SetLineColor(kOrange);
        banana1fit[i]->Draw("same");
        banana2fit[i]->Draw("same");
        banana3fit[i]->Draw("same");
        banana4fit[i]->Draw("same");
        banana5fit[i]->Draw("same");
        banana6fit[i]->Draw("same");
        //banana7fit[i]->Draw("same");
        TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
        leg->AddEntry(banana1fit[i],"banana1fit","l");
        leg->AddEntry(banana2fit[i],"banana2fit","l");
        leg->AddEntry(banana3fit[i],"banana3fit","l");
        leg->AddEntry(banana4fit[i],"banana4fit","l");
        leg->AddEntry(banana5fit[i],"banana5fit","l");
        leg->AddEntry(banana6fit[i],"banana6fit","l");
        //leg->AddEntry(banana7fit[i],"banana7fit","l");
        leg->Draw();
        c->SaveAs(Form("%s/bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), i));
        c->Write(Form("bananas_2Dfitfunction_L%d", i));
        delete c;
    }
    dir_aligned->cd();

    // ________________________________________________________________
    Double_t minEnergy = 100;
    Double_t maxEnergy = 8000.0;
    Double_t energyStep = 1.0;

    // Iterate over the energy values and compute the precomputed values
    for (int i = 0; i < NumDetectors; i++)
    {
        printf("Calculating precomputed values for detector %d...\n", i);
        for (Double_t e = minEnergy; e <= maxEnergy; e += energyStep)
        {
            if ((e >= 2000) && (e <= 8000))
            {
                precomputedValues1[std::make_pair(i, e)] = banana1fit[i]->Eval(e);
                precomputedValues4[std::make_pair(i, e)] = banana4fit[i]->Eval(e);
                precomputedValues2[std::make_pair(i, e)] = banana2fit[i]->Eval(e);
                
            }
            else if (e <= 8000)
            {
                precomputedValues3[std::make_pair(i, e)] = banana3fit[i]->Eval(e);
            }
            else if (e < 2000)
            {
                precomputedValues6[std::make_pair(i, e)] = banana6fit[i]->Eval(e);

                if (e < 1000) 
                {
                    precomputedValues5[std::make_pair(i, e)] = banana5fit[i]->Eval(e);
                }
            }
        }
    }
    
    // ________________________________________________________________

    // Loop through the entries to update timeSlowCorrected[i] and fill offsetTree
    for (int i = 0; i < NumDetectors; i++) 
    {
        printf("Calculating the corrected time offset for detector %d...\n", i);
        
        for (unsigned long long entry = threadStart; entry < threadEnd; entry++) 
        {
            offsetTree->GetEntry(entry);     
            Double_t t = timeOffsetBeforeCorrection[i];

            if (t == 0) continue;

            Double_t e = static_cast<Double_t>(static_cast<int>(energySlowCp[i]));
            Double_t timeOffset = 0.0;

            if ((t <= 264) && (e >= 2000)) 
            {
                timeOffset = precomputedValues1[std::make_pair(i, e)];
                // printf("Detector %d: Time Offset function 1 = %f\n", i, timeOffset);
            } 
            else if ((t > 264) && (t <= 275) && (e >= 1000)) 
            {
                timeOffset = precomputedValues2[std::make_pair(i, e)];
                // printf("Detector %d: Time Offset function 2 = %f\n", i, timeOffset);
            }
            else if ((t > 275) && (t <= 300)) 
            {
                timeOffset = precomputedValues3[std::make_pair(i, e)];
                // printf("Detector %d: Time Offset function 3 = %f\n", i, timeOffset);
            } 
            else if ((t > 300) && (t <= 370) && (e >= 2000)) 
            {
                timeOffset = precomputedValues4[std::make_pair(i, e)];
                // printf("Detector %d: Time Offset function 4 = %f\n", i, timeOffset);
            }
            else if ((e < 1000) && (t <= 340)) 
            {
                timeOffset = precomputedValues5[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
            } 
            else if ((e < 2000) && (t > 340)) 
            {
                    timeOffset = precomputedValues6[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 6 = %f\n", i, timeOffset);
            }
            else continue;

            // Calculate timeSlowCorrected and fill correctedTree using timeOffset
            if (timeOffset != 0.0) 
            {
                timeSlowCorrected[i] = timeSlowCp[i] - timeOffset * 10;
                correctedTree->Fill();
            }
        }
    }

    // calculate the time difference between the first timeSlowCorrected[0] >0 and the last timeSlowCorrected[0]>0 in correctedTree
    Double_t fTime = 0;
    Double_t lTime = 0;


    for (int i = 0;  i < 1; i++)
    {
        for (unsigned long long entry = 0; entry < threadEnd; entry++) 
        {
            correctedTree->GetEntry(entry);
            if (timeSlowCorrected[i] != 0)
            {
                fTime = timeSlowCorrected[i];
                break;
            }
        }

        // Find the last non-zero entry in the timeSlowCorrected[0] branch of the tree
        for (unsigned long long entry = threadEnd - 1; entry >= 0; entry--) 
        {
            correctedTree->GetEntry(entry);
            if (timeSlowCorrected[i] != 0)
            {
                lTime = timeSlowCorrected[i];
                break;
            }
        }
    }

    Double_t timeDuration = (lTime - fTime)*1e-9/60;

    printf("First timeSlowCorrected[0] = %f\n", fTime);
    printf("Last timeSlowCorrected[0] = %f\n", lTime);
    printf("----- Time duration of the run = %f minutes (CHECK IF NECESSARY TO TIME SORT DATA) -----\n", timeDuration);

    for(int i = 0; i < NumDetectors; i++) 
    {

        // CAN COMMENT THIS OUT IF NOT LOOKING TO SEE THE LINEARITY OF THE ENERGY VS TIME OFFSET HISTOGRAM
        printf("Calculating corrected time offset for to plot hist for detector %d...\n", i);
        for (unsigned long long entry = threadStart; entry < threadEnd; entry++) 
        {
            correctedTree->GetEntry(entry);
            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextEnergySlow = 0;
            Double_t nextTimeSlow = 0;
            Double_t deltaTime = 0;

            for (unsigned long long nextEntry = entry + 1; nextEntry < threadEnd; nextEntry++) 
            {
                correctedTree->GetEntry(nextEntry);            
                if (timeSlowCorrected[i] == 0) continue;
                nextEnergySlow = energySlowCp[i];
                nextTimeSlow = timeSlowCorrected[i];
                
                //printf("Detector %d: Fast Time = %f, Slow Time = %f\n", i, refTimeFast, nextTimeSlow);
                //printf("Detector %d: Corrected Time Offset = %f\n", i, deltaTime);
                deltaTime = (nextTimeSlow - refTimeFast)/10;
                if (deltaTime > 600.0 || deltaTime < -600.0) break; // Stop searching if time difference is greater than 600 ns or less than -600 ns
                //printf("Detector %d: Corrected Time Offset = %f\n", i, deltaTime);
                timeOffsetCorrectedHists[i]->Fill(deltaTime);
                timeEnergyOffsetCorrectedHists[i]->Fill(deltaTime, nextEnergySlow);
            }
        }
    
        timeOffsetCorrectedHists[i]->Write();
        timeEnergyOffsetCorrectedHists[i]->Write();
        
        TCanvas* c2 = new TCanvas(Form("c2_%d", i), Form("Uncorrected time offset vs energy: L%d", i), 800, 600);
        slowfasttimediff_energy_hists[i]->SetStats(0);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetRangeUser(0, 600);
        slowfasttimediff_energy_hists[i]->SetMarkerStyle(20);
        slowfasttimediff_energy_hists[i]->SetMarkerSize(0.8);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_energy_hists[i]->GetYaxis()->SetTitle("Energy (keV)");
        slowfasttimediff_energy_hists[i]->SetMaximum(1200);
        slowfasttimediff_energy_hists[i]->Draw("COLZ");
        c2->SaveAs(Form("%s/uncorrected_time_offset_vs_energy_hist_L%d.root", dirPathplots.c_str(), i));
        c2->Write(Form("uncorrected_time_offset_vs_energy_hist_L%d", i));
        delete c2;

        TCanvas* c3 = new TCanvas(Form("c3_%d", i), Form("Corrected time offset: L%d", i), 800, 600);
        //timeOffsetCorrectedHists[i]->SetStats(0);
        timeOffsetCorrectedHists[i]->SetFillColor(kBlue);
        timeOffsetCorrectedHists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        timeOffsetCorrectedHists[i]->GetYaxis()->SetTitle("Counts/ns");
        timeOffsetCorrectedHists[i]->Draw();
        c3->SaveAs(Form("%s/corrected_time_offset_banana_hist_L%d.root", dirPathplots.c_str(), i));
        c3->Write(Form("corrected_time_offset_banana_hist_L%d", i));

        TCanvas* c4 = new TCanvas(Form("c4_%d", i), Form("Corrected time offset vs energy: L%d", i), 800, 600);
        //timeEnergyOffsetCorrectedHists[i]->SetStats(0);
        timeEnergyOffsetCorrectedHists[i]->SetMarkerStyle(20);
        timeEnergyOffsetCorrectedHists[i]->SetMarkerSize(0.8);
        timeEnergyOffsetCorrectedHists[i]->GetYaxis()->SetTitle("Energy (keV)");
        timeEnergyOffsetCorrectedHists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        timeEnergyOffsetCorrectedHists[i]->SetMaximum(1200);
        timeEnergyOffsetCorrectedHists[i]->Draw("COLZ");
        c4->SaveAs(Form("%s/corrected_energy_vs_time_offset_banana_hist_L%d.root", dirPathplots.c_str(), i));
        c4->Write(Form("corrected_energy_vs_time_offset_banana_hist_L%d", i));

        delete c3;
        delete c4;        
    }

    // Write the histograms to the output file
    for (int i = 0; i < NumDetectors; i++) 
    {
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
    for (int i = 0; i < NumDetectors +2; i++) 
    {
        printf("Aligning data for detector %d...\n", i);
        for (unsigned long long entry = threadStart; entry < threadEnd; entry++) 
        {
            correctedTree->GetEntry(entry);

            if (i == RFIndex) 
            {
                detectorID = RFIndex;
                timeGlobal[0] = timeFastCp[i];
                alignedTimeFast[0] = timeFastCp[i];
                alignedEnergyFast[0] = 0;
                alignedTimeSlow[0] = 0;
                alignedEnergySlow[0] = 0;
                alignedTree->Fill();
                // printf("RF Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i == POLARISIndex) 
            {
                detectorID = POLARISIndex;
                timeGlobal[0] = fastTimePOLARISCp[0];
                alignedTimeFast[0] = fastTimePOLARISCp[0];
                alignedEnergyFast[0] = fastEnergyPOLARISCp[0];
                alignedTimeSlow[0] = 0;
                alignedEnergySlow[0] = 0;
                alignedTree->Fill();
                // printf("POLARIS Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i < RFIndex)
            {
                detectorID = i;
                // printf("Detector %d: timeFastCp = %f | , timeSlowCorrected = %f\n", i, timeFastCp[i], timeSlowCorrected[i]);
                if ((timeFastCp[i] == 0) && (timeSlowCorrected[i] != 0)) 
                {
                    timeGlobal[0] = timeSlowCorrected[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrected[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                    alignedTree->Fill();
                    // printf("Loop1 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrected = %f, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrected[i] == 0)) 
                {
                    timeGlobal[0] = timeFastCp[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrected[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                    alignedTree->Fill();
                    // printf("Loop2 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrected = %f\n, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrected[i] != 0)) 
                {
                    // The ratio in the print below should be 1
                    //printf("Error: Both fast and slow times are non-zero for detector %d\n, timeFastCp = %f | , timeSlowCorrected + timeOffsets = %f\n, ratio of timeFastCp to timeSlowCp + timeOffsets = %f\n", i, timeFastCp[i], timeSlowCp[i] + timeOffsets[i], (timeFastCp[i]/(timeSlowCp[i] + timeOffsets[i])));
                    timeGlobal[0] = timeFastCp[i];
                    alignedTimeFast[0] = timeFastCp[i];
                    alignedEnergyFast[0] = 0;
                    alignedTimeSlow[0] = timeSlowCorrected[i];
                    alignedEnergySlow[0] = energySlowCp[i];
                    alignedTree->Fill();
                }
                else 
                {
                    // The slow energy in the print below should be zero
                    //printf("Error: Both fast and slow times are zero for detector %d\n, Slow Energy = %f\n", i, energySlowCp[i]);
                    continue;
                }
            }
        }
    }

    alignedTree->Write();
}

void align_offset(TFile *outputFile, std::string dirPathplots, TTree* tree, TDirectory *dir_banana, TDirectory *dir_aligned, std::vector<Double_t> mean_fasttimediff, bool& tensOfNanoSeconds, bool& hundredsOfNanoSeconds, unsigned long long numEntries) {
    printf("Aligning offsets...\n");
    
    dir_aligned->cd();
    // Set branch addresses for each detector
    for (int i = 0; i < NumDetectors; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
    tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
    tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

    printf("Cloning the TTree...\n");
    // Set branch addresses for each detector
    for (int i = 0; i < NumDetectors; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
    tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
    tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

    TTree* copyTree = new TTree("LaBrDataCopy", "LaBrDataCopy");
    for (int i = 0; i < NumDetectors; i++) 
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
        for (int i = 0; i < NumDetectors; i++) 
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
    for (int i = 0; i < NumDetectors; i++) 
    {
        slowfasttimediff_hists[i] = new TH1D(Form("slowfasttimediff_hists%d", i), Form("Uncorrected time offset: L%d", i),  600, 0, 600);
        
        // Full
        slowfasttimediff_energy_hists[i] = new TH2D(Form("slowfasttimediff_energy_hists%d", i), Form("Uncorrected time offset vs energy: L%d", i),  600, 0, 600, 8000, 0, 8000);
        energy_slowfasttimediff_hists[i] = new TH2D(Form("energy_slowfasttimediff_hists%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 600, 0, 600);

        // banana 1
        slowfasttimediff_energy_hists1[i] = new TH2D(Form("slowfasttimediff_energy_hists1_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  14, 250, 264, 6000, 2000, 8000);
        energy_slowfasttimediff_hists1[i] = new TH2D(Form("energy_slowfasttimediff_hists1_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  6000, 2000, 8000, 14, 250, 264);

        // banana 2
        slowfasttimediff_energy_hists2[i] = new TH2D(Form("slowfasttimediff_energy_hists2_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  17, 250, 272, 7000, 1000, 8000);
        energy_slowfasttimediff_hists2[i] = new TH2D(Form("energy_slowfasttimediff_hists2_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  7000, 1000, 8000, 17, 255, 272);

        // banana 3
        slowfasttimediff_energy_hists3[i] = new TH2D(Form("slowfasttimediff_energy_hists3_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  28, 272, 300, 7000, 1000, 8000);
        energy_slowfasttimediff_hists3[i] = new TH2D(Form("energy_slowfasttimediff_hists3_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  7000, 1000, 8000, 28, 272, 300);

        // banana 4
        slowfasttimediff_energy_hists4[i] = new TH2D(Form("slowfasttimediff_energy_hists4_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  70, 300, 370, 6000, 2000, 8000);
        energy_slowfasttimediff_hists4[i] = new TH2D(Form("energy_slowfasttimediff_hists4_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  6000, 2000, 8000, 70, 300, 370);

        // banana 5
        slowfasttimediff_energy_hists5[i] = new TH2D(Form("slowfasttimediff_energy_hists5_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  70, 270, 340, 1000, 0, 1000);
        energy_slowfasttimediff_hists5[i] = new TH2D(Form("energy_slowfasttimediff_hists5_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  1000, 0, 1000, 70, 270, 340);

        // banana 6
        slowfasttimediff_energy_hists6[i] = new TH2D(Form("slowfasttimediff_energy_hists6_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  280, 320, 600, 2000, 0, 2000);
        energy_slowfasttimediff_hists6[i] = new TH2D(Form("energy_slowfasttimediff_hists6_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 280, 320, 600);

        // banana 7
        // slowfasttimediff_energy_hists7[i] = new TH2D(Form("slowfasttimediff_energy_hists7_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  50, 370, 420, 300, 700, 1000);
        // energy_slowfasttimediff_hists7[i] = new TH2D(Form("energy_slowfasttimediff_hists7_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  300, 700, 1000, 50, 370, 420);

        timeOffsetCorrectedHists[i] = new TH1D(Form("timeOffsetCorrectedHists%d", i), Form("Corrected time offset: L%d", i),  1200, -600, 600);
        timeEnergyOffsetCorrectedHists[i] = new TH2D(Form("timeEnergyOffsetCorrectedHists%d", i), Form("Corrected time offset vs energy: L%d", i),  1200, -600, 600, 8000, 0, 8000);
    }
    
    //________________________________________________________________________________________________________________________

    const int numThreads = 1; // You can adjust the number of threads according to your system
    std::vector<std::thread> threads;

    // Split the range of entries into equal parts for each thread
    const unsigned long long entriesPerThread = (numEntries) / numThreads;
    unsigned long long threadStart = 0;
    unsigned long long endEntry = numEntries;
    unsigned long long threadEnd = threadStart + entriesPerThread;

    // Create and start the threads
    for (int i = 0; i < numThreads; i++) 
    {
        printf("Thread %d: threadStart = %lld, threadEnd = %lld\n", i, threadStart, threadEnd);
        try 
        {
            threads.emplace_back(process_entries_range, outputFile, dirPathplots, copyTree, dir_banana, dir_aligned, threadStart, threadEnd);
        } 
        catch (const std::system_error& e) 
        {
            // Handle the exception (e.g., log an error message or take appropriate action)
        }
        // Update the entry ranges for the next thread
        threadStart = threadEnd;
        threadEnd = i == numThreads - 2 ? endEntry : threadEnd + entriesPerThread; 
        

    }

    // Wait for all threads to finish
    for (int i = 0; i < numThreads; i++) 
    {
        try 
        {
            threads[i].join();
        } 
        catch (const std::system_error& e) 
        {
            // Handle the exception (e.g., log an error message or take appropriate action)
        }
    }

    printf("All threads have finished processing.\n");
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
    unsigned long long numEntries = tree->GetEntries();
    numEntries = numEntries/10;

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
