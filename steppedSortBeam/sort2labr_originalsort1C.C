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

// END BANANA FIT FUNCTIONS //

// banana1 good fit
TF1* banana1(TFile *outputFile, int i, TH2D* slowfasttimediff_energy_hists1, TH2D* energy_slowfasttimediff_hists1, std::string dirPathplots)
{
    std::cout <<"\nFitting banana1 (0-8000 keV) for CFD - detector L \n" << i << std::endl;

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
    gr1->GetYaxis()->SetRangeUser(0, 500);
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
    std::cout <<"Fitting banana2 (0-8000 keV) for CFD - detector L \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists2->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists2->GetXaxis()->GetNbins();

    int slicewidth = 50, bin = 1, binmax = 1;

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

    int slicewidth = 100, bin = 1, binmax = 1;

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

    //TF1* fit3 = new TF1("fit3", " pol0(0)", 0, 8000); 
    TF1* fit3 = new TF1("fit3", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)", 500, 8000); 
    fit3->SetParameters(5.63514e+00, 1.52326e-06, 3.13394e+00, -1.07626e-03);
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
    TF1* fit5 = new TF1("fit5", "pol0(0)", 0, 1000);
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
    std::cout <<"Fitting banana6 (0-500 keV) for detector Leading Edge \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists6->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists6->GetXaxis()->GetNbins();

    int slicewidth = 50, bin = 1, binmax = 1;

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
    gr6->GetXaxis()->SetRangeUser(0, 500);

    TF1* fit6 = new TF1("fit6", "gaus(0)", 0, 500 ); 
    fit6->SetParameters(4.35397e+02, -2.67761e+02, 1.43843e+03);
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
    std::cout <<"Fitting banana7 (0-1000 keV) for detector Leading Edge \n" << i << std::endl;

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
    gr7->GetYaxis()->SetRangeUser(0, 600);
    gr7->GetXaxis()->SetRangeUser(0, 1000);
    //TF1* fit7 = new TF1("fit7", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)+exp([p4]+[p5]*x)+exp([p6]+[p7]*x)", 0, 1000);
    TF1* fit7 = new TF1("fit7", "pol0(0)", 0, 1000);
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
    std::cout <<"Fitting banana8 (500-4500 keV) for detector Leading Edge \n" << i << std::endl;

    int yBins = slowfasttimediff_energy_hists8->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists8->GetXaxis()->GetNbins();

    int slicewidth = 120, bin = 1, binmax = 1;

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

        slice_charge.push_back((bin+binmax)/2.+500);
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
    gr8->GetYaxis()->SetRangeUser(0, 600);
    gr8->GetXaxis()->SetRangeUser(500, 4300);
    TF1* fit8 = new TF1("fit8", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)", 500, 4300 ); 
    fit8->SetParameters(5.84260e+00, -2.37579e-06, 4.21317e+00, -1.31619e-03);
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
        slowfasttimediff_energy_hists1[i] = new TH2D(Form("slowfasttimediff_energy_hists1_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  14, 250, 264, 8000, 0, 8000);
        energy_slowfasttimediff_hists1[i] = new TH2D(Form("energy_slowfasttimediff_hists1_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 14, 250, 264);

        // banana 2
        slowfasttimediff_energy_hists2[i] = new TH2D(Form("slowfasttimediff_energy_hists2_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  17, 250, 272, 7000, 1000, 8000);
        energy_slowfasttimediff_hists2[i] = new TH2D(Form("energy_slowfasttimediff_hists2_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  7000, 1000, 8000, 17, 255, 272);

        // banana 3
        slowfasttimediff_energy_hists3[i] = new TH2D(Form("slowfasttimediff_energy_hists3_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  28, 272, 300, 8000, 0, 8000);
        energy_slowfasttimediff_hists3[i] = new TH2D(Form("energy_slowfasttimediff_hists3_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  8000, 0, 8000, 28, 272, 300);

        // banana 4
        slowfasttimediff_energy_hists4[i] = new TH2D(Form("slowfasttimediff_energy_hists4_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  70, 300, 370, 6000, 2000, 8000);
        energy_slowfasttimediff_hists4[i] = new TH2D(Form("energy_slowfasttimediff_hists4_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  6000, 2000, 8000, 70, 300, 370);

        // banana 5
        slowfasttimediff_energy_hists5[i] = new TH2D(Form("slowfasttimediff_energy_hists5_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  20, 325, 345, 1000, 0, 1000);
        energy_slowfasttimediff_hists5[i] = new TH2D(Form("energy_slowfasttimediff_hists5_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  1000, 0, 1000, 20, 325, 345);

        // banana 6
        slowfasttimediff_energy_hists6[i] = new TH2D(Form("slowfasttimediff_energy_hists6_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  280, 320, 600, 500, 0, 500);
        energy_slowfasttimediff_hists6[i] = new TH2D(Form("energy_slowfasttimediff_hists6_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  500, 0, 500, 280, 320, 600);

        // banana 7
        slowfasttimediff_energy_hists7[i] = new TH2D(Form("slowfasttimediff_energy_hists7_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  25, 300, 325, 1000, 0, 1000);
        energy_slowfasttimediff_hists7[i] = new TH2D(Form("energy_slowfasttimediff_hists7_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  1000, 0, 1000, 25, 300, 325);

        // banana 8
        slowfasttimediff_energy_hists8[i] = new TH2D(Form("slowfasttimediff_energy_hists8_%d", i), Form("Uncorrected time offset vs energy: L%d", i),  280, 320, 600, 4000, 500, 4500);
        energy_slowfasttimediff_hists8[i] = new TH2D(Form("energy_slowfasttimediff_hists8_%d", i), Form("Uncorrected energy vs time offset: L%d", i),  4000, 500, 4500, 280, 320, 600);

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

                deltaTime = (nextTimeSlow - refTimeFast)/10; 
                energy = energySlowCp[i];
                if ((deltaTime > 600.0) || (deltaTime < -600.0)) break;
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
        banana1fit[i]->Draw("same");
        banana2fit[i]->Draw("same");
        banana3fit[i]->Draw("same");
        banana4fit[i]->Draw("same");
        banana5fit[i]->Draw("same");
        banana6fit[i]->Draw("same");
        banana7fit[i]->Draw("same");
        banana8fit[i]->Draw("same");
        TLegend *leg1 = new TLegend(0.1,0.7,0.48,0.9);
        leg1->AddEntry(banana1fit[i],"banana1fit","l");
        leg1->AddEntry(banana2fit[i],"banana2fit","l");
        leg1->AddEntry(banana3fit[i],"banana3fit","l");
        leg1->AddEntry(banana4fit[i],"banana4fit","l");
        leg1->AddEntry(banana5fit[i],"banana5fit","l");
        leg1->AddEntry(banana6fit[i],"banana6fit","l");
        leg1->AddEntry(banana7fit[i],"banana7fit","l");
        leg1->AddEntry(banana8fit[i],"banana8fit","l");
        leg1->Draw();
        c2->SaveAs(Form("%s/bananas_2Dfitfunction_L%d.root", dirPathplots.c_str(), i));
        c2->Write(Form("bananas_2Dfitfunction_L%d", i));
        delete c2;
    
        dir_aligned->cd();

        // ________________________________________________________________
        Double_t minEnergy = 1;
        Double_t maxEnergy = 8000.0;
        Double_t energyStep = 1.0;

        // Iterate over the energy values and compute the precomputed values

        printf("\n\nCalculating precomputed values for detector %d...\n", i);
        for (Double_t e = minEnergy; e <= maxEnergy; e += energyStep)
        {
            if ((e >= 2000) && (e <= 8000))
            {
                precomputedValues4[std::make_pair(i, e)] = banana4fit[i]->Eval(e);
                
            }
            precomputedValues1[std::make_pair(i, e)] = banana1fit[i]->Eval(e);
            precomputedValues2[std::make_pair(i, e)] = banana2fit[i]->Eval(e);
            precomputedValues3[std::make_pair(i, e)] = banana3fit[i]->Eval(e);

            if ((e>=500)&&(e<=4500))
            {
                precomputedValues8[std::make_pair(i, e)] = banana8fit[i]->Eval(e);
            }

            if (e < 2500)
            {
                precomputedValues5[std::make_pair(i, e)] = banana5fit[i]->Eval(e);
                precomputedValues7[std::make_pair(i, e)] = banana7fit[i]->Eval(e);
                if (e < 500)
                {
                    precomputedValues6[std::make_pair(i, e)] = banana6fit[i]->Eval(e);
                }
            }
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

            if ((t <= 264))
            {
                timeOffset = precomputedValues1[std::make_pair(i, e)];
                // printf("Detector %d: Time Offset function 1 = %f\n", i, timeOffset);
            } 
            else if ((t > 264) && (t <= 275) && (e>=1000)) 
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
            else if ((e <= 1000) && (t >= 325) && (t <= 345)) 
            {
                timeOffset = precomputedValues5[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
            } 
            else if ((e <= 1000) && (t >= 300) && (t <= 345)) 
            {
                timeOffset = precomputedValues7[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 5 = %f\n", i, timeOffset);
            } 
            else if ((e < 500) && (t > 340)) 
            {
                    timeOffset = precomputedValues6[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 6 = %f\n", i, timeOffset);
            }
            else if ((e >= 500)&& (e<=4500) && (t > 340)) 
            {
                    timeOffset = precomputedValues8[std::make_pair(i, e)];
                    // printf("Detector %d: Time Offset function 6 = %f\n", i, timeOffset);
            }

            // Calculate timeSlowCorrected and fill correctedTree using timeOffset
            //if (timeOffset != 0.0) 
            {
                timeSlowCorrected[i] = timeSlowCp[i] - timeOffset * 10;
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
                
                deltaTime = (nextTimeSlow - refTimeFast)/10;
                // printf("Detector %d: Time Offset Forward = %f\n", i, deltaTime2);
                if ((deltaTime > 600.0) || (deltaTime < -600.0)) break;
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
    unsigned long long numEntries = (tree->GetEntries())/10;

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
