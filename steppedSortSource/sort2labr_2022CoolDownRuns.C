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
#include <Math/WrappedTF1.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <ROOT/RDataFrame.hxx> 
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

std::vector<TH1D*> slowfasttimediff_hists(NumDetectors);
std::vector<TH1D*> energy_hists(NumDetectors);
std::vector<TH1D*> timeOffsetCorrectedHists(NumDetectors);
std::vector<TH2D*> slowfasttimediff_energy_hists(NumDetectors);
std::vector<TH2D*> energy_slowfasttimediff_hists(NumDetectors);  
std::vector<TH2D*> timeEnergyOffsetCorrectedHists(NumDetectors);

std::vector<Double_t> mean_fasttimediff(NumDetectors-1);
std::vector<Double_t> sigma_fasttimediff(NumDetectors-1);
bool tensOfNanoSeconds = false;
bool hundredsOfNanoSeconds = false;

TF1* banana(int i, TH2D* slowfasttimediff_energy_hists, TH2D* energy_slowfasttimediff_hists, std::string dirPathplots);
TF1* bananafit[NumDetectors];

TF1* banana(int i, TH2D* slowfasttimediff_energy_hists, TH2D* energy_slowfasttimediff_hists, std::string dirPathplots)
{
    std::cout <<"\n Fitting banana for detector " << i << std::endl;
    int yBins = slowfasttimediff_energy_hists->GetYaxis()->GetNbins();
    int xBins = slowfasttimediff_energy_hists->GetXaxis()->GetNbins();

    int slicewidth = 50, bin = 1, binmax = 1;

    std::vector<Double_t> slice_charge1(yBins);
    std::vector<Double_t> slice_charge_err1(yBins);
    std::vector<Double_t> slice_time1(yBins);
    std::vector<Double_t> slice_time_err1(yBins);

    TCanvas * p = new TCanvas("p", "slices", 800, 600);
    p->Divide(5,5);
    int l = 3;
    p->cd(1);
    slowfasttimediff_energy_hists->Draw();
    p->cd(2);
    energy_slowfasttimediff_hists->Draw();

    yBins = slowfasttimediff_energy_hists->GetYaxis()->GetNbins();

    while (bin < yBins)
    {
        //if (bin > 1000) slicewidth = 30;
        binmax = bin + slicewidth;
        if (binmax > yBins) binmax = yBins;

        TH1D* ProjX = slowfasttimediff_energy_hists->ProjectionX(Form("ProjX%d", bin), bin, binmax);        
        if (ProjX->GetEntries() == 0) 
        {
            delete ProjX;
            bin = binmax;
            continue;
        }
        TF1* fitSlice = new TF1("fitSlice", "gaus(0)",  50, 200);

        p->cd(l);
        ProjX->Fit("fitSlice", "Q");
        ProjX->Draw("same");
        fitSlice->Draw("same");

        slice_charge1.push_back((bin+binmax)/2.);
        slice_charge_err1.push_back(slicewidth/2.);
        slice_time1.push_back(fitSlice->GetParameter(1));
        slice_time_err1.push_back(fitSlice->GetParameter(2));

        l++;
        bin = binmax;
    }

    p -> SaveAs(Form("%s/bananaSlices_expo1_L%d.root", dirPathplots.c_str(), i));
    delete p;

    printf("\n Fitting bananas for detector %d...\n", i);

    slowfasttimediff_energy_hists->GetXaxis()->SetRange( 50, 2000);
    energy_slowfasttimediff_hists->GetYaxis()->SetRange( 50, 2000);

    TCanvas* cb = new TCanvas();
	// TGraphErrors* gr1 = new TGraphErrors(slice_charge1.size(),&slice_charge1[0], &slice_time1[0], &slice_charge_err1[0], &slice_time_err1[0]);
	TGraphErrors* gr1 = new TGraphErrors(slice_charge1.size(), &slice_charge1[0], &slice_time1[0], &slice_charge_err1[0], &slice_time_err1[0]);
    gr1 -> SetTitle(Form("Banana Fit L%d", i));
    gr1 -> SetStats(0);
	gr1 -> GetYaxis()->SetTitle("Mean Time Difference (ns)");
	gr1 -> GetXaxis()->SetTitle("Mean Energy (keV)");
    gr1 -> SetMarkerStyle(20);
    gr1 -> SetMarkerSize(0.8);
    gr1 -> SetMarkerColor(kBlue);
    gr1->GetYaxis()->SetRangeUser( 50, 2000);
    gr1->GetXaxis()->SetRangeUser( 0,2000);
    TF1* fit1 = new TF1("fit1", "exp([p0]+[p1]*x)+exp([p2]+[p3]*x)+exp([p4]+[p5]*x)+exp([p6]+[p7]*x)", 0, 2000);
    fit1->SetParameters(5.02954e+00, -1.06957e-04, 3.91213e+00, -2.90836e-03, 3.92668e+00, -2.83541e-03, 2.74410e+00, -2.85524e-03);

    gr1->Draw("AL*");
    gr1->Fit("fit1", "");
    fit1->Draw("same");
    cb -> SaveAs(Form("%s/bananaFit_L%d.root", dirPathplots.c_str(), i));
    cb->Write();

    delete gr1;
    printf("\nClearing variables for detector...\n");
    slice_charge1.clear();
    slice_charge_err1.clear();
    slice_time1.clear();
    slice_time_err1.clear();

    printf("\nReturning banana fit for detector %d...\n", i);
    return fit1;
}

void runduration(TTree *tree)
{
    printf("\nCalculating run duration...\n");

    for (int i = 0; i < NumDetectors; i++) 
    {
        tree->SetBranchAddress(Form("slowECalibL%d", i), &energySlow[i]);
        tree->SetBranchAddress(Form("timeSL%d", i), &timeSlow[i]);
        tree->SetBranchAddress(Form("timeFL%d", i), &timeFast[i]);
    }
    tree->SetBranchAddress("timeRF", &timeFast[RFIndex]);
    tree->SetBranchAddress("fastEPOLARIS", &fastEnergyPOLARIS[0]);
    tree->SetBranchAddress("fastTPOLARIS", &fastTimePOLARIS[0]);

    // Find the first non-zero entry in the timeFL0 branch of the tree
    uint64_t numEntries = (tree->GetEntries());
    Double_t firstTime;
    Double_t lastTime;

    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);
        if (timeFast[0] != 0) 
        {
            firstTime = timeFast[0];
            break;
        }
    }

    // Find the last non-zero entry in the timeFL0 branch of the tree
    for (uint64_t entry = numEntries - 1; entry >= 0; entry--) 
    {
        tree->GetEntry(entry);
        if (timeFast[0] != 0) 
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
    printf("\nTime duration of the run = %f minutes\n", timeDuration);
}

std::vector<Double_t> fasttime_calibration_constants(TTree *tree, std::string dirPathplots, TDirectory *dir_fasttimecalib, bool beam)
{
    if (beam)
    {
        printf("\n Beam run: Calculating inter-detector fast time calibration constants (used in next function)...\n");
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

        uint64_t numEntries = (tree->GetEntries());

        for (uint64_t entry = 0; entry < numEntries; entry++) 
        {
            tree->GetEntry(entry);
            if (timeFast[0] == 0) continue;
            Double_t timeFastDet0 = 0;  // Initialize
            timeFastDet0 = timeFast[0];

            for (int i = 1; i < NumDetectors; i++) 
            {
                tree->GetEntry(entry + 1); // Get the next entry for timeFastDetX

                if (timeFast[i] != 0) 
                {
                    Double_t timeFastDetX = 0;  // Initialize
                    
                    timeFastDetX = timeFast[i];
                    Double_t deltaTime = (timeFastDetX - timeFastDet0);

                    //printf("\n timeFastL0: %f, timeFastL%d: %f\n", timeFastDet0, i, timeFastDetX);
                    //printf("\nDetector %d - 0: Fast time diff = %f\n", i, deltaTime);

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
            printf("Detector %d: Fast time diff sigma = %f\n\n", i, sigma_fasttimediff[i]);
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
        printf("\n Source run: not calculating inter-detector fast time calibration constants ...\n");
        for (int i = 1; i < NumDetectors; i++) 
        {
            mean_fasttimediff[i] = 0;
            sigma_fasttimediff[i] = 0;
        }
    }
    
    return mean_fasttimediff;
}

void align_offset(TTree *tree, std::string dirPathplots,TDirectory *dir_banana, TDirectory *dir_aligned, std::vector<Double_t> mean_fasttimediff)
{

    printf("\n Cloning the TTree...\n");
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

    // Get number of entries in the tree
    uint64_t numEntries = (tree->GetEntries());

    // variables for copytree data
    std::vector<Double_t> energySlowCp(NumDetectors);
    std::vector<Double_t> timeSlowCp(NumDetectors);
    std::vector<Double_t> timeSlowCorrected(NumDetectors);
    std::vector<Double_t> timeFastCp(NumDetectors+1); 
    std::vector<Double_t> fastEnergyPOLARISCp(1);
    std::vector<Double_t> fastTimePOLARISCp(1);

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

    double slowTimeScaleFactor = 1.0;
    if (tensOfNanoSeconds) slowTimeScaleFactor = 0.1;
    else if (hundredsOfNanoSeconds) slowTimeScaleFactor = 0.01;

    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        tree->GetEntry(entry);

        // Apply the scaling factor for all detectors
        for (int i = 0; i < NumDetectors; i++) 
        {
            energySlowCp[i] = energySlow[i];
            timeSlowCp[i] = timeSlow[i] * slowTimeScaleFactor;
        }

        timeFastCp[0] = timeFast[0] * slowTimeScaleFactor;
        for (int i = 1; i <= 3; i++) 
        {
            timeFastCp[i] = (timeFast[i] - mean_fasttimediff[i]) * slowTimeScaleFactor;
        }

        timeFastCp[RFIndex] = timeFast[RFIndex] * slowTimeScaleFactor;
        fastTimePOLARISCp[0] = fastTimePOLARIS[0] * slowTimeScaleFactor;
        fastEnergyPOLARISCp[0] = fastEnergyPOLARIS[0];

        copyTree->Fill();
    }

    // Initialise the histograms
    for (int i = 0; i < NumDetectors; i++) 
    {
        slowfasttimediff_hists[i] = new TH1D(Form("slowfasttimediff_hists%d", i), Form("Uncorrected time offset: L%d", i),  2000, 0, 2000);
        energy_hists[i] = new TH1D(Form("energy_hists%d", i), Form("Uncorrected energy: L%d", i),  2000, 0, 2000);
        slowfasttimediff_energy_hists[i] = new TH2D(Form("slowfasttimediff_energy_hists%d", i), Form("Uncorrected time offset vs energy: L%d", i),  2000, 0, 2000, 2000, 0, 2000);
        energy_slowfasttimediff_hists[i] = new TH2D(Form("energy_slowfasttimediff_hists%d", i), Form("Uncorrected energy vs time offset: L%d", i),  2000, 0, 2000, 2000, 0, 2000);
        timeOffsetCorrectedHists[i] = new TH1D(Form("timeOffsetCorrectedHists%d", i), Form("Corrected time offset: L%d", i),  1500, -750, 750);
        timeEnergyOffsetCorrectedHists[i] = new TH2D(Form("timeEnergyOffsetCorrectedHists%d", i), Form("Corrected time offset vs energy: L%d", i),  1500, -750, 750, 2000, 0, 2000);
    }
    
    // ________________________________________________________________
    for (int i = 0; i < NumDetectors; i++) 
    {
        printf("\nCalculating time offset for detector %d...\n", i);
        for (uint64_t entry = 0; entry < numEntries; entry++) 
        {
            copyTree->GetEntry(entry);

            //if (timeFastCp[i] == 0) continue;
            if (timeFast[i] < 0) printf("\n NEG timeFastL%d: %f -- This should *never* happen\n", i, timeFast[i]);

            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;  // Initialize nextTimeSlow
            Double_t deltaTime = 0;     // Initialize deltaTime

            for (uint64_t nextEntry = entry + 1; nextEntry < numEntries; nextEntry++) 
            {
                copyTree->GetEntry(nextEntry);                
                //if (timeSlowCp[i] == 0) continue;
                nextTimeSlow = timeSlowCp[i];
                deltaTime = nextTimeSlow - refTimeFast;
                if (deltaTime > 10000.0) break; // Stop searching if time difference exceeds 1 microsecond
                // printf("\nDetector %d: Time Offset Forward = %f\n", i, deltaTime);
                slowfasttimediff_hists[i]->Fill(deltaTime);
                slowfasttimediff_energy_hists[i]->Fill(deltaTime, energySlowCp[i]);
                energy_slowfasttimediff_hists[i]->Fill(energySlowCp[i], deltaTime);
                energy_hists[i]->Fill(energySlowCp[i]);
                //if ((i == 2) && (energySlowCp[i] > 1) && (energySlowCp[i] < 600)) printf("\nDetector %d: Time Offset Forward = %f for energy = %f\n", i, deltaTime, energySlowCp[i]);
                {
                    timeOffsetCorrectedHists[i]->Fill(deltaTime);
                    timeEnergyOffsetCorrectedHists[i]->Fill(deltaTime, energySlowCp[i]);
                }
            }
        }
        
        TCanvas* c1 = new TCanvas(Form("c1_%d", i), Form("Uncorrected time offset: L%d", i), 800, 600);
        slowfasttimediff_hists[i]->GetXaxis()->SetRangeUser(0, 2000);
        slowfasttimediff_hists[i]->SetStats(0);
        slowfasttimediff_hists[i]->SetFillColor(kBlue);
        slowfasttimediff_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_hists[i]->GetYaxis()->SetTitle("Counts/ns");
        slowfasttimediff_hists[i]->Draw();
        c1->SaveAs(Form("%s/time_offset_hist_L%d.root", dirPathplots.c_str(), i));
        c1->Write(Form("time_offset_hist_L%d", i));
        delete c1;
        energy_slowfasttimediff_hists[i]->GetYaxis()->SetRangeUser(0, 1500);

        dir_banana->cd();
        bananafit[i] = banana(i, slowfasttimediff_energy_hists[i], energy_slowfasttimediff_hists[i], dirPathplots);
        
        TCanvas * c = new TCanvas("c", "banana", 800, 600);
        energy_slowfasttimediff_hists[i]->Draw();
        bananafit[i]->Draw("same");
        c->SaveAs(Form("%s/banana_2Dfitfunction_L%d.root", dirPathplots.c_str(), i));
        c->Write(Form("banana_2Dfitfunction_L%d", i));
        delete c;
    }

    dir_aligned->cd();
    
    TTree* correctedTree = new TTree("LaBrDataCorrected", "LaBrDataCorrected"); // I had to create another tree as you cannot read in entries and write to the same tree :( 
    for (int i = 0; i < NumDetectors; i++) 
    {
        correctedTree->Branch(Form("slowECalibL%d", i), &energySlowCp[i], Form("slowECalibL%d/D", i));
        correctedTree->Branch(Form("timeSL%d", i), &timeSlowCp[i], Form("timeSL%d/D", i));
        correctedTree->Branch(Form("timeFL%d", i), &timeFastCp[i], Form("timeFL%d/D", i));
        correctedTree->Branch(Form("timeSCorrectedL%d", i), &timeSlowCorrected[i], Form("timeSCorrectedL%d/D", i)); //new Branch
    }
    correctedTree->Branch("timeRF", &timeFastCp[RFIndex], "timeRF/D");
    correctedTree->Branch("fastEPOLARIS", &fastEnergyPOLARISCp[0], "fastEPOLARIS/D");
    correctedTree->Branch("fastTPOLARIS", &fastTimePOLARISCp[0], "fastTPOLARIS/D");
    
    const int energyRangeMin = 300;
    const int energyRangeMax = 2000;
    const int numEnergyBins = energyRangeMax - energyRangeMin + 1;

    // Create arrays to store energy values and time offsets for each detector
    std::vector<std::vector<double>> energyValues(NumDetectors, std::vector<double>(numEnergyBins));
    std::vector<std::vector<double>> timeOffsets(NumDetectors, std::vector<double>(numEnergyBins));

    // Loop over a range of energies and calculate time offsets for each detector
    for (int energy = energyRangeMin; energy <= energyRangeMax; energy++) 
    {
        for (int i = 0; i < NumDetectors; i++) {
            // Calculate time offset using the bananafit function for each detector
            int energyIndex = energy - energyRangeMin;  // Index in the arrays
            energyValues[i][energyIndex] = energy;
            timeOffsets[i][energyIndex] = bananafit[i]->Eval(static_cast<double>(energy));
        }
    }

    // Loop through the entries to update timeSlowCorrected[i] and fill correctedTree
    for (uint64_t entry = 0; entry < numEntries; entry++) 
    {
        for (int i = 0; i < NumDetectors; i++) 
        {
            copyTree->GetEntry(entry);
            //printf("\nDetector %d: Slow Time = %f, Slow Energy = %f\n", i, timeSlowCp[i], energySlowCp[i]);
            if ((timeSlowCp[i] <=0) || (energySlowCp[i] <=170))continue;
             // Look up the time offset from the array based on energy
            int energyIndex = static_cast<int>(energySlowCp[i]) - energyRangeMin;
            if (energyIndex >= 0 && energyIndex < numEnergyBins) 
            {
                timeOffsets[i][energyIndex]; 
                timeSlowCorrected[i] = timeSlowCp[i] - timeOffsets[i][energyIndex];
            }
            else 
            {
                timeSlowCorrected[i] = 0;
            }

            //printf("\nDetector %d: Time Corrected = %f\n", i, timeSlowCorrected[i]);
            correctedTree->Fill();
        }   
    }


    for(int i = 0; i < NumDetectors; i++) 
    {
        printf("\nCalculating corrected time offset -> only for the purpose of plotting the histograms for detector %d...\n", i);
        for (uint64_t entry = 0; entry < numEntries; entry++) 
        {
            correctedTree->GetEntry(entry);
            if (timeFastCp[i] == 0) continue;
            Double_t refTimeFast = timeFastCp[i];
            Double_t nextTimeSlow = 0;  // Initialize nextTimeSlow
            Double_t deltaTime = 0; // Initialize deltaTime
            Double_t nextEnergySlow = 0; // Initialize nextEnergySlow

            for (uint64_t nextEntry = entry + 1; nextEntry < numEntries; nextEntry++) 
            {
                correctedTree->GetEntry(nextEntry);            
                if (energySlowCp[i] == 0) continue;
                nextEnergySlow = energySlowCp[i];
                nextTimeSlow = timeSlowCorrected[i];
                //printf("\nDetector %d: Fast Time = %f, Slow Time = %f\n", i, refTimeFast, nextTimeSlow);
                deltaTime = (nextTimeSlow - refTimeFast);
                if (deltaTime > 100.0) break; // Stop searching if time difference exceeds 1 microsecond
                //printf("\nDetector %d: Corrected Time Offset = %f\n", i, deltaTime);
                timeOffsetCorrectedHists[i]->Fill(deltaTime);
                timeEnergyOffsetCorrectedHists[i]->Fill(deltaTime, nextEnergySlow);
            }
        }
    }

    for (int i = 0; i < NumDetectors; i++) 
    {
        timeOffsetCorrectedHists[i]->Write();
        timeEnergyOffsetCorrectedHists[i]->Write();
        
        TCanvas* c2 = new TCanvas(Form("c2_%d", i), Form("Uncorrected time offset vs energy: L%d", i), 800, 600);
        slowfasttimediff_energy_hists[i]->SetStats(0);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetRangeUser(0, 2000);
        slowfasttimediff_energy_hists[i]->SetMarkerStyle(20);
        slowfasttimediff_energy_hists[i]->SetMarkerSize(0.8);
        slowfasttimediff_energy_hists[i]->GetXaxis()->SetTitle("Time Slow - Time Fast (ns)");
        slowfasttimediff_energy_hists[i]->GetYaxis()->SetTitle("Energy (keV)");
        //slowfasttimediff_energy_hists[i]->SetMaximum(1200);
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
        //timeEnergyOffsetCorrectedHists[i]->SetMaximum(1200);
        timeEnergyOffsetCorrectedHists[i]->Draw("COLZ");
        c4->SaveAs(Form("%s/corrected_energy_vs_time_offset_banana_hist_L%d.root", dirPathplots.c_str(), i));
        c4->Write(Form("corrected_energy_vs_time_offset_banana_hist_L%d", i));

        delete c3;
        delete c4;  

        slowfasttimediff_hists[i]->Write();
        energy_hists[i]->Write();
        slowfasttimediff_energy_hists[i]->Write();
        energy_slowfasttimediff_hists[i]->Write();      
    }

    // ________________________________________________________________
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

    // Vector to store entry data
    std::vector<EntryData> entryDataVec;

    // Align the slow time and energy data for each detector based on the calculated offsets
    for (int i = 0; i < NumDetectors +2; i++) 
    {
        printf("\nAligning data for detector %d...\n", i);
        for (uint64_t entry = 0; entry < numEntries; entry++) 
        {
            correctedTree->GetEntry(entry);

            if (i == RFIndex) 
            {
                detectorID = RFIndex;
                timeGlobal = timeFastCp[i];
                alignedTimeFast = timeFastCp[i];
                alignedEnergyFast = 0;
                alignedTimeSlow = 0;
                alignedEnergySlow = 0;
                // printf("\nRF Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i == POLARISIndex) 
            {
                detectorID = POLARISIndex;
                timeGlobal = fastTimePOLARISCp[0];
                alignedTimeFast = fastTimePOLARISCp[0];
                alignedEnergyFast = fastEnergyPOLARISCp[0];
                alignedTimeSlow = 0;
                alignedEnergySlow = 0;
                // printf("\nPOLARIS Time entry %lld: | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCp = %f, energyFast = %f | , energySlowCp = %f\n", entry, timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
            }
            else if (i < RFIndex)
            {
                detectorID = i;
                // printf("\nDetector %d: timeFastCp = %f | , timeSlowCorrected = %f\n", i, timeFastCp[i], timeSlowCorrected[i]);
                if ((timeFastCp[i] == 0) && (timeSlowCorrected[i] != 0)) 
                {
                    timeGlobal = timeSlowCorrected[i];
                    alignedTimeFast = timeFastCp[i];
                    alignedEnergyFast = 0;
                    alignedTimeSlow = timeSlowCorrected[i];
                    alignedEnergySlow = energySlowCp[i];
                    //printf("Loop1 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrected = %f, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrected[i] == 0)) 
                {
                    timeGlobal = timeFastCp[i];
                    alignedTimeFast = timeFastCp[i];
                    alignedEnergyFast = 0;
                    alignedTimeSlow = timeSlowCorrected[i];
                    alignedEnergySlow = energySlowCp[i];
                    // printf("Loop2 | , timeGlobal = %f | , timeFastCp = %f | , timeSlowCorrected = %f\n, energyFast = %f | , energySlowCp = %f\n", timeGlobal[0], alignedTimeFast[0], alignedTimeSlow[0], alignedEnergyFast[0], alignedEnergySlow[0]);
                }
                else if ((timeFastCp[i] != 0) && (timeSlowCorrected[i] != 0)) 
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

    // Sort the entryDataVec based on timeGlobal in ascending order
    std::sort(entryDataVec.begin(), entryDataVec.end(), [](const EntryData& a, const EntryData& b) {return a.timeGlobal < b.timeGlobal;});

    // Set the branches for the aligned data in the new copyTree
    TTree *alignedTree1 = new TTree("AlignedData", "AlignedData");
    alignedTree1->Branch("detectorID", &detectorID, "detectorID/I");
    alignedTree1->Branch("globalTime", &timeGlobal, "globalTime/D");
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
        timeGlobal = entryData.timeGlobal;
        alignedTimeFast = entryData.alignedTimeFast;
        alignedEnergyFast = entryData.alignedEnergyFast;
        alignedTimeSlow = entryData.alignedTimeSlow;
        alignedEnergySlow = entryData.alignedEnergySlow;

        alignedTree1->Fill();
    }

    // TTree *alignedTree = new TTree("OrderedData", "OrderedData");
    // alignedTree->Branch("detectorID", &detectorID, "detectorID/I");
    // alignedTree->Branch("globalTime", &timeGlobal[0], "globalTime/D");
    // alignedTree->Branch("timeF", &alignedTimeFast[0], "timeF/D");
    // alignedTree->Branch("energyF", &alignedEnergyFast[0], "energyF/D");
    // alignedTree->Branch("timeS", &alignedTimeSlow[0], "timeS/D");
    // alignedTree->Branch("energyS", &alignedEnergySlow[0], "energyS/D");

    // // Copy data from the original tree to the new tree while sorting by globalTime
    // alignedTree1->BuildIndex("globalTime");
    // TTreeIndex *index = (TTreeIndex*)alignedTree1->GetTreeIndex();
    // uint64_t nEntries = (alignedTree1->GetEntries());

    // printf("\n Sorting the data by globalTime...\n")
    // for (uint64_t i = 0; i < nEntries; ++i) {
    //     uint64_t originalEntry = index->GetIndex()[i];
    //     alignedTree1->GetEntry(originalEntry);
    //     alignedTree->Fill();
    // }

    uint64_t nEntries = (alignedTree1->GetEntries());

    uint64_t lastEntry;
    uint64_t firstEntry;
    Double_t lastTime;
    Double_t firstTime;

    for (uint64_t entry = 0; entry < nEntries; entry++) 
    {
        alignedTree1->GetEntry(entry);
        if (timeGlobal == 0 )continue;
        firstTime = timeGlobal;
        break;
    }
   
    lastEntry = nEntries - 1;
    for (uint64_t entry = lastEntry; entry >= 0; entry--) 
    {
        alignedTree1->GetEntry(entry);
        lastTime = timeGlobal;
        break;
    }
    
    alignedTree1->Write();

    Double_t runLength = (lastTime - firstTime)*1e-9/60;
    printf("\n First entry timeGlobal = %f\n", firstTime);
    printf("\n Last entry timeGlobal = %f\n", lastTime);
    printf("\nRun Length (According to timeGlobal variable. Used to check if it's in sequential time) = %f\n", runLength);
}

// ________________________________________ Main function ________________________________________
int main(int argc, char* argv[], char* argv2[])
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
    std::string dataType = argv[2];

    // the below bool is used to determine whether the data is beam or angle data
    // this is used to determine if the fast time detector intercalibration should be applied 
    // it is not applied for source data
    bool beam = false;
    if (dataType == "beam") 
    {
        beam = true;
    }
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

    // Reserve memory for the vectors
    energySlow.reserve(NumDetectors);
    timeSlow.reserve(NumDetectors);
    timeFast.reserve(NumDetectors + 1);
    fastEnergyPOLARIS.reserve(1);
    fastTimePOLARIS.reserve(1);

    // ________________________________________________________________________________________________________________________

    // Set the number of threads to use
    int NumThreads = 4;
    ROOT::EnableImplicitMT(NumThreads);

    // ________________________________________________________________________________________________________________________
    
    // Open the ROOT file
    TFile*file = new TFile(inputFileName.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file->Get("LaBrData"));
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");

    if (!file->IsOpen()) 
    {
        std::cout << "Error: File not open" << std::endl;
        return 1;
    }

    if (!outputFile->IsOpen()) 
    {
        std::cout << "Error: Output file not open" << std::endl;
        return 1;
    }

    TDirectory *dir_fasttimecalib = outputFile->mkdir("dir_fasttimecalib");
    TDirectory *dir_banana = outputFile->mkdir("dir_banana");
    TDirectory *dir_aligned = outputFile->mkdir("dir_aligned");

    // ________________________________________________________________________________________________________________________
    // FUNCTIONS

    runduration(tree);
    fasttime_calibration_constants(tree, dirPathplots, dir_fasttimecalib, beam);
    align_offset(tree, dirPathplots, dir_banana, dir_aligned, mean_fasttimediff);

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
