/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms and TTrees*/

/* Tested up to ROOT 6.268/08 */

/* Changes to implement 4 detectors and generate event data */
/* based on time window */

/* Version 1:  (c) S. Hart (01/06/2023) */

/* This code reads in the R63_sorted.root file produced by the time-ordering code sort3labr.c */
/* The data is ordered by ascending time. This is done by restructuring the data to have the format detectorID (0,1,2,3,RF,POLARIS),timeF, slowT, slowE */
/* The time offset between the fast and slow times is calculated and the slow time & energy is then shifted by this value to align the data. */
/* The code then starts at the first timestamp and takes the time difference between this and the next timestamp, incrementing until the time difference exceeds some "window". */
/* When this window is exceeded, all of the data that falls within the window is one event. */
/* The cyclotron beam RF gate is then applied to in-sync events to reduce the statictics of out of sync data events.*/

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort4labr.C -o exe4 `root-config --cflags --libs` -lSpectrum 
     ./exe4 ~/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63_sorted.root 
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

const int NumDetectors = 4;
const int RFIndex = 4;
const int POLARISIndex = 5;

time_t start, end;

void event_windowing(TFile *outputFile, TDirectory *dir_rf, TDirectory * dir_sorted, TTree* tree, std::string dirPathplots)
{
    dir_sorted->cd();
    TTree *sortedTree = tree->CloneTree();
    sortedTree->Write();
    printf("\n Tree copied into directory...\n");

    dir_rf->cd();
    printf("\n Windowing data into events...\n");
    int64_t item;

    // Initialise the sortedTree variables
    Int_t sortedDetectorID;
    Double_t sortedGlobalTime;
    Double_t sortedTimeF;
    Double_t sortedEnergyF;
    Double_t sortedTimeS;
    Double_t sortedEnergyS;

    tree->SetBranchAddress("detectorID", &sortedDetectorID);
    tree->SetBranchAddress("globalTime", &sortedGlobalTime);
    tree->SetBranchAddress("timeF", &sortedTimeF);
    tree->SetBranchAddress("energyF", &sortedEnergyF);
    tree->SetBranchAddress("timeS", &sortedTimeS);
    tree->SetBranchAddress("energyS", &sortedEnergyS);
    nEntries = tree->GetEntries();

    // Declare vectors to store the event data
    std::vector<double> FT[6], SE[4], ST[4], FE[1];
    double FTL[5] = {0,0,0,0,0}, SEL[4] = {0,0,0,0}, STL[4] = {0,0,0,0} , POLFT[1] = {0}, POLFE[1] = {0};
    int fastcount[4]={0,0,0,0}, slowcount[4]={0,0,0,0}, fastslowcount[4]={0,0,0,0};
    double tdF;

    TTree *EventTree = new TTree("EventTree", "EventTree");
    EventTree->Branch("RFT", &FTL[4], "RFT/D");
    EventTree->Branch("FT0", &FTL[0], "FT0/D");
    EventTree->Branch("SE0", &SEL[0], "SE0/D");
    EventTree->Branch("ST0", &STL[0], "ST0/D");
    EventTree->Branch("FT1", &FTL[1], "FT1/D");
    EventTree->Branch("SE1", &SEL[1], "SE1/D");
    EventTree->Branch("ST1", &STL[1], "ST1/D");
    EventTree->Branch("FT2", &FTL[2], "FT2/D");
    EventTree->Branch("SE2", &SEL[2], "SE2/D");
    EventTree->Branch("ST2", &STL[2], "ST2/D");
    EventTree->Branch("FT3", &FTL[3], "FT3/D");
    EventTree->Branch("SE3", &SEL[3], "SE3/D");
    EventTree->Branch("ST3", &STL[3], "ST3/D");
    EventTree->Branch("POLFT", &POLFT[0], "POLFT/D");
    EventTree->Branch("POLFE", &POLFE[0], "POLFE/D");

    TH1F** eventCounts = new TH1F*[4];
    TH1D** slowE = new TH1D*[4];
    TH1D** insyncE = new TH1D*[4];
    TH1D** fastTD4 = new TH1D*[4];
    TH2D** slowEfastTD4 = new TH2D*[4];
    TH1D** oneMeVgateTD4 = new TH1D*[4];
    TH1D** threeMeVgateTD4 = new TH1D*[4];
    TH1D** fiveMeVgateTD4 = new TH1D*[4];
    TH2D** slowEfastTD4_1MeV = new TH2D*[4];
    TH2D** slowEfastTD4_3MeV = new TH2D*[4];
    TH2D** slowEfastTD4_5MeV = new TH2D*[4];

    printf("\n Making RF histograms...\n");
    for (int j = 0; j <= 3; j++)
    {
        eventCounts[j] = new TH1F(TString::Format("Event Counts L%2d", j), TString::Format("Event Counts for L%2d", j), 3, 0.5, 3.5);
        slowE[j] = new TH1D(TString::Format("slowE L%2d", j), TString::Format("Slow Signal Energy Spectrum for L%2d", j), 8000, 0, 8000);
        insyncE[j] = new TH1D(TString::Format("insyncE L%2d", j), TString::Format("In Sync Signal Energy Spectrum for L%2d", j), 8000, 0, 8000);
        fastTD4[j] = new TH1D(TString::Format("fastTD4 L%2d", j), TString::Format("Fast Time Difference for L%2d", j), 1600, -800, 800);
        slowEfastTD4[j] = new TH2D(TString::Format("slowEfastTD4 L%2d", j), TString::Format("Slow Signal Energy Vs Fast Time Difference for L%2d", j), 1600, -800, 800, 8000, 0, 8000);
        oneMeVgateTD4[j] = new TH1D(TString::Format("oneMeVgateTD4 L%2d", j), TString::Format("Fast Time Difference for RF-L%2d, slow energy gate [1000;2000)keV]", j), 1600, -800, 800);
        threeMeVgateTD4[j] = new TH1D(TString::Format("threeMeVgateTD4 L%2d", j), TString::Format("Fast Time Difference for RF-L%2d, slow energy gate [3000;4000)keV]", j), 1600, -800, 800);
        fiveMeVgateTD4[j] = new TH1D(TString::Format("fiveMeVgateTD4 L%2d", j), TString::Format("Fast Time Difference for RF-L%2d, slow energy gate [5000;6000)keV]", j), 1600, -800, 800);
        slowEfastTD4_1MeV[j] = new TH2D(TString::Format("slowEfastTD4_1MeV L%2d", j), TString::Format("Slow Signal Energy Vs Fast Time Difference for L%2d, slow energy gate [1000;2000)keV]", j), 1600, -800, 800, 1000, 1000, 2000);
        slowEfastTD4_3MeV[j] = new TH2D(TString::Format("slowEfastTD4_3MeV L%2d", j), TString::Format("Slow Signal Energy Vs Fast Time Difference for L%2d, slow energy gate [3000;4000)keV]", j), 1600, -800, 800, 1000, 3000, 4000);
        slowEfastTD4_5MeV[j] = new TH2D(TString::Format("slowEfastTD4_5MeV L%2d", j), TString::Format("Slow Signal Energy Vs Fast Time Difference for L%2d, slow energy gate [5000;6000)keV]", j), 1600, -800, 800, 1000, 5000, 6000);
    }
    printf("\n RF histograms made...\n");
    Long64_t eventloop = 0;
    int64_t TSdiff, TS, TSinit, TSfirst, eventWindow = 1000, counter=0;

    // Loop over the entries in the aligned tree
    for (Long64_t entry = 0; entry < nEntries-1; entry++) 
    {
        // printf("\nEvent %lld: | , Entry %lld: | , TSdiff = %lld | , TS = %lld | , TSinit = %lld\n", eventloop, entry, TSdiff, TS, TSinit);
        tree->GetEntry(entry);
        TS = timeGlobal[0];

        if (counter == 0) 
        {
            TSfirst = TS;
            counter = 1;
        }

        if (item == 0) TSinit = TS;

        TSdiff = (TS - TSinit)/10;

        if (TSdiff <= eventWindow) 
        {
            //printf("\nEvent %lld: | , TSdiff = %lld | , detectorID = %d", eventloop, TSdiff, detectorID);
            if (detectorID >= 0 && detectorID <= 3)
            {
                if (FT[detectorID].size() == 0)
                {
                    FT[detectorID].push_back(alignedTimeFast[0]);
                    SE[detectorID].push_back(alignedEnergySlow[0]);
                    ST[detectorID].push_back(alignedTimeSlow[0]);
                }
            }
            else if (detectorID == 4) 
            {
                if (FT[detectorID].size() == 0)
                {
                    FT[detectorID].push_back(alignedTimeFast[0]);
                }
            }
            else if (detectorID == 5) 
            {
                if (FT[detectorID].size() == 0)
                {
                    FT[detectorID].push_back(alignedTimeFast[0]);
                    FE[detectorID].push_back(alignedEnergyFast[0]);
                }
            }
            if (item == 0) TSinit=TS;
            item++;
            //printf("\nItem %lld: | , TSdiff = %lld | , detectorID = %d, FT[detectorID].size() = %lu\n, ST[detectorID].size() = %lu\n, SE[detectorID].size() = %lu\n", item, TSdiff, detectorID, FT[detectorID].size(), ST[detectorID].size(), SE[detectorID].size());
        }

        if (TSdiff > eventWindow) 
        {
            //printf("\nEnd of event %lld\n", eventloop);
            for (int j = 0; j < 5; j++) 
            {
                if (FT[j].size() == 0) FT[j].push_back(0); // if there is no data for a detector, push back a 0
                if (SE[j].size() == 0) SE[j].push_back(0);
                if (ST[j].size() == 0) ST[j].push_back(0);
            }
            if (FE[0].size() == 0) FE[0].push_back(0); 
            if (FT[5].size() == 0) FT[5].push_back(0);

            // Fill the vectors with the data for each detector
            FTL[0] = FT[0][0];
            FTL[1] = FT[1][0];
            FTL[2] = FT[2][0];
            FTL[3] = FT[3][0];
            FTL[4] = FT[4][0];
            SEL[0] = SE[0][0];
            SEL[1] = SE[1][0];
            SEL[2] = SE[2][0];
            SEL[3] = SE[3][0];
            STL[0] = ST[0][0];
            STL[1] = ST[1][0];
            STL[2] = ST[2][0];
            STL[3] = ST[3][0];
            POLFT[0] = FT[5][0];
            POLFE[0] = FE[0][0];
            //printf("\nEvent %lld: | , FTL[0] = %f | , FTL[1] = %f | , FTL[2] = %f | , FTL[3] = %f | , FTL[4] = %f | , SEL[0] = %f | , SEL[1] = %f | , SEL[2] = %f | , SEL[3] = %f | , STL[0] = %f | , STL[1] = %f | , STL[2] = %f | , STL[3] = %f | , POLFT[0] = %f | , POLFE[0] = %f\n", eventloop, FTL[0], FTL[1], FTL[2], FTL[3], FTL[4], SEL[0], SEL[1], SEL[2], SEL[3], STL[0], STL[1], STL[2], STL[3], POLFT[0], POLFE[0]);

            EventTree->Fill();

            for (int j = 0; j <= 3; j++)
            {
                std::cout << "L" << j << "| energySlow[j] " << SEL[j] << "| timeFast[j] " << FTL[j] << "| timeSlow[j] " << STL[j] << std::endl;
                if (SEL[j] >= 1000)
                {
                    if (FTL[j] != 0 && STL[j] == 0) fastcount[j]++;
                    else if (FTL[j] == 0 && STL[j] != 0) slowcount[j]++;
                    else if (FTL[j] != 0 && STL[j] != 0) fastslowcount[j]++;

                    eventCounts[j]->Fill(1,fastcount[j]);
                    eventCounts[j]->Fill(2,slowcount[j]);
                    eventCounts[j]->Fill(3,fastslowcount[j]);

                    eventCounts[j]->GetXaxis()->SetBinLabel(1,"Only Fast Time");
                    eventCounts[j]->GetXaxis()->SetBinLabel(2,"Only Slow Time");
                    eventCounts[j]->GetXaxis()->SetBinLabel(3,"Coinc Fast+Slow Time");
                }                      
                slowE[j]->Fill(SEL[j]);
            }

            for (int j = 0; j <= 4; j++) 
            { 
                for (int k = 0; k <= 3; k++) 
                {
                    if (j != k ) 
                    {
                        tdF=(FTL[j]-FTL[k]); // Measured in 1e13s which is 100ns
                        tdF=(tdF/10);
                        if (j == 4) 
                        {
                            slowEfastTD4[k]->Fill(tdF, SEL[k]);
                            fastTD4[k]->Fill(tdF);
                            if ((SEL[k]>=1000) && (SEL[k]<2000))
                            {
                                oneMeVgateTD4[k]->Fill(tdF);
                                slowEfastTD4_1MeV[k]->Fill(tdF, SEL[k]);
                            }
                            else if ((SEL[k]>=3000) && (SEL[k]<4000)) 
                            {
                                threeMeVgateTD4[k]->Fill(tdF);
                                slowEfastTD4_3MeV[k]->Fill(tdF, SEL[k]);
                            }
                            else if ((SEL[k]>=5000) && (SEL[k]<6000)) 
                            {
                                fiveMeVgateTD4[k]->Fill(tdF);
                                slowEfastTD4_5MeV[k]->Fill(tdF, SEL[k]);
                            }
                            else continue;
                        }

                    }
                }                        
            }
            
            // Clear the vectors
            for (int i = 0; i < 5; i++) 
            {
                FT[i].clear();
                SE[i].clear();
                ST[i].clear();
            }
            FT[5].clear();
            FE[0].clear();
            //Reinitialise arrays
            for (int i = 0; i <= 4; i++) 
            {
                FTL[i] = 0;
                SEL[i] = 0;
                STL[i] = 0;
            }
            POLFT[0] = 0;
            POLFE[0] = 0;

            delete[] eventCounts;
            item = 0;
            eventloop++;
        }

        TSinit=TS;

    }

    for (int j = 0; j <= 3; j++)
    {
        slowE[j]->Write();
        fastTD4[j]->Write();
        slowEfastTD4[j]->Write();
        eventCounts[j]->Write();
        oneMeVgateTD4[j]->Write();
        threeMeVgateTD4[j]->Write();
        fiveMeVgateTD4[j]->Write();
        slowEfastTD4_1MeV[j]->Write();
        slowEfastTD4_3MeV[j]->Write();
        slowEfastTD4_5MeV[j]->Write();
    }
    
    for (int j = 0; j <= 3; j++)
    {
        // Plot time difference histogram
        TCanvas* c1 = new TCanvas(Form("c1_%d", j), Form("Time Difference: RF -  Detector %d", j), 800, 600);
        fastTD4[j]->SetStats(0);
        fastTD4[j]->SetFillColor(kBlue);
        fastTD4[j]->Draw();
        c1->SaveAs(Form("%s/time_diff_hist_L%d.root", dirPathplots.c_str(), j));
        c1->Write(Form("time_diff_hist_L%d", j));

        // Plot time difference vs alignedEnergySlow histogram
        TCanvas* c2 = new TCanvas(Form("c2_%d", j), Form("Time Difference RF - Detector vs Aligned Energy: Detector %d", j), 800, 600);
        slowEfastTD4[j]->SetStats(0);
        slowEfastTD4[j]->SetMarkerStyle(20);
        slowEfastTD4[j]->SetMarkerSize(0.8);
        slowEfastTD4[j]->GetXaxis()->SetTitle("Time Difference (Fast) RF - LaBr_3:Ce (ns)");
        slowEfastTD4[j]->GetYaxis()->SetTitle("Aligned Energy (Slow) (keV)");        
        slowEfastTD4[j]->Draw("COLZ");
        c2->SaveAs(Form("%s/time_diff_vs_energy_hist_L%d.root", dirPathplots.c_str(), j));
        c2->Write(Form("time_diff_vs_energy_hist_L%d", j));

        TCanvas* c3 = new TCanvas(Form("c3_%d", j), Form("Event Counts: Detector %d", j), 800, 600);
        eventCounts[j]->SetStats(0);
        eventCounts[j]->SetFillColor(kBlue);
        eventCounts[j]->Draw("BAR");
        c3->SaveAs(Form("%s/event_counts_L%d.root", dirPathplots.c_str(), j));
        c3->Write(Form("event_counts_L%d", j));

        TCanvas* c4 = new TCanvas(Form("c4_%d", j), Form("Time Difference (Fast) RF - LaBr_3:Ce (ns): Detector %d vs Slow Energy (gated [1000;2000)keV)", j), 800, 600);
        slowEfastTD4_1MeV[j]->SetStats(0);
        slowEfastTD4_1MeV[j]->SetMarkerStyle(20);
        slowEfastTD4_1MeV[j]->SetMarkerSize(0.8);   
        slowEfastTD4_1MeV[j]->GetXaxis()->SetTitle("Time Difference (Fast) RF - LaBr_3:Ce (ns)");
        slowEfastTD4_1MeV[j]->GetYaxis()->SetTitle("Aligned Energy (Slow) (keV)");
        slowEfastTD4_1MeV[j]->Draw("COLZ");
        c4->SaveAs(Form("%s/time_diff_vs_energy_hist_L%d_1MeV.root", dirPathplots.c_str(), j));
        c4->Write(Form("time_diff_vs_energy_hist_L%d_1MeV", j));

        TCanvas* c5 = new TCanvas(Form("c5_%d", j), Form("Time Difference (Fast) RF - LaBr_3:Ce (ns): Detector %d vs Slow Energy (gated [3000;4000)keV)", j), 800, 600);
        slowEfastTD4_3MeV[j]->SetStats(0);
        slowEfastTD4_3MeV[j]->SetMarkerStyle(20);
        slowEfastTD4_3MeV[j]->SetMarkerSize(0.8);
        slowEfastTD4_3MeV[j]->GetXaxis()->SetTitle("Time Difference (Fast) RF - LaBr_3:Ce (ns)");
        slowEfastTD4_3MeV[j]->GetYaxis()->SetTitle("Aligned Energy (Slow) (keV)");
        slowEfastTD4_3MeV[j]->Draw("COLZ");
        c5->SaveAs(Form("%s/time_diff_vs_energy_hist_L%d_3MeV.root", dirPathplots.c_str(), j));
        c5->Write(Form("time_diff_vs_energy_hist_L%d_3MeV", j));

        TCanvas* c6 = new TCanvas(Form("c6_%d", j), Form("Time Difference (Fast) RF - LaBr_3:Ce (ns): Detector %d vs Slow Energy (gated [5000;6000)keV)", j), 800, 600);
        slowEfastTD4_5MeV[j]->SetStats(0);
        slowEfastTD4_5MeV[j]->SetMarkerStyle(20);
        slowEfastTD4_5MeV[j]->SetMarkerSize(0.8);
        slowEfastTD4_5MeV[j]->GetXaxis()->SetTitle("Time Difference (Fast) RF - LaBr_3:Ce (ns)");
        slowEfastTD4_5MeV[j]->GetYaxis()->SetTitle("Aligned Energy (Slow) (keV)");
        slowEfastTD4_5MeV[j]->Draw("COLZ");
        c6->SaveAs(Form("%s/time_diff_vs_energy_hist_L%d_5MeV.root", dirPathplots.c_str(), j));
        c6->Write(Form("time_diff_vs_energy_hist_L%d_5MeV", j));

        delete c1;
        delete c2;
        delete c3;
        delete c4;
        delete c5;
        delete c6;
    }
    EventTree->Write();

    // Calculating total time taken by the program.
    double TimeDiff = (double) (TS - TSfirst)*(1e-8)/60.0;
    printf("The first time stamp is (e-8 s): %" PRIu64 "\n", TSfirst);
    printf("The last time stamp is (e-8 s): %" PRIu64 "\n", TS);
    printf("The time difference is (min): %f\n", TimeDiff);
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
    std::string runStr = inputFileName.substr(inputFileName.find_last_of("/") + 2);
    std::string dirPath = inputFileName.substr(0, inputFileName.find_last_of("/") + 1);
    int run = std::atoi(runStr.c_str());
    std::string outputFileName = dirPath + "R" + std::to_string(run) + "_final.root";
    std::cout << "Run: " << run << std::endl;
    std::cout << "Directory Path: " << dirPath << std::endl;
    std::cout << "Output File Name: " << outputFileName << std::endl;

    // ________________________________________________________________________________________________________________________
    // Open the ROOT file
    TFile*file = new TFile(inputFileName.c_str(), "OPEN");
    if (!file->IsOpen()) {
        std::cout << "Error: File not open" << std::endl;
        return 1;
    }    

    TTree* tree = static_cast<TTree*>(file->Get("SortedData"));
    if (!tree) {
        std::cout << "Error: Tree not found" << std::endl;
        return 1;
    }

    std::string dirPathplots = dirPath + "plots/sort4labr/";
    // check to see if dirPath contains a folder called plots. if not, create it
    if (dirPath.find("plots/sort4labr/") == std::string::npos)
    {
        std::string dirPath1 = dirPath + "plots/";
        std::string dirPathplots = dirPath1 + "sort4labr/";
        std::string command = "mkdir " + dirPath1;
        std::system(command.c_str());
        command = "mkdir " + dirPathplots;
        system(command.c_str());
    }

    // Open the output ROOT file
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cout << "Error: Output file not open" << std::endl;
        return 1;
    }
    TDirectory *dir_rf = outputFile->mkdir("dir_rf");
    TDirectory *dir_sorted = outputFile->mkdir("dir_sorted");

    // ________________________________________________________________________________________________________________________

    event_windowing(outputFile, dir_rf, dir_sorted, tree, dirPathplots);

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
