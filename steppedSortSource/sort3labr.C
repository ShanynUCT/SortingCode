/* Advanced sorting code for TDR format data */
/* Produces ROOT histograms and TTrees*/
/* Tested up to ROOT 6.268/08 */
/* Version 1:  (c) S. Hart (19/06/2023) */

/* This code reads in the RXX_aligned.root file produced by sort2labr.c */

/* 
    The tree has the following stucture:   
    TTree *tree = new TTree("AlignedData", "AlignedData");
    tree->Branch("detectorID", &detectorID, "detectorID/I");
    tree->Branch("globalTime", &timeGlobal, "globalTime/D");
    tree->Branch("timeF", &alignedTimeFast, "timeF/D");
    tree->Branch("energyF", &alignedEnergyFast, "energyF/D");
    tree->Branch("timeS", &alignedTimeSlow, "timeS/D");
    tree->Branch("energyS", &alignedEnergySlow, "energyS/D");
*/

/* The tree of the ROOT file contains information about the fast time, slow time & slow energy for the detectors, */
/* as well as the fast time of the RF (detectorID 4) & the fast energy and time of the POLARIS sync pulse (detectorID 5). */

/* The purpose of this code is to order the tree by ascending globalTime. It produces a root file called R63_sorted.root for event sorting. */

/*   COMPILE & RUN: 
     clear && g++ -std=c++0x -O3 sort3labr.C -o exe3 `root-config --cflags --libs` -lSpectrum 
     ./exe3 ~/Documents/PhD/exp/2022/220615/analysis/angle/water/run27/R63_aligned.root 
*/


#include <TFile.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TTreeIndex.h>

#include <time.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstdlib>
#include <inttypes.h> 
#include <stdexcept>
#include <deque>
#include <tuple>
#include <algorithm> 

const int NumDetectors = 4;
const int RFIndex = 4;
const int POLARISIndex = 5;

time_t start, end;

void timesort(TFile *outputFile, TTree* tree, std::string dirPath, std::string outputFileName)
{
    printf("\n Time sorting data...\n");

    // Initialise the tree variables
    Int_t detectorID;
    Double_t timeGlobal;
    Double_t alignedTimeFast;
    Double_t alignedEnergyFast;
    Double_t alignedTimeSlow;
    Double_t alignedEnergySlow;

    // Access the tree variables
    tree->SetBranchAddress("detectorID", &detectorID);
    tree->SetBranchAddress("globalTime", &timeGlobal);
    tree->SetBranchAddress("timeF", &alignedTimeFast);
    tree->SetBranchAddress("energyF", &alignedEnergyFast);
    tree->SetBranchAddress("timeS", &alignedTimeSlow);
    tree->SetBranchAddress("energyS", &alignedEnergySlow);

    // Initialise the sortedTree variables
    Int_t sortedDetectorID;
    Double_t sortedGlobalTime;
    Double_t sortedTimeF;
    Double_t sortedEnergyF;
    Double_t sortedTimeS;
    Double_t sortedEnergyS;

    // Create a new tree to store the sorted data
    TTree* sortedTree = new TTree("SortedData", "SortedData");
    sortedTree->Branch("detectorID", &sortedDetectorID, "detectorID/I");
    sortedTree->Branch("globalTime", &sortedGlobalTime, "globalTime/D");
    sortedTree->Branch("timeF", &sortedTimeF, "timeF/D");
    sortedTree->Branch("energyF", &sortedEnergyF, "energyF/D");
    sortedTree->Branch("timeS", &sortedTimeS, "timeS/D");
    sortedTree->Branch("energyS", &sortedEnergyS, "energyS/D");

    // Copy data from the original tree to the new tree while sorting by globalTime
    tree->BuildIndex("globalTime");
    TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();
    Long64_t nEntries = (tree->GetEntries());

    printf("\n Sorting the data by globalTime...\n");
    for (Long64_t i = 0; i < nEntries; ++i) {
        Long64_t originalEntry = index->GetIndex()[i];
        tree->GetEntry(originalEntry);
        sortedTree->Fill();
    }

    // Write the sortedTree to the output ROOT file
    sortedTree->Write();

    // Find the first non-zero entry
    double firstTime = 0.0;
    for (int entryIndex = 0; entryIndex < sortedTree->GetEntries(); ++entryIndex) {
        sortedTree->GetEntry(entryIndex);
        if ((sortedGlobalTime >0)) {
            firstTime = sortedGlobalTime;
            break;
        }
    }

    // Find the last globalTime entry
    sortedTree->GetEntry(sortedTree->GetEntries() - 1);
    double lastTime = sortedGlobalTime;

    // Calculate and print run duration
    double runDurationMinutes = (lastTime - firstTime)*1e-9/ 60.0;
    printf("First time: %.2f \n", firstTime/60.0);
    printf("Last time: %.2f \n", lastTime/60.0);
    printf("Run duration: %.2f min\n", runDurationMinutes);

    
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
    std::string outputFileName = dirPath + "R" + std::to_string(run) + "_sorted.root";
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

    // Get the offsets directory
    TDirectory *dir_aligned = (TDirectory*)file->Get("dir_aligned");
    if (!dir_aligned) {
        std::cout << "Error: Directory not found" << std::endl;
        return 1;
    }
    // print which directory is being accessed
    std::cout << "Directory: " << dir_aligned->GetName() << std::endl;

    // Get the list of trees in the input file
    // TList* treeList = dir_aligneds->GetListOfKeys();
    // treeList->Print();

    TTree* tree = static_cast<TTree*>(dir_aligned->Get("AlignedData"));
    if (!tree) {
        std::cout << "Error: Tree not found" << std::endl;
        return 1;
    }

    // Open the output ROOT file
    TFile *outputFile = new TFile(outputFileName.c_str(), "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cout << "Error: Output file not open" << std::endl;
        return 1;
    }

    // ________________________________________________________________________________________________________________________
    // Call function to time sort the data
    timesort(outputFile, tree, dirPath, outputFileName);    
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
