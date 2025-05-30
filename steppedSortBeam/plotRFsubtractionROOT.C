

#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <iostream>
#include <string>
// Open the ROOT files

void plotRFsubtractionROOT() 
{
    const char* dir1 = "/Users/shanyn/Documents/PhD/Codes/SortingCode/steppedSortBeam/";

    // LaBr3Ce spectrum 
    std::string file1 = std::string(dir1) + "insync_outsync_prompt0.root";
    const char* canvasName1 = "c1"; // Replace with the actual canvas name in the file1
    const char* histName1 = "In Sync L 0"; // Efficiency
    const char* histName2 = "Out Sync L 0"; // Efficiency

    TFile *f1 = TFile::Open(file1.c_str());

    if (!f1 ) {
        std::cerr << "Failed to open files" << std::endl;
        return;
    }

    // Retrieve the canvases
    TCanvas *c1 = (TCanvas*)f1->Get(canvasName1);

    if (!c1) {
        std::cerr << "Failed to get canvases" << std::endl;
        return;
    }

    // Retrieve the histograms from the canvases
    TH1D *h1 = (TH1D*)c1->GetPrimitive(histName1); // Histogram from first canvas
    TH1D *h2 = (TH1D*)c1->GetPrimitive(histName2); // Histogram from second canvas

    if (!h1 || !h2) {
        std::cerr << "Failed to get histograms" << std::endl;
        return;
    }

    if (h1->GetNbinsX() != h2->GetNbinsX()) {
    std::cerr << "Histograms have different binning." << std::endl;
    return;
    }

    if (h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() || 
    h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax()) {
    std::cerr << "Histograms have different ranges." << std::endl;
    return;
    }

    for (int i = 1; i <= h1->GetNbinsX(); i++) 
   {
        double binContent = h1->GetBinContent(i);
        double newBinContent = binContent*2.5;
        h1->SetBinContent(i, newBinContent);
    } 

    /* const char* dir1 = "/Users/shanyn/Documents/PhD/Codes/SortingCode/steppedSortBeam/";

    // LaBr3Ce spectrum 
    std::string file1 = std::string(dir1) + "InSync.root";
    std::string file2 = std::string(dir1) + "OutSync.root";
    const char* canvasName1 = "Canvas_1"; // Replace with the actual canvas name in the file1
    const char* canvasName2 = "Canvas_1"; // Replace with the actual canvas name in the file1
    const char* histName1 = "Slow_Energy_L00"; // Efficiency
    const char* histName2 = "Slow_Energy_L00"; // Efficiency 

    TFile *f1 = TFile::Open(file1.c_str());
    TFile *f2 = TFile::Open(file2.c_str());

    if (!f1 ) {
        std::cerr << "Failed to open files" << std::endl;
        return;
    }

    // Retrieve the canvases
    TCanvas *c1 = (TCanvas*)f1->Get(canvasName1);
    //TCanvas *c3 = (TCanvas*)f2->Get(canvasName2);

    if (!c1) {
        std::cerr << "Failed to get canvases" << std::endl;
        return;
    }

    // Retrieve the histograms from the canvases
    TH1D *h1 = (TH1D*)c1->GetPrimitive(histName1); // Histogram from first canvas
     TH1D *h2 = (TH1D*)c1->GetPrimitive(histName2); // Histogram from first canvas
    //TH1D *h2 = (TH1D*)c3->GetPrimitive(histName2); // Histogram from second canvas

    if (!h1 || !h2) {
        std::cerr << "Failed to get histograms" << std::endl;
        return;
    }

    if (h1->GetNbinsX() != h2->GetNbinsX()) {
    std::cerr << "Histograms have different binning." << std::endl;
    return;
    }

    if (h1->GetXaxis()->GetXmin() != h2->GetXaxis()->GetXmin() || 
    h1->GetXaxis()->GetXmax() != h2->GetXaxis()->GetXmax()) {
    std::cerr << "Histograms have different ranges." << std::endl;
    return;
    }

     for (int i = 500; i <=560; i++) 
   {
        double binContent = h2->GetBinContent(i);
        double newBinContent = binContent/3;
        h2->SetBinContent(i, newBinContent);
    } 

    for (int i = 1; i <= h2->GetNbinsX(); i++) 
   {
        double binContent = h2->GetBinContent(i);
        double newBinContent = binContent;
        h2->SetBinContent(i, newBinContent);
    }  */


    //create a new histogram called prompt that is the h2-h1
    TH1D *prompt = (TH1D*)h2->Clone("prompt");
    prompt->Add(h1, -1);

/*     for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double binCenter = h1->GetBinCenter(i);
        double newBinCenter = binCenter + gRandom->Gaus(-0.1,0.1);
        h1->SetBinContent(newBinCenter, h1->GetBinContent(i));

        binCenter = h2->GetBinCenter(i);
        newBinCenter = binCenter + gRandom->Gaus(-0.1,0.1);
        h2->SetBinContent(newBinCenter, h2->GetBinContent(i));

        binCenter = prompt->GetBinCenter(i);
        newBinCenter = binCenter + gRandom->Gaus(-0.1,0.1);
        prompt->SetBinContent(newBinCenter, prompt->GetBinContent(i));
    }   */


    // plot them on a new canvas
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    prompt->SetLineColor(kBlack);
    h2->GetXaxis()->SetTitle("Energy (keV)");
    h2->GetYaxis()->SetTitle("Counts (1/keV)"); 
    h2->SetTitle("");
    h2->GetXaxis()->SetTitleSize(0.06);
    h2->GetYaxis()->SetTitleSize(0.06);
    h2->GetXaxis()->SetTitleOffset(0.8);
    h2->GetYaxis()->SetTitleOffset(0.6);
    h2->GetXaxis()->SetTitleFont(22);
    h2->GetYaxis()->SetTitleFont(22);
    h2->GetXaxis()->SetLabelFont(132);
    h2->GetYaxis()->SetLabelFont(132);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetTitleSize(0.06);
    h2->GetYaxis()->SetTitleSize(0.06);
    h2->GetXaxis()->SetRangeUser(0, 8000);
    h2->SetStats(0);
    h1->SetStats(0);
    prompt->SetStats(0);
    h2->Draw();
    h1->Draw("same");
    prompt->Draw("same");

    TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->AddEntry(h1, "In Sync", "l");
    legend->AddEntry(h2, "Out Sync", "l");
    legend->AddEntry(prompt, "Prompt", "l");
    legend->Draw();

    c2->SaveAs("RFSUBRACTION_Run27_water.root");

    std::ofstream file;
    file.open("insync_R27.txt");
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        file << h1->GetBinCenter(i) << " " << h1->GetBinContent(i) << std::endl;
    }
    file.close();

    file.open("outsync_R27.txt");
    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        file << h2->GetBinCenter(i) << " " << h2->GetBinContent(i) << std::endl;
    }
    file.close();

    file.open("prompt_R27.txt");
    for (int i = 1; i <= prompt->GetNbinsX(); i++) {
        file << prompt->GetBinCenter(i) << " " << prompt->GetBinContent(i) << std::endl;
    }
    file.close();
}

