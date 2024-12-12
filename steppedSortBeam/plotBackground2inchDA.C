

#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TLegend.h>
#include <iostream>
#include <string>
// Open the ROOT files

void plotBackground2inchDA() 
{
    const char* dir1 = "/Users/shanyn/Documents/PhD/PANGoLINS/Measurements/DetectorAssemblies/LaBr3/sn230824-05/13082024/backgroundCompare2inchDA/";

    /* // LaBr3Ce spectrum 
    std::string file1 = std::string(dir1) + "EnergyCalibDA_R35.root";
    std::string file2 = std::string(dir1) + "LaBr3Energy2inch_run6.root";
    const char* canvasName1 = "Canvas_1"; // Replace with the actual canvas name in the file1
    const char* canvasName2 = "Canvas_1_n2"; // Replace with the actual canvas name in the file1
    const char* histName1 = "Slow_Energy_Calib_L 0"; // Efficiency
    const char* histName2 = "Slow_Energy_Calib_L 0"; // Efficiency

    TFile *f1 = TFile::Open(file1.c_str());
    TFile *f2 = TFile::Open(file2.c_str());

    if (!f1 ) {
        std::cerr << "Failed to open files" << std::endl;
        return;
    }

    // Retrieve the canvases
    TCanvas *c1 = (TCanvas*)f1->Get(canvasName1);
    TCanvas *c4 = (TCanvas*)f2->Get(canvasName2);

    if (!c1 || !c4) {
        std::cerr << "Failed to get canvases" << std::endl;
        return;
    }

    // Retrieve the histograms from the canvases
    TH1D *h1 = (TH1D*)c1->GetPrimitive(histName1); // Histogram from first canvas
    TH1D *h2 = (TH1D*)c4->GetPrimitive(histName2); // Histogram from second canvas

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

    // write out the histograms as text files
    std::ofstream file;
    file.open("EnergyCalibDA_R35.txt");
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        file << h1->GetBinCenter(i) << " " << h1->GetBinContent(i) << std::endl;
    }
    file.close();

    file.open("LaBr3Energy2inch_run6.txt");
    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        file << h2->GetBinCenter(i) << " " << h2->GetBinContent(i) << std::endl;
    }
    file.close();
 */
   /*  for (int i = 1; i <= h1->GetNbinsX(); i++) 
   {
        double binContent = h1->GetBinContent(i);
        double newBinContent = binContent*2.5;
        h1->SetBinContent(i, newBinContent);
    }  */

   // read in the text files and plot them
    TH1D *h1 = new TH1D("h1", "h1", 1600, 0, 1600);
    TH1D *h2 = new TH1D("h2", "h2", 1600, 0, 1600);

    std::ifstream file;
    file.open("EnergyCalibDA_R35.txt");
    double energy, counts;
    while (file >> energy >> counts) {
        h1->Fill(energy, counts);
    }
    file.close();

    file.open("LaBr3Energy2inch_run6.txt");
    while (file >> energy >> counts) {
        h2->Fill(energy, counts);
    }
    file.close();


    // plot them on a new canvas
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h2->GetXaxis()->SetTitle("Energy (keV)");
    h2->GetYaxis()->SetTitle("Counts (1/keV)"); 
    h2->SetTitle("");
    h2->GetXaxis()->SetTitleSize(0.06);
    h2->GetYaxis()->SetTitleSize(0.06);
    h2->GetXaxis()->SetTitleOffset(0.8);
    h2->GetYaxis()->SetTitleOffset(0.8);
    h2->GetXaxis()->SetTitleFont(22);
    h2->GetYaxis()->SetTitleFont(22);
    h2->GetXaxis()->SetLabelFont(132);
    h2->GetYaxis()->SetLabelFont(132);
    h2->GetXaxis()->SetLabelSize(0.05);
    h2->GetYaxis()->SetLabelSize(0.05);
    h2->GetXaxis()->SetTitleSize(0.06);
    h2->GetYaxis()->SetTitleSize(0.06);
    h2->GetXaxis()->SetRangeUser(0, 1600);
    h2->SetStats(0);
    h1->SetStats(0);
    h2->Draw();
    h1->Draw("same");

    TLegend *legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->AddEntry(h1, "Detector Assembly", "l");
    legend->AddEntry(h2, " 2''#\times2'' Detector ", "l");
    legend->Draw();

    c2->SaveAs("CompareInternalRadoiActivity.root");
}

