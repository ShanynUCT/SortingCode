#include<iostream> //INPUT OUTPUT ON SCRREN
#include<string>  //DYNAMIC ARRAY
#include<vector>  //DYNAMIC ARRAYS
#include<string>  //DYNAMIC ARRAYS OF CHAR
#include<dirent.h> //DIRECTORY
#include<fstream> //INPUT OUTPUT FILES
#include<cmath> //MATH FUNCTIONS

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TPad.h"
#include "TChain.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TRandom3.h"
using namespace std; //USE STANDARD NAMESPACE
#include "labrsortincode_withcalib_latest.h"

int main (int argc, char* argv[])
{
    TFile *g = new TFile(argv[1], "UPDATE");
    TTree *LaBr0Data = (TTree*)g->Get("LaBr0Data");

 /******************************************************************************/
 /*******************             RF CORRECTION             ********************/
 /******************************************************************************/

    TF1* sbanano = 0;
    TCanvas *c = new TCanvas();
    int RF_mask = 1, Out_RF_mask = 0;
    double TOF_Sbanano = 0, EnergyNoRF, EnergyRF, slowLaBr0calib1stFullRange, slowtimeLaBr0, RFsignaltime;
    Long64_t mintime, maxtime, bins, entriesn, minRFtime, maxRFtime, binsRF;
    
    entriesn = LaBr0Data->GetEntries();
    entriesn = entriesn;

    LaBr0Data->SetBranchAddress("EnergySlow_LaBr_Det0_1stOrder_Calib_FullRange", &slowLaBr0calib1stFullRange);
    LaBr0Data->SetBranchAddress("TimingSlow_LaBr_Det0", &slowtimeLaBr0); // time is actually in 0.1 ns units
    LaBr0Data->SetBranchAddress("TimingSlow_RF", &RFsignaltime);

    TBranch* Energy_L0_Minus_RF = LaBr0Data->Branch("PromptEnergy", &EnergyNoRF, "EnergyNoRF/I");
    TBranch* RF_Energy = LaBr0Data->Branch("OutOfSyncEnergyRF", &EnergyRF, "EnergyRF/I");

    // the mintime is the first slowtimeLaBr0 value > 1e13 for which there is a non zero slowLaBr0calib1stFullRange value 
    for (Long64_t i = 0; i < entriesn; i++)
    {
        LaBr0Data->GetEntry(i);
        if (slowLaBr0calib1stFullRange > 1 && slowtimeLaBr0 > 1e13)
        {
            mintime = slowtimeLaBr0;
            minRFtime = RFsignaltime;
            break;
        }
    }

    for (Long64_t i = entriesn-1; i > 0; i--)
    {
        LaBr0Data->GetEntry(i);
        maxtime = slowtimeLaBr0;
        maxRFtime = RFsignaltime;
        break;
    }

    bins = (maxRFtime - mintime)*1e-4; // 1e-4 is the conversion from seconds to tens of microseconds 

    binsRF = bins; // 1e-4 is the conversion from seconds to tens of microseconds

    cout << "Min LaBr3 time: " << mintime << " s"<< " \nMax RF time: " << maxtime << " s"<< "\nTime duration (min LaBr3 - Max RF): " << ((maxtime-mintime)*1e-10)/60 << " minutes" << "\nbins: " << bins <<endl;

    TH2D *h_tdc_qdc = new TH2D("h_tdc_qdc", "htdcqdc", 8000, 0, 8000, 8000, 0, 8000); // time in us bins
    TH2D *h_qdc_tdc = new TH2D("h_qdc_tdc", "hqdctdc", 8000, 0, 10000, 8000, 0, 8000);
    TH1D *h_tdc = new TH1D("h_tdc", "htdc", 50000, 0, 50000);
    
    TH1D *RF_tdc = new TH1D("RF_tdc", "RF_tdc", 5000000, mintime , mintime + 5e7);
    TH1  *LaBr0Time = new TH1D("LaBr0Time", "LaBr0Time", 5000000, mintime, mintime + 5e7);
 
    TH1D *h_qdc0_0_Sync = new TH1D("h_qdc0_0_Sync", "In-sync gamma energy calibration", 8000, 0, 8000);
    TH1D *h_qdc0_0_outSync = new TH1D("h_qdc0_0_outSync", "out-sync gamma energy calibration", 8000, 0, 8000);

    for (Long64_t i=0; i < entriesn; i++)
    {
        LaBr0Data->GetEntry(i);
        if (RFsignaltime <= maxtime && RFsignaltime >= mintime)
        {
        RF_tdc->Fill(RFsignaltime/10); // RF signal time in ns units
        float TOF = slowtimeLaBr0 - RFsignaltime;
        //cout << "TOF: " << TOF << endl;
        // take absolute value TOF
        TOF = abs(TOF);
        h_tdc_qdc->Fill(slowLaBr0calib1stFullRange, TOF);
        h_qdc_tdc->Fill(TOF, slowLaBr0calib1stFullRange);
        h_tdc->Fill(TOF);
        LaBr0Time->Fill(slowtimeLaBr0/10); // LaBr0 time in ns units
        }
    }
    // plot LaBr0 time and RF_tdc on the same axes
    c = new TCanvas("cx", "cx", 1000, 800);
    c-> cd();
    LaBr0Time->Draw();
    LaBr0Time->SetLineColor(kRed);
    LaBr0Time->GetXaxis()->SetTitle("Time (ns)");
    RF_tdc->Draw("same");
    RF_tdc->SetLineColor(kBlue);
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(LaBr0Time,"LaBr0Time","l");
    leg->AddEntry(RF_tdc,"RF_tdc","l");
    leg->Draw();
    c->Write();
    c->SaveAs("LaBr0Time_RF_tdc.png");
    c->Close();

    RF_tdc->Write();
    h_tdc->Write();
    
    h_qdc_tdc->Write();
    h_tdc_qdc->Write();

    // find the peak of the RF_tdc histogram
    sbanano = Sbanano(h_qdc_tdc, h_tdc_qdc); //fit for sbanano function TF1, see .h file
    TCanvas* plot = new TCanvas();
    h_tdc_qdc->Draw();
    sbanano->Draw("same");
    plot->Write("Sbanano_2D_fitfunction");

    
    for (Long64_t i=0; i < entriesn; i++)
    {
        LaBr0Data->GetEntry(i);
        if (RFsignaltime <= maxtime && RFsignaltime >= mintime)
        {
            float TOF = slowtimeLaBr0 - RFsignaltime;
            TOF = abs(TOF);
            //cout << "TOF: " << TOF << endl;
            TOF_Sbanano = (double)TOF - sbanano->Eval(slowLaBr0calib1stFullRange); // Evaluate the time of flight minus the sbanano function which tells us the time of flight at which the RF signal is at its peak
            TOF_Sbanano = abs(TOF_Sbanano);
            //cout << "Evaluating TOF_Sbanano: " << TOF_Sbanano << endl;

            //I label each event if in-sync or out-of-sync with RF
            RF_mask = 0;
            Out_RF_mask =0;
        
            //first loop for time > 60 ns is the RF period
            for(double t =0; t<30; t+=6.09) 
            {
                if(abs(TOF_Sbanano-t) < 0.15) // if the time of flight is within 1.5 ns of the RF period, then it is in sync with the RF
                {
                    RF_mask = 1;
                    //cout << "abs(TOF_Sbanano-t): " << abs(TOF_Sbanano-t) << endl;
                    break;
                }

                if(abs(TOF_Sbanano-(t+3)) < 0.15) // if the time of flight is within 1.5 ns of the RF period + 5 ns, then it is out of sync with the RF
                {
                    Out_RF_mask =1;
                    //cout << "abs(TOF_Sbanano-(t+3): " << abs(TOF_Sbanano-(t+3)) << endl;
                    break;
                }
            }
        

            if(RF_mask==1)
            {
                h_qdc0_0_Sync->Fill(slowLaBr0calib1stFullRange);
                EnergyNoRF = slowLaBr0calib1stFullRange;
                LaBr0Data -> Fill();
            }
            
            if(Out_RF_mask==1)
            {
                h_qdc0_0_outSync->Fill(slowLaBr0calib1stFullRange);
                EnergyRF = slowLaBr0calib1stFullRange;
                LaBr0Data -> Fill();
            }
        }
    }

    h_qdc0_0_Sync -> Write();
    h_qdc0_0_outSync -> Write(); 
    
    LaBr0Data->Write();

    //SUBTRACTION 
    //ONLY EVENTS PG = SYNC - OutSync
    TH1D* h_qdc0_0_Prompt = (TH1D*)h_qdc0_0_Sync->Clone("h_qdc0_0_Prompt");
    h_qdc0_0_Prompt -> Add(h_qdc0_0_outSync, -1);

    //PROMPT GAMMA SPECTRA
    h_qdc0_0_Prompt -> Scale(1., "width");
    h_qdc0_0_Prompt->Write();

    delete g; // automatically deletes all "owned" histograms, too
    return 0;
}
