#include<iostream> //INPUT OUTPUT ON SCRREN
#include<string>  //DYNAMIC ARRA#include<iostream> //INPUT OUTPUT ON SCRREN
#include<vector>  //DYNAMIC ARRAYS
#include<string>  //DYNAMIC ARRAYS OF CHAR
#include<dirent.h>
#include<fstream> //INPUT OUTPUT FILES
#include<cmath>
#include <stdio.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <inttypes.h>
#include <initializer_list>
#include <iomanip>
#include <sstream>

// ROOT libraries
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
#include "TROOT.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TH3D.h"
#include "TObject.h"
#include "TString.h"
#include "TRandom.h"
#include "TH2F.h"

using namespace std;

// ********************
struct Calibration
{
	vector<double> intercalib01 = {0., 1.}; //y=a+b*x
	vector<double> intercalib02 = {0., 1.}; //y=a+b*x
	vector<double> enecalib = {0., 1., 0.}; //it can be a quadratic, now it's linear

	double QadraticCalib(double ch, double a0=0., double a1=1., double a2=0.)
	{
		return a0 + a1*ch + a2*ch*ch;
	}
};


//********************
// FUNCTION DEFINITION
//********************
vector<double> eneCalib(TH1* hist);
void ChannelToEneHisto(TH1D* h_ch, TH1D* h_ene, Calibration calib);
void Transform(TH1* hist,Double_t a0,Double_t a1,Double_t scale,Double_t dscale);


//*** Transform histograms, used to calibrate the data (I apply a LINEAR calibration)*******************************/
void Transform(TH1* hist,Double_t a0=0,	Double_t a1=1,Double_t scale=0,Double_t dscale=0) {
	TAxis *xaxis;
	Int_t nbins;
	Double_t xmin,xmax;

	nbins = hist->GetNbinsX();
	if(a0 != 0.0 || a1 != 1.0) {
		xaxis = hist->GetXaxis();
		xmin = a0 + a1*xaxis->GetXmin();
		xmax = a0 + a1*xaxis->GetXmax();
		xaxis->Set(nbins,xmin,xmax);
	}
	cout << "xmin = " << xmin << " xmax = " << xmax << endl;
	/*if(scale != 0.0) {
	  Int_t sw2 = hist->GetSumw2N();
	  if(sw2 == 0) hist->Sumw2();
	  if(dscale != 0.0) {
	  Double_t ds = dscale/scale;
	  ds = ds*ds;
	  for(Int_t i=1; i<=hist->GetNbinsX(); i++) {
	  Double_t err,val;
	//err = hist->GetBinError(i);
	val = hist->GetBinContent(i);
	//err = sqrt(err*err+val*val*ds);
	//hist->SetBinError(i,err);
	}
	}
	hist->Scale(scale);
	}*/
}

//*****************************
// GET ENERGY CALIBRATION CURVE
//***************************** 
vector<double> eneCalib(TH1* hist)
{
	cout << "\n******************\n ENERGY CALIBRATION \n*******************\n";
	TCanvas* can = new TCanvas();
	can -> Divide(2,1);

	can -> cd(1);
	hist -> Draw();
	int maxpeaks = 20;
	TSpectrum* s = new TSpectrum();
	int nfound = s->Search(hist, 2, "nobackground new", 0.005);
	//cout << "Number of peaks found: " << nfound << endl;
	
	//fit gauss+bkg for energy calib
	vector<double> xpeaks6fit(6, 0.);
	vector<double> ene4fit = {121.8, 224.7,344.3,778.9,964.1, 1112}; //peaks in keV // 160.56, 327.,460.,549.,589.,1036.,1281., 1725, 1871
	
	double *xpeaks;
	xpeaks = s->GetPositionX();
	for (int p=0; p<nfound; p++) 
	{
		double xp = xpeaks[p];
		int bin = hist->GetXaxis()->FindBin(xp);
		double yp = hist->GetBinContent(bin);
        
		double xStart = 150, xEnd = 170;
		if(xp > xStart && xp < xEnd) //found by graphs
		{
			cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 121.8 keV" <<endl;
			TF1* fit_peak1 = new TF1("fit_peak1", "gaus(0)+pol1(3)", xp-50, xp+50);
			fit_peak1 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak1", "Q");
			fit_peak1 -> Draw("same");
			xpeaks6fit[0] = fit_peak1->GetParameter(1);

		}

		xStart = 317; xEnd = 337;
		if(xp > xStart && xp < xEnd)
		{
			cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 224.7 keV" <<endl;
			TF1* fit_peak2 = new TF1("fit_peak2", "gaus(0)+pol1(3)", xp-30, xp+30);
			fit_peak2 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak2", "Q");
			fit_peak2 -> Draw("same");
			xpeaks6fit[1] = fit_peak2->GetParameter(1);
		}

		xStart = 450; xEnd = 470;
		if(xp > xStart && xp < xEnd)
		{
			cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 344.3 keV" <<endl;
			TF1* fit_peak3 = new TF1("fit_peak3", "gaus(0)+pol1(3)", xp-30, xp+30);
			fit_peak3 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak3", "Q");
			fit_peak3 -> Draw("same");
			xpeaks6fit[2] = fit_peak3->GetParameter(1);
		}

		xStart = 539; 	xEnd = 559;
		if(xp > xStart && xp < xEnd)
		{
			cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 778.9 keV" <<endl;
			TF1* fit_peak4 = new TF1("fit_peak4", "gaus(0)+pol1(3)", xp-30, xp+30);
			fit_peak4 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak4", "Q");
			fit_peak4 -> Draw("same");
			xpeaks6fit[3] = fit_peak4->GetParameter(1);
		}

        xStart = 1715; 	xEnd = 1735;
		if(xp > xStart && xp < xEnd)
		{
			//cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 964.1 keV" <<endl;
			TF1* fit_peak5 = new TF1("fit_peak5", "gaus(0)+pol1(3)", xp-30, xp+30);
			fit_peak5 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak5", "Q");
			fit_peak5 -> Draw("same");
			xpeaks6fit[4] = fit_peak5->GetParameter(1);
		}

        xStart = 1861; 	xEnd = 1881;
		if(xp > xStart && xp < xEnd)
		{
			cout << "PEAK FOUD AT CHANNEL: " << xp <<" --> 1112 keV" <<endl;
			TF1* fit_peak6 = new TF1("fit_peak6", "gaus(0)+pol1(3)", xp-30, xp+30);
			fit_peak6 -> SetParameters(500, xp, 100, 500, -1);
			hist -> GetXaxis()->SetRangeUser(xStart, xEnd);
			hist ->Fit("fit_peak6", "Q");
			fit_peak6 -> Draw("same");
			xpeaks6fit[5] = fit_peak6->GetParameter(1);
		}

	}

	can -> cd(2);
	sort(xpeaks6fit.begin(), xpeaks6fit.end());
	TF1* ecalib = new TF1("ecalib", "pol1", 800, 2500);
	TGraph* gg = new TGraph(ene4fit.size(), &xpeaks6fit[0], &ene4fit[0]);
	gg->GetXaxis()->SetTitle("ADC channel");
	gg->GetYaxis()->SetTitle("Energy [keV]");
	gg->Draw("A*");
	gg->Fit("ecalib");

	vector<double> result = {ecalib->GetParameter(0), ecalib->GetParameter(1)};


	can -> Write("eneCalib");
	can -> Close();
	//cout << "\n******************\n";

	return result;

}

//*******************************************************
// MAKE HISTOGRAM WITH DEFINED BINWIDTH WITHOUT ARTIFACTS
//*******************************************************
void ChannelToEneHisto(TH1D* h_ch, TH1D* h_ene, Calibration c)
{
	TRandom3* randy = new TRandom3();

	int N = h_ch->GetXaxis()->GetNbins();
	double bw = h_ch->GetBinWidth(1);
	for(int i=1; i<=N; i++)
	{
		int n = h_ch->GetBinContent(i);
		double ch = h_ch -> GetBinCenter(i);

		for(int j=1; j<=n; j++)
		{
			double rndUnif = randy->Rndm()*bw-(bw/2.); //generate uniform rnd number iside binwidth
			double rndchannel = ch+rndUnif; //sum to bincenter
			double enechannel = c.QadraticCalib(rndchannel, c.enecalib[0], c.enecalib[1], 0.); //get energy calibration

			h_ene ->Fill(enechannel); //fill energy calibrated histo
		}
	}


}
