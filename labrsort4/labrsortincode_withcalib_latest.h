#include<iostream> //INPUT OUTPUT ON SCRREN
#include<string>  //DYNAMIC ARRA#include<iostream> //INPUT OUTPUT ON SCRREN
#include<vector>  //DYNAMIC ARRAYS
#include<string>  //DYNAMIC ARRAYS OF CHAR
#include<dirent.h>
#include<fstream> //INPUT OUTPUT FILES
#include<cmath>

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
using namespace std;


TF1* Sbanano(TH2D* LaBr_Vs_TOF, TH2D* TOF_Vs_LaBr);
//**************************************
// GET FUNCTION FOR TIME WALK CORRECTION
//**************************************
TF1* Sbanano(TH2D* h2_qdc_tdc, TH2D* h2_tdc_qdc)
{
	cout <<"\n*****************\n SBANANO \n*******************\n";
	h2_qdc_tdc->GetXaxis()->SetRangeUser(5900, 6620);
	int yBins = h2_qdc_tdc->GetYaxis()->GetNbins();
	int xBins = h2_qdc_tdc->GetXaxis()->GetNbins();

	int slicewidth = 100, bin =1, binmax = 1;

	vector<double> slice_charge, slice_charge_err, slice_time, slice_time_err;

	TCanvas* c = new TCanvas();
	c->Divide(5,5);
	int i=3;
	c->cd(1);
	h2_qdc_tdc->Draw();
	c->cd(2);
	h2_tdc_qdc->Draw();


	bin = 0; //bin start
	yBins = 8000; //bin stop
	while(bin < yBins)
	{
		binmax = bin+slicewidth; //bin stop

		if(binmax > yBins) 
			break;

		TH1D* ProjX = h2_qdc_tdc->ProjectionX(Form("h%d",bin), bin, binmax); //PROJECTION ON X AXIS where FORM is the name of the histogram and the x axis is the time
		TF1* fitSlice = new TF1("fitSlice", "gaus(0)", 60, 120); //
		fitSlice ->SetParameter(1,ProjX->GetXaxis()->GetBinCenter(ProjX->GetMaximumBin()));
		fitSlice ->SetParameter(2,160);

		c->cd(i);
		ProjX -> Fit("fitSlice" , "Q");
		ProjX->Draw();
		fitSlice->Draw("same");

		slice_charge.push_back((bin+binmax)/2.); //The energy of the bin is stored in the vector slice_charge.
		slice_charge_err.push_back(slicewidth/2.); //The error on the energy of the bin is stored in the vector slice_charge_err.
		slice_time.push_back(fitSlice ->GetParameter(1)); //The mean of the Gaussian fit is stored in the vector slice_time.
		slice_time_err.push_back(fitSlice ->GetParameter(2)); //The error on the mean of the Gaussian fit is stored in the vector slice_time_err.

		i++; // counter
		bin = binmax; //bin start

	}
	c->Write("Sbanano_Slices");

	TCanvas* c1 = new TCanvas();
	TGraphErrors* gr = new TGraphErrors(slice_charge.size(),&slice_charge[0], &slice_time[0], &slice_charge_err[0], &slice_time_err[0]); //The vectors slice_charge, slice_charge_err, slice_time and slice_time_err are used to create a TGraphErrors object.
	gr -> GetYaxis()->SetTitle("Mean TOF");
	gr -> GetXaxis()->SetTitle("Mean QDC");
	TF1* fit = new TF1("fit", "expo(0)+expo(2)+expo(4)+expo(6)", 0, 8000); //The fit function is a sum of four exponentials. It is used to fit the TGraphErrors object to obtain the time walk correction function.
	//fit -> SetParameters(6.83841e+00, -1.26326e-06, 7.43936e+00, -7.26286e-03, 7.40885e+00, -7.30679e-03, 3.56114e+00, -1.49377e-03); // The initial values of the parameters of the fit function are set.
	//fit -> SetParameters(6.84204e+00, -1.36413e-06, 7.48302e+00, -6.38097e-03, 3.38462e+00, -1.40599e-03, -2.79675e-01, -6.42847e-05); // The initial values of the parameters of the fit function are set.
	gr->Fit("fit", ""); //The TGraphErrors object is fitted with the fit function.
	//fit -> SetParameters(3.98202, -0.00254556, 3.1248, -3.12981e-5, 2.56719, 4.5255e-5, 1.81332, -0.0001458);
	gr->Draw("AL*");
	fit->Draw("same");

	c1->Write("FitSbanano");

	return fit;

}