#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TText.h"

// Function to calculate the Voigt profile
double voigt_profile(double x, double sigma, double gamma) {
    // Lorentzian part
    double lorentz = gamma / (M_PI * (x * x + gamma * gamma));
    
    // Gaussian part
    double gauss = (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * (x / sigma) * (x / sigma));
    
    // Voigt profile
    return gauss * lorentz;
}

// Wrapper function for ROOT fitting
double voigt_function(double *x, double *par) {
    // par[0]: Amplitude
    // par[1]: Mean
    // par[2]: Sigma (Gaussian width)
    // par[3]: Gamma (Lorentzian width)
    
    double xx = x[0] - par[1];
    return par[0] * voigt_profile(xx, par[2], par[3]);
}

void fit_voigt() {
    // Access the histogram by name
    TH1 *h = (TH1*)gDirectory->Get("timediffcoinc");
    if (!h) {
        std::cerr << "Histogram 'timediffcoinc' not found!" << std::endl;
        return;
    }

    // Create a TF1 with the Voigt function
    TF1 *voigt = new TF1("voigt", voigt_function, -2, 2, 4);

    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

    voigt -> SetParameters(359,0.0,0.09,0.08);
    // force parameters 
    voigt->SetParLimits(0, 0, 1000);
    voigt->SetParNames("Amplitude", "Mean", "Sigma", "Gamma");

    // Fit the histogram
    h->Fit("voigt");

    // Retrieve the fit parameters
    double amplitude = voigt->GetParameter(0);
    double mean = voigt->GetParameter(1);
    double sigma = voigt->GetParameter(2);
    double gamma = voigt->GetParameter(3);
    double chi2 = voigt->GetChisquare();
    double ndf = voigt->GetNDF();

    // Print fit parameters and Chi-squared/ndf to terminal
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Sigma: " << sigma << std::endl;
    std::cout << "Gamma: " << gamma << std::endl;
    std::cout << "Chi2/ndf: " << chi2 << "/" << ndf << std::endl;

    // Draw the histogram and the fit
    TCanvas *c = new TCanvas("c", "Voigt Fit", 800, 600);
    h->Draw();
    voigt->Draw("same");

    // Add text with Chi2/ndf and Sigma on the plot
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.65, 0.65, Form("#chi^{2}/ndf = %.2f / %d", chi2, int(ndf)));
    latex.DrawLatex(0.65, 0.60, Form("#sigma = %.4f", sigma));
}
