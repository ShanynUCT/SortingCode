#include <cmath>
#include "TF1.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"

// Lorentzian function
double lorentzian(double *x, double *par) {
    // par[0]: Amplitude
    // par[1]: Mean
    // par[2]: Gamma (Lorentzian width)

    double xx = x[0] - par[1];
    return par[0] / (1 + (xx * xx) / (par[2] * par[2]));
}

void fit_lorentzian() {
    // Access the histogram by name
    TH1 *h = (TH1*)gDirectory->Get("timediffcoinc");
    if (!h) {
        std::cerr << "Histogram 'timediffcoinc' not found!" << std::endl;
        return;
    }

    // Create a TF1 with the Lorentzian function
    TF1 *lorentz = new TF1("lorentz", lorentzian, -0.5, 0.5, 3);

    // Set initial parameters
    lorentz->SetParameters(1000.79, -0.005, 0.09);
    lorentz->SetParNames("Amplitude", "Mean", "Gamma");

    // Use Minuit2 with Migrad
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");

    // Fit the histogram
    h->Fit("lorentz", "R");

    // Retrieve the fit parameters
    double amplitude = lorentz->GetParameter(0);
    double mean = lorentz->GetParameter(1);
    double gamma = lorentz->GetParameter(2);
    double chi2 = lorentz->GetChisquare();
    double ndf = lorentz->GetNDF();

    // Print fit parameters and Chi-squared/ndf to terminal
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Gamma: " << gamma << std::endl;
    std::cout << "Chi2/ndf: " << chi2 << "/" << ndf << std::endl;

    // Draw the histogram and the fit
    TCanvas *c = new TCanvas("c", "Lorentzian Fit", 800, 600);
    h->Draw();
    lorentz->Draw("same");

    // Add text with Chi2/ndf on the plot
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.65, 0.65, Form("#chi^{2}/ndf = %.2f / %d", chi2, int(ndf)));
    latex.DrawLatex(0.65, 0.60, Form("#Gamma = %.2f", gamma));
}
