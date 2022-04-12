#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <math.h>
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include <vector>
#include <map>

// ==========================
// Available methods
// --------------------------
void     PolarisationHeMinuit2D    ();
// ==========================
// Global data
// --------------------------
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsY;
std::vector< Double_t > values;
std::vector< Double_t > errors;
Double_t GlobalChi = 0;
//______________________________________________
void BeautifyPad(){
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
}
//______________________________________________
void BeautifyHisto(TH1* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");
}
//______________________________________________
void BeautifyHisto2D(TH2* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void PolarisationHeMinuit2D(){
  TFile* f[24];
  TH1F*  h[24];
  TH1F*  generated[24];
  Double_t chi2_all[24];
  Double_t chi2_total = 0.;
  for (Int_t i = 5; i < 19; i++) {
    f[i]         = new TFile(Form("SignalExtractionCoarse/Unfolding/UnfoldHeV2_%d.root", i));
    h[i]         = (TH1F*) f[i]->Get(Form("histo_%d", i)); // actual
  }



  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();









  // =========================================
  // Filling the 2D histogram
  // -----------------------------------------
  TH2F* histo2 = new TH2F("histo2", "histo2", 6, 0., 2.*TMath::Pi(), 24, -1., 1.);
  Double_t PhiCenters[6];
  Double_t Spacing = TMath::Pi()/6.;
  for (size_t i = 0; i < 6; i++) {
    PhiCenters[i] = (2.*((Double_t) i) + 1.) *Spacing;
    cout << "Phi centres = " << PhiCenters[i] << endl;

  }
  Double_t CosThetaCenters[24];
  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
    cout << "CosTheta centres = " << CosThetaCenters[i] << endl;

  }
  /// reset data structure
  coordsX = std::vector<Double_t>();
  coordsY = std::vector<Double_t>();
  values  = std::vector<Double_t>();
  errors  = std::vector<Double_t>();
  Int_t nBins = 1;
  for (Int_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
    nBins = h[iCosThetaBins]->GetNbinsX();
    for (Int_t iPhiBins = 0; iPhiBins < nBins; iPhiBins++) {
      coordsX.push_back( TMath::Pi()/((Double_t) nBins) + iPhiBins*2.*TMath::Pi()/((Double_t) nBins) );
      coordsY.push_back( CosThetaCenters[iCosThetaBins] );
      values.push_back(  h[iCosThetaBins]->GetBinContent(iPhiBins+1)        );
      errors.push_back(  h[iCosThetaBins]->GetBinError(iPhiBins+1)        );

    }
  }
  Double_t NormalisationSum      = 0.;
  Double_t NormalisationSumError = 0.;
  for (Int_t i = 0; i < coordsX.size(); i++) {
    NormalisationSum      += values[i];
    NormalisationSumError += errors[i];
  }


  /// fill data structure
  Double_t Mn = 6.;
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
    for (size_t iPhiBins = 0; iPhiBins < 6; iPhiBins++) {

      Int_t binx = -1;
      binx = iPhiBins+1;
      Int_t biny2 = histo2->GetYaxis()->FindBin(CosThetaCenters[iCosThetaBins]);
      Int_t binx2 = histo2->GetXaxis()->FindBin(PhiCenters[iPhiBins]);

      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx)*Mn );
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx)*Mn );


    }
  }










  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  BeautifyHisto2D(histo2);
  histo2->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );









    const Int_t NRGBs = 6;
    const Int_t NCont = 999;

    Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };


    // TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    // gStyle->SetNumberContours(NCont);

    gStyle->SetOptStat(0);

    //here the actually interesting code starts
    const Double_t min = histo2->GetMaximum()*0.3;
    const Double_t max = histo2->GetMaximum()*1.3;

    const Int_t nLevels = 999;
    Double_t levels[nLevels];


    for(int i = 1; i < nLevels; i++) {
      levels[i] = min + (max - min) / (nLevels - 1) * (i);
    }
    levels[0] = 0.01;

    // histo2->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    histo2->DrawClone("col");// draw "axes", "contents", "statistics box"
    histo2->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
    histo2->Draw("colz same"); // draw the "color palette"



  TLatex* latex = new TLatex();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.12,0.94,"This thesis, Data, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.8,0.94,"Helicity");
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.6,0.79,"Raw yields" );
  gPad->SaveAs("SignalExtractionCoarse/Fitting/2Dmaps-uncorrected-data.pdf", "recreate");


  new TCanvas;
  histo2->Draw("surf same");
  latex->DrawLatex(0.12,0.94,"This thesis, Data, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.8,0.94,"Helicity");
  latex->SetTextSize(0.04);
  latex->DrawLatex(0.6,0.79,"Raw yields" );
  gPad->SaveAs("SignalExtractionCoarse/Fitting/3D-uncorrected-data.pdf", "recreate");



}
