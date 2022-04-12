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
Double_t calcChi2      (TH1* h1, TH1* h2, int nbins);
Double_t calcChi2v2    (TH1* h1, TH1* h2, int nbins);
void     BeautifyPad   ();
void     BeautifyHisto (TH1* histogram);
void     ToyData       (Int_t Stripe = 14);
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
//_____________________________________________________________________________
Double_t calcChi2(TH1* h1, TH1* h2, Int_t nbins) {
    // Int_t nb = h1->GetNbinsX();
    Int_t nb = nbins;
    // std::cout << "chi2 Nbins:" << nb << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nb+1; i++) {
        Double_t b1  = h1->GetBinContent(i);
        Double_t Eb1 = h1->GetBinError(i);
        Double_t Eb2 = h2->GetBinError(i);
        Double_t b2  = h2->GetBinContent(i);
        // cout << " partial = " << (b1-b2)*(b1-b2)/Eb1 << endl;
        chi2 += (b1-b2)*(b1-b2)/Eb1;
        // chi2 += (b1-b1)*(b2-b1)/Eb2;
    }
    return chi2;
}
//_____________________________________________________________________________
Double_t calcChi2v2(TH1* h1, TH1* h2, Int_t nbins) {
    // Int_t nb = h1->GetNbinsX();
    Int_t nb = nbins;
    // std::cout << "chi2 Nbins:" << nb << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nb+1; i++) {
        Double_t b1  = h1->GetBinContent(i);
        Double_t Eb1 = h1->GetBinError(i);
        Double_t Eb2 = h2->GetBinError(i);
        Double_t b2  = h2->GetBinContent(i);
        // cout << " partial = " << (b1-b2)*(b1-b2)/b1 << endl;
        chi2 += (b1-b2)*(b1-b2)/b2;
        // chi2 += (b1-b1)*(b2-b1)/Eb2;
    }
    return chi2;
}
//_____________________________________________________________________________
void UnfoldModulation(){

  // ===================
  // Retrieve REAL data
  // -------------------
  TFile* fData[24];
  TH1F* hdata[24];
  TH1F* hgen[24];
  for (Int_t i = 5; i < 20; i++) {
    fData[i]  = new TFile(Form("ScanModulation/LPhi005/ModulationYields/MonteCarloYieldsHe_phimodulation_%d.root", i));
    hdata[i]  = (TH1F*)fData[i]->Get(Form("h_%d",  i));
    hgen[i]   = (TH1F*)fData[i]->Get(Form("hg_%d", i));
  }
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight or similar
  // -------------------
  TFile* f[24];
  RooUnfoldResponse*  h[24];
  TH1F* hdata_formatted[24];
  Int_t M = 6;
  for (Int_t i = 5; i < 19; i++) {
    f[i]    = new TFile(Form("ScanModulation/LPhi005/UnfoldingModulation/UnfoldedClosureHe_phimodulation_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f[i]->Get("response;1");
    hdata_formatted[i] = new TH1F(Form("hdata_formatted_%d", i), Form("hdata_formatted_%d", i),  M, 0, 2.*TMath::Pi());
    for (Int_t j = 1; j < M+1; j++) {
      hdata_formatted[i]->SetBinContent(j, hdata[i] ->GetBinContent(j));
      hdata_formatted[i]->SetBinError(  j, hdata[i] ->GetBinError(j));
    }
  }




  // ===================
  // Unfolding
  // -------------------
  RooUnfoldBayes unfold[24];
  TH1D* hunfold[24];
  TH1D* relativeuncertainties[24];
  Int_t N = 15;
  // -------------------
  // How many iterations?
  // From the analysis
  // of the correlation
  // coefficients
  // -------------------
  for (Int_t i = 5; i < 19; i++) {
    unfold[i]  = RooUnfoldBayes(h[i], hdata_formatted[i], N);
    hunfold[i] = (TH1D*) unfold[i].Hreco();
    relativeuncertainties[i] = new TH1D(Form("relativeuncertainties_%d", i), Form("relativeuncertainties_%d", i),  M, 0, 2.*TMath::Pi());
    for (Int_t j = 1; j < M+1; j++) {
      relativeuncertainties[i]->SetBinContent(j, hunfold[i]->GetBinError(j)/hunfold[i]->GetBinContent(j));
    }
  }


  // -------------------
  // Draw unfolded data
  // -------------------
  new TCanvas;
  BeautifyPad();
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  for (Int_t i = 5; i < 19; i++) {
    BeautifyHisto(hunfold[i]);
    hunfold[i]->GetXaxis()->SetTitle("#varphi");
    hunfold[i]->GetYaxis()->SetTitle(Form("Unfolded dN/d#varphi  [1/%f]", hunfold[i]->GetXaxis()->GetBinWidth(1)));
    hunfold[i]->GetYaxis()->SetRangeUser(0.,hunfold[i]->GetMaximum()*(2.));
    hunfold[i]->Draw();
    latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    latex5->DrawLatex(0.31,0.82,Form("Iterations %d", N));
    latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", i));
    gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/unfolded-data-%d.pdf", i), "recreate");

  }
  // -------------------
  // Draw relative
  // uncertainties
  // -------------------
  new TCanvas;
  BeautifyPad();
  for (Int_t i = 5; i < 19; i++) {
    BeautifyHisto(relativeuncertainties[i]);
    relativeuncertainties[i]->GetXaxis()->SetTitle("#varphi");
    // relativeuncertainties[i]->GetYaxis()->SetTitle(Form("Relative uncertainties [a.u.]", relativeuncertainties[i]->GetXaxis()->GetBinWidth(1)));
    relativeuncertainties[i]->GetYaxis()->SetTitle("Relative uncertainties [a.u.]");
    relativeuncertainties[i]->GetYaxis()->SetRangeUser(0.,relativeuncertainties[i]->GetMaximum()*(1.5));
    relativeuncertainties[i]->Draw();
    latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    latex5->DrawLatex(0.31,0.82,Form("Iterations %d", N));
    latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", i));
    gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/unfolded-data-relative-uncertainties-%d.pdf", i), "recreate");

  }
















  new TCanvas;
  BeautifyPad();
  TH1F* UnfoldGenChiPlot = new TH1F("UnfoldGenChiPlot", "UnfoldGenChiPlot", 100, -0.5, 99.5);
  BeautifyHisto(UnfoldGenChiPlot);
  Double_t ChiSquareTestRef  = 0;
  Double_t ChiSquareTestRef2 = 0;
  Int_t    Counter6          = 1;
  for (Int_t i = 5; i < 19; i++) {
    ChiSquareTestRef  += calcChi2v2(hgen[i], hunfold[i], 6);
    ChiSquareTestRef2 += calcChi2(hunfold[i],hgen[i],    6);
    for (size_t j = 1; j < 7; j++) {
      UnfoldGenChiPlot->SetBinContent(Counter6, (hunfold[i]->GetBinContent(j)-hgen[i]->GetBinContent(j))/hgen[i]->GetBinContent(j));
      Counter6++;
    }
  }
  UnfoldGenChiPlot->GetXaxis()->SetTitle("Bin number");
  UnfoldGenChiPlot->GetYaxis()->SetTitle("(Unfolded DATA - GEN)/GEN [%]");
  UnfoldGenChiPlot->GetYaxis()->SetRangeUser(-1.,1.);
  UnfoldGenChiPlot->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.82,Form("Iterations %d", N));
  latex5->DrawLatex(0.2,0.72,Form("#chi^{2}_{1} = #sum (Unfolded DATA - GEN)^{2}/GEN = %0.02f",            ChiSquareTestRef));
  latex5->DrawLatex(0.2,0.62,Form("#chi^{2}_{2} = #sum (Unfolded DATA - GEN)^{2}/#DeltaUnfolded = %0.02f", ChiSquareTestRef2));
  gPad->SaveAs("ScanModulation/LPhi005/ConvergenceModulation/unfolded-vs-gen-modulation.pdf", "recreate");






}
