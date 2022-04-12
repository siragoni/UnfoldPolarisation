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
void UnfoldClosure(Int_t Iterations = 1){

  // ===================
  // Retrieve REAL data
  // -------------------
  TFile* fData[24];
  TH1F* hdata[24];
  TH1F* hgen[24];
  for (Int_t i = 5; i < 20; i++) {
    // fData[i]  = new TFile(Form("SignalExtractionCoarse/RawYieldsHeV3_%d.root", i));
    fData[i]  = new TFile(Form("SignalExtractionCoarse/ClosureYields/MonteCarloYieldsHe_%d.root", i));
    hdata[i]  = (TH1F*)fData[i]->Get(Form("h_%d",  i));
    hgen[i]  = (TH1F*)fData[i]->Get(Form("hg_%d",  i));
  }
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight or similar
  // -------------------
  TFile* f[24];
  RooUnfoldResponse*  h[24];
  TH1F* hdata_formatted[24];
  TH1F* hdata_formatted_gen[24];
  Int_t M = 6;
  for (Int_t i = 5; i < 19; i++) {
    f[i]    = new TFile(Form("SignalExtractionCoarse/Unfolding/UnfoldHeV2_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f[i]->Get("response;1");
    hdata_formatted[i]     = new TH1F(Form("hdata_formatted_%d", i),     Form("hdata_formatted_%d", i),      M, 0, 2.*TMath::Pi());
    hdata_formatted_gen[i] = new TH1F(Form("hdata_formatted_gen_%d", i), Form("hdata_formatted_gen_%d", i),  M, 0, 2.*TMath::Pi());
    for (Int_t j = 1; j < M+1; j++) {
      hdata_formatted[i]->SetBinContent(    j, hdata[i] ->GetBinContent(j));
      hdata_formatted[i]->SetBinError(      j, hdata[i] ->GetBinError(j));
      hdata_formatted_gen[i]->SetBinContent(j, hgen[i]  ->GetBinContent(j));
      hdata_formatted_gen[i]->SetBinError(  j, hgen[i]  ->GetBinError(j));
    }
  }




  // ===================
  // Unfolding
  // -------------------
  RooUnfoldBayes unfold[24];
  TH1D* hunfold[24];
  TH1D* refolded[24];
  TH1D* relativeuncertainties[24];
  Int_t N = Iterations;
  Int_t Counter = 1;
  TH1F* refoldedrecoH     = new TH1F("refoldedrecoH","refoldedrecoH", 300, -0.5, 299.5);
  TH1F* unfoldedgenPullH  = new TH1F("unfoldedgenPullH","unfoldedgenPullH", 100, -4000, 4000);
  TH1F* unfoldedgenResH   = new TH1F("unfoldedgenResH","unfoldedgenResH", 100, -40, 40);
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
    refolded[i] = (TH1D*) h[i]->ApplyToTruth(hunfold[i]);
    for (Int_t j = 1; j < M+1; j++) {
      refoldedrecoH->SetBinContent(Counter, (refolded[i]->GetBinContent(j) - hdata_formatted[i]->GetBinContent(j))/hdata_formatted[i]->GetBinContent(j) );
      Counter++;
      unfoldedgenPullH->Fill((hunfold[i]->GetBinContent(j)-hdata_formatted_gen[i]->GetBinContent(j))/hunfold[i]->GetBinError(j));
      unfoldedgenResH ->Fill(hunfold[i]->GetBinContent(j)-hdata_formatted_gen[i]->GetBinContent(j));
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
    gPad->SaveAs(Form("SignalExtractionCoarse/ConvergenceClosure/unfolded-data-%d.pdf", i), "recreate");

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
    gPad->SaveAs(Form("SignalExtractionCoarse/ConvergenceClosure/unfolded-data-relative-uncertainties-%d.pdf", i), "recreate");

  }
  // -------------------
  // Draw refolded
  // reco ratio
  // distribution
  // -------------------
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(refoldedrecoH);
  refoldedrecoH->GetXaxis()->SetTitle("Bin number");
  // refoldedrecoH->GetYaxis()->SetTitle(Form("Relative uncertainties [a.u.]", relativeuncertainties[i]->GetXaxis()->GetBinWidth(1)));
  refoldedrecoH->GetYaxis()->SetTitle("(Refolded - REC)/REC [%]");
  refoldedrecoH->GetYaxis()->SetRangeUser(-1.e-15,1.e-15);
  refoldedrecoH->GetXaxis()->SetRangeUser(-0.5, 85.5);
  refoldedrecoH->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.78,"1 iteration");}
  else {latex5->DrawLatex(0.31,0.62,Form("%d iterations", Iterations));}
  // latex5->DrawLatex(0.2,0.80,Form("#chi^{2}_{1} = #sum (Refolded - REC)^{2}/REC = %0.02f",            ChiSquareTestRef));
  // latex5->DrawLatex(0.2,0.70,Form("#chi^{2}_{2} = #sum (Refolded - REC)^{2}/#DeltaRefolded = %0.02f", ChiSquareTestRef2));
  gPad->SaveAs(Form("SignalExtractionCoarse/ConvergenceClosure/refolded-iterations=%d.pdf", Iterations), "recreate");
















  // -------------------
  // Draw pull
  // closure
  // distribution
  // -------------------
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(unfoldedgenPullH);
  unfoldedgenPullH->GetXaxis()->SetTitle("Pulls");
  // unfoldedgenPullH->GetYaxis()->SetTitle(Form("Relative uncertainties [a.u.]", relativeuncertainties[i]->GetXaxis()->GetBinWidth(1)));
  unfoldedgenPullH->GetYaxis()->SetTitle("Counts [a.u.]");
  // unfoldedgenPullH->GetYaxis()->SetRangeUser(-1.e-15,1.e-15);
  // unfoldedgenPullH->GetXaxis()->SetRangeUser(-0.5, 85.5);
  unfoldedgenPullH->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.78,"1 iteration");}
  else {latex5->DrawLatex(0.31,0.62,Form("%d iterations", Iterations));}
  latex5->DrawLatex(0.2,0.80,"Pulls = #frac{Unfolded - GEN}{#DeltaUnfolded}");
  // latex5->DrawLatex(0.2,0.80,Form("#chi^{2}_{1} = #sum (Refolded - REC)^{2}/REC = %0.02f",            ChiSquareTestRef));
  // latex5->DrawLatex(0.2,0.70,Form("#chi^{2}_{2} = #sum (Refolded - REC)^{2}/#DeltaRefolded = %0.02f", ChiSquareTestRef2));
  gPad->SaveAs(Form("SignalExtractionCoarse/ConvergenceClosure/pulls-iterations=%d.pdf", Iterations), "recreate");





  // -------------------
  // Draw residual
  // closure
  // distribution
  // -------------------
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(unfoldedgenResH);
  unfoldedgenResH->GetXaxis()->SetTitle("Residuals");
  // unfoldedgenResH->GetYaxis()->SetTitle(Form("Relative uncertainties [a.u.]", relativeuncertainties[i]->GetXaxis()->GetBinWidth(1)));
  unfoldedgenResH->GetYaxis()->SetTitle("Counts [a.u.]");
  // unfoldedgenResH->GetYaxis()->SetRangeUser(-1.e-15,1.e-15);
  // unfoldedgenResH->GetXaxis()->SetRangeUser(-0.5, 85.5);
  unfoldedgenResH->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.78,"1 iteration");}
  else {latex5->DrawLatex(0.31,0.62,Form("%d iterations", Iterations));}
  latex5->DrawLatex(0.2,0.80,"Residuals = Unfolded - GEN");
  // latex5->DrawLatex(0.2,0.80,Form("#chi^{2}_{1} = #sum (Refolded - REC)^{2}/REC = %0.02f",            ChiSquareTestRef));
  // latex5->DrawLatex(0.2,0.70,Form("#chi^{2}_{2} = #sum (Refolded - REC)^{2}/#DeltaRefolded = %0.02f", ChiSquareTestRef2));
  gPad->SaveAs(Form("SignalExtractionCoarse/ConvergenceClosure/residuals-iterations=%d.pdf", Iterations), "recreate");

}
