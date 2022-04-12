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
void ToyModulation(Int_t Stripe = 14){

  // ===================
  // Retrieve REAL data
  // -------------------
  TFile* fData[24];
  TH1F* hdata[24];
  TH1F* htruth[24];
  for (Int_t i = 5; i < 20; i++) {
    fData[i]  = new TFile(Form("ScanModulation/LPhi005/ModulationYields/MonteCarloYieldsHe_phimodulation_%d.root", i));
    hdata[i]  = (TH1F*)fData[i]->Get(Form("h_%d",  i));
    htruth[i] = (TH1F*)fData[i]->Get(Form("hg_%d", i));
  }
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight or similar
  // -------------------
  TFile* f[24];
  RooUnfoldResponse*  h[24];
  for (Int_t i = 5; i < 20; i++) {
    f[i]    = new TFile(Form("ScanModulation/LPhi005/UnfoldingModulation/UnfoldedClosureHe_phimodulation_differentprior_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f[i]->Get("response;1");
  }
  Int_t M = 6;
  // if( Stripe == 0  || Stripe == 1  || Stripe == 2  || Stripe == 3  ||
  //     Stripe == 4  || Stripe == 23 || Stripe == 22 || Stripe == 21 ||
  //     Stripe == 20 || Stripe == 19 )
  // {
  //   M = 1;
  // } else if ( Stripe == 5  || Stripe == 6  || Stripe == 18  || Stripe == 17 ) {
  //   M = 6;
  // } else if ( Stripe == 7  || Stripe == 8  || Stripe == 16  || Stripe == 15 ) {
  //   M = 12;
  // } else if ( Stripe == 9  || Stripe == 10  || Stripe == 11  ||
  //             Stripe == 12 || Stripe == 13  || Stripe == 14 ) {
  //   M = 24;
  // }
  TH1F* hdata_formatted  = new TH1F("hdata_formatted", "hdata_formatted",  M, 0, 2.*TMath::Pi());
  TH1F* htruth_formatted = new TH1F("htruth_formatted","htruth_formatted", M, 0, 2.*TMath::Pi());
  for (Int_t i = 1; i < M+1; i++) {
    hdata_formatted ->SetBinContent(i, hdata[Stripe] ->GetBinContent(i));
    hdata_formatted ->SetBinError(  i, hdata[Stripe] ->GetBinError(i));
    htruth_formatted->SetBinContent(i, htruth[Stripe]->GetBinContent(i));
    htruth_formatted->SetBinError(  i, htruth[Stripe]->GetBinError(i));
  }

  // ===================
  // Unfolding
  // -------------------
  RooUnfoldBayes unfold[400];
  TH1D* hunfold[400];
  Double_t chi2_errorsfromunfold[1000];
  Double_t chi2_errorsfromgen[1000];
  TH1F* chiUnfoldH  = new TH1F("chiUnfoldH", "chiUnfoldH",  2000, -0.5, 1999.5 );
  TH1F* chiGenH     = new TH1F("chiGenH",    "chiGenH",     2000, -0.5, 1999.5 );
  for (Int_t i = 1; i < 401; i++) {
    cout << "i = " << i << endl;
    unfold[i-1]                = RooUnfoldBayes(h[Stripe], hdata_formatted, i*5);
    hunfold[i-1]               = (TH1D*) unfold[i-1].Hreco();
    chi2_errorsfromunfold[i-1] = calcChi2  (hunfold[i-1], htruth_formatted, M);
    chi2_errorsfromgen[i-1]    = calcChi2v2(hunfold[i-1], htruth_formatted, M);
    // ===================
    // Filling Chi2 plots
    // -------------------
    chiUnfoldH->SetBinContent( i*5, chi2_errorsfromunfold[i-1]/80. );
    chiUnfoldH->SetBinError(   i*5, 0. );
    chiGenH   ->SetBinContent( i*5, chi2_errorsfromgen[i-1]/80. );
    chiGenH   ->SetBinError(   i*5, 0. );
  }




  new TCanvas;
  BeautifyPad();
  BeautifyHisto(chiUnfoldH);
  gPad->SetLogy();
  chiUnfoldH->GetXaxis()->SetTitle("Iterations");
  chiUnfoldH->GetYaxis()->SetTitle("(Unfolded - GEN)/(#DeltaUnfolded NDF)  [a.u.]");
  chiUnfoldH->GetYaxis()->SetRangeUser(0.1,chiUnfoldH->GetMaximum()*(10.));
  chiUnfoldH->GetXaxis()->SetRangeUser(-0.5, 1999.5);
  chiUnfoldH->Draw();
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis");
  latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", Stripe));
  gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/chi2-unfold-errors-modulation-different-prior-%d.pdf", Stripe), "recreate");
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(chiGenH);
  gPad->SetLogy();
  chiGenH->GetXaxis()->SetTitle("Iterations");
  chiGenH->GetYaxis()->SetTitle("(Unfolded - GEN)/(GEN NDF)  [a.u.]");
  chiGenH->GetYaxis()->SetRangeUser(0.1,chiGenH->GetMaximum()*(10.));
  chiGenH->GetXaxis()->SetRangeUser(-0.5, 1999.5);
  chiGenH->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.74,"This thesis");
  latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", Stripe));
  gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/chi2-generated-errors-modulation-different-prior-%d.pdf", Stripe), "recreate");









  // ===================================
  // Covariance Matrices and average Rho
  // -----------------------------------
  TMatrixD** Covariance;
  Covariance = new TMatrixD*[24];
  Covariance[0]  = new TMatrixD(6,6);
  Covariance[1]  = new TMatrixD(6,6);
  Covariance[2]  = new TMatrixD(6,6);
  Covariance[3]  = new TMatrixD(6,6);
  Covariance[4]  = new TMatrixD(6,6);
  Covariance[5]  = new TMatrixD(6,6);
  Covariance[6]  = new TMatrixD(6,6);
  Covariance[7]  = new TMatrixD(6,6);
  Covariance[8]  = new TMatrixD(6,6);
  Covariance[9]  = new TMatrixD(6,6);
  Covariance[10] = new TMatrixD(6,6);
  Covariance[11] = new TMatrixD(6,6);
  Covariance[12] = new TMatrixD(6,6);
  Covariance[13] = new TMatrixD(6,6);
  Covariance[14] = new TMatrixD(6,6);
  Covariance[15] = new TMatrixD(6,6);
  Covariance[16] = new TMatrixD(6,6);
  Covariance[17] = new TMatrixD(6,6);
  Covariance[18] = new TMatrixD(6,6);
  Covariance[19] = new TMatrixD(6,6);
  Covariance[20] = new TMatrixD(6,6);
  Covariance[21] = new TMatrixD(6,6);
  Covariance[22] = new TMatrixD(6,6);
  Covariance[23] = new TMatrixD(6,6);
  TMatrixD** CovarianceCopy;
  CovarianceCopy = new TMatrixD*[24];
  CovarianceCopy[0]  = new TMatrixD(6,6);
  CovarianceCopy[1]  = new TMatrixD(6,6);
  CovarianceCopy[2]  = new TMatrixD(6,6);
  CovarianceCopy[3]  = new TMatrixD(6,6);
  CovarianceCopy[4]  = new TMatrixD(6,6);
  CovarianceCopy[5]  = new TMatrixD(6,6);
  CovarianceCopy[6]  = new TMatrixD(6,6);
  CovarianceCopy[7]  = new TMatrixD(6,6);
  CovarianceCopy[8]  = new TMatrixD(6,6);
  CovarianceCopy[9]  = new TMatrixD(6,6);
  CovarianceCopy[10] = new TMatrixD(6,6);
  CovarianceCopy[11] = new TMatrixD(6,6);
  CovarianceCopy[12] = new TMatrixD(6,6);
  CovarianceCopy[13] = new TMatrixD(6,6);
  CovarianceCopy[14] = new TMatrixD(6,6);
  CovarianceCopy[15] = new TMatrixD(6,6);
  CovarianceCopy[16] = new TMatrixD(6,6);
  CovarianceCopy[17] = new TMatrixD(6,6);
  CovarianceCopy[18] = new TMatrixD(6,6);
  CovarianceCopy[19] = new TMatrixD(6,6);
  CovarianceCopy[20] = new TMatrixD(6,6);
  CovarianceCopy[21] = new TMatrixD(6,6);
  CovarianceCopy[22] = new TMatrixD(6,6);
  CovarianceCopy[23] = new TMatrixD(6,6);
  Double_t Determinant[24];
  Double_t rho[400];
  TH1F* averagerhoH = new TH1F("averagerhoH", "averagerhoH", 2000, -0.5, 1999.5 );
  for (Int_t i = 1; i < 401; i++) {
    *Covariance[Stripe]     = unfold[i-1].Ereco(RooUnfold::kCovariance);
    *CovarianceCopy[Stripe] = *Covariance[Stripe];
    CovarianceCopy[Stripe]->InvertFast(&Determinant[Stripe]);
    for (Int_t j = 0; j < M; j++) {
      rho[i-1] += TMath::Sqrt(1. - 1./( (*Covariance[Stripe])(j,j) * (*CovarianceCopy[Stripe])(j,j) ) );
    }
    rho[i-1] /= ((Double_t)M);
    // =========================
    // Filling average Rho plots
    // -------------------------
    averagerhoH->SetBinContent( i*5, rho[i-1] );
    averagerhoH->SetBinError(   i*5, 0. );
  }
  new TCanvas;
  BeautifyPad();
  BeautifyHisto(averagerhoH);
  // gPad->SetLogy();
  averagerhoH->GetXaxis()->SetTitle("Iterations");
  averagerhoH->GetYaxis()->SetTitle("<#rho_{i}> [a.u.]");
  averagerhoH->GetYaxis()->SetRangeUser(averagerhoH->GetBinContent(900)*0.99,1.001);
  if( Stripe == 5  ) averagerhoH->GetYaxis()->SetRangeUser(0.650,1.000);
  if( Stripe == 6  ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 7  ) averagerhoH->GetYaxis()->SetRangeUser(0.600,0.990);
  if( Stripe == 8  ) averagerhoH->GetYaxis()->SetRangeUser(0.600,0.990);
  if( Stripe == 9  ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 10 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 11 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 12 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 13 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 14 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 15 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,0.990);
  if( Stripe == 16 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,0.990);
  if( Stripe == 17 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  if( Stripe == 18 ) averagerhoH->GetYaxis()->SetRangeUser(0.600,1.000);
  averagerhoH->GetXaxis()->SetRangeUser(-0.5, 1999.5);
  averagerhoH->Draw();
  latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", Stripe));
  gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/average-rho-data-prior-%d.pdf", Stripe), "recreate");
  averagerhoH->GetXaxis()->SetRangeUser(-0.5, 49.5);
  gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/average-rho-data-prior-zoom-%d.pdf", Stripe), "recreate");
  averagerhoH->GetXaxis()->SetRangeUser(-0.5, 1999.5);
  gPad->SaveAs(Form("ScanModulation/LPhi005/ConvergenceModulation/average-rho-data-prior-%d.pdf", Stripe), "recreate");





















}
