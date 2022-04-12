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
void UnfoldData(){

  // ===================
  // Retrieve REAL data
  // -------------------
  TFile* fData[24];
  TH1F* hdata[24];
  for (Int_t i = 5; i < 20; i++) {
    fData[i]  = new TFile(Form("SignalExtraction/RawYieldsHeV3_%d.root", i));
    hdata[i]  = (TH1F*)fData[i]->Get(Form("h_%d",  i));
  }
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight or similar
  // -------------------
  TFile* f[24];
  RooUnfoldResponse*  h[24];
  TH1F* hdata_formatted[24];
  Int_t M = 1;
  for (Int_t i = 5; i < 20; i++) {
    f[i]    = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f[i]->Get("response;1");
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      M = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      M = 6;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      M = 12;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      M = 24;
    }
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
  Int_t N = 1;
  // -------------------
  // How many iterations?
  // From the analysis
  // of the correlation
  // coefficients
  // -------------------
  for (Int_t i = 5; i < 19; i++) {
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 15;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 200;
    } else if ( i == 11 || i == 13 || i == 14 ) {
      N = 700;
    } else if ( i == 9  || i == 10  ) {
      N = 1200;
    } else if ( i == 12 ) {
      N = 1000;
    }
    unfold[i]  = RooUnfoldBayes(h[i], hdata_formatted[i], N);
    hunfold[i] = (TH1D*) unfold[i].Hreco();
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      M = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      M = 6;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      M = 12;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      M = 24;
    }
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
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 15;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 200;
    } else if ( i == 11 || i == 13 || i == 14 ) {
      N = 700;
    } else if ( i == 9  || i == 10  ) {
      N = 1200;
    } else if ( i == 12 ) {
      N = 1000;
    }
    latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    latex5->DrawLatex(0.31,0.82,Form("Iterations %d", N));
    latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", i));
    gPad->SaveAs(Form("Unfolding2D/Convergence/unfolded-data-%d.pdf", i), "recreate");

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
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 15;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 200;
    } else if ( i == 11 || i == 13 || i == 14 ) {
      N = 700;
    } else if ( i == 9  || i == 10  ) {
      N = 1200;
    } else if ( i == 12 ) {
      N = 1000;
    }
    latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
    latex5->DrawLatex(0.31,0.82,Form("Iterations %d", N));
    latex5->DrawLatex(0.31,0.74,Form("This thesis, Slice in cos#theta = %d", i));
    gPad->SaveAs(Form("Unfolding2D/Convergence/unfolded-data-relative-uncertainties-%d.pdf", i), "recreate");

  }






















}
