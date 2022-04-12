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




//______________________________________________
void BeautifyPad(){
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(11111);
}
//______________________________________________
void BeautifyHisto(TGraphErrors* histogram){
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
  histogram->SetMarkerSize(1.5);
  histogram->SetMarkerStyle(kCircle);
  histogram->SetMarkerColor(kBlue);
  histogram->Draw("");
}
//______________________________________________
void results(){

  const Int_t n = 5;
  Double_t x[n]  = {0.9806, 0.4928, 0.1978, 0.0988, 0.018};
  Double_t ex[n] = {0.0005, 0.0008, 0.0009, 0.0009, 0.001};
  Double_t y[n]  = {0.999,  0.501,  0.200,  0.099,  0.01 };
  Double_t ey[n] = {0.002,  0.003,  0.004,  0.004,  0.004};
  auto gr = new TGraphErrors(n,x,y,ex,ey);
  TCanvas* results = new TCanvas("results", "results", 900, 800);
  BeautifyPad();
  BeautifyHisto(gr);
  gr->SetTitle("Scan modulation in #varphi");
  // gr->SetMarkerColor(4);
  // gr->SetMarkerStyle(21);
  gr->GetXaxis()->SetTitle("#lambda_{#varphi}^{GENERATED}");
  gr->GetYaxis()->SetTitle("#lambda_{#varphi}^{UNFOLDED}");
  gr->Draw("ALP");
  TF1* fitfunc = new TF1("fitfunc", "[0] + [1]*x", -1, 10);
  gr->Fit("fitfunc");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.2,0.80,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.68,Form("Slope     = %.3f #pm %.3f", fitfunc->GetParameter(1),  fitfunc->GetParError(1) ));
  latex10->DrawLatex(0.2,0.62,Form("Intercept = %.3f #pm %.3f", fitfunc->GetParameter(0),  fitfunc->GetParError(0) ));
  gPad->SaveAs("ScanModulation/results.pdf", "recreate");

}
