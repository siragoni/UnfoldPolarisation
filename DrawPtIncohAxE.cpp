#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>


#include "TH2.h"


//_____________________________________________________________________________
/* -
 * - Original macro:
 * - https://root.cern/doc/v610/graphShade_8C.html
 */
void DrawPtIncohAxE()
{
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetLeftMargin(0.2);
  // // gPad->SetLeftMargin(0.12);
  // // gPad->SetRightMargin(0.1);
  // gPad->SetRightMargin(0.01);
  // gPad->SetBottomMargin(0.1);
  // gPad->SetLogy();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);



  Double_t x[13]    = { 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  Double_t y[13]    = { 8.60678e-03, 2.49939e-02, 4.12754e-02, 8.18689e-02,1.63746e-01 ,2.45026e-01,3.26109e-01,4.06911e-01,4.87838e-01,5.68418e-01,6.48417e-01,7.28044e-01 ,8.07328e-01 };
  Double_t yErr[13] = { 2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04  };

  TGraphErrors *UnlikeSignDimuon = new TGraphErrors(13,x,y, yErr);
  UnlikeSignDimuon->GetHistogram()->SetMaximum(1.);   // along
  UnlikeSignDimuon->GetHistogram()->SetMinimum(0.);  //   Y
  // UnlikeSignDimuon->Rebin(5);
  // UnlikeSignDimuonMC->Rebin(5);
  // UnlikeSignDimuon->Divide(UnlikeSignDimuonMC);
  UnlikeSignDimuon->SetTitle("#lambda_{#varphi} from Toy Model, unfolded by STARlight");
  // UnlikeSignDimuon->Rebin(5);
  // UnlikeSignDimuon->Rebin(2);
  UnlikeSignDimuon->GetXaxis()->SetTitleOffset(1.15);
  // UnlikeSignDimuon->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuon->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuon->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuon->GetYaxis()->SetLabelFont(42);

  UnlikeSignDimuon->GetXaxis()->SetTitle("#lambda_{#varphi}");
  UnlikeSignDimuon->GetYaxis()->SetTitle("Value");
  // UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.,0.2);
  // UnlikeSignDimuon->GetXaxis()->SetRangeUser(0.,1.8);
  UnlikeSignDimuon->SetLineWidth(5);
  UnlikeSignDimuon->SetLineColor(2);
  // UnlikeSignDimuon->SetFillColor(kRed-3);
  // UnlikeSignDimuon->SetFillStyle(1001);
  // UnlikeSignDimuon->Draw("ep");
  UnlikeSignDimuon->Draw("");








  Double_t x2[4]    = { 0.7, 0.8, 0.9, 1.0 };
  Double_t y2[4]    = { 6.06515e-01, 6.91999e-01, 7.77218e-01, 8.62064e-01 };
  Double_t yErr2[4] = { 2.e-04,2.e-04,2.e-04,2.e-04  };

  TGraphErrors *UnlikeSignDimuon2 = new TGraphErrors(4,x2,y2, yErr2);
  UnlikeSignDimuon2->GetHistogram()->SetMaximum(1.);   // along
  UnlikeSignDimuon2->GetHistogram()->SetMinimum(0.);  //   Y
  // UnlikeSignDimuon2->SetTitle("#lambda_{#varphi} from Toy Model, unfolded by STARlight");
  UnlikeSignDimuon2->GetXaxis()->SetTitleOffset(1.15);
  // UnlikeSignDimuon->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon2->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon2->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon2->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon2->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon2->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon2->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuon2->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuon2->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuon2->GetYaxis()->SetLabelFont(42);
  UnlikeSignDimuon2->GetXaxis()->SetTitle("#lambda_{#varphi}");
  UnlikeSignDimuon2->GetYaxis()->SetTitle("Value");
  UnlikeSignDimuon2->SetLineWidth(5);
  UnlikeSignDimuon2->SetLineColor(3);




  Double_t x22[4]    = { 0.7, 0.8, 0.9, 1.0 };
  Double_t y22[4]    = { 6.58900e-01, 7.52308e-01, 8.45479e-01, 9.38363e-01 };
  Double_t yErr22[4] = { 2.e-04,2.e-04,2.e-04,2.e-04  };

  TGraphErrors *UnlikeSignDimuon22 = new TGraphErrors(4,x22,y22, yErr22);
  UnlikeSignDimuon22->GetHistogram()->SetMaximum(1.);   // along
  UnlikeSignDimuon22->GetHistogram()->SetMinimum(0.);  //   Y
  // UnlikeSignDimuon22->SetTitle("#lambda_{#varphi} from Toy Model, unfolded by STARlight");
  UnlikeSignDimuon22->GetXaxis()->SetTitleOffset(1.15);
  // UnlikeSignDimuon2->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon22->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon22->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon22->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon22->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon22->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon22->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuon22->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuon22->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuon22->GetYaxis()->SetLabelFont(42);
  UnlikeSignDimuon22->GetXaxis()->SetTitle("#lambda_{#varphi}");
  UnlikeSignDimuon22->GetYaxis()->SetTitle("Value");
  UnlikeSignDimuon22->SetLineWidth(5);
  UnlikeSignDimuon22->SetLineColor(9);






  Double_t x222[11]    = { 0.00,         0.03,        0.05,        0.1,        0.3,        0.5,        0.6,        0.7,        0.8,        0.9,        1.0 };
  Double_t y222[11]    = { -2.01537e-03, 2.25089e-02, 3.87694e-02, 7.93042e-02,2.42238e-01,4.03985e-01,4.84846e-01,5.65347e-01,6.45271e-01,7.24822e-01,8.04012e-01 };
  Double_t yErr222[11] = { 2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04,2.e-04  };

  TGraphErrors *UnlikeSignDimuon222 = new TGraphErrors(11,x222,y222, yErr222);
  UnlikeSignDimuon222->GetHistogram()->SetMaximum(1.);   // along
  UnlikeSignDimuon222->GetHistogram()->SetMinimum(0.);  //   Y
  // UnlikeSignDimuon222->SetTitle("#lambda_{#varphi} from Toy Model, unfolded by STARlight");
  UnlikeSignDimuon222->GetXaxis()->SetTitleOffset(1.15);
  // UnlikeSignDimuon2->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon222->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon222->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon222->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon222->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon222->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon222->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuon222->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuon222->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuon222->GetYaxis()->SetLabelFont(42);
  UnlikeSignDimuon222->GetXaxis()->SetTitle("#lambda_{#varphi}");
  UnlikeSignDimuon222->GetYaxis()->SetTitle("Value");
  UnlikeSignDimuon222->SetLineWidth(5);
  UnlikeSignDimuon222->SetLineColor(9);



  TMultiGraph *mg = new TMultiGraph();
  mg->Add(UnlikeSignDimuon222);
  mg->Add(UnlikeSignDimuon22);
  mg->Add(UnlikeSignDimuon2);
  mg->Add(UnlikeSignDimuon);
  mg->SetMaximum(1.6);
  mg->SetMinimum(0.0000001);
  mg->GetXaxis()->SetLimits(0.0000001,1.6);
  mg->GetXaxis()->SetTitle("Generated #lambda_{#varphi}");
  mg->GetYaxis()->SetTitle("Unfolded  #lambda_{#varphi}");
  mg->Draw("AP");


  // TLatex* latex5 = new TLatex();
  // latex5->SetTextSize(0.045);
  // latex5->SetTextFont(42);
  // latex5->SetTextAlign(11);
  // latex5->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");



  TLegend *leg_pt = new TLegend(0.5,0.45,0.85,0.79);
  leg_pt->SetFillStyle(0);
  leg_pt->SetBorderSize(0);
  leg_pt->SetTextSize(0.042);
  leg_pt->AddEntry(UnlikeSignDimuon222,"Fit with #cos(2#varphi) and #sin(#varphi + #phi)", "LP");
  // leg_pt->AddEntry(UnlikeSignDimuon22,"Uniform flag rest > 0.7", "LP");
  // leg_pt->AddEntry(UnlikeSignDimuon2, "Uniform flag rest > 0.3", "LP");
  leg_pt->AddEntry(UnlikeSignDimuon,  "Fit with only #cos(2#varphi)", "LP");
  leg_pt->Draw();

  gPad->SaveAs("lambda.pdf", "recreate");




}
