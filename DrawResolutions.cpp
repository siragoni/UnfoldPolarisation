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
void DrawResolutions()
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



  // TFile* file = new TFile("MCtrainResults/2019-09-17/kIncohJpsiToMu/AnalysisResults.root");
  TFile* file = new TFile("SavingFileData.root");
  // TDirectory* dir;
  // dir = file->GetDirectory("MyTask");
  // TList* listings;
  // dir->GetObject("MyOutputContainer", listings);

  TH1F* UnlikeSignDimuon    = (TH1F*)file->Get("CosThetaRecMinusGenH");
  TH1F* UnlikeSignDimuonMC  = (TH1F*)file->Get("PhiRecMinusGenH");
  // UnlikeSignDimuon->Rebin(5);
  // UnlikeSignDimuonMC->Rebin(5);
  // UnlikeSignDimuon->Divide(UnlikeSignDimuonMC);
  UnlikeSignDimuon->SetTitle("");
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

  UnlikeSignDimuon->GetXaxis()->SetTitle("cos#theta_{rec} - cos#theta_{gen}");
  UnlikeSignDimuon->GetYaxis()->SetTitle("Counts [a.u.]");
  // UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.,0.2);
  UnlikeSignDimuon->GetXaxis()->SetRangeUser(-2. , 2.);
  UnlikeSignDimuon->SetLineWidth(5);
  UnlikeSignDimuon->SetLineColor(2);
  // UnlikeSignDimuon->SetFillColor(kRed-3);
  // UnlikeSignDimuon->SetFillStyle(1001);
  // UnlikeSignDimuon->Draw("ep");
  UnlikeSignDimuon->Draw("");




  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");

  gPad->SaveAs("costhetaresolution.pdf", "recreate");











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

  UnlikeSignDimuonMC->SetTitle("");
  UnlikeSignDimuonMC->GetXaxis()->SetTitleOffset(1.15);
  UnlikeSignDimuonMC->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuonMC->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuonMC->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuonMC->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuonMC->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuonMC->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuonMC->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuonMC->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuonMC->GetYaxis()->SetLabelFont(42);

  UnlikeSignDimuonMC->GetXaxis()->SetTitle("#varphi_{rec} - #varphi_{gen}");
  UnlikeSignDimuonMC->GetYaxis()->SetTitle("Counts [a.u.]");
  // UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.,0.2);
  UnlikeSignDimuonMC->GetXaxis()->SetRangeUser(-2.*TMath::Pi() , 2.*TMath::Pi());
  UnlikeSignDimuonMC->SetLineWidth(5);
  UnlikeSignDimuonMC->SetLineColor(2);
  // UnlikeSignDimuon->SetFillColor(kRed-3);
  // UnlikeSignDimuon->SetFillStyle(1001);
  // UnlikeSignDimuon->Draw("ep");
  UnlikeSignDimuonMC->Draw("");




  TLatex* latex55 = new TLatex();
  latex55->SetTextSize(0.045);
  latex55->SetTextFont(42);
  latex55->SetTextAlign(11);
  latex55->SetNDC();
  latex55->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");

  gPad->SaveAs("phiresolution.pdf", "recreate");



}
