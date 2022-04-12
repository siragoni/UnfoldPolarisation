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
  histogram->SetMarkerSize(1);
  histogram->SetMarkerStyle(kOpenCircle);
  histogram->SetMarkerColor(2);
  histogram->Draw("");
}
//______________________________________________
void DrawTriggerSysHe(){
  TH1F* histoLTheta = new TH1F("TriggerLTheta", "TriggerLTheta", 2000, 0., 2.);
  histoLTheta->SetBinContent( histoLTheta->FindBin(0.85), 0.804  );
  histoLTheta->SetBinContent( histoLTheta->FindBin(0.90), 0.81   );
  histoLTheta->SetBinContent( histoLTheta->FindBin(0.95), 0.845  );
  histoLTheta->SetBinContent( histoLTheta->FindBin(1.00), 0.949  );
  histoLTheta->SetBinContent( histoLTheta->FindBin(1.05), 1.093  );
  histoLTheta->SetBinContent( histoLTheta->FindBin(1.10), 1.136  );
  histoLTheta->SetBinContent( histoLTheta->FindBin(1.15), 1.137  );

  histoLTheta->SetBinError( histoLTheta->FindBin(0.85), 0.272  );
  histoLTheta->SetBinError( histoLTheta->FindBin(0.90), 0.273  );
  histoLTheta->SetBinError( histoLTheta->FindBin(0.95), 0.270  );
  histoLTheta->SetBinError( histoLTheta->FindBin(1.00), 0.271  );
  histoLTheta->SetBinError( histoLTheta->FindBin(1.05), 0.273  );
  histoLTheta->SetBinError( histoLTheta->FindBin(1.10), 0.266  );
  histoLTheta->SetBinError( histoLTheta->FindBin(1.15), 0.253  );

  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLTheta);
  histoLTheta->GetXaxis()->SetTitle("Single muon p_{T} [GeV/#it{c}]");
  histoLTheta->GetYaxis()->SetTitle("#lambda_{#theta}");
  histoLTheta->GetXaxis()->SetRangeUser(0.75,1.2);
  histoLTheta->GetYaxis()->SetRangeUser(0., 1.6);
  histoLTheta->Draw("*");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Helicity");
  gPad->SaveAs("TriggerCoarse/trigger-ltheta-he.pdf", "recreate");
















  TH1F* histoLPhi = new TH1F("TriggerLPhi", "TriggerLPhi", 2000, 0., 2.);
  histoLPhi->SetBinContent( histoLTheta->FindBin(0.85), 0.028  );
  histoLPhi->SetBinContent( histoLTheta->FindBin(0.90), 0.028   );
  histoLPhi->SetBinContent( histoLTheta->FindBin(0.95), 0.028  );
  histoLPhi->SetBinContent( histoLTheta->FindBin(1.00), 0.032  );
  histoLPhi->SetBinContent( histoLTheta->FindBin(1.05), 0.047  );
  histoLPhi->SetBinContent( histoLTheta->FindBin(1.10), 0.062  );
  histoLPhi->SetBinContent( histoLTheta->FindBin(1.15), 0.069  );

  histoLPhi->SetBinError( histoLTheta->FindBin(0.85), 0.027  );
  histoLPhi->SetBinError( histoLTheta->FindBin(0.90), 0.027  );
  histoLPhi->SetBinError( histoLTheta->FindBin(0.95), 0.027  );
  histoLPhi->SetBinError( histoLTheta->FindBin(1.00), 0.027  );
  histoLPhi->SetBinError( histoLTheta->FindBin(1.05), 0.027  );
  histoLPhi->SetBinError( histoLTheta->FindBin(1.10), 0.026  );
  histoLPhi->SetBinError( histoLTheta->FindBin(1.15), 0.025  );

  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLPhi);
  histoLPhi->GetXaxis()->SetTitle("Single muon p_{T} [GeV/#it{c}]");
  histoLPhi->GetYaxis()->SetTitle("#lambda_{#varphi}");
  histoLPhi->GetXaxis()->SetRangeUser(0.75,1.2);
  histoLPhi->GetYaxis()->SetRangeUser(0., 0.2);
  histoLPhi->Draw("*");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Helicity");
  gPad->SaveAs("TriggerCoarse/trigger-lphi-he.pdf", "recreate");












  TH1F* histoLThetaPhi = new TH1F("TriggerLThetaPhi", "TriggerLThetaPhi", 2000, 0., 2.);
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(0.85), 0.210  );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(0.90), 0.210   );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(0.95), 0.199  );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(1.00), 0.196  );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(1.05), 0.183  );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(1.10), 0.186  );
  histoLThetaPhi->SetBinContent( histoLTheta->FindBin(1.15), 0.185  );

  histoLThetaPhi->SetBinError( histoLTheta->FindBin(0.85), 0.058  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(0.90), 0.058  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(0.95), 0.058  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(1.00), 0.058  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(1.05), 0.058  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(1.10), 0.055  );
  histoLThetaPhi->SetBinError( histoLTheta->FindBin(1.15), 0.053  );

  new TCanvas;
  BeautifyPad();
  BeautifyHisto(histoLThetaPhi);
  histoLThetaPhi->GetXaxis()->SetTitle("Single muon p_{T} [GeV/#it{c}]");
  histoLThetaPhi->GetYaxis()->SetTitle("#lambda_{#theta#varphi}");
  histoLThetaPhi->GetXaxis()->SetRangeUser(0.75,1.2);
  histoLThetaPhi->GetYaxis()->SetRangeUser(0., 0.4);
  histoLThetaPhi->Draw("*");
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex10->DrawLatex(0.2,0.80,"Helicity");
  gPad->SaveAs("TriggerCoarse/trigger-lthetaphi-he.pdf", "recreate");

}
