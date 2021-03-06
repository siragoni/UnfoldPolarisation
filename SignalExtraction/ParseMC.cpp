#include <TROOT.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TBox.h>
#include <TWbox.h>
#include <TLine.h>
#include <TGraphAsymmErrors.h>
#include <TMarker.h>
#include <stdio.h>
#include <iostream>
// #include "RooUnfoldResponse.h"
// #include <RooUnfold.h>


void ParseMC(){

  // TFile* file   = new TFile("AnalysisResultsCoh.root", "READ");
  TFile* file   = new TFile("AnalysisResultsCohWithSinglePt.root", "READ");
  TTree* ReconTree = (TTree*)file->Get("MyTask/fOutputTree");
  TTree* GenerTree = (TTree*)file->Get("MyTask/fOutputTreeMC");



  Double_t pTrec;
  Double_t pTgen;
  ReconTree->SetBranchAddress("fTrkTrkPt",  &pTrec);
  GenerTree->SetBranchAddress("fMCTrkTrkPt",&pTgen);
  Double_t Mrec;
  Double_t Mgen;
  ReconTree->SetBranchAddress("fTrkTrkM",  &Mrec);
  GenerTree->SetBranchAddress("fMCTrkTrkM",&Mgen);
  Double_t Yrec;
  Double_t Ygen;
  ReconTree->SetBranchAddress("fTrkTrkY",  &Yrec);
  GenerTree->SetBranchAddress("fMCTrkTrkY",&Ygen);
  Double_t CosThetaHErec;
  Double_t CosThetaHEgen;
  Double_t CosThetaHErec2;
  Double_t CosThetaHEgen2;
  ReconTree->SetBranchAddress("fCosThetaHE",  &CosThetaHErec);
  GenerTree->SetBranchAddress("fMCCosThetaHE",&CosThetaHEgen);
  Double_t CosThetaCSrec;
  Double_t CosThetaCSgen;
  ReconTree->SetBranchAddress("fCosThetaCS",  &CosThetaCSrec);
  GenerTree->SetBranchAddress("fMCCosThetaCS",&CosThetaCSgen);
  Double_t PhiHErec;
  Double_t PhiHEgen;
  Double_t PhiHErec2;
  Double_t PhiHEgen2;
  ReconTree->SetBranchAddress("fPhiHE",  &PhiHErec);
  GenerTree->SetBranchAddress("fMCPhiHE",&PhiHEgen);
  Double_t PhiCSrec;
  Double_t PhiCSgen;
  ReconTree->SetBranchAddress("fPhiCS",  &PhiCSrec);
  GenerTree->SetBranchAddress("fMCPhiCS",&PhiCSgen);


  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


  TH1F *PtRecH = new TH1F("PtRecH", "PtRecH", 4000, 0, 20);
  TH1F *PtGenH = new TH1F("PtGenH", "PtGenH", 4000, 0, 20);

  TH1F *PhiRecMinusGenH      = new TH1F("PhiRecMinusGenH",      "PhiRecMinusGenH",      50, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH1F *CosThetaRecMinusGenH = new TH1F("CosThetaRecMinusGenH", "CosThetaRecMinusGenH", 50, -2., 2.);


  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 25, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 25, -1., 1.);
  // TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      25, 0., 2.*TMath::Pi());
  // TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      25, 0., 2.*TMath::Pi());
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      250, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      250, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());



  TH1F* InvMassH_binzero[1];  // -1 + 1*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binzero[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binzero_%d", iCosThetaBins),
                // Form("InvMassH_binzero_%d", iCosThetaBins),
                Form("InvMassH_0_%d", iCosThetaBins),
                Form("InvMassH_0_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binone[1];  // -1 + 2*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binone[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binone_%d", iCosThetaBins),
                // Form("InvMassH_binone_%d", iCosThetaBins),
                Form("InvMassH_1_%d", iCosThetaBins),
                Form("InvMassH_1_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwo[1];  // -1 + 3*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwo[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwo_%d", iCosThetaBins),
                // Form("InvMassH_bintwo_%d", iCosThetaBins),
                Form("InvMassH_2_%d", iCosThetaBins),
                Form("InvMassH_2_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binthree[1];  // -1 + 4*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binthree[iCosThetaBins] = new TH1F(
                // Form("InvMassH_3_%d", iCosThetaBins),
                // Form("InvMassH_3_%d", iCosThetaBins),
                Form("InvMassH_3_%d", iCosThetaBins),
                Form("InvMassH_3_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfour[1];  // -1 + 5*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binfour[iCosThetaBins] = new TH1F(
                // Form("InvMassH_4_%d", iCosThetaBins),
                // Form("InvMassH_4_%d", iCosThetaBins),
                Form("InvMassH_4_%d", iCosThetaBins),
                Form("InvMassH_4_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfive[6];  // -1 + 6*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binfive[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfive_%d", iCosThetaBins),
                // Form("InvMassH_binfive_%d", iCosThetaBins),
                Form("InvMassH_5_%d", iCosThetaBins),
                Form("InvMassH_5_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binsix[6];  // -1 + 7*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binsix[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binsix_%d", iCosThetaBins),
                // Form("InvMassH_binsix_%d", iCosThetaBins),
                Form("InvMassH_6_%d", iCosThetaBins),
                Form("InvMassH_6_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binseven[12];  // -1 + 8*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binseven[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binseven_%d", iCosThetaBins),
                // Form("InvMassH_binseven_%d", iCosThetaBins),
                Form("InvMassH_7_%d", iCosThetaBins),
                Form("InvMassH_7_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineight[12];  // -1 + 9*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_bineight[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineight_%d", iCosThetaBins),
                // Form("InvMassH_bineight_%d", iCosThetaBins),
                Form("InvMassH_8_%d", iCosThetaBins),
                Form("InvMassH_8_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binnine[24];  // -1 + 10*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binnine[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binnine_%d", iCosThetaBins),
                // Form("InvMassH_binnine_%d", iCosThetaBins),
                Form("InvMassH_9_%d", iCosThetaBins),
                Form("InvMassH_9_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binten[24];  // -1 + 11*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binten[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binten_%d", iCosThetaBins),
                // Form("InvMassH_binten_%d", iCosThetaBins),
                Form("InvMassH_10_%d", iCosThetaBins),
                Form("InvMassH_10_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineleven[24];  // -1 + 12*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_bineleven[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineleven_%d", iCosThetaBins),
                // Form("InvMassH_bineleven_%d", iCosThetaBins),
                Form("InvMassH_11_%d", iCosThetaBins),
                Form("InvMassH_11_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwelve[24];  // -1 + 13*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_bintwelve[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwelve_%d", iCosThetaBins),
                // Form("InvMassH_bintwelve_%d", iCosThetaBins),
                Form("InvMassH_12_%d", iCosThetaBins),
                Form("InvMassH_12_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binthirteen[24];  // -1 + 14*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binthirteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binthirteen_%d", iCosThetaBins),
                // Form("InvMassH_binthirteen_%d", iCosThetaBins),
                Form("InvMassH_13_%d", iCosThetaBins),
                Form("InvMassH_13_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfourteen[24];  // -1 + 15*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binfourteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfourteen_%d", iCosThetaBins),
                // Form("InvMassH_binfourteen_%d", iCosThetaBins),
                Form("InvMassH_14_%d", iCosThetaBins),
                Form("InvMassH_14_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfifteen[12];  // -1 + 16*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binfifteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfifteen_%d", iCosThetaBins),
                // Form("InvMassH_binfifteen_%d", iCosThetaBins),
                Form("InvMassH_15_%d", iCosThetaBins),
                Form("InvMassH_15_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binsixteen[12];  // -1 + 17*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binsixteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binsixteen_%d", iCosThetaBins),
                // Form("InvMassH_binsixteen_%d", iCosThetaBins),
                Form("InvMassH_16_%d", iCosThetaBins),
                Form("InvMassH_16_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binseventeen[6];  // -1 + 18*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binseventeen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binseventeen_%d", iCosThetaBins),
                // Form("InvMassH_binseventeen_%d", iCosThetaBins),
                Form("InvMassH_17_%d", iCosThetaBins),
                Form("InvMassH_17_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineighteen[6];  // -1 + 19*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_bineighteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineighteen_%d", iCosThetaBins),
                // Form("InvMassH_bineighteen_%d", iCosThetaBins),
                Form("InvMassH_18_%d", iCosThetaBins),
                Form("InvMassH_18_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binnineteen[1];  // -1 + 20*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binnineteen[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binnineteen_%d", iCosThetaBins),
                // Form("InvMassH_binnineteen_%d", iCosThetaBins),
                Form("InvMassH_19_%d", iCosThetaBins),
                Form("InvMassH_19_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwenty[1];  // -1 + 21*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwenty[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwenty_%d", iCosThetaBins),
                // Form("InvMassH_bintwenty_%d", iCosThetaBins),
                Form("InvMassH_20_%d", iCosThetaBins),
                Form("InvMassH_20_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentyone[1];  // -1 + 22*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentyone[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentyone_%d", iCosThetaBins),
                // Form("InvMassH_bintwentyone_%d", iCosThetaBins),
                Form("InvMassH_21_%d", iCosThetaBins),
                Form("InvMassH_21_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentytwo[1];  // -1 + 23*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentytwo[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentytwo_%d", iCosThetaBins),
                // Form("InvMassH_bintwentytwo_%d", iCosThetaBins),
                Form("InvMassH_22_%d", iCosThetaBins),
                Form("InvMassH_22_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentythree[1];  // -1 + 24*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentythree[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentythree_%d", iCosThetaBins),
                // Form("InvMassH_bintwentythree_%d", iCosThetaBins),
                Form("InvMassH_23_%d", iCosThetaBins),
                Form("InvMassH_23_%d", iCosThetaBins),
                2000, 0, 20
                );
  }



  TH1F* InvMassH_closure[24];  // -1 + 24*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_closure[iCosThetaBins] = new TH1F(
                // Form("InvMassH_closure%d", iCosThetaBins),
                // Form("InvMassH_closure%d", iCosThetaBins),
                Form("InvMassH_closure_%d", iCosThetaBins),
                Form("InvMassH_closure_%d", iCosThetaBins),
                2000, 0, 20
                );
  }











  // Gen level plots
  TH1F* InvMassH_binzerog[1];  // -1 + 1*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binzerog[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binzero_%d", iCosThetaBins),
                // Form("InvMassH_binzero_%d", iCosThetaBins),
                Form("InvMassHg_0_%d", iCosThetaBins),
                Form("InvMassHg_0_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binoneg[1];  // -1 + 2*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binoneg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binone_%d", iCosThetaBins),
                // Form("InvMassH_binone_%d", iCosThetaBins),
                Form("InvMassHg_1_%d", iCosThetaBins),
                Form("InvMassHg_1_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwog[1];  // -1 + 3*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwog[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwo_%d", iCosThetaBins),
                // Form("InvMassH_bintwo_%d", iCosThetaBins),
                Form("InvMassHg_2_%d", iCosThetaBins),
                Form("InvMassHg_2_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binthreeg[1];  // -1 + 4*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binthreeg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_3_%d", iCosThetaBins),
                // Form("InvMassH_3_%d", iCosThetaBins),
                Form("InvMassHg_3_%d", iCosThetaBins),
                Form("InvMassHg_3_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfourg[1];  // -1 + 5*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binfourg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_4_%d", iCosThetaBins),
                // Form("InvMassH_4_%d", iCosThetaBins),
                Form("InvMassHg_4_%d", iCosThetaBins),
                Form("InvMassHg_4_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfiveg[6];  // -1 + 6*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binfiveg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfive_%d", iCosThetaBins),
                // Form("InvMassH_binfive_%d", iCosThetaBins),
                Form("InvMassHg_5_%d", iCosThetaBins),
                Form("InvMassHg_5_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binsixg[6];  // -1 + 7*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binsixg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binsix_%d", iCosThetaBins),
                // Form("InvMassH_binsix_%d", iCosThetaBins),
                Form("InvMassHg_6_%d", iCosThetaBins),
                Form("InvMassHg_6_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binseveng[12];  // -1 + 8*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binseveng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binseven_%d", iCosThetaBins),
                // Form("InvMassH_binseven_%d", iCosThetaBins),
                Form("InvMassHg_7_%d", iCosThetaBins),
                Form("InvMassHg_7_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineightg[12];  // -1 + 9*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_bineightg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineight_%d", iCosThetaBins),
                // Form("InvMassH_bineight_%d", iCosThetaBins),
                Form("InvMassHg_8_%d", iCosThetaBins),
                Form("InvMassHg_8_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binnineg[24];  // -1 + 10*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binnineg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binnine_%d", iCosThetaBins),
                // Form("InvMassH_binnine_%d", iCosThetaBins),
                Form("InvMassHg_9_%d", iCosThetaBins),
                Form("InvMassHg_9_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binteng[24];  // -1 + 11*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binteng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binten_%d", iCosThetaBins),
                // Form("InvMassH_binten_%d", iCosThetaBins),
                Form("InvMassHg_10_%d", iCosThetaBins),
                Form("InvMassHg_10_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineleveng[24];  // -1 + 12*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_bineleveng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineleven_%d", iCosThetaBins),
                // Form("InvMassH_bineleven_%d", iCosThetaBins),
                Form("InvMassHg_11_%d", iCosThetaBins),
                Form("InvMassHg_11_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwelveg[24];  // -1 + 13*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_bintwelveg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwelve_%d", iCosThetaBins),
                // Form("InvMassH_bintwelve_%d", iCosThetaBins),
                Form("InvMassHg_12_%d", iCosThetaBins),
                Form("InvMassHg_12_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binthirteeng[24];  // -1 + 14*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binthirteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binthirteen_%d", iCosThetaBins),
                // Form("InvMassH_binthirteen_%d", iCosThetaBins),
                Form("InvMassHg_13_%d", iCosThetaBins),
                Form("InvMassHg_13_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfourteeng[24];  // -1 + 15*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 24; iCosThetaBins++ ){
    InvMassH_binfourteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfourteen_%d", iCosThetaBins),
                // Form("InvMassH_binfourteen_%d", iCosThetaBins),
                Form("InvMassHg_14_%d", iCosThetaBins),
                Form("InvMassHg_14_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfifteeng[12];  // -1 + 16*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binfifteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfifteen_%d", iCosThetaBins),
                // Form("InvMassH_binfifteen_%d", iCosThetaBins),
                Form("InvMassHg_15_%d", iCosThetaBins),
                Form("InvMassHg_15_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binsixteeng[12];  // -1 + 17*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 12; iCosThetaBins++ ){
    InvMassH_binsixteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binsixteen_%d", iCosThetaBins),
                // Form("InvMassH_binsixteen_%d", iCosThetaBins),
                Form("InvMassHg_16_%d", iCosThetaBins),
                Form("InvMassHg_16_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binseventeeng[6];  // -1 + 18*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_binseventeeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binseventeen_%d", iCosThetaBins),
                // Form("InvMassH_binseventeen_%d", iCosThetaBins),
                Form("InvMassHg_17_%d", iCosThetaBins),
                Form("InvMassHg_17_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bineighteeng[6];  // -1 + 19*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 6; iCosThetaBins++ ){
    InvMassH_bineighteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bineighteen_%d", iCosThetaBins),
                // Form("InvMassH_bineighteen_%d", iCosThetaBins),
                Form("InvMassHg_18_%d", iCosThetaBins),
                Form("InvMassHg_18_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binnineteeng[1];  // -1 + 20*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binnineteeng[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binnineteen_%d", iCosThetaBins),
                // Form("InvMassH_binnineteen_%d", iCosThetaBins),
                Form("InvMassHg_19_%d", iCosThetaBins),
                Form("InvMassHg_19_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentyg[1];  // -1 + 21*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentyg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwenty_%d", iCosThetaBins),
                // Form("InvMassH_bintwenty_%d", iCosThetaBins),
                Form("InvMassHg_20_%d", iCosThetaBins),
                Form("InvMassHg_20_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentyoneg[1];  // -1 + 22*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentyoneg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentyone_%d", iCosThetaBins),
                // Form("InvMassH_bintwentyone_%d", iCosThetaBins),
                Form("InvMassHg_21_%d", iCosThetaBins),
                Form("InvMassHg_21_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentytwog[1];  // -1 + 23*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentytwog[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentytwo_%d", iCosThetaBins),
                // Form("InvMassH_bintwentytwo_%d", iCosThetaBins),
                Form("InvMassHg_22_%d", iCosThetaBins),
                Form("InvMassHg_22_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_bintwentythreeg[1];  // -1 + 24*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_bintwentythreeg[iCosThetaBins] = new TH1F(
                // Form("InvMassH_bintwentythree_%d", iCosThetaBins),
                // Form("InvMassH_bintwentythree_%d", iCosThetaBins),
                Form("InvMassHg_23_%d", iCosThetaBins),
                Form("InvMassHg_23_%d", iCosThetaBins),
                2000, 0, 20
                );
  }




  TH2F* generatedlevel = new TH2F("generatedlevel","generatedlevel", 1000, 0, 2.*TMath::Pi(), 1000, -1, 1);
  TH2F* reconlevel     = new TH2F("reconlevel",    "reconlevel",     1000, 0, 2.*TMath::Pi(), 1000, -1, 1);







  Double_t controlFlag  = 0;
  Double_t controlFlag2 = 0;
  Double_t controlFlag3 = 0;

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);

      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {

        CosThetaHEgen2 = -1.*CosThetaHEgen;
        if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
        PhiHEgen2 = PhiHEgen+TMath::Pi();
        generatedlevel->Fill(PhiHEgen2, CosThetaHEgen2);















        controlFlag3 = 0;
        // if (        CosThetaHErec < (-1. + 1.*(0.08+0.01/3.)) ){
        if (        CosThetaHEgen2 < (-1. + 1.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_binzerog[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 2.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 2.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_binoneg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 3.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 3.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_bintwog[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 4.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 4.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_binthreeg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 5.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 5.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_binfourg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 6.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 6.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 6; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
              InvMassH_binfiveg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 7.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 7.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 6; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
              InvMassH_binsixg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 8.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 8.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 12; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
              InvMassH_binseveng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 9.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 9.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 12; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
              InvMassH_bineightg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 10.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 10.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_binnineg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 11.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 11.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_binteng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 12.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 12.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_bineleveng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 13.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 13.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_bintwelveg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 14.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 14.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_binthirteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 15.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 15.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 24; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
              InvMassH_binfourteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 16.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 16.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 12; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
              InvMassH_binfifteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 17.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 17.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 12; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
              InvMassH_binsixteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 18.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 18.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 6; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
              InvMassH_binseventeeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 19.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 19.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 6; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
              InvMassH_bineighteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 20.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 20.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_binnineteeng[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 21.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 21.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_bintwentyg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 22.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 22.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_bintwentyoneg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 23.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 23.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_bintwentytwog[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        // } else if ( CosThetaHErec < (-1. + 24.*(0.08+0.01/3.)) ){
        } else if ( CosThetaHEgen2 < (-1. + 24.*(0.08+0.01/3.)) ){
          for(Int_t i = 0; i < 1; i++) {
            if( controlFlag3 == 1) break;
            if( (PhiHEgen2) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
              InvMassH_bintwentythreeg[i]->Fill(Mrec);
              controlFlag3 = 1;
            }
          }
        }












          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {

                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  reconlevel->Fill(PhiHErec, CosThetaHErec);


                  controlFlag2 = 0;
                  if (        CosThetaHEgen2 > -1. && CosThetaHEgen2 < 1. ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag2 == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_closure[i]->Fill(Mrec);
                        controlFlag2 = 1;
                      }
                    }
                  }



                  controlFlag = 0;
                  // if (        CosThetaHErec < (-1. + 1.*(0.08+0.01/3.)) ){
                  if (        CosThetaHEgen2 < (-1. + 1.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binzero[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 2.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 2.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binone[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 3.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 3.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwo[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 4.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 4.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binthree[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 5.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 5.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binfour[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 6.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 6.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binfive[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 7.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 7.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binsix[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 8.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 8.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binseven[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 9.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 9.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_bineight[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 10.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 10.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binnine[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 11.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 11.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binten[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 12.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 12.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_bineleven[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 13.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 13.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_bintwelve[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 14.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 14.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binthirteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 15.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 15.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binfourteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 16.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 16.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binfifteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 17.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 17.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binsixteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 18.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 18.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binseventeen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 19.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 19.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_bineighteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 20.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 20.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binnineteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 21.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 21.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwenty[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 22.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 22.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentyone[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 23.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 23.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentytwo[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  // } else if ( CosThetaHErec < (-1. + 24.*(0.08+0.01/3.)) ){
                  } else if ( CosThetaHEgen2 < (-1. + 24.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentythree[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  }




          }
      }
    }
  }













  TFile *SavingFile = new TFile("SignalExtraction/CohHe.root", "RECREATE");
  InvMassH_binzero[0]->Write();
  InvMassH_binone[0]->Write();
  InvMassH_bintwo[0]->Write();
  InvMassH_binthree[0]->Write();
  InvMassH_binfour[0]->Write();
  InvMassH_binnineteen[0]->Write();
  InvMassH_bintwenty[0]->Write();
  InvMassH_bintwentyone[0]->Write();
  InvMassH_bintwentytwo[0]->Write();
  InvMassH_bintwentythree[0]->Write();
  // InvMassH_bintwentyfour[0]->Write();
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binfive[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binsix[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binseven[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_bineight[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binnine[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binten[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_bineleven[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_bintwelve[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binthirteen[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binfourteen[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binfifteen[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binsixteen[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binseventeen[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_bineighteen[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_closure[i]->Write();
  }
  InvMassH_binzerog[0]->Write();
  InvMassH_binoneg[0]->Write();
  InvMassH_bintwog[0]->Write();
  InvMassH_binthreeg[0]->Write();
  InvMassH_binfourg[0]->Write();
  InvMassH_binnineteeng[0]->Write();
  InvMassH_bintwentyg[0]->Write();
  InvMassH_bintwentyoneg[0]->Write();
  InvMassH_bintwentytwog[0]->Write();
  InvMassH_bintwentythreeg[0]->Write();
  // InvMassH_bintwentyfour[0]->Write();
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binfiveg[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binsixg[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binseveng[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_bineightg[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binnineg[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binteng[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_bineleveng[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_bintwelveg[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binthirteeng[i]->Write();
  }
  for(Int_t i = 0; i < 24; i++ ){
    InvMassH_binfourteeng[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binfifteeng[i]->Write();
  }
  for(Int_t i = 0; i < 12; i++ ){
    InvMassH_binsixteeng[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_binseventeeng[i]->Write();
  }
  for(Int_t i = 0; i < 6; i++ ){
    InvMassH_bineighteeng[i]->Write();
  }
  generatedlevel->Write();
  reconlevel->Write();
  SavingFile       ->Close();














  Int_t Counter[24] = {1,1,1,1,1, 6,6, 12,12, 24,24,24,24,24,24, 12,12, 6,6, 1,1,1,1,1 };

  Double_t CosThetaCenters[24];
  Double_t PhiCenters[24];

  PhiCenters[0] = PhiCenters[1] = PhiCenters[2]  = PhiCenters[3]  = PhiCenters[4] = PhiCenters[23] = PhiCenters[22] = PhiCenters[21] = PhiCenters[20] = PhiCenters[19] = TMath::Pi();
  Double_t Spacing1 = TMath::Pi()/12.;
  PhiCenters[5] = PhiCenters[6] = PhiCenters[18] = PhiCenters[17] = Spacing1;
  Double_t Spacing2 = TMath::Pi()/24.;
  PhiCenters[7] = PhiCenters[8] = PhiCenters[16] = PhiCenters[15] = Spacing2;
  Double_t Spacing3 = TMath::Pi()/48.;
  PhiCenters[9] = PhiCenters[10] = PhiCenters[11] = PhiCenters[12] = PhiCenters[13] = PhiCenters[14] = Spacing3;


  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
  }

  TH1F* RawYields[24];
  for (Int_t i = 0; i < 24; i++) {
    RawYields[i] = new TH1F(Form("h_%d", i), Form("h_%d", i), 100, -0.5, 99.5);
  }
  TH1F* RawYieldsg[24];
  for (Int_t i = 0; i < 24; i++) {
    RawYieldsg[i] = new TH1F(Form("hg_%d", i), Form("hg_%d", i), 100, -0.5, 99.5);
  }
  TFile* parsedMC = new TFile("SignalExtraction/CohHe.root");

  Double_t JPsiPeakValue    = 0;
  Double_t JPsiPeakValueErr = 0;
  Double_t BkgValue         = 0;
  Double_t BkgValueError    = 0;
  Double_t JPsiPeakValue2    = 0;
  Double_t JPsiPeakValueErr2 = 0;
  Double_t BkgValue2         = 0;
  Double_t BkgValueError2    = 0;
  TH1F* fCohJpsiToMu        = 0x0;
  TH1F* fCohJpsiToMu2        = 0x0;
  for (Int_t iCosThetaBins = 4; iCosThetaBins < 20; iCosThetaBins++) {

    for (Int_t iPhiBins = 0; iPhiBins < Counter[iCosThetaBins]; iPhiBins++) {

      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;
      JPsiPeakValue2    = 0;
      JPsiPeakValueErr2 = 0;
      BkgValue2         = 0;
      BkgValueError2    = 0;

      fCohJpsiToMu     = (TH1F*)parsedMC->Get( Form( "InvMassH_%d_%d", iCosThetaBins, iPhiBins ) );
      JPsiPeakValue    = (Double_t) fCohJpsiToMu->GetEntries();
      JPsiPeakValueErr = TMath::Sqrt((Double_t) fCohJpsiToMu->GetEntries());
      fCohJpsiToMu2     = (TH1F*)parsedMC->Get( Form( "InvMassHg_%d_%d", iCosThetaBins, iPhiBins ) );
      JPsiPeakValue2    = (Double_t) fCohJpsiToMu2->GetEntries();
      JPsiPeakValueErr2 = TMath::Sqrt((Double_t) fCohJpsiToMu2->GetEntries());

      // RawYields[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue   /((PhiCenters[iCosThetaBins]*2.)*(0.08+0.01/3.)));
      // RawYields[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr/((PhiCenters[iCosThetaBins]*2.)*(0.08+0.01/3.)));
      RawYields[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue   );
      RawYields[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr);
      RawYieldsg[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue2   );
      RawYieldsg[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr2);
      // RawYields[iCosThetaBins]->SetBinContent(       iPhiBins+1,   JPsiPeakValue   );
      // RawYields[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr);

    }
    TFile f(Form("SignalExtraction/MonteCarloYieldsHe_%d.root", iCosThetaBins),   "recreate");
    RawYields[iCosThetaBins]->Write();
    RawYieldsg[iCosThetaBins]->Write();
    f.Close();

  }


  TH1F* RawYieldsSimpleClosure = new TH1F("h", "h", 100, -0.5, 99.5);
  for (Int_t iPhiBins = 0; iPhiBins < 24; iPhiBins++) {

    JPsiPeakValue    = 0;
    JPsiPeakValueErr = 0;
    BkgValue         = 0;
    BkgValueError    = 0;

    fCohJpsiToMu     = (TH1F*)parsedMC->Get( Form( "InvMassH_closure_%d", iPhiBins ) );
    JPsiPeakValue    = (Double_t) fCohJpsiToMu->GetEntries();
    JPsiPeakValueErr = TMath::Sqrt((Double_t) fCohJpsiToMu->GetEntries());

    // RawYieldsSimpleClosure->Fill(       iPhiBins,   JPsiPeakValue   /((PhiCenters[7]*2.)*(1.2)));
    // RawYieldsSimpleClosure->SetBinError(iPhiBins+1, JPsiPeakValueErr/((PhiCenters[7]*2.)*(1.2)));
    RawYieldsSimpleClosure->Fill(       iPhiBins,   JPsiPeakValue   );
    RawYieldsSimpleClosure->SetBinError(iPhiBins+1, JPsiPeakValueErr);

  }

  TFile* SimpleClosure = new TFile("SignalExtraction/SimpleClosure.root", "RECREATE");
  RawYieldsSimpleClosure->Write();
  SimpleClosure->Close();


}
