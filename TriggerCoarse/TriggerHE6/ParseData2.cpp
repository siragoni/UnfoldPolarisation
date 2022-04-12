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


void ParseData(Int_t TriggerFlag = 0){

  TFile* file   = new TFile("AnalysisResultsLHC18qr15o_28032021.root", "READ");
  TTree* ReconTree = (TTree*)file->Get("MyTask/fOutputTree");
  // TTree* GenerTree = (TTree*)file->Get("MyTask/fOutputTreeMC");



  Double_t pTrec;
  // Double_t pTgen;
  ReconTree->SetBranchAddress("fTrkTrkPt",  &pTrec);
  // GenerTree->SetBranchAddress("fMCTrkTrkPt",&pTgen);
  Double_t Mrec;
  // Double_t Mgen;
  ReconTree->SetBranchAddress("fTrkTrkM",  &Mrec);
  // GenerTree->SetBranchAddress("fMCTrkTrkM",&Mgen);
  Double_t Yrec;
  // Double_t Ygen;
  ReconTree->SetBranchAddress("fTrkTrkY",  &Yrec);
  // GenerTree->SetBranchAddress("fMCTrkTrkY",&Ygen);
  Double_t CosThetaHErec;
  // Double_t CosThetaHEgen;
  Double_t CosThetaHErec2;
  // Double_t CosThetaHEgen2;
  ReconTree->SetBranchAddress("fCosThetaHE",  &CosThetaHErec);
  // GenerTree->SetBranchAddress("fMCCosThetaHE",&CosThetaHEgen);
  Double_t CosThetaCSrec;
  // Double_t CosThetaCSgen;
  ReconTree->SetBranchAddress("fCosThetaCS",  &CosThetaCSrec);
  // GenerTree->SetBranchAddress("fMCCosThetaCS",&CosThetaCSgen);
  Double_t PhiHErec;
  // Double_t PhiHEgen;
  Double_t PhiHErec2;
  // Double_t PhiHEgen2;
  ReconTree->SetBranchAddress("fPhiHE",  &PhiHErec);
  // GenerTree->SetBranchAddress("fMCPhiHE",&PhiHEgen);
  Double_t PhiCSrec;
  // Double_t PhiCSgen;
  ReconTree->SetBranchAddress("fPhiCS",  &PhiCSrec);
  // GenerTree->SetBranchAddress("fMCPhiCS",&PhiCSgen);
  Int_t fV0ADecision;
  ReconTree->SetBranchAddress("fV0ADecision",  &fV0ADecision);
  Int_t fV0CDecision;
  ReconTree->SetBranchAddress("fV0CDecision",  &fV0CDecision);
  Int_t fADADecision;
  ReconTree->SetBranchAddress("fADADecision",  &fADADecision);
  Int_t fADCDecision;
  ReconTree->SetBranchAddress("fADCDecision",  &fADCDecision);
  Int_t fV0TotalNCells;
  ReconTree->SetBranchAddress("fV0TotalNCells",  &fV0TotalNCells);
  Double_t pTsingle1;
  Double_t pTsingle2;
  ReconTree->SetBranchAddress("fTrkPt1",  &pTsingle1);
  ReconTree->SetBranchAddress("fTrkPt2",  &pTsingle2);

  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  // Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


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





  TH1F*** InvMassH;
  InvMassH = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassH[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j] = new TH1F(Form("InvMassH_%d_%d", i, j),Form("InvMassH_%d_%d", i, j),2000, 0, 20);
    }
  }


  TH1F* FullInvMassH = new TH1F("FullInvMassH","FullInvMassH",2000, 0, 20);


  Double_t controlFlag  = 0;
  Double_t controlFlag2 = 0;

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);


          if( (Mrec > 2.85 && Mrec < 3.35) && (Yrec > -4. && Yrec < -2.5) && (CosThetaHErec > -0.2 && CosThetaHErec < -0.12) ) {
            PtRecH->Fill(pTrec);
          }
          if( (fV0TotalNCells > 2) ||
            (fADADecision != 0)  ||
            (fADCDecision != 0)  ||
            (fV0ADecision != 0)  ||
            !(fV0CDecision == 0 || fV0CDecision == 1) )
          {
            continue;
          }

          Double_t TriggerSelectionFlag = 1;
          if(        TriggerFlag == 1 ) {
            if( (pTsingle1 > 0.85) && (pTsingle2 > 0.85) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 2 ) {
            if( (pTsingle1 > 0.90) && (pTsingle2 > 0.90) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 3 ) {
            if( (pTsingle1 > 0.95) && (pTsingle2 > 0.95) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 4 ) {
            if( (pTsingle1 > 1.00) && (pTsingle2 > 1.00) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 5 ) {
            if( (pTsingle1 > 1.05) && (pTsingle2 > 1.05) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 6 ) {
            if( (pTsingle1 > 1.10) && (pTsingle2 > 1.10) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else if( TriggerFlag == 7 ) {
            if( (pTsingle1 > 1.15) && (pTsingle2 > 1.15) ){
              TriggerSelectionFlag = 1;
            } else {
              TriggerSelectionFlag = 0;
            }
          } else {
            TriggerSelectionFlag = 1;
          }

          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5)  && (TriggerSelectionFlag > 0.5) ) {




                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();

                  //================================
                  // Translating Phi and CosTheta
                  // into bin numbers
                  //--------------------------------
                  Double_t TraslatedCosThetaGen = 0.5*(CosThetaHErec + 1.)*24.;
                  Double_t iCosThetaBins2 = -1;
                  Double_t RemainderCosTheta = 0.;
                  RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
                  Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;
                  //--------------------------------
                  // Binning in phi depends
                  // on CosTheta
                  //--------------------------------
                  Double_t M = 6.;
                  Double_t TraslatedPhiGen = (PhiHErec)*M/(2.*TMath::Pi());
                  Double_t iPhiBins2 = -1;
                  Double_t RemainderPhi = 0.;
                  RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
                  Int_t iPhiBins = (Int_t)  iPhiBins2;
                  //--------------------------------------
                  InvMassH[iCosThetaBins][iPhiBins]->Fill(Mrec);



                  if(iCosThetaBins > 5 && iCosThetaBins < 18) FullInvMassH->Fill(Mrec);

    }
  }













  TFile *SavingFile = new TFile("TriggerCoarse/TriggerHE6/CohHe_data2.root", "RECREATE");
  for (Int_t i = 0; i < 24; i++) {
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j]->Write();
    }
  }
  FullInvMassH     ->Write();
  SavingFile       ->Close();

}
