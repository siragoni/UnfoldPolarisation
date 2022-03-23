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


void ParseMC(){

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



  TH1F*** InvMassH;
  InvMassH = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassH[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j] = new TH1F(Form("InvMassH_%d_%d", i, j),Form("InvMassH_%d_%d", i, j),2000, 0, 20);
    }
  }

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);

      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {

        CosThetaHEgen2 = -1.*CosThetaHEgen;
        if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
        PhiHEgen2 = PhiHEgen;


          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {

                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();


                  //================================
                  // Translating Phi and CosTheta
                  // into bin numbers
                  //--------------------------------
                  Double_t TraslatedCosThetaGen = 0.5*(CosThetaHEgen2 + 1.)*24.;
                  Double_t iCosThetaBins2 = -1;
                  Double_t RemainderCosTheta = 0.;
                  RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
                  Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;
                  //--------------------------------
                  // Binning in phi depends
                  // on CosTheta
                  //--------------------------------
                  Double_t M = 6.;
                  // Double_t TraslatedPhiGen = (PhiHEgen2+TMath::Pi())*M/(2.*TMath::Pi());
                  Double_t TraslatedPhiGen = (PhiHErec)*M/(2.*TMath::Pi());
                  Double_t iPhiBins2 = -1;
                  Double_t RemainderPhi = 0.;
                  RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
                  Int_t iPhiBins = (Int_t)  iPhiBins2;
                  //--------------------------------------
                  InvMassH[iCosThetaBins][iPhiBins]->Fill(Mrec);

          }
      }
    }
  }




  TFile *SavingFile = new TFile("SignalExtractionCoarse/CohHe.root", "RECREATE");
  for (Int_t i = 0; i < 24; i++) {
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j]->Write();
    }
  }
  SavingFile->Close();
















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
  TFile* parsedMC = new TFile("SignalExtraction/CohHeSmart.root");

  Double_t JPsiPeakValue    = 0;
  Double_t JPsiPeakValueErr = 0;
  Double_t BkgValue         = 0;
  Double_t BkgValueError    = 0;
  TH1F* fCohJpsiToMu        = 0x0;
  for (Int_t iCosThetaBins5 = 4; iCosThetaBins5 < 20; iCosThetaBins5++) {

    for (Int_t iPhiBins5 = 0; iPhiBins5 < Counter[iCosThetaBins5]; iPhiBins5++) {

      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;

      fCohJpsiToMu     = (TH1F*)parsedMC->Get( Form( "InvMassH_%d_%d", iCosThetaBins5, iPhiBins5 ) );
      JPsiPeakValue    = (Double_t) fCohJpsiToMu->GetEntries();
      JPsiPeakValueErr = TMath::Sqrt((Double_t) fCohJpsiToMu->GetEntries());

      // RawYields[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue   /((PhiCenters[iCosThetaBins]*2.)*(0.08+0.01/3.)));
      // RawYields[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr/((PhiCenters[iCosThetaBins]*2.)*(0.08+0.01/3.)));
      RawYields[iCosThetaBins5]->Fill(       iPhiBins5,   JPsiPeakValue   );
      RawYields[iCosThetaBins5]->SetBinError(iPhiBins5+1, JPsiPeakValueErr);

    }
    TFile f(Form("SignalExtraction/MonteCarloYieldsHe_%d.root", iCosThetaBins5),   "recreate");
    RawYields[iCosThetaBins5]->Write();
    f.Close();

  }





}
