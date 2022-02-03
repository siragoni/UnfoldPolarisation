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

  TFile* file   = new TFile("AnalysisResultsLHC18qr15o_20210117.root", "READ");
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
                // Form("InvMassH_binthree_%d", iCosThetaBins),
                // Form("InvMassH_binthree_%d", iCosThetaBins),
                Form("InvMassH_3_%d", iCosThetaBins),
                Form("InvMassH_3_%d", iCosThetaBins),
                2000, 0, 20
                );
  }

  TH1F* InvMassH_binfour[1];  // -1 + 5*(0.08+0.01/3.)
  for(Int_t iCosThetaBins = 0; iCosThetaBins < 1; iCosThetaBins++ ){
    InvMassH_binfour[iCosThetaBins] = new TH1F(
                // Form("InvMassH_binfour_%d", iCosThetaBins),
                // Form("InvMassH_binfour_%d", iCosThetaBins),
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



  Double_t controlFlag  = 0;
  Double_t controlFlag2 = 0;

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    // GenerTree->GetEntry(i);

      // if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {
      //
      //   CosThetaHEgen2 = -1.*CosThetaHEgen;
      //   if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
      //   PhiHEgen2 = PhiHEgen;


          if( (Mrec > 2.85 && Mrec < 3.35) && (Yrec > -4. && Yrec < -2.5) && (CosThetaHErec > -0.2 && CosThetaHErec < -0.12) ) {
            PtRecH->Fill(pTrec);
          }

          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {

                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();

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
                  if (        CosThetaHErec < (-1. + 1.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binzero[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 2.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binone[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 3.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwo[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 4.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binthree[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 5.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binfour[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 6.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binfive[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 7.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binsix[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 8.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binseven[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 9.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_bineight[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 10.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binnine[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 11.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binten[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 12.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_bineleven[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 13.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_bintwelve[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 14.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binthirteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 15.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 24; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/24. ) {
                        InvMassH_binfourteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 16.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binfifteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 17.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 12; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/12. ) {
                        InvMassH_binsixteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 18.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_binseventeen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 19.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 6; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/6. ) {
                        InvMassH_bineighteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 20.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_binnineteen[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 21.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwenty[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 22.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentyone[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 23.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentytwo[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  } else if ( CosThetaHErec < (-1. + 24.*(0.08+0.01/3.)) ){
                    for(Int_t i = 0; i < 1; i++) {
                      if( controlFlag == 1) break;
                      if( (PhiHErec) < ((Double_t)i + 1.)*2.*TMath::Pi()/1. ) {
                        InvMassH_bintwentythree[i]->Fill(Mrec);
                        controlFlag = 1;
                      }
                    }
                  }




      //     }
      // }
    }
  }













  TFile *SavingFile = new TFile("SignalExtraction/CohHe_data.root", "RECREATE");
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
  PtRecH->Write();

  SavingFile       ->Close();

}
