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



void AxE(Int_t TriggerFlag = 0){

  // #ifdef __CINT__
  // // if (!TClass::GetDict("RooUnfold")) gSystem->Load("../RooUnfold/libRooUnfold");
  //   gSystem->Load("../RooUnfold/libRooUnfold");
  // #endif

  TFile* file   = new TFile("AnalysisResultsCohWithSinglePt.root", "READ");
  TTree* ReconTree = (TTree*)file->Get("MyTask/fOutputTree");
  TTree* GenerTree = (TTree*)file->Get("MyTask/fOutputTreeMC");


  Double_t pTsingle1;
  Double_t pTsingle2;
  ReconTree->SetBranchAddress("fTrkPt1",  &pTsingle1);
  ReconTree->SetBranchAddress("fTrkPt2",  &pTsingle2);



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
  Double_t TildePhiHEposrec;
  Double_t TildePhiHEposgen;
  Double_t TildePhiHEposrec2;
  Double_t TildePhiHEposgen2;
  ReconTree->SetBranchAddress("fTildePhiHEpos",  &TildePhiHEposrec);
  GenerTree->SetBranchAddress("fMCTildePhiHEpos",&TildePhiHEposgen);
  Double_t TildePhiCSposrec;
  Double_t TildePhiCSposgen;
  ReconTree->SetBranchAddress("fTildePhiCSpos",  &TildePhiCSposrec);
  GenerTree->SetBranchAddress("fMCTildePhiCSpos",&TildePhiCSposgen);
  Double_t TildePhiHEnegrec;
  Double_t TildePhiHEneggen;
  Double_t TildePhiHEnegrec2;
  Double_t TildePhiHEneggen2;
  ReconTree->SetBranchAddress("fTildePhiHEneg",  &TildePhiHEnegrec);
  GenerTree->SetBranchAddress("fMCTildePhiHEneg",&TildePhiHEneggen);
  Double_t TildePhiCSnegrec;
  Double_t TildePhiCSneggen;
  ReconTree->SetBranchAddress("fTildePhiCSneg",  &TildePhiCSnegrec);
  GenerTree->SetBranchAddress("fMCTildePhiCSneg",&TildePhiCSneggen);


  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


  TH1F *PtRecH = new TH1F("PtRecH", "PtRecH", 4000, 0, 20);
  TH1F *PtGenH = new TH1F("PtGenH", "PtGenH", 4000, 0, 20);

  TH1F *PhiRecMinusGenH      = new TH1F("PhiRecMinusGenH",      "PhiRecMinusGenH",      50, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH1F *CosThetaRecMinusGenH = new TH1F("CosThetaRecMinusGenH", "CosThetaRecMinusGenH", 50, -2., 2.);


  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 25, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 25, -1., 1.);
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      25, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      25, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());

  Int_t N = 1;

  Double_t Costhetahelp = -999.;
  Double_t phihelp      = -999.;
  Double_t Costhetahelp2 = -999.;
  Double_t phihelp2      = -999.;
  Double_t TriggerSelectionFlag = 1;





  TH1F* PhiGenBinsCosTheta[24];
  TH1F* PhiRecBinsCosTheta[24];
  for (Int_t CosThetaFlag = 0; CosThetaFlag < 24; CosThetaFlag++) {
    if( CosThetaFlag == 0  || CosThetaFlag == 1  || CosThetaFlag == 2  || CosThetaFlag == 3  ||
        CosThetaFlag == 4  || CosThetaFlag == 23 || CosThetaFlag == 22 || CosThetaFlag == 21 ||
        CosThetaFlag == 20 || CosThetaFlag == 19 )
    {
      N = 1;
    } else if ( CosThetaFlag == 5  || CosThetaFlag == 6  || CosThetaFlag == 18  || CosThetaFlag == 17 ) {
      N = 6;
    } else if ( CosThetaFlag == 7  || CosThetaFlag == 8  || CosThetaFlag == 16  || CosThetaFlag == 15 ) {
      N = 12;
    } else if ( CosThetaFlag == 9  || CosThetaFlag == 10  || CosThetaFlag == 11  ||
                CosThetaFlag == 12 || CosThetaFlag == 13  || CosThetaFlag == 14 ) {
      N = 24;
    }

    PhiGenBinsCosTheta[CosThetaFlag] = new TH1F(Form("PhiGenBinsCosTheta_%d", CosThetaFlag), Form("PhiGenBinsCosTheta_%d", CosThetaFlag), N, 0, 2.*TMath::Pi() );
    PhiRecBinsCosTheta[CosThetaFlag] = new TH1F(Form("PhiRecBinsCosTheta_%d", CosThetaFlag), Form("PhiRecBinsCosTheta_%d", CosThetaFlag), N, 0, 2.*TMath::Pi() );

  }





  TH2F* phithetarecon = new TH2F("phithetarecon","phithetarecon", 100,-1,1, 100, -2.*TMath::Pi(), 2.*TMath::Pi());
  TH2F* phithetagener = new TH2F("phithetagener","phithetagener", 100,-1,1, 100, -2.*TMath::Pi(), 2.*TMath::Pi());

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);
    Costhetahelp = -999.;
    phihelp = -999.;

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




      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) && (TriggerSelectionFlag > 0.5) ) {
          PhiGenH     ->Fill( (PhiHEgen+TMath::Pi()) );
          CosThetaGenH->Fill( CosThetaHEgen          );
          CosThetaHEgen2 = -1.*CosThetaHEgen;
          PhiHEgen2 = PhiHEgen;

          Costhetahelp = CosThetaHEgen2;
          Costhetahelp2 = CosThetaHErec - CosThetaHEgen2;
          phihelp2 = PhiHErec - PhiHEgen2;


          phithetagener->Fill(CosThetaHEgen, PhiHEgen+TMath::Pi());


            if(  ( CosThetaHEgen > -1.                      && CosThetaHEgen < (-1. + 1.*(0.08+0.01/3.) ) ) ) {PhiGenBinsCosTheta[0]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 1.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 2.*(0.08+0.01/3.) )) ) {PhiGenBinsCosTheta[1]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 2.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 3.*(0.08+0.01/3.) )) ) {PhiGenBinsCosTheta[2]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 3.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 4.*(0.08+0.01/3.) )) ) {PhiGenBinsCosTheta[3]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 4.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 5.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[4]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 5.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 6.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[5]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 6.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 7.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[6]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 7.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 8.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[7]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 8.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 9.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[8]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 9.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 10.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[9]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 10.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 11.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[10]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 11.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 12.*(0.08+0.01/3.) )) ) {PhiGenBinsCosTheta[11]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 12.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 13.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[12]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 13.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 14.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[13]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 14.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 15.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[14]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 15.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 16.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[15]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 16.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 17.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[16]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 17.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 18.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[17]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 18.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 19.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[18]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 19.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 20.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[19]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 20.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 21.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[20]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 21.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 22.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[21]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 22.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 23.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[22]->Fill(PhiHEgen+TMath::Pi());}
            if(  ( CosThetaHEgen > (-1. + 23.*(0.08+0.01/3.)) && CosThetaHEgen < (-1. + 24.*(0.08+0.01/3.) )) )  {PhiGenBinsCosTheta[23]->Fill(PhiHEgen+TMath::Pi());}










          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {



              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  // responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  PhiRecH->Fill( PhiHErec );
                  phithetarecon->Fill(CosThetaHErec, PhiHErec);

                  if(  ( CosThetaHErec > -1.                        && CosThetaHErec < (-1. + 1.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[0]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 1.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 2.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[1]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 2.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 3.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[2]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 3.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 4.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[3]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 4.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 5.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[4]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 5.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 6.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[5]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 6.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 7.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[6]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 7.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 8.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[7]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 8.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 9.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[8]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 9.*(0.08+0.01/3.))  && CosThetaHErec < (-1. + 10.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[9]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 10.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 11.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[10]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 11.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 12.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[11]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 12.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 13.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[12]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 13.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 14.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[13]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 14.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 15.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[14]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 15.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 16.*(0.08+0.01/3.) ))){PhiRecBinsCosTheta[15]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 16.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 17.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[16]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 17.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 18.*(0.08+0.01/3.) )) ) {PhiRecBinsCosTheta[17]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 18.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 19.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[18]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 19.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 20.*(0.08+0.01/3.) )) ) {PhiRecBinsCosTheta[19]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 20.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 21.*(0.08+0.01/3.) )) ){PhiRecBinsCosTheta[20]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 21.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 22.*(0.08+0.01/3.) )) ) {PhiRecBinsCosTheta[21]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 22.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 23.*(0.08+0.01/3.) )) ) {PhiRecBinsCosTheta[22]->Fill(PhiHErec);}
                  if(  ( CosThetaHErec > (-1. + 23.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 24.*(0.08+0.01/3.) )) ) {PhiRecBinsCosTheta[23]->Fill(PhiHErec);}







              }

          } else {

          }

      }

  }




  TFile* fileDataRaw[24];
  TH1F *Data[24];
  TH1F *Data2[24];
  for(Int_t i = 4; i < 20; i++){
    fileDataRaw[i] = new TFile(Form("SignalExtraction/RawYieldsHe_%d.root", i));
    Data[i] = (TH1F*)fileDataRaw[i]->Get(Form("h_%d", i));
  }

  TH1F* histo[24];
  for(Int_t i = 4; i < 20; i++){
    histo[i] = new TH1F(Form("histo_%d", i), Form("histo_%d", i), N, 0, 2.*TMath::Pi() );
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 6;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 12;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      N = 24;
    }

    for (Int_t j = 0; j < N; j++) {
      histo[i]->SetBinContent(j+1, Data[i]->GetBinContent(j+1));
      histo[i]->SetBinError(  j+1, Data[i]->GetBinError(j+1));
    }
  }






  TH1F* acceptance[24];
  for(Int_t i = 4; i < 20; i++){
    cout << "i" << i << endl;
    acceptance[i] = (TH1F*) PhiRecBinsCosTheta[i]->Clone(Form("acceptance_%d", i));
    acceptance[i]->Divide(PhiGenBinsCosTheta[i]);
    Data2[i] = (TH1F*) histo[i]->Clone(Form("AxEcorrected_%d", i));
    // Data2[i]->Divide(acceptance[i]);

  }



  for(Int_t i = 4; i < 20; i++){
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 6;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 12;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      N = 24;
    }

    for (Int_t j = 0; j < N; j++) {
      Data2[i]->SetBinContent(j+1, histo[i]->GetBinContent(j+1)/acceptance[i]->GetBinContent(j+1));
      Data2[i]->SetBinError(j+1, histo[i]->GetBinError(j+1)+acceptance[i]->GetBinError(j+1));
    }
  }




  TH2F* acceptance2D = (TH2F*)phithetarecon->Clone("acc2D");
  acceptance2D->Divide(phithetagener);

  TFile* SavingFile[24];
  for(Int_t i = 4; i < 20; i++){
    SavingFile[i] = new TFile(Form("Unfolding2D/AxE_%d.root", i), "RECREATE");
    PhiGenBinsCosTheta[i]->Write();
    PhiRecBinsCosTheta[i]->Write();
    acceptance[i]->Write();
    Data2[i]->Write();
    acceptance2D->Write();
    phithetarecon->Write();
    phithetagener->Write();
    SavingFile[i]->Close();
  }




}
