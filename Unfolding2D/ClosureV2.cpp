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
#include "TMatrixD.h"
// #include "RooUnfoldResponse.h"
// #include <RooUnfold.h>



void UnfoldHe2D(Int_t Iterations = 1){

  #ifdef __CINT__
  // if (!TClass::GetDict("RooUnfold")) gSystem->Load("../RooUnfold/libRooUnfold");
    gSystem->Load("../RooUnfold/libRooUnfold");
  #endif

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
  Double_t CosThetaHEgen3;
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
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      24, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      24, 0., 2.*TMath::Pi());
  TH1F *PhiGenH2      = new TH1F("PhiGenH2",      "PhiGenH2",      24, 0., 2.*TMath::Pi());
  TH1F *PhiGenBinsH[24];
  TH1F *PhiRecBinsH[24];

  Int_t N = 1;
  RooUnfoldResponse responsePhi[24];
  for (Int_t i = 0; i < 24; i++) {
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
    responsePhi[i] = RooUnfoldResponse(N,  0., 2.*TMath::Pi());
    PhiGenBinsH[i] = new TH1F(Form("PhiGenBinsH_%d", i), Form("PhiGenBinsH_%d", i), N, 0., 2.*TMath::Pi());
    PhiRecBinsH[i] = new TH1F(Form("PhiRecBinsH_%d", i), Form("PhiRecBinsH_%d", i), N, 0., 2.*TMath::Pi());
  }








  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);



      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {
          PhiGenH     ->Fill( (PhiHEgen+TMath::Pi()) );
          CosThetaGenH->Fill( CosThetaHEgen          );
          CosThetaHEgen2 = -1.*CosThetaHEgen;
          PhiHEgen2 = PhiHEgen;

          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {


              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ) {
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  PhiRecH->Fill( PhiHErec );


                  for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                    if(  ( CosThetaHEgen2 > (-1. + ( (Double_t) iCounter)*(0.08+0.01/3.)) &&
                           CosThetaHEgen2 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                    {
                      responsePhi[iCounter].Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                      PhiGenBinsH[iCounter]->Fill( (PhiHEgen2+TMath::Pi()) );
                      PhiRecBinsH[iCounter]->Fill( (PhiHErec) );
                    }

                  }


              }

          } else {


                  for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                    if(  ( CosThetaHEgen2 > (-1. + ((Double_t) iCounter)*(0.08+0.01/3.)) &&
                           CosThetaHEgen2 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                    {
                      responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                      PhiGenBinsH[iCounter]->Fill( (PhiHEgen2+TMath::Pi()) );

                    }

                  }


          }

      }

  }





  TFile *SavingFile[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    SavingFile[iC]  = new TFile(Form("RomanExercise2_%d.root", iC), "RECREATE");
    PhiRecH   ->Write();
    PhiGenH   ->Write();
    PhiGenBinsH[iC]->Write();
    PhiRecBinsH[iC]->Write();
    responsePhi[iC].Write();
    PhiGenBinsH[iC]->Write();
    SavingFile[iC]->Close();
  }

}
