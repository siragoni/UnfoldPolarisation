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



void UnfoldHe2D(Int_t TriggerFlag = 0){

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
  TH1F *PhiGenBinsH[24];

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
  }


  Double_t Costhetahelp = -999.;
  Double_t phihelp      = -999.;
  Double_t Costhetahelp2 = -999.;
  Double_t phihelp2      = -999.;
  Double_t TriggerSelectionFlag = 1;


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
          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {


                // if(  ( CosThetaHErec > -1.                      && CosThetaHErec < (-1. + 1.*(0.08+0.01/3.) ) ) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 1.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 2.*(0.08+0.01/3.) )) ){continue;}
                // if(  ( CosThetaHErec > (-1. + 2.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 3.*(0.08+0.01/3.) )) ){continue;}
                // if(  ( CosThetaHErec > (-1. + 3.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 4.*(0.08+0.01/3.) )) ){continue;}
                // if(  ( CosThetaHErec > (-1. + 4.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 5.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 5.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 6.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 6.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 7.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 7.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 8.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 8.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 9.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 9.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 10.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 10.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 11.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 11.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 12.*(0.08+0.01/3.) )) ){continue;}
                // if(  ( CosThetaHErec > (-1. + 12.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 13.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 13.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 14.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 14.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 15.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 15.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 16.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 16.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 17.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 17.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 18.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 18.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 19.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 19.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 20.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 20.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 21.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 21.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 22.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 22.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 23.*(0.08+0.01/3.) )) ) {continue;}
                // if(  ( CosThetaHErec > (-1. + 23.*(0.08+0.01/3.)) && CosThetaHErec < (-1. + 24.*(0.08+0.01/3.) )) ) {continue;}

              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  // responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  PhiRecH->Fill( PhiHErec );


                  // for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                  //   if(  ( CosThetaHErec > (-1. + ((Double_t) iCounter)*(0.08+0.01/3.)) && CosThetaHErec < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                  //   {
                  //     responsePhi[iCounter].Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  //   } else {
                  //     responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                  //   }
                  //
                  // }
                  for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                    if(  ( CosThetaHEgen2 > (-1. + ((Double_t) iCounter)*(0.08+0.01/3.)) && CosThetaHEgen2 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                    {
                      responsePhi[iCounter].Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                      PhiGenBinsH[iCounter]->Fill( (PhiHEgen2+TMath::Pi()) );
                    // } else {
                    //   responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                    }

                  }







              }

          } else {
                  // for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                  //   responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                  // }
                  for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                    if(  ( CosThetaHEgen2 > (-1. + ((Double_t) iCounter)*(0.08+0.01/3.)) && CosThetaHEgen2 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                    {
                      responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                      PhiGenBinsH[iCounter]->Fill( (PhiHEgen2+TMath::Pi()) );

                    }

                  }


          }

      }

  }




  // TFile* fileDataRaw = 0x0;
  // if(        TriggerFlag == 1 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 2 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 3 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 4 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 5 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 6 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else if( TriggerFlag == 7 ) {
  //   fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d_%d.root", CosThetaFlag, TriggerFlag));
  // } else {
    // fileDataRaw = new TFile(Form("SignalExtraction/RawYieldsHe_%d.root", CosThetaFlag));
  // }
  // TH1F *Data = (TH1F*)fileDataRaw->Get(Form("h_%d", CosThetaFlag));

  TFile* fileDataRaw[24];
  TH1F *Data[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    fileDataRaw[iC] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_%d.root", iC));
    // fileDataRaw[iC] = new TFile(Form("SignalExtraction/RawYieldsHe_%d.root", iC));
    Data[iC]        = (TH1F*)fileDataRaw[iC]->Get(Form("h_%d", iC));
  }


  TH1F* histo[24];
  for (Int_t i = 4; i < 20; i++) {
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

    histo[i] = new TH1F(Form("histo_%d", i), Form("histo_%d", i), N, 0., 2.*TMath::Pi());
    for(Int_t j = 1; j < 30; j++){
      histo[i]->SetBinContent(j, Data[i]->GetBinContent(j));
      histo[i]->SetBinError(  j, Data[i]->GetBinError(j)  );
    }
  }



  // RooUnfoldBayes    unfold  (&responsePhi, Data, 1);
  // RooUnfoldBayes    unfold1 (&responsePhi, Data, 2);
  // RooUnfoldBayes    unfold2 (&responsePhi, Data, 3);
  // RooUnfoldBayes    unfold3 (&responsePhi, Data, 4);
  // // RooUnfoldSvd      unfold (&response, hist_measured, kterm);
  // // RooUnfoldBinByBin unfold (&response, hist_measured);
  RooUnfoldBayes    unfold[24];
  RooUnfoldBayes    unfold1[24];
  RooUnfoldBayes    unfold2[24];
  RooUnfoldBayes    unfold3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    // unfold[iC]  = RooUnfoldBayes(&responsePhi[iC], Data[iC], 1);
    // unfold1[iC] = RooUnfoldBayes(&responsePhi[iC], Data[iC], 2);
    // unfold2[iC] = RooUnfoldBayes(&responsePhi[iC], Data[iC], 3);
    // unfold3[iC] = RooUnfoldBayes(&responsePhi[iC], Data[iC], 4);
    unfold[iC]  = RooUnfoldBayes(&responsePhi[iC], histo[iC], 1);
    unfold1[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 2);
    unfold2[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 3);
    unfold3[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 4);
  }





  // TH1D* unfolded  = (TH1D*) unfold.Hreco();
  // TH1D* unfolded1 = (TH1D*) unfold1.Hreco();
  // TH1D* unfolded2 = (TH1D*) unfold2.Hreco();
  // TH1D* unfolded3 = (TH1D*) unfold3.Hreco();
  TH1D* unfolded[24];
  TH1D* unfolded1[24];
  TH1D* unfolded2[24];
  TH1D* unfolded3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    unfolded[iC]  = (TH1D*) unfold[iC].Hreco();
    unfolded1[iC] = (TH1D*) unfold1[iC].Hreco();
    unfolded2[iC] = (TH1D*) unfold2[iC].Hreco();
    unfolded3[iC] = (TH1D*) unfold3[iC].Hreco();
  }










  // TH1F* histo = new TH1F("histo", "histo", 100, -0.5, 99.5);
  // for(Int_t i = 1; i < 30; i++){
  //   histo->SetBinContent(i, unfolded->GetBinContent(i));
  //   histo->SetBinError(i, unfolded->GetBinError(i));
  // }
  // TH1F* histo[24];
  // for (Int_t iC = 4; iC < 20; iC++) {
  //   histo[iC] = new TH1F(Form("histo_%d", iC), Form("histo_%d", iC), 100, -0.5, 99.5);
  //   for(Int_t i = 1; i < 30; i++){
  //     histo[iC]->SetBinContent(i, unfolded[iC]->GetBinContent(i));
  //     histo[iC]->SetBinError(  i, unfolded[iC]->GetBinError(i)  );
  //   }
  // }



  TH1F* histo2[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    histo2[iC] = new TH1F(Form("histo2_%d", iC), Form("histo2_%d", iC), 100, -0.5, 99.5);
    for(Int_t i = 1; i < 30; i++){
      histo2[iC]->SetBinContent(i, unfolded[iC]->GetBinContent(i));
      histo2[iC]->SetBinError(  i, unfolded[iC]->GetBinError(i)  );
    }
  }


  TH1F* histo3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    histo3[iC] = new TH1F(Form("histo3_%d", iC), Form("histo3_%d", iC), 100, -0.5, 99.5);
    if( iC == 0  || iC == 1  || iC == 2  || iC == 3  ||
        iC == 4  || iC == 23 || iC == 22 || iC == 21 ||
        iC == 20 || iC == 19 )
    {
      N = 1;
    } else if ( iC == 5  || iC == 6  || iC == 18  || iC == 17 ) {
      N = 6;
    } else if ( iC == 7  || iC == 8  || iC == 16  || iC == 15 ) {
      N = 12;
    } else if ( iC == 9  || iC == 10  || iC == 11  ||
                iC == 12 || iC == 13  || iC == 14 ) {
      N = 24;
    }

    for(Int_t i = 1; i < 30; i++){
      histo3[iC]->SetBinContent(i, PhiGenBinsH[iC]->GetBinContent(i)*((Double_t) N));
      histo3[iC]->SetBinError(  i, PhiGenBinsH[iC]->GetBinError(i)*((Double_t) N)  );
    }
  }



  TFile *SavingFile[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    SavingFile[iC]  = new TFile(Form("Unfolding2D/UnfoldedClosureHe_%d.root", iC), "RECREATE");
    PhiRecH   ->Write();
    PhiGenH   ->Write();
    responsePhi[iC].Write();
    unfolded[iC] ->Write();
    unfolded1[iC]->Write();
    unfolded2[iC]->Write();
    unfolded3[iC]->Write();
    histo[iC]->Write();
    histo2[iC]->Write();
    histo3[iC]->Write();
    PhiGenBinsH[iC]->Write();
    SavingFile[iC]->Close();
  }


}
