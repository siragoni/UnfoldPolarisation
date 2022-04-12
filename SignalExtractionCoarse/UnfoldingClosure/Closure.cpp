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



void UnfoldHe2D(Int_t TriggerFlag = 0, Int_t Iterations = 1){

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
  // TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      25, 0., 2.*TMath::Pi());
  // TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      25, 0., 2.*TMath::Pi());
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      24, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      24, 0., 2.*TMath::Pi());
  TH1F *PhiGenH2      = new TH1F("PhiGenH2",      "PhiGenH2",      24, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());
  TH1F *PhiGenBinsH[24];
  TH1F *PhiRecBinsH[24];

  Int_t N = 6;
  RooUnfoldResponse responsePhi[24];
  for (Int_t i = 0; i < 24; i++) {
    responsePhi[i] = RooUnfoldResponse(N,  0., 2.*TMath::Pi());
    PhiGenBinsH[i] = new TH1F(Form("PhiGenBinsH_%d", i), Form("PhiGenBinsH_%d", i), N, 0., 2.*TMath::Pi());
    PhiRecBinsH[i] = new TH1F(Form("PhiRecBinsH_%d", i), Form("PhiRecBinsH_%d", i), N, 0., 2.*TMath::Pi());
  }



  RooUnfoldResponse responsePhiSimple = RooUnfoldResponse(24,  0., 2.*TMath::Pi());



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
    CosThetaHEgen3 = -1.*CosThetaHEgen;
    // for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
    //   if(  ( CosThetaHEgen3 > (-1. + ( (Double_t) iCounter)*(0.08+0.01/3.)) &&
    //          CosThetaHEgen3 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
    //   {
    //     PhiGenBinsH[iCounter]->Fill( (PhiHEgen+TMath::Pi()) );
    //   }
    //
    // }



      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) && (TriggerSelectionFlag > 0.5) ) {
          PhiGenH     ->Fill( (PhiHEgen+TMath::Pi()) );
          CosThetaGenH->Fill( CosThetaHEgen          );
          CosThetaHEgen2 = -1.*CosThetaHEgen;
          PhiHEgen2 = PhiHEgen;

          Costhetahelp = CosThetaHEgen2;
          Costhetahelp2 = CosThetaHErec - CosThetaHEgen2;
          phihelp2 = PhiHErec - PhiHEgen2;
          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {


              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ) {
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  // responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  PhiRecH->Fill( PhiHErec );

                  if(  ( CosThetaHEgen2 > (-0.6) && CosThetaHEgen2 < (0.6) ) ) {
                    responsePhiSimple.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                    PhiGenH2     ->Fill( (PhiHEgen2+TMath::Pi()) );

                  }

                  for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                    if(  ( CosThetaHEgen2 > (-1. + ( (Double_t) iCounter)*(0.08+0.01/3.)) &&
                           CosThetaHEgen2 < (-1. + (((Double_t) iCounter)+1.)*(0.08+0.01/3.) )) )
                    {
                      responsePhi[iCounter].Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                      PhiGenBinsH[iCounter]->Fill( (PhiHEgen2+TMath::Pi()) );
                      PhiRecBinsH[iCounter]->Fill( (PhiHErec) );
                    // } else {
                    //   responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                    }

                  }


              }

          } else {
                  // for (Int_t iCounter = 0; iCounter < 24; iCounter++) {
                  //   responsePhi[iCounter].Miss((PhiHEgen2+TMath::Pi()));
                  // }
                  if(  ( CosThetaHEgen2 > (-0.6) && CosThetaHEgen2 < (0.6) ) ) {
                    responsePhiSimple.Miss((PhiHEgen2+TMath::Pi()));
                    PhiGenH2     ->Fill( (PhiHEgen2+TMath::Pi()) );

                  }

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
    fileDataRaw[iC] = new TFile(Form("SignalExtractionCoarse/ClosureYields/MonteCarloYieldsHe_%d.root", iC));
    // fileDataRaw[iC] = new TFile(Form("SignalExtraction/RawYieldsHe_%d.root", iC));
    Data[iC]        = (TH1F*)fileDataRaw[iC]->Get(Form("h_%d", iC));
  }


  TFile* fileDataRawSimple = new TFile("SignalExtraction/SimpleClosure.root");
  TH1F *DataSimple = (TH1F*)fileDataRawSimple->Get("h");



  TH1F* histo[24];
  for (Int_t i = 4; i < 20; i++) {
    histo[i] = new TH1F(Form("histo_%d", i), Form("histo_%d", i), 6, 0., 2.*TMath::Pi());
    for(Int_t j = 1; j < 30; j++){
      histo[i]->SetBinContent(j, Data[i]->GetBinContent(j));
      histo[i]->SetBinError(  j, Data[i]->GetBinError(j)  );
    }
  }


  TH1F* histoSimple = new TH1F("histoSimple", "histoSimple", 24, 0., 2.*TMath::Pi());
  for(Int_t j = 1; j < 30; j++){
    histoSimple->SetBinContent(j, DataSimple->GetBinContent(j));
    histoSimple->SetBinError(  j, DataSimple->GetBinError(j)  );
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
    unfold[iC]  = RooUnfoldBayes(&responsePhi[iC], histo[iC], Iterations);
    unfold1[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 4);
    unfold2[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 10);
    unfold3[iC] = RooUnfoldBayes(&responsePhi[iC], histo[iC], 20);
    // unfold[iC]  = RooUnfoldBayes(&responsePhi[iC], PhiRecBinsH[iC], Iterations);
    // unfold1[iC] = RooUnfoldBayes(&responsePhi[iC], PhiRecBinsH[iC], 4);
    // unfold2[iC] = RooUnfoldBayes(&responsePhi[iC], PhiRecBinsH[iC], 10);
    // unfold3[iC] = RooUnfoldBayes(&responsePhi[iC], PhiRecBinsH[iC], 20);
  }



  RooUnfoldBayes    unfoldSimple;
  RooUnfoldBayes    unfoldSimple1;
  RooUnfoldBayes    unfoldSimple2;
  RooUnfoldBayes    unfoldSimple3;
  unfoldSimple  = RooUnfoldBayes(&responsePhiSimple, histoSimple, Iterations);
  unfoldSimple1 = RooUnfoldBayes(&responsePhiSimple, histoSimple, 4);
  unfoldSimple2 = RooUnfoldBayes(&responsePhiSimple, histoSimple, 10);
  unfoldSimple3 = RooUnfoldBayes(&responsePhiSimple, histoSimple, 20);





  // TH1D* unfolded  = (TH1D*) unfold.Hreco();
  // TH1D* unfolded1 = (TH1D*) unfold1.Hreco();
  // TH1D* unfolded2 = (TH1D*) unfold2.Hreco();
  // TH1D* unfolded3 = (TH1D*) unfold3.Hreco();
  TH1D* unfolded[24];
  TH1D* unfolded1[24];
  TH1D* unfolded2[24];
  TH1D* unfolded3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    // unfolded[iC]  = (TH1D*) unfold[iC].Hreco(RooUnfold::kCovariance);
    // unfolded1[iC] = (TH1D*) unfold1[iC].Hreco(RooUnfold::kCovariance);
    // unfolded2[iC] = (TH1D*) unfold2[iC].Hreco(RooUnfold::kCovariance);
    // unfolded3[iC] = (TH1D*) unfold3[iC].Hreco(RooUnfold::kCovariance);
    unfolded[iC]  = (TH1D*) unfold[iC].Hreco();
    unfolded1[iC] = (TH1D*) unfold1[iC].Hreco();
    unfolded2[iC] = (TH1D*) unfold2[iC].Hreco();
    unfolded3[iC] = (TH1D*) unfold3[iC].Hreco();
  }

  TH1D* refolded[24];
  TH1D* refolded1[24];
  TH1D* refolded2[24];
  TH1D* refolded3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    refolded[iC]  = (TH1D*) responsePhi[iC].ApplyToTruth(unfolded[iC]);
    refolded1[iC] = (TH1D*) responsePhi[iC].ApplyToTruth(unfolded1[iC]);
    refolded2[iC] = (TH1D*) responsePhi[iC].ApplyToTruth(unfolded2[iC]);
    refolded3[iC] = (TH1D*) responsePhi[iC].ApplyToTruth(unfolded3[iC]);
  }



  TH1D* unfoldedSimple  = (TH1D*) unfoldSimple.Hreco();
  TH1D* unfoldedSimple1 = (TH1D*) unfoldSimple1.Hreco();
  TH1D* unfoldedSimple2 = (TH1D*) unfoldSimple2.Hreco();
  TH1D* unfoldedSimple3 = (TH1D*) unfoldSimple3.Hreco();




  TH1D* refoldedSimple  = (TH1D*) responsePhiSimple.ApplyToTruth(unfoldedSimple);
  TH1D* refoldedSimple1 = (TH1D*) responsePhiSimple.ApplyToTruth(unfoldedSimple1);
  TH1D* refoldedSimple2 = (TH1D*) responsePhiSimple.ApplyToTruth(unfoldedSimple2);
  TH1D* refoldedSimple3 = (TH1D*) responsePhiSimple.ApplyToTruth(unfoldedSimple3);



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

  TH1F* histoSimple2 = new TH1F("histoSimple2", "histoSimple2", 100, -0.5, 99.5);
    for(Int_t i = 1; i < 30; i++){
      histoSimple2->SetBinContent(i, unfoldedSimple->GetBinContent(i));
      histoSimple2->SetBinError(  i, unfoldedSimple->GetBinError(i)  );
    }




  TH1F* histo3[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    histo3[iC] = new TH1F(Form("histo3_%d", iC), Form("histo3_%d", iC), 100, -0.5, 99.5);

    for(Int_t i = 1; i < 30; i++){
      histo3[iC]->SetBinContent(i, PhiGenBinsH[iC]->GetBinContent(i)*(6.0));
      histo3[iC]->SetBinError(  i, PhiGenBinsH[iC]->GetBinError(i)*(6.0)  );
    }
  }






  TFile *SavingFileSimpleUnfold;
    SavingFileSimpleUnfold  = new TFile("Unfolding2D/UnfoldedClosureHeSimple.root", "RECREATE");
    responsePhiSimple.Write();
    unfoldedSimple->Write();
    unfoldedSimple1->Write();
    unfoldedSimple2->Write();
    unfoldedSimple3->Write();
    SavingFileSimpleUnfold->Close();



    new TCanvas;
    unfoldedSimple->SetLineColor(kRed);
    unfoldedSimple1->SetLineColor(kBlue);
    unfoldedSimple2->SetLineColor(kYellow);
    unfoldedSimple3->SetLineColor(kGreen);
    unfoldedSimple->Draw();
    unfoldedSimple1->Draw("same");
    unfoldedSimple2->Draw("same");
    unfoldedSimple3->Draw("same");
    PhiGenH2->SetLineColor(kBlack);
    PhiGenH2->SetLineWidth(4);
    PhiGenH2->Draw("same");
    new TCanvas;
    PhiGenH2->SetLineColor(kBlack);
    PhiGenH2->SetLineWidth(4);
    PhiGenH2->Draw("same");




    new TCanvas;
    refoldedSimple->Divide(histoSimple);
    refoldedSimple1->Divide(histoSimple);
    refoldedSimple2->Divide(histoSimple);
    refoldedSimple3->Divide(histoSimple);
    refoldedSimple->SetLineColor(kRed);
    refoldedSimple1->SetLineColor(kBlue);
    refoldedSimple2->SetLineColor(kYellow);
    refoldedSimple3->SetLineColor(kGreen);
    histoSimple->SetLineColor(kBlack);
    histoSimple->SetLineWidth(3);
    refoldedSimple->Draw("same");
    refoldedSimple1->Draw("same");
    refoldedSimple2->Draw("same");
    refoldedSimple3->Draw("same");
    histoSimple->Draw("same");





    TLatex* latex = new TLatex();
    latex->SetTextSize(0.055);
    latex->SetTextFont(42);
    latex->SetTextAlign(11);
    latex->SetNDC();
    latex->DrawLatex(0.12,0.94,"This thesis");




    TLegend *leg_pt[24];
    for (Int_t iC = 4; iC < 20; iC++) {
    new TCanvas;
    gPad->SetMargin(0.13,0.10,0.12,0.10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    unfolded[iC]->SetLineColor(kRed);
    unfolded[iC]->SetTitle("Generated MC vs unfolded reconstructed MC");
    // TLatex* latex = new TLatex();
    // latex->SetTextSize(0.055);
    // latex->SetTextFont(42);
    // latex->SetTextAlign(11);
    // latex->SetNDC();
    // latex->DrawLatex(0.12,0.94,"This thesis");
    // latex->SetTextSize(0.055);
    // latex->DrawLatex(0.7,0.94,"Helicity");

    unfolded[iC]->GetXaxis()->SetTitleOffset(1.15);
    unfolded[iC]->GetYaxis()->SetTitleOffset(1.25);
    unfolded[iC]->GetXaxis()->SetTitleSize(0.045);
    unfolded[iC]->GetYaxis()->SetTitleSize(0.045);
    unfolded[iC]->GetXaxis()->SetLabelSize(0.045);
    unfolded[iC]->GetYaxis()->SetLabelSize(0.045);
    unfolded[iC]->GetXaxis()->SetTitleFont(42);
    unfolded[iC]->GetYaxis()->SetTitleFont(42);
    unfolded[iC]->GetXaxis()->SetLabelFont(42);
    unfolded[iC]->GetYaxis()->SetLabelFont(42);

    unfolded[iC]->GetXaxis()->SetTitle("#varphi");
    unfolded[iC]->GetYaxis()->SetTitle("Unfolded MC");
    unfolded[iC]->GetYaxis()->SetRangeUser(unfolded[iC]->GetMinimum()*0.9, unfolded[iC]->GetMaximum()*1.3);
    // unfolded[iC]->SetLineWidth(5);

    unfolded1[iC]->SetLineColor(kBlue);
    unfolded2[iC]->SetLineColor(kYellow);
    unfolded3[iC]->SetLineColor(kGreen);
    unfolded[iC]->Draw();
    unfolded1[iC]->Draw("same");
    unfolded2[iC]->Draw("same");
    unfolded3[iC]->Draw("same");
    PhiGenBinsH[iC]->SetLineColor(kBlack);
    PhiGenBinsH[iC]->SetLineWidth(4);
    PhiGenBinsH[iC]->Draw("same");
    // TLatex* latex = new TLatex();
    // latex->SetTextSize(0.055);
    // latex->SetTextFont(42);
    // latex->SetTextAlign(11);
    // latex->SetNDC();
    latex->DrawLatex(0.2,0.74,"This thesis");
    leg_pt[iC] = new TLegend(0.5,0.65,0.85,0.79);
    leg_pt[iC]->SetFillStyle(0);
    leg_pt[iC]->SetBorderSize(0);
    leg_pt[iC]->SetTextSize(0.042);
    leg_pt[iC]->AddEntry(unfolded[iC],"1 iteration", "LP");
    leg_pt[iC]->AddEntry(unfolded1[iC],"4 iteration", "LP");
    leg_pt[iC]->AddEntry(unfolded2[iC],"10 iteration", "LP");
    leg_pt[iC]->AddEntry(unfolded3[iC],"20 iteration", "LP");
    leg_pt[iC]->Draw();
    gPad->SaveAs(Form("Unfolding2D/closureplots/unfolded_%d.pdf", iC), "recreate");


    }


    TH1F* histo6 = new TH1F("histo6", "histo6", 1000, -0.5, 999.5);
    Int_t Counter  = 1;
    Int_t Counter5 = 1;
    TH1F* histo9 = new TH1F("histo9", "histo9", 1000, -0.5, 999.5);
    TLegend *leg_pt2[24];
    for (Int_t iC = 5; iC < 19; iC++) {
    new TCanvas;
    gPad->SetMargin(0.13,0.10,0.12,0.10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gStyle->SetOptStat(0);

    refolded[iC]->SetLineColor(kRed);
    refolded[iC]->SetTitle("Ratio of refolded MC to reconstructed MC");
    refolded[iC]->GetXaxis()->SetTitleOffset(1.15);
    refolded[iC]->GetYaxis()->SetTitleOffset(1.35);
    refolded[iC]->GetXaxis()->SetTitleSize(0.045);
    refolded[iC]->GetYaxis()->SetTitleSize(0.045);
    refolded[iC]->GetXaxis()->SetLabelSize(0.045);
    refolded[iC]->GetYaxis()->SetLabelSize(0.045);
    refolded[iC]->GetXaxis()->SetTitleFont(42);
    refolded[iC]->GetYaxis()->SetTitleFont(42);
    refolded[iC]->GetXaxis()->SetLabelFont(42);
    refolded[iC]->GetYaxis()->SetLabelFont(42);

    refolded[iC]->GetXaxis()->SetTitle("#varphi");
    refolded[iC]->GetYaxis()->SetTitle("Refolded ratio to reconstructed MC");
    refolded[iC]->GetYaxis()->SetRangeUser(0.975, 1.025);
    // unfolded[iC]->SetLineWidth(5);
    for (Int_t j = 1; j < 30; j++) {
      if(refolded[iC]->GetBinContent(j) > 0.000000001){
        Double_t avalue = refolded[iC]->GetBinContent(j);
        Double_t bvalue = histo[iC]   ->GetBinContent(j);
        Double_t cvalue = unfolded[iC]->GetBinContent(j);
        Double_t aerror = refolded[iC]->GetBinError(j);
        Double_t berror = histo[iC]   ->GetBinError(j);
        Double_t cerror = unfolded[iC]->GetBinError(j);
        cout << "====================" << j << endl;
        cout << "avalue = " << avalue << endl;
        cout << "bvalue = " << bvalue << endl;
        cout << "cvalue = " << cvalue << endl;
        cout << "aerror = " << aerror << endl;
        cout << "berror = " << berror << endl;
        cout << "cerror = " << cerror << endl;
        // if (aerror == 0){
        // if (aerror == 0 || (aerror < (cerror * avalue/cvalue))){
          // if (berror == 0) cout << "PROBLEM2" << endl;
          // if (cerror == 0) cout << "PROBLEM3" << endl;
          aerror = cerror * avalue/cvalue;
        // }
        // histo9->SetBinContent(Counter5, (refolded[iC]->GetBinContent(j) - histo[iC]->GetBinContent(j))*(refolded[iC]->GetBinContent(j)-histo[iC]->GetBinContent(j))/(refolded[iC]->GetBinError(j)));
        histo9->SetBinContent(Counter5, (avalue - bvalue)*(avalue - bvalue)/(aerror));
        Counter5++;
        cout << "AAAAAA    " << j << endl;
      }
    }
    for (Int_t j = 1; j < 7; j++) {
      if(refolded[iC]->GetBinContent(j) > 0.000000001){
        // histo6->SetBinContent(Counter, (refolded[iC]->GetBinContent(j)-1.)*(refolded[iC]->GetBinContent(j)-1.)*(histo[iC]->GetBinContent(j)));
        histo6->SetBinContent(Counter, (refolded[iC]->GetBinContent(j)-histo[iC]->GetBinContent(j))*(refolded[iC]->GetBinContent(j)-histo[iC]->GetBinContent(j))/(histo[iC]->GetBinContent(j)));
        Counter++;
      }
    }

    refolded1[iC]->SetLineColor(kBlue);
    refolded2[iC]->SetLineColor(kYellow);
    refolded3[iC]->SetLineColor(kGreen);
    refolded[iC]->Divide(histo[iC]);
    refolded1[iC]->Divide(histo[iC]);
    refolded2[iC]->Divide(histo[iC]);
    refolded3[iC]->Divide(histo[iC]);
    refolded[iC]->Draw();
    refolded1[iC]->Draw("same");
    refolded2[iC]->Draw("same");
    refolded3[iC]->Draw("same");
    // TLatex* latex = new TLatex();
    // latex->SetTextSize(0.055);
    // latex->SetTextFont(42);
    // latex->SetTextAlign(11);
    // latex->SetNDC();
    latex->DrawLatex(0.2,0.74,"This thesis");


    leg_pt2[iC] = new TLegend(0.5,0.65,0.85,0.79);
    leg_pt2[iC]->SetFillStyle(0);
    leg_pt2[iC]->SetBorderSize(0);
    leg_pt2[iC]->SetTextSize(0.042);
    leg_pt2[iC]->AddEntry(refolded[iC],"1 iteration", "LP");
    leg_pt2[iC]->AddEntry(refolded1[iC],"4 iteration", "LP");
    leg_pt2[iC]->AddEntry(refolded2[iC],"10 iteration", "LP");
    leg_pt2[iC]->AddEntry(refolded3[iC],"20 iteration", "LP");
    leg_pt2[iC]->Draw();
    gPad->SaveAs(Form("Unfolding2D/closureplots/refolded_%d.pdf", iC), "recreate");

    }


    Double_t M = 6.;
    TH1F* histo4[24];
    for (Int_t iC = 4; iC < 20; iC++) {
      histo4[iC] = new TH1F(Form("histo4_%d", iC), Form("histo4_%d", iC), 100, -0.5, 99.5);
      for(Int_t i = 1; i < 30; i++){
        histo4[iC]->SetBinContent(i, histo2[iC]->GetBinContent(i)*((Double_t) M));
        histo4[iC]->SetBinError(  i, histo2[iC]->GetBinError(i)*((Double_t) M));
      }
    }



    TH1F* histo5 = new TH1F("histo5", "histo5", 1000, -0.5, 999.5);
    for (Int_t iC = 5; iC < 19; iC++) {
      for (Int_t j = 0; j < 6; j++) {
        histo5->SetBinContent( iC*30+j, refolded[iC]->GetBinContent( j+1 ) );
        histo5->SetBinError(   iC*30+j, refolded[iC]->GetBinError( j+1 ) );
      }

    }




    TFile *SavingFile[24];
    for (Int_t iC = 4; iC < 20; iC++) {
      SavingFile[iC]  = new TFile(Form("SignalExtractionCoarse/UnfoldingClosure/UnfoldedClosureHe_%d.root", iC), "RECREATE");
      PhiRecH   ->Write();
      PhiGenH   ->Write();
      CosThetaGenH   ->Write();
      PhiGenBinsH[iC]->Write();
      PhiRecBinsH[iC]->Write();
      responsePhi[iC].Write();
      unfold[iC].Write();
      unfolded[iC] ->Write();
      unfolded1[iC]->Write();
      unfolded2[iC]->Write();
      unfolded3[iC]->Write();
      histo[iC]->Write();
      histo2[iC]->Write();
      histo3[iC]->Write();
      histo4[iC]->Write();
      histo5->Write();
      histo6->Write();
      histo9->Write();
      PhiGenBinsH[iC]->Write();
      SavingFile[iC]->Close();
    }





}
