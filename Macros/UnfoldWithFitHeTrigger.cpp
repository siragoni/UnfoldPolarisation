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

Double_t CosTheta();
Double_t Phi();
Double_t PhiV2();
Double_t DummyPhi();
Double_t MixWithPositiveCosTheta();
Double_t MixWithNegativeCosTheta();
Double_t TildePhi();
Double_t TildePhiV2();
Double_t SimultaneousFitComplete();
Double_t SimultaneousFit();
Double_t SimultaneousFitLastHopeComplete();
Double_t SimultaneousFitLastHope();
Double_t DummyFitCosTheta();
Double_t DummyFitPhi();
void     FcnForMinimisation();
void     FcnForMinimisationV2();
void     FcnForMinimisationV3();
void     PolarisationHeMinuit1D(TH1* CorrectedCosTheta, TH1* CorrectedPhi, TH1* CorrectedTildePhi, Int_t FitRangeMode = 0, Int_t SignalRangeSelectionMode = 0 );


void UnfoldWithFit(Int_t RangeFlag = 0, Int_t TriggerFlag = 0){

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



  RooUnfoldResponse responsePhi      (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseTildePhi (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseCosTheta (25, -1., 1.);


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



      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5)  ) {
          PhiGenH     ->Fill( (PhiHEgen+TMath::Pi()) );
          CosThetaGenH->Fill( CosThetaHEgen          );
          CosThetaHEgen2 = -1.*CosThetaHEgen;
          PhiHEgen2 = PhiHEgen;
          /* -
           * - Retranslating the distributions...
           */
           // TildePhiHEposgen2  = PhiHEgen2 - 0.25 * TMath::Pi() ;
           // TildePhiHEneggen2  = PhiHEgen2 - 0.75 * TMath::Pi() ;
           TildePhiHEposgen2  = PhiHEgen2 - 1.25 * TMath::Pi() ;
           TildePhiHEneggen2  = PhiHEgen2 - 1.75 * TMath::Pi() ;
          if( TildePhiHEposgen2 < 0. ) {
            TildePhiHEposgen2 += 2. * TMath::Pi();
          }
          if( TildePhiHEposgen2 < 0. ) {
            TildePhiHEposgen2 += 2. * TMath::Pi();
          }
          if( TildePhiHEneggen2 < 0. ) {
            TildePhiHEneggen2 += 2. * TMath::Pi();
          }
          if( TildePhiHEneggen2 < 0. ) {
            TildePhiHEneggen2 += 2. * TMath::Pi();
          }

          // TildePhiHEposgen2 = TildePhiHEposgen;
          // TildePhiHEneggen2 = TildePhiHEneggen;
          Costhetahelp = CosThetaHEgen2;
          Costhetahelp2 = CosThetaHErec - CosThetaHEgen2;
          phihelp2 = PhiHErec - PhiHEgen2;
          if ( CosThetaHEgen2 > 0. ){
              TildePhiGenH->Fill( TildePhiHEposgen2 );
          } else {
              TildePhiGenH->Fill( TildePhiHEneggen2 );
          }
          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {
              // if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
              //     responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen2+TMath::Pi()) );
              //     responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen2+TMath::Pi()) );
              //     PhiRecH->Fill( (PhiHErec+TMath::Pi()) );
              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  // responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen2+TMath::Pi()) );
                  responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  // PhiRecH->Fill( (PhiHErec+TMath::Pi()) );
                  PhiRecH->Fill( PhiHErec );

                  // PhiRecMinusGenH->Fill(PhiHErec-phihelp);
                  PhiRecMinusGenH->Fill(PhiHErec - (PhiHEgen2+TMath::Pi()));
                  // cout << "Costhetahelp2 = " << Costhetahelp2 << endl;
                  // cout << "phihelp2      = " << phihelp2 << endl;



                  /*
                   * TildePhi response matrix
                  **/
                  TildePhiHEposrec2  = PhiHErec - 0.25 * TMath::Pi() ;
                  TildePhiHEnegrec2  = PhiHErec - 0.75 * TMath::Pi() ;
                  // TildePhiHEposrec2  = PhiHErec - 1.*TMath::Pi() - 0.25 * TMath::Pi() ;
                  // TildePhiHEnegrec2  = PhiHErec - 1.*TMath::Pi() - 0.75 * TMath::Pi() ;
                  if( TildePhiHEposrec2 < 0. ) {
                    TildePhiHEposrec2 += 2. * TMath::Pi();
                  }
                  if( TildePhiHEnegrec2 < 0. ) {
                    TildePhiHEnegrec2 += 2. * TMath::Pi();
                  }
                  if ( CosThetaHErec > 0. ){
                      // if( TildePhiHEposrec2 > 0. && TildePhiHEposrec2 < 2.*TMath::Pi() ){
                          responseTildePhi.Fill( TildePhiHEposrec2, TildePhiHEposgen2 );
                          TildePhiRecH->Fill( TildePhiHEposrec2 );
                      // }
                  } else {
                      // if( TildePhiHEnegrec2 > 0. && TildePhiHEnegrec2 < 2.*TMath::Pi() ){
                          responseTildePhi.Fill( TildePhiHEnegrec2, TildePhiHEneggen2 );
                          TildePhiRecH->Fill( TildePhiHEnegrec2 );
                      // }
                  }



              }
              // if ( CosThetaHErec > 0. ){
              //     if( TildePhiHEposrec > 0. && TildePhiHEposrec < 2.*TMath::Pi() ){
              //         responseTildePhi.Fill( TildePhiHEposrec, TildePhiHEposgen2 );
              //         TildePhiRecH->Fill( TildePhiHEposrec );
              //     }
              // } else {
              //     if( TildePhiHEnegrec > 0. && TildePhiHEnegrec < 2.*TMath::Pi() ){
              //         responseTildePhi.Fill( TildePhiHEnegrec, TildePhiHEneggen2 );
              //         TildePhiRecH->Fill( TildePhiHEnegrec );
              //     }
              // }
              if( CosThetaHErec > -1. && CosThetaHErec < 1. ){
                  if(TriggerSelectionFlag > 0.5){
                    responseCosTheta.Fill( CosThetaHErec, CosThetaHEgen2 );
                  } else {
                    responseCosTheta.Miss( CosThetaHEgen2 );
                  }
                  CosThetaRecH->Fill( CosThetaHErec );
                  // CosThetaRecMinusGenH->Fill(CosThetaHErec-Costhetahelp);
                  CosThetaRecMinusGenH->Fill(Costhetahelp2);
                  // cout << "Costhetahelp2 = " << Costhetahelp2 << endl;
                  // cout << "phihelp2      = " << phihelp2 << endl;
              }

          } else {
                  responsePhi.Miss( (PhiHEgen2+TMath::Pi()) );
                  if ( CosThetaHEgen2 > 0. ){
                      responseTildePhi.Miss( TildePhiHEposgen2 );
                  } else {
                      responseTildePhi.Miss( TildePhiHEneggen2 );
                  }
                  responseCosTheta.Miss( CosThetaHEgen2 );

          }

      }

  }




  // TFile* fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
  // TH1F *DataTildePhi         = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  // TFile* fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
  // TH1F *DataPhi              = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  // TFile* fileDataRawCosTheta = new TFile("RealData/CosThetaHeFrame.root");
  // TH1F *DataCosTheta         = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  TFile* fileDataRawTildePhi = 0x0;
  // TH1F *DataTildePhi         = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  TFile* fileDataRawPhi      = 0x0;
  // TH1F *DataPhi              = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  TFile* fileDataRawCosTheta = 0x0;
  // TH1F *DataCosTheta         = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");
  if(        TriggerFlag == 1 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger1/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 2 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger2/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 3 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger3/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 4 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger4/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 5 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger5/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 6 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger6/CosThetaHE/CosThetaHeFrame.root");
  } else if( TriggerFlag == 7 ) {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("TriggerUncert/PolTrigger7/CosThetaHE/CosThetaHeFrame.root");
  } else {
    fileDataRawTildePhi = new TFile("RealData/TildePhiHeFrameV2.root");
    fileDataRawPhi      = new TFile("RealData/PhiHeFrameV2.root");
    fileDataRawCosTheta = new TFile("RealData/CosThetaHeFrame.root");
  }
  TH1F *DataTildePhi         = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  TH1F *DataPhi              = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  TH1F *DataCosTheta         = (TH1F*)fileDataRawCosTheta->Get("CosThetaAfterSignalExtractionErrorsH");


  RooUnfoldBayes    unfold  (&responsePhi, DataPhi, 1);
  RooUnfoldBayes    unfold1 (&responsePhi, DataPhi, 2);
  RooUnfoldBayes    unfold2 (&responsePhi, DataPhi, 3);
  RooUnfoldBayes    unfold3 (&responsePhi, DataPhi, 4);
  // RooUnfoldSvd      unfold (&response, hist_measured, kterm);
  // RooUnfoldBinByBin unfold (&response, hist_measured);


  RooUnfoldBayes    unfoldCosTheta  (&responseCosTheta, DataCosTheta, 1);
  RooUnfoldBayes    unfoldTildePhi1  (&responseTildePhi, DataTildePhi, 1);
  RooUnfoldBayes    unfoldTildePhi2  (&responseTildePhi, DataTildePhi, 2);
  RooUnfoldBayes    unfoldTildePhi3  (&responseTildePhi, DataTildePhi, 3);
  RooUnfoldBayes    unfoldTildePhi4  (&responseTildePhi, DataTildePhi, 4);
  // RooUnfoldBayes    unfoldTildePhi1  (&responseTildePhi, TildePhiRecH, 1);
  // RooUnfoldBayes    unfoldTildePhi2  (&responseTildePhi, TildePhiRecH, 2);
  // RooUnfoldBayes    unfoldTildePhi3  (&responseTildePhi, TildePhiRecH, 3);
  // RooUnfoldBayes    unfoldTildePhi4  (&responseTildePhi, TildePhiRecH, 4);



  TH1D* unfolded  = (TH1D*) unfold.Hreco();
  TH1D* unfolded1 = (TH1D*) unfold1.Hreco();
  TH1D* unfolded2 = (TH1D*) unfold2.Hreco();
  TH1D* unfolded3 = (TH1D*) unfold3.Hreco();



  TH1D* unfoldedCosTheta   = (TH1D*) unfoldCosTheta.Hreco();
  TH1D* unfoldedTildePhi1  = (TH1D*) unfoldTildePhi1.Hreco();
  TH1D* unfoldedTildePhi2  = (TH1D*) unfoldTildePhi2.Hreco();
  TH1D* unfoldedTildePhi3  = (TH1D*) unfoldTildePhi3.Hreco();
  TH1D* unfoldedTildePhi4  = (TH1D*) unfoldTildePhi4.Hreco();


  new TCanvas;
  unfolded ->SetLineColor(kRed);
  unfolded1->SetLineColor(kBlue);
  unfolded2->SetLineColor(kYellow);
  unfolded3->SetLineColor(kGreen);
  PhiGenH  ->SetLineColor(kBlack);
  unfolded ->SetLineWidth(3);
  unfolded1->SetLineWidth(3);
  unfolded2->SetLineWidth(3);
  unfolded3->SetLineWidth(3);
  PhiGenH  ->SetLineWidth(3);
  unfolded ->GetYaxis()->SetRangeUser(0, 1.e+6);
  unfolded ->Draw("same");
  unfolded1->Draw("same");
  unfolded2->Draw("same");
  unfolded3->Draw("same");
  PhiGenH  ->Draw("same");


  new TCanvas;
  unfoldedCosTheta->SetLineColor(kRed);
  CosThetaGenH    ->SetLineColor(kBlue);
  unfoldedCosTheta->SetLineWidth(3);
  CosThetaGenH    ->SetLineWidth(3);
  unfoldedCosTheta->GetYaxis()->SetRangeUser(0, 1.e+6);
  unfoldedCosTheta->Draw("same");
  CosThetaGenH    ->Draw("same");


  new TCanvas;
  unfoldedTildePhi1->SetLineColor(kRed);
  unfoldedTildePhi2->SetLineColor(kGreen);
  unfoldedTildePhi3->SetLineColor(kYellow);
  unfoldedTildePhi4->SetLineColor(kMagenta);
  TildePhiGenH     ->SetLineColor(kBlue);
  unfoldedTildePhi1->SetLineWidth(3);
  unfoldedTildePhi2->SetLineWidth(3);
  unfoldedTildePhi3->SetLineWidth(3);
  unfoldedTildePhi4->SetLineWidth(3);
  TildePhiGenH     ->SetLineWidth(3);
  unfoldedTildePhi1->GetYaxis()->SetRangeUser(0, 1.e+6);
  unfoldedTildePhi1->Draw("same");
  unfoldedTildePhi2->Draw("same");
  unfoldedTildePhi3->Draw("same");
  unfoldedTildePhi4->Draw("same");
  TildePhiGenH     ->Draw("same");





  /**
    * Applying the response matrix
    * to a few distributions just to see
    * what happens afterwards.
    */
  // TH1F* DeltaFunctionH = new TH1F("DeltaFunctionH", "DeltaFunctionH", 25, 0., 2.*TMath::Pi());
  // DeltaFunctionH->Fill(1.,10000);
  // new TCanvas;
  // DeltaFunctionH->Draw();
  // new TCanvas;
  // TH1F* FoldedTruthH = (TH1F*) responsePhi.ApplyToTruth(DeltaFunctionH);
  // FoldedTruthH->Draw();
  //
  // // TF1 *fa1 = new TF1("fa1","sin(x)/x",0,2.*TMath::Pi());
  // auto form1 = new TFormula("form1","abs(sin(x)/x)");
  // auto sqroot = new TF1("sqroot","x*gaus(0) + [3]*form1",0,2.*TMath::Pi());
  // sqroot->SetParameters(10,4,1,20);
  // TH1F* SinFunctionH = new TH1F("SinFunctionH", "SinFunctionH", 25, 0., 2.*TMath::Pi());
  // SinFunctionH->FillRandom("sqroot",20000);
  // new TCanvas;
  // SinFunctionH->Draw();
  // new TCanvas;
  // TH1F* FoldedTruthH2 = (TH1F*) responsePhi.ApplyToTruth(SinFunctionH);
  // FoldedTruthH2->Draw();



  // PolarisationHeMinuit1D( unfoldedCosTheta, unfolded, unfoldedTildePhi1);
  PolarisationHeMinuit1D( unfoldedCosTheta, unfolded, unfoldedTildePhi1, RangeFlag, 0);



  TFile *SavingFile = new TFile("SavingFileData.root", "RECREATE");
  PhiRecH          ->Write();
  TildePhiRecH->Write();
  PhiGenH          ->Write();
  PhiRecMinusGenH->Write();
  CosThetaRecMinusGenH->Write();
  responsePhi.Write();
  responseCosTheta.Write();
  responseTildePhi.Write();
  unfolded         ->Write();
  unfoldedCosTheta ->Write();
  // unfoldedTildePhi1->Write();
  // unfoldedTildePhi2->Write();
  // unfoldedTildePhi3->Write();
  // unfoldedTildePhi4->Write();
  SavingFile       ->Close();

}
//______________________________________________
Int_t switchFlag = 0;

Double_t ReducedChiSquare = 0;


//_____________________________________________________________________________
/* - Coding in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
//______________________________________________
/* -
 * - CosTheta's distribution
 */
Double_t CosTheta(Double_t *x, Double_t *par) {
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t returnValue     = 1. + par[0] * CosSquaredTheta;
  // returnValue              = par[1] * returnValue / ( 3. + par[0] );
  returnValue              = par[4] * 3. * returnValue / ( 3. + par[0] );
  return   returnValue;
}
//______________________________________________
/* -
 * - Phi's distribution
 */
Double_t Phi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2. * x[1] );
  // Double_t returnValue = par[2] * ( 1. + 2. * par[3] * CosOfTwoPhi / ( 3. + par[0] ) );
  Double_t returnValue = par[4] * ( 1. + 2. * par[3] * CosOfTwoPhi / ( 3. + par[0] ) );
  return   returnValue;
}
//______________________________________________
/* -
 * - Phi's distribution but taken along the same axis as CosTheta
 */
Double_t PhiV2(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2. * x[0] );
  // Double_t returnValue = par[2] * ( 1. + 2. * par[3] * CosOfTwoPhi / ( 3. + par[0] ) );
  Double_t returnValue = par[4] * ( 1. + 2. * par[3] * CosOfTwoPhi / ( 3. + par[0] ) );
  return   returnValue;
}
//______________________________________________
Double_t DummyPhi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi = TMath::Cos( 2. * x[1] );
  Double_t returnValue = par[0] * ( 1. + 2. * par[1] * CosOfTwoPhi / ( 3. + 1.13220 ) );
  return   returnValue;
}
//______________________________________________
Double_t MixWithPositiveCosTheta(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2. * x[1] - 0.50 * TMath::Pi() );
  Double_t returnValue    = par[4] * ( 1. + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3. + par[0] ) );
  return   returnValue;
}
//______________________________________________
Double_t MixWithNegativeCosTheta(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2. * x[1] - 1.50 * TMath::Pi() );
  Double_t returnValue    = par[4] * ( 1. + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3. + par[0] ) );
  return   returnValue;
}
//______________________________________________
/* -
 * - TildePhi's distribution.
 */
Double_t TildePhi(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2. * x[2] );
  Double_t returnValue    = par[4] * ( 1. + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3. + par[0] ) );
  return   returnValue;
}
//______________________________________________
/* -
 * - TildePhi's distribution but along the same axis as x[0]
 */
Double_t TildePhiV2(Double_t *x, Double_t *par) {
  Double_t CosTwoTildePhi = TMath::Cos( 2. * x[0] );
  Double_t returnValue    = par[4] * ( 1. + TMath::Sqrt(2) * par[5] * CosTwoTildePhi / ( 3. + par[0] ) );
  return   returnValue;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFitComplete(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  if ( x[0] < 0 ) {
    sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par ) + MixWithNegativeCosTheta( x, par );
  } else {
    sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par ) + MixWithPositiveCosTheta( x, par );
  }
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFit(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = CosTheta( x, par ) + Phi( x, par ) + TildePhi( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFitLastHopeComplete(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  if        ( x[0] < 2*TMath::Pi() ){
    sumOfTheSubFits = CosTheta( x, par );
  } else if ( x[0] < 6*TMath::Pi() ){
    sumOfTheSubFits = PhiV2( x, par );
  } else                            {
    sumOfTheSubFits = TildePhiV2( x, par );
  }
  return sumOfTheSubFits;
}

//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t SimultaneousFitLastHope(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = ( x[0] < 2*TMath::Pi() ) ? CosTheta( x, par ) : PhiV2( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Dummy fit CosTheta only.
 * -
 */
Double_t DummyFitCosTheta(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = CosTheta( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Simultaneous fit complete.
 * -
 */
Double_t DummyFitPhi(double *x, double *par) {
  Double_t sumOfTheSubFits = 0;
  sumOfTheSubFits = DummyPhi( x, par );
  return sumOfTheSubFits;
}
//_____________________________________________________________________________
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coords;
std::vector< Double_t > values;
std::vector< Double_t > errors;

// void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
void FcnForMinimisation(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
  // cout << "HI" << flush << endl;
  Int_t n = coords.size();
  // cout << "n is " << n << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[3];
  for ( Int_t i = 0; i < n; ++i ) {
    if        ( i < 15 ) {
      x[0] = coords[i];
      x[1] = 0;
      x[2] = 0;
    } else if ( i < 40 ) {
      x[0] = 0;
      x[1] = coords[i];
      x[2] = 0;
    } else               {
      x[0] = 0;
      x[1] = 0;
      x[2] = coords[i];
    }

    // cout << "HI2" << flush << endl;
    if ( values[i] != 0 ) {
      tmp = ( values[i] - SimultaneousFit( x, p ) ) / errors[i];
      // tmp = ( values[i] - DummyFitCosTheta( x, p ) ) / errors[i];
      // tmp = ( values[i] - DummyFitPhi( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  f = chi2;
  cout << "ChiSquared = " << chi2 << endl;
}
//_____________________________________________________________________________
Int_t SignalRangeModeFromBash = 0;
Int_t Counter = 0;
void FcnForMinimisationV2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
  Int_t n = coords.size();
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  // cout << "SignalRangeModeFromBash = " << SignalRangeModeFromBash << endl;
  // cout << "Counter = " << Counter << endl;
  for ( Int_t i = 0; i < n; ++i ) {
    if        ( i < 15 ) {
    // if        ( i < Counter ) {
      x[0] = coords[i];
      x[1] = 0;
    } else if ( i < 40 ) {
    // } else if ( i < 25 + Counter ) {
      // x[0] = coords[i] + 4*TMath::Pi();
      x[0] = coords[i] + 4*3.14;
      x[1] = 0;
    } else               {
      // x[0] = coords[i] + 8*TMath::Pi();
      x[0] = coords[i] + 8*3.14;
      x[1] = 0;
    }
    if ( values[i] != 0 ) {
      tmp = ( values[i] - SimultaneousFitLastHopeComplete( x, p ) ) / errors[i];
      // tmp = ( values[i] - SimultaneousFit( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  f = chi2;
  cout << "ChiSquared = " << chi2 << endl;
}
//_____________________________________________________________________________
void FcnForMinimisationV3(Int_t &npar, Double_t *gin, Double_t &f, Double_t *p, Int_t iflag)
{
  Int_t n = coords.size();
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  // cout << "SignalRangeModeFromBash = " << SignalRangeModeFromBash << endl;
  // cout << "Counter = " << Counter << endl;
  for ( Int_t i = 0; i < n; ++i ) {
    if        ( i < 25 ) {
      // x[0] = coords[i] + 4*TMath::Pi();
      x[0] = coords[i] + 4*3.14;
      x[1] = 0;
    } else if ( i < 50 ) {
      // x[0] = coords[i] + 8*TMath::Pi();
      x[0] = coords[i] + 8*3.14;
      x[1] = 0;
    } else {
      x[0] = coords[i];
      x[1] = 0;
    }
    if ( values[i] != 0 ) {
      tmp = ( values[i] - SimultaneousFitLastHopeComplete( x, p ) ) / errors[i];
      // tmp = ( values[i] - SimultaneousFit( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  f = chi2;
  cout << "ChiSquared = " << chi2 << endl;
  ReducedChiSquare = chi2;

}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void PolarisationHeMinuit1D( TH1* CorrectedCosTheta, TH1* CorrectedPhi, TH1* CorrectedTildePhi, Int_t FitRangeMode = 0, Int_t SignalRangeSelectionMode = 0 ){
  // #ifdef __CINT__
  // // if (!TClass::GetDict("RooUnfold")) gSystem->Load("../RooUnfold/libRooUnfold");
  //   gSystem->Load("../RooUnfold/libRooUnfold");
  // #endif
  //
  // SignalRangeModeFromBash = SignalRangeSelectionMode;
  // TDatime d;
  // // TFile* file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // TFile* file1D = 0x0;
  // if        ( SignalRangeSelectionMode == 0 ) {
  //   file1D = new TFile("SavingFileData.root");
  //   // file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_long.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else if ( SignalRangeSelectionMode == 1 ) {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_long.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else if ( SignalRangeSelectionMode == 2 ) {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else if ( SignalRangeSelectionMode == 3 ) {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else if ( SignalRangeSelectionMode == 4 ) {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else if ( SignalRangeSelectionMode == 5 ) {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // } else {
  //   file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  // }
  // TH1D* CorrectedCosTheta = (TH1D*) file1D->Get("response;5");
  // TH1D* CorrectedPhi      = (TH1D*) file1D->Get("response;4");
  // TH1D* CorrectedTildePhi = (TH1D*) file1D->Get("response;6");


  Double_t CosThetaLowLimit   = -1;
  Double_t CosThetaUpperLimit = +1;
  Double_t PhiLowLimit        = -3.14;
  Double_t PhiUpperLimit      = +3.14;


  // TF2 * helicitySimultaneously1d = new TF2( "helicitySimultaneously1d",
  //                                           my2Dfunc,xlow2,xup2,ylow2,yup2, 10);


  Int_t nBinsCosTheta = CorrectedCosTheta->GetNbinsX();
  Int_t nBinsPhi      = CorrectedPhi     ->GetNbinsX();
  Int_t nBinsTildePhi = CorrectedTildePhi->GetNbinsX();


  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure

  // for (Int_t ix = 6; ix <= nBinsCosTheta-5; ++ix) {
  //   if        ( FitRangeMode == 1 ) {
  //     if ( (ix == 6) || (ix == nBinsCosTheta-5) )  continue;
  //   } else if ( FitRangeMode == 2 ) {
  //     if ( (ix == 6) || (ix == 7) || (ix == nBinsCosTheta-6) || (ix == nBinsCosTheta-5) )  continue;
  //   } else {
  //   }
  //   Counter+=1;
  //   coords.push_back( CorrectedCosTheta->GetXaxis()->GetBinCenter(ix) );
  //   values.push_back( CorrectedCosTheta->GetBinContent(ix)            );
  //   errors.push_back( CorrectedCosTheta->GetBinError(ix)              );
  // }
  for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
    coords.push_back( CorrectedPhi     ->GetXaxis()->GetBinCenter(iy) );
    values.push_back( CorrectedPhi     ->GetBinContent(iy)            );
    errors.push_back( CorrectedPhi     ->GetBinError(iy)              );
  }
  for (Int_t iy = 1; iy <= nBinsTildePhi; ++iy) {
    coords.push_back( CorrectedTildePhi->GetXaxis()->GetBinCenter(iy) );
    values.push_back( CorrectedTildePhi->GetBinContent(iy)            );
    errors.push_back( CorrectedTildePhi->GetBinError(iy)              );
  }

  for (Int_t ix = 6; ix <= nBinsCosTheta-5; ++ix) {
    if        ( FitRangeMode == 1 ) {
      if ( (ix == 6) || (ix == nBinsCosTheta-5) )  continue;
    } else if ( FitRangeMode == 2 ) {
      if ( (ix == 6) || (ix == 7) || (ix == nBinsCosTheta-6) || (ix == nBinsCosTheta-5) )  continue;
    } else if ( FitRangeMode == 3 ) {
      if ( ix == 6 )  continue;
    } else if ( FitRangeMode == 4 ) {
      if ( ix == nBinsCosTheta-5 )  continue;
    } else if ( FitRangeMode == 5 ) {
      if ( (ix == 6) || (ix == 7)|| (ix == 8)|| (ix == 9) || (ix == nBinsCosTheta-6) || (ix == nBinsCosTheta-5) )  continue;
    } else {
    }
    Counter+=1;
    coords.push_back( CorrectedCosTheta->GetXaxis()->GetBinCenter(ix) );
    values.push_back( CorrectedCosTheta->GetBinContent(ix)            );
    errors.push_back( CorrectedCosTheta->GetBinError(ix)              );
  }


  for( Int_t i = 0; i < 75; i++ ){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }


  TMinuit *gMinuit = new TMinuit(6);
  // gMinuit->SetFCN(FcnForMinimisation);
  // gMinuit->SetFCN(FcnForMinimisationV2);
  gMinuit->SetFCN(FcnForMinimisationV3);
  gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.01, -2, 2);
  gMinuit->DefineParameter(1, "NormalTheta", 0, 0,  -1., 3.2e+04);
  gMinuit->DefineParameter(2, "NormalisPhi",      0, 0,  0, 8700);
  gMinuit->DefineParameter(3, "LambdaPhi",           0, 0.01,    -2, 2   );
  // gMinuit->DefineParameter(3, "LambdaPhi",           0.01, 0.0,    -2, 2   );
  // gMinuit->DefineParameter(4, "NormalisTildePhi", 8300, 100,  7700, 8700);
  gMinuit->DefineParameter(4, "NormalisTildePhi", 7700, 100,  6000, 8700);
  gMinuit->DefineParameter(5, "LambdaThetaPhi",      0, 0.01,    -2, 2   );
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,    LambdaPhi,    NormalTheta,    NormalisPhi,    LambdaThetaPhi,    NormalisTildePhi;
  Double_t LambdaThetaErr, LambdaPhiErr, NormalThetaErr, NormalisPhiErr, LambdaThetaPhiErr, NormalisTildePhiErr;
  gMinuit->GetParameter(0, LambdaTheta,      LambdaThetaErr     );
  gMinuit->GetParameter(1, NormalTheta,      NormalThetaErr     );
  gMinuit->GetParameter(2, NormalisPhi,      NormalisPhiErr     );
  gMinuit->GetParameter(3, LambdaPhi,        LambdaPhiErr       );
  gMinuit->GetParameter(4, NormalisTildePhi, NormalisTildePhiErr);
  gMinuit->GetParameter(5, LambdaThetaPhi,   LambdaThetaPhiErr  );
  printf("LambdaTheta     : %+.7f +- %.7f\n",LambdaTheta,     LambdaThetaErr     );
  printf("NormalTheta     : %+.7f +- %.7f\n",NormalTheta,     NormalThetaErr     );
  printf("NormalisPhi     : %+.7f +- %.7f\n",NormalisPhi,     NormalisPhiErr     );
  printf("LambdaPhi       : %+.7f +- %.7f\n",LambdaPhi,       LambdaPhiErr       );
  printf("NormalisTildePhi: %+.7f +- %.7f\n",NormalisTildePhi,NormalisTildePhiErr);
  printf("LambdaThetaPhi  : %+.7f +- %.7f\n",LambdaThetaPhi,  LambdaThetaPhiErr  );
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);
  gMinuit->mnmatu(1);


  gStyle->SetOptStat(0);

  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -0.6 ,0.6 );
  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedCosTheta->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedCosTheta->GetYaxis()->SetTitleOffset(1.25);
  CorrectedCosTheta->GetYaxis()->SetTitleOffset(1.55);
  CorrectedCosTheta->GetXaxis()->SetTitleSize(0.045);
  CorrectedCosTheta->GetYaxis()->SetTitleSize(0.045);
  CorrectedCosTheta->GetXaxis()->SetLabelSize(0.045);
  CorrectedCosTheta->GetYaxis()->SetLabelSize(0.045);
  CorrectedCosTheta->GetXaxis()->SetTitleFont(42);
  CorrectedCosTheta->GetYaxis()->SetTitleFont(42);
  CorrectedCosTheta->GetXaxis()->SetLabelFont(42);
  CorrectedCosTheta->GetYaxis()->SetLabelFont(42);
  CorrectedCosTheta->GetXaxis()->SetNdivisions(408);
  CorrectedCosTheta->GetYaxis()->SetRangeUser(0., CorrectedCosTheta->GetMaximum()*0.016);
  // CorrectedCosTheta->GetXaxis()->SetRangeUser(2, 6);
  CorrectedCosTheta->SetTitle(  Form(  ";cos(#theta); ACCxEFF Corrected Counts / %.3f",
                           CorrectedCosTheta->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedCosTheta->Draw();
  TLatex* latex = new TLatex();
  latex->SetTextSize(0.05);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.045);
  latex->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  latex->DrawLatex(0.55,0.78,"Simultaneous Minuit Fit");
  latex->DrawLatex(0.55,0.70,Form("#lambda_{#theta} = %.3f #pm %.3f", LambdaTheta, LambdaThetaErr));
  // latex->DrawLatex(0.55,0.62,Form("#tilde{#chi} = %.3f #pm %.3f", LambdaTheta, LambdaThetaErr));
  latex->DrawLatex(0.55,0.18,Form(   "#tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     ReducedChiSquare,
                                     65 - gMinuit->GetNumFreePars(),
                                     ReducedChiSquare/((Double_t)  (65 - gMinuit->GetNumFreePars()))
                                     )
                                    );
  Model->SetParameter( 0, LambdaTheta );
  Model->SetParameter( 1, NormalisTildePhi*3. );
  // Model->SetParameter( 1, NormalTheta );
  Model->SetNpx(500);
  Model->Draw("same");
  if ( SignalRangeSelectionMode == 0 || FitRangeMode == 0 ) gPad->SaveAs("CosThetaHeMinuit.png", "recreate");
  gPad->SaveAs(Form("CosThetaHeMinuit_SigEx_%d_FitRange_%d_HE.png", SignalRangeSelectionMode, FitRangeMode), "recreate");


  // TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", 0., 2.*TMath::Pi() );
  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedPhi->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedPhi->GetYaxis()->SetTitleOffset(1.25);
  CorrectedPhi->GetYaxis()->SetTitleOffset(1.55);
  CorrectedPhi->GetXaxis()->SetTitleSize(0.045);
  CorrectedPhi->GetYaxis()->SetTitleSize(0.045);
  CorrectedPhi->GetXaxis()->SetLabelSize(0.045);
  CorrectedPhi->GetYaxis()->SetLabelSize(0.045);
  CorrectedPhi->GetXaxis()->SetTitleFont(42);
  CorrectedPhi->GetYaxis()->SetTitleFont(42);
  CorrectedPhi->GetXaxis()->SetLabelFont(42);
  CorrectedPhi->GetYaxis()->SetLabelFont(42);
  CorrectedPhi->GetXaxis()->SetNdivisions(408);
  CorrectedPhi->GetYaxis()->SetRangeUser(0., CorrectedPhi->GetMaximum()*0.016);
  // CorrectedPhi->GetXaxis()->SetRangeUser(2, 6);
  CorrectedPhi->SetTitle(  Form(  ";#phi; ACCxEFF Corrected Counts / %.3f",
                           CorrectedPhi->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedPhi->Draw();
  TLatex* latex2 = new TLatex();
  latex2->SetTextSize(0.05);
  latex2->SetTextFont(42);
  latex2->SetTextAlign(11);
  latex2->SetNDC();
  latex2->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex2->SetTextSize(0.045);
  latex2->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  latex2->DrawLatex(0.55,0.78,"Simultaneous Minuit Fit");
  latex2->DrawLatex(0.55,0.70,Form("#lambda_{#phi} = %.3f #pm %.3f",   LambdaPhi,   LambdaPhiErr));
  latex2->DrawLatex(0.55,0.62,Form("#lambda_{#theta} = %.3f #pm %.3f", LambdaTheta, LambdaThetaErr));
  // latex2->DrawLatex(0.55,0.62,Form("#tilde{#chi} = %.3f #pm %.3f", LambdaTheta, LambdaThetaErr));
  latex2->DrawLatex(0.55,0.18,Form(   "#tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     ReducedChiSquare,
                                     65 - gMinuit->GetNumFreePars(),
                                     ReducedChiSquare/((Double_t)  (65 - gMinuit->GetNumFreePars()))
                                     )
                                    );
  Model2->SetParameter( 0, LambdaTheta );
  Model2->SetParameter( 2, LambdaPhi );
  // Model2->SetParameter( 1, NormalisPhi );
  Model2->SetParameter( 1, NormalisTildePhi );
  Model2->SetNpx(500);
  Model2->Draw("same");
  if ( SignalRangeSelectionMode == 0 || FitRangeMode == 0 ) gPad->SaveAs("PhiHeMinuit.png", "recreate");
  gPad->SaveAs(Form("PhiHeMinuit_SigEx_%d_FitRange_%d_HE.png", SignalRangeSelectionMode, FitRangeMode), "recreate");

  TF1* Model3 = new TF1("Model3", "[1]*(1+TMath::Sqrt(2)*[2]*cos(2*x)/(3+[0]))", 0 ,6.2 );
  new TCanvas;
  gPad->SetMargin(0.13,0.01,0.12,0.01);
  CorrectedTildePhi->GetXaxis()->SetTitleOffset(1.15);
  // CorrectedTildePhi->GetYaxis()->SetTitleOffset(1.25);
  CorrectedTildePhi->GetYaxis()->SetTitleOffset(1.55);
  CorrectedTildePhi->GetXaxis()->SetTitleSize(0.045);
  CorrectedTildePhi->GetYaxis()->SetTitleSize(0.045);
  CorrectedTildePhi->GetXaxis()->SetLabelSize(0.045);
  CorrectedTildePhi->GetYaxis()->SetLabelSize(0.045);
  CorrectedTildePhi->GetXaxis()->SetTitleFont(42);
  CorrectedTildePhi->GetYaxis()->SetTitleFont(42);
  CorrectedTildePhi->GetXaxis()->SetLabelFont(42);
  CorrectedTildePhi->GetYaxis()->SetLabelFont(42);
  CorrectedTildePhi->GetXaxis()->SetNdivisions(408);
  CorrectedTildePhi->GetYaxis()->SetRangeUser(0., CorrectedPhi->GetMaximum()*4);
  // CorrectedTildePhi->GetXaxis()->SetRangeUser(2, 6);
  CorrectedTildePhi->SetTitle(  Form(  ";#tilde{#phi}; ACCxEFF Corrected Counts / %.3f",
                           CorrectedTildePhi->GetXaxis()->GetBinWidth(1)  )  );
  CorrectedTildePhi->Draw();
  TLatex* latex3 = new TLatex();
  latex3->SetTextSize(0.05);
  latex3->SetTextFont(42);
  latex3->SetTextAlign(11);
  latex3->SetNDC();
  latex3->DrawLatex(0.17,0.94,"ALICE Performance, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex3->SetTextSize(0.045);
  latex3->DrawLatex(0.55,0.84,"UPC, Run 2 dataset, HE");
  latex3->DrawLatex(0.55,0.78,"Simultaneous Minuit Fit");
  // latex3->DrawLatex(0.55,0.70,Form("#lambda_{#phi} = %.3f #pm %.3f",       LambdaPhi,      LambdaPhiErr));
  latex3->DrawLatex(0.55,0.70,Form("#lambda_{#theta} = %.3f #pm %.3f",     LambdaTheta,    LambdaThetaErr));
  latex3->DrawLatex(0.55,0.62,Form("#lambda_{#theta#phi} = %.3f #pm %.3f", LambdaThetaPhi, LambdaThetaPhiErr));
  // latex3->DrawLatex(0.55,0.62,Form("#tilde{#chi} = %.3f #pm %.3f", LambdaTheta, LambdaThetaErr));
  latex3->DrawLatex(0.55,0.18,Form(   "#tilde{#chi}^{2} = %.2f / %.2d = %.2f  ",
                                     ReducedChiSquare,
                                     65 - gMinuit->GetNumFreePars(),
                                     ReducedChiSquare/((Double_t)  (65 - gMinuit->GetNumFreePars()))
                                     )
                                    );
  Model3->SetParameter( 0, LambdaTheta      );
  Model3->SetParameter( 2, LambdaThetaPhi   );
  Model3->SetParameter( 1, NormalisTildePhi );
  Model3->SetNpx(500);
  Model3->Draw("same");
  if ( SignalRangeSelectionMode == 0 || FitRangeMode == 0 ) gPad->SaveAs("TildePhiHeMinuit.png", "recreate");
  gPad->SaveAs(Form("TildePhiHeMinuit_SigEx_%d_FitRange_%d_HE.png", SignalRangeSelectionMode, FitRangeMode), "recreate");


}
