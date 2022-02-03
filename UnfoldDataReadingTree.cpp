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


void UnfoldPhiReadingTree(){

  #ifdef __CINT__
  // if (!TClass::GetDict("RooUnfold")) gSystem->Load("../RooUnfold/libRooUnfold");
    gSystem->Load("../RooUnfold/libRooUnfold");
  #endif

  TFile* file   = new TFile("AnalysisResultsCoh.root", "READ");
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
  ReconTree->SetBranchAddress("fCosThetaHE",  &CosThetaHErec);
  GenerTree->SetBranchAddress("fMCCosThetaHE",&CosThetaHEgen);
  Double_t CosThetaCSrec;
  Double_t CosThetaCSgen;
  ReconTree->SetBranchAddress("fCosThetaCS",  &CosThetaCSrec);
  GenerTree->SetBranchAddress("fMCCosThetaCS",&CosThetaCSgen);
  Double_t PhiHErec;
  Double_t PhiHEgen;
  ReconTree->SetBranchAddress("fPhiHE",  &PhiHErec);
  GenerTree->SetBranchAddress("fMCPhiHE",&PhiHEgen);
  Double_t PhiCSrec;
  Double_t PhiCSgen;
  ReconTree->SetBranchAddress("fPhiCS",  &PhiCSrec);
  GenerTree->SetBranchAddress("fMCPhiCS",&PhiCSgen);
  Double_t TildePhiHEposrec;
  Double_t TildePhiHEposgen;
  ReconTree->SetBranchAddress("fTildePhiHEpos",  &TildePhiHEposrec);
  GenerTree->SetBranchAddress("fMCTildePhiHEpos",&TildePhiHEposgen);
  Double_t TildePhiCSposrec;
  Double_t TildePhiCSposgen;
  ReconTree->SetBranchAddress("fTildePhiCSpos",  &TildePhiCSposrec);
  GenerTree->SetBranchAddress("fMCTildePhiCSpos",&TildePhiCSposgen);
  Double_t TildePhiHEnegrec;
  Double_t TildePhiHEneggen;
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

  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 25, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 25, -1., 1.);
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      25, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      25, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());



  RooUnfoldResponse responsePhi      (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseTildePhi (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseCosTheta (25, -1., 1.);


  RooUnfoldResponse responsePhiBins1 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins2 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins3 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins4 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins5 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins6 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins7 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins8 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins9 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins10 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins11 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins12 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins13 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins14 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins15 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins16 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins17 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins18 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins19 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins20 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins21 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins22 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins23 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins24 (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhiBins25 (25,  0., 2.*TMath::Pi());


  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);
      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {
          PhiGenH     ->Fill( (PhiHEgen+TMath::Pi()) );
          CosThetaGenH->Fill( CosThetaHEgen          );
          if ( CosThetaHEgen > 0. ){
              TildePhiGenH->Fill( TildePhiHEposgen );
          } else {
              TildePhiGenH->Fill( TildePhiHEneggen );
          }
          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {
              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
                  responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  PhiRecH->Fill( (PhiHErec+TMath::Pi()) );
                  if (        CosThetaHEgen > -1. && CosThetaHEgen < -0.92 ) {
                    responsePhiBins1.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.92 && CosThetaHEgen < -0.84 ) {
                    responsePhiBins2.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.84 && CosThetaHEgen < -0.76 ) {
                    responsePhiBins3.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.76 && CosThetaHEgen < -0.68 ) {
                    responsePhiBins4.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.68 && CosThetaHEgen < -0.60 ) {
                    responsePhiBins5.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.60 && CosThetaHEgen < -0.52 ) {
                    responsePhiBins6.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.52 && CosThetaHEgen < -0.44 ) {
                    responsePhiBins7.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.44 && CosThetaHEgen < -0.36 ) {
                    responsePhiBins8.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.36 && CosThetaHEgen < -0.28 ) {
                    responsePhiBins9.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.28 && CosThetaHEgen < -0.20 ) {
                    responsePhiBins10.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.20 && CosThetaHEgen < -0.12 ) {
                    responsePhiBins11.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.12 && CosThetaHEgen < -0.04 ) {
                    responsePhiBins12.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen > -0.04 && CosThetaHEgen <  0.04 ) {
                    responsePhiBins13.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.04 && CosThetaHEgen <  0.12 ) {
                    responsePhiBins14.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.12 && CosThetaHEgen <  0.20 ) {
                    responsePhiBins15.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.20 && CosThetaHEgen <  0.28 ) {
                    responsePhiBins16.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.28 && CosThetaHEgen <  0.36 ) {
                    responsePhiBins17.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.36 && CosThetaHEgen <  0.44 ) {
                    responsePhiBins18.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.44 && CosThetaHEgen <  0.52 ) {
                    responsePhiBins19.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.52 && CosThetaHEgen <  0.60 ) {
                    responsePhiBins20.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.60 && CosThetaHEgen <  0.68 ) {
                    responsePhiBins21.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.68 && CosThetaHEgen <  0.76 ) {
                    responsePhiBins22.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.76 && CosThetaHEgen <  0.84 ) {
                    responsePhiBins23.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.84 && CosThetaHEgen <  0.92 ) {
                    responsePhiBins24.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  } else if ( CosThetaHEgen >  0.92 && CosThetaHEgen <  1.00 ) {
                    responsePhiBins25.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  }
              }
              if ( CosThetaHErec > 0. ){
                  if( TildePhiHEposrec > 0. && TildePhiHEposrec < 2.*TMath::Pi() ){
                      responseTildePhi.Fill( TildePhiHEposrec, TildePhiHEposgen );
                      TildePhiRecH->Fill( TildePhiHEposrec );
                  }
              } else {
                  if( TildePhiHEnegrec > 0. && TildePhiHEnegrec < 2.*TMath::Pi() ){
                      responseTildePhi.Fill( TildePhiHEnegrec, TildePhiHEneggen );
                      TildePhiRecH->Fill( TildePhiHEnegrec );
                  }
              }
              if( CosThetaHErec > -1. && CosThetaHErec < 1. ){
                  responseCosTheta.Fill( CosThetaHErec, CosThetaHEgen );
                  CosThetaRecH->Fill( CosThetaHErec );
              }

          } else {
                  responsePhi.Miss( (PhiHEgen+TMath::Pi()) );
                  if ( CosThetaHErec > 0. ){
                      responseTildePhi.Miss( TildePhiHEposgen );
                  } else {
                      responseTildePhi.Miss( TildePhiHEneggen );
                  }
                  responseCosTheta.Miss( CosThetaHEgen );

          }

      }

  }




  TFile* fileDataRawTildePhi = new TFile("../MyUPC/pngResults/2021-09-21/TildePhiHEv2/TildePhiHeFrameV2.root");
  TH1F *DataTildePhi         = (TH1F*)fileDataRawTildePhi->Get("TildePhiAfterSignalExtractionErrorsH");
  TFile* fileDataRawPhi      = new TFile("../MyUPC/pngResults/2021-09-21/PhiHEv2/PhiHeFrameV2.root");
  TH1F *DataPhi              = (TH1F*)fileDataRawPhi->Get("PhiAfterSignalExtractionErrorsH");
  TFile* fileDataRawCosTheta = new TFile("../MyUPC/pngResults/2021-09-21/CosThetaHE/CosThetaHeFrame.root");
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


  TFile *SavingFile = new TFile("SavingFileData.root", "RECREATE");
  PhiRecH          ->Write();
  PhiGenH          ->Write();
  responsePhi.Write();
  responsePhiBins1.Write();
  responsePhiBins2.Write();
  responsePhiBins3.Write();
  responsePhiBins4.Write();
  responsePhiBins5.Write();
  responsePhiBins6.Write();
  responsePhiBins7.Write();
  responsePhiBins8.Write();
  responsePhiBins9.Write();
  responsePhiBins10.Write();
  responsePhiBins11.Write();
  responsePhiBins12.Write();
  responsePhiBins13.Write();
  responsePhiBins14.Write();
  responsePhiBins15.Write();
  responsePhiBins16.Write();
  responsePhiBins17.Write();
  responsePhiBins18.Write();
  responsePhiBins19.Write();
  responsePhiBins20.Write();
  responsePhiBins21.Write();
  responsePhiBins22.Write();
  responsePhiBins23.Write();
  responsePhiBins24.Write();
  responsePhiBins25.Write();
  responseCosTheta.Write();
  responseTildePhi.Write();
  unfolded         ->Write();
  unfoldedCosTheta ->Write();
  unfoldedTildePhi1->Write();
  unfoldedTildePhi2->Write();
  unfoldedTildePhi3->Write();
  unfoldedTildePhi4->Write();
  SavingFile       ->Close();

}
