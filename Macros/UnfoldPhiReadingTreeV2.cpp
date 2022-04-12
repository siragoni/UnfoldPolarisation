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


void UnfoldPhiReadingTree(Int_t Iterations = 1){

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


  TH1F *CosThetaRecH = new TH1F("CosThetaRecH", "CosThetaRecH", 24, -1., 1.);
  TH1F *CosThetaGenH = new TH1F("CosThetaGenH", "CosThetaGenH", 24, -1., 1.);
  // TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      25, 0., 2.*TMath::Pi());
  // TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      25, 0., 2.*TMath::Pi());
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      250, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      250, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());



  // RooUnfoldResponse responsePhi      (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responsePhi      (250,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseTildePhi (25,  0., 2.*TMath::Pi());
  RooUnfoldResponse responseCosTheta (25, -1., 1.);


  Double_t Flag = -999.;
  Double_t Flag2 = -999.;
  TRandom3 *r3=new TRandom3();

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);
    Flag = -999.;
    Flag2 = -999.;
    // CosThetaHEgen2 = -1.*CosThetaHEgen;
    // if((CosThetaHEgen < 1.) || (CosThetaHEgen > -1.)){
      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {
        Flag2 = r3->Uniform(0,1);

        CosThetaHEgen2 = -1.*CosThetaHEgen;
        // if((CosThetaHEgen2 < 0.3) && (CosThetaHEgen2 > -0.3)){
        if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
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

        if ( CosThetaHEgen2 > 0. ){
            TildePhiGenH->Fill( TildePhiHEposgen2 );
        } else {
            TildePhiGenH->Fill( TildePhiHEneggen2 );
        }

        // if ( PhiGenH->GetBinContent(PhiGenH->FindBin(PhiHEgen2+TMath::Pi())) < (14500.*25./(2.*TMath::Pi()))*(1.+0.5*1.*TMath::Cos(2.*(PhiHEgen2+TMath::Pi()) ) ) ){
        if ( PhiGenH->GetBinContent(PhiGenH->FindBin(PhiHEgen2+TMath::Pi())) < (1000.*25./(2.*TMath::Pi()))*(1.+0.5*1.*TMath::Cos(2.*(PhiHEgen2+TMath::Pi()) ) ) ){

          PhiGenH     ->Fill( (PhiHEgen2+TMath::Pi()) );
          // PhiGenH     ->Fill( PhiHEgen2 );
          CosThetaGenH->Fill( CosThetaHEgen2          );
          Flag = 1;
        }

          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) ) {

              if( (PhiHErec+TMath::Pi()) > 0. && (PhiHErec+TMath::Pi()) < 2.*TMath::Pi() ){
                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();
                  // Flag2 = r3->Uniform(0,1);
                  if (Flag > 0.){
                    // responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                    responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  // } else {
                  //   if (Flag2 > 0.0){
                  //     // responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                  //     responsePhi.Fill( PhiHErec, (PhiHEgen2+TMath::Pi()) );
                  //   }
                  }

                  if ( Flag > 0. ){
                    PhiRecH->Fill( PhiHErec );
                  }
              }
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


              if( CosThetaHErec > -1. && CosThetaHErec < 1. ){
                  responseCosTheta.Fill( CosThetaHErec, CosThetaHEgen2 );
                  CosThetaRecH->Fill( CosThetaHErec );
                  CosThetaRecMinusGenH->Fill(CosThetaHErec-CosThetaHEgen2);

              }
              // }
          } else {
                if ( Flag > 0. ){
                  responsePhi.Miss( (PhiHEgen2+TMath::Pi()) );
                // } else {
                //   if (Flag2 > 0.0){
                //     // responsePhi.Fill( (PhiHErec+TMath::Pi()), (PhiHEgen+TMath::Pi()) );
                //     responsePhi.Miss( (PhiHEgen2+TMath::Pi()) );
                //   }
                }
                  if ( CosThetaHErec > 0. ){
                      responseTildePhi.Miss( TildePhiHEposgen2 );
                  } else {
                      responseTildePhi.Miss( TildePhiHEneggen2 );
                  }
                  responseCosTheta.Miss( CosThetaHEgen2 );

          }
        // }
      }
    }
  }



  RooUnfoldBayes    unfold  (&responsePhi, PhiRecH, Iterations);
  RooUnfoldBayes    unfold1 (&responsePhi, PhiRecH, 5);
  RooUnfoldBayes    unfold2 (&responsePhi, PhiRecH, 10);
  RooUnfoldBayes    unfold3 (&responsePhi, PhiRecH, 20);
  // RooUnfoldSvd      unfold (&response, hist_measured, kterm);
  // RooUnfoldBinByBin unfold (&response, hist_measured);


  RooUnfoldBayes    unfoldCosTheta  (&responseCosTheta, CosThetaRecH, 1);
  RooUnfoldBayes    unfoldTildePhi1  (&responseTildePhi, TildePhiRecH, 1);
  RooUnfoldBayes    unfoldTildePhi2  (&responseTildePhi, TildePhiRecH, 2);
  RooUnfoldBayes    unfoldTildePhi3  (&responseTildePhi, TildePhiRecH, 3);
  RooUnfoldBayes    unfoldTildePhi4  (&responseTildePhi, TildePhiRecH, 4);



  TH1D* unfolded  = (TH1D*) unfold.Hreco();
  TH1D* unfolded1 = (TH1D*) unfold1.Hreco();
  TH1D* unfolded2 = (TH1D*) unfold2.Hreco();
  TH1D* unfolded3 = (TH1D*) unfold3.Hreco();


  new TCanvas;
  // TF1* fit = new TF1("fit", "[0]*(1+0.5*[1]*cos(2.*x))", 0, 2.*TMath::Pi());
  TF1* fit = new TF1("fit", "[0]*(1.+0.5*[1]*cos(2.*x) -[3]*sin(x+[2]) )", 0, 2.*TMath::Pi());
  // TF1* fit = new TF1("fit", "[0]*(1-[1]*sin(x+[2]))", 0, 2.*TMath::Pi());
  fit->SetParameter(0,400);
  fit->SetParameter(1,0.5);
  fit->FixParameter(3,4.29585e-02);
  fit->FixParameter(2,6.09946);
  // fit->SetParameter(2,0.000000001);
  unfolded1->Draw();
  unfolded1->Fit("fit");


  TF1* fit2 = new TF1("fit2", "[0]*(1.+0.5*[1]*cos(2.*x) )", 0, 2.*TMath::Pi());
  fit2->SetParameter(0,400);
  fit2->SetParameter(1,0.5);
  // unfolded->Fit("fit2");
  // PhiGenH->Fit("fit2");

  // new TCanvas;
  // gPad->SetMargin(0.13,0.10,0.12,0.10);
  // gPad->SetTickx(1);
  // gPad->SetTicky(1);
  // gPad->SetGridx();
  // gPad->SetGridy();
  // gStyle->SetOptStat(0);
  // unfolded->SetTitle("");
  // unfolded->GetXaxis()->SetTitleOffset(1.15);
  // unfolded->GetYaxis()->SetTitleOffset(1.25);
  // unfolded->GetXaxis()->SetTitleSize(0.045);
  // unfolded->GetYaxis()->SetTitleSize(0.045);
  // unfolded->GetXaxis()->SetLabelSize(0.045);
  // unfolded->GetYaxis()->SetLabelSize(0.045);
  // unfolded->GetXaxis()->SetTitleFont(42);
  // unfolded->GetYaxis()->SetTitleFont(42);
  // unfolded->GetXaxis()->SetLabelFont(42);
  // unfolded->GetYaxis()->SetLabelFont(42);
  // unfolded->GetXaxis()->SetTitle("#varphi");
  // unfolded->GetYaxis()->SetTitle("Counts [a.u.]");
  // // unfolded->GetYaxis()->SetRangeUser(0.,0.2);
  // unfolded->GetXaxis()->SetRangeUser(0.*TMath::Pi() , 2.*TMath::Pi());
  // unfolded->SetLineWidth(5);
  // unfolded->SetLineColor(2);
  // // unfolded->Scale(fit2->GetParameter(0)/fit->GetParameter(0));
  // unfolded->Draw("same");
  // PhiGenH ->Draw("same");
  // TLegend *leg_pt = new TLegend(0.5,0.45,0.85,0.79);
  // leg_pt->SetFillStyle(0);
  // leg_pt->SetBorderSize(0);
  // leg_pt->SetTextSize(0.042);
  // leg_pt->AddEntry(unfolded,"Unfolded distribution", "LP");
  // leg_pt->AddEntry(PhiGenH, "Generated modulation", "LP");
  // leg_pt->Draw();
  // gPad->SaveAs("genvsfolded.pdf", "recreate");




  // TH1D* unfoldedCosTheta   = (TH1D*) unfoldCosTheta.Hreco();
  // TH1D* unfoldedTildePhi1  = (TH1D*) unfoldTildePhi1.Hreco();
  // TH1D* unfoldedTildePhi2  = (TH1D*) unfoldTildePhi2.Hreco();
  // TH1D* unfoldedTildePhi3  = (TH1D*) unfoldTildePhi3.Hreco();
  // TH1D* unfoldedTildePhi4  = (TH1D*) unfoldTildePhi4.Hreco();


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


  // new TCanvas;
  // unfoldedCosTheta->SetLineColor(kRed);
  // CosThetaGenH    ->SetLineColor(kBlue);
  // unfoldedCosTheta->SetLineWidth(3);
  // CosThetaGenH    ->SetLineWidth(3);
  // unfoldedCosTheta->GetYaxis()->SetRangeUser(0, 1.e+6);
  // unfoldedCosTheta->Draw("same");
  // CosThetaGenH    ->Draw("same");


  // new TCanvas;
  // unfoldedTildePhi1->SetLineColor(kRed);
  // unfoldedTildePhi2->SetLineColor(kGreen);
  // unfoldedTildePhi3->SetLineColor(kYellow);
  // unfoldedTildePhi4->SetLineColor(kMagenta);
  // TildePhiGenH     ->SetLineColor(kBlue);
  // unfoldedTildePhi1->SetLineWidth(3);
  // unfoldedTildePhi2->SetLineWidth(3);
  // unfoldedTildePhi3->SetLineWidth(3);
  // unfoldedTildePhi4->SetLineWidth(3);
  // TildePhiGenH     ->SetLineWidth(3);
  // unfoldedTildePhi1->GetYaxis()->SetRangeUser(0, 1.e+6);
  // unfoldedTildePhi1->Draw("same");
  // unfoldedTildePhi2->Draw("same");
  // unfoldedTildePhi3->Draw("same");
  // unfoldedTildePhi4->Draw("same");
  // TildePhiGenH     ->Draw("same");


  TFile *SavingFile = new TFile("SavingFile.root", "RECREATE");
  PhiRecH          ->Write();
  PhiGenH          ->Write();
  CosThetaRecH          ->Write();
  CosThetaGenH          ->Write();
  PhiRecMinusGenH->Write();
  CosThetaRecMinusGenH->Write();
  responsePhi.Write();
  responseCosTheta.Write();
  responseTildePhi.Write();
  unfolded         ->Write();
  // unfoldedCosTheta ->Write();
  // unfoldedTildePhi1->Write();
  // unfoldedTildePhi2->Write();
  // unfoldedTildePhi3->Write();
  // unfoldedTildePhi4->Write();
  SavingFile       ->Close();

}
