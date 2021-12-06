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


void PtAxE(){


  TFile* file   = new TFile("AnalysisResultsMC_incoh_16r_2.root", "READ");
  TTree* ReconTree = (TTree*)file->Get("MyTask/fOutputTree");
  TTree* GenerTree = (TTree*)file->Get("MyTask/fOutputTreeMC");

  // create a pointer to an event object. This will be used
  // to read the branch values.
  // Event *event      = new Event();
  // Event *eventGener = new Event();


  Double_t pTrec;
  Double_t pTgen;
  ReconTree->SetBranchAddress("fTrkTrkPt",  &pTrec);
  GenerTree->SetBranchAddress("fMCTrkTrkPt",&pTgen);


  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


  TH1F *PtRecH = new TH1F("PtRecH", "PtRecH", 4000, 0, 20);
  TH1F *PtGenH = new TH1F("PtGenH", "PtGenH", 4000, 0, 20);

  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    PtRecH   ->Fill(pTrec);
  }


  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesGen; i++) {
    GenerTree->GetEntry(i);
    PtGenH   ->Fill(pTgen);
  }



  new TCanvas;
  PtRecH->Draw();
  new TCanvas;
  PtGenH->Draw();


  PtRecH->Sumw2();
  PtGenH->Sumw2();
  PtRecH->Rebin(10);
  PtGenH->Rebin(10);



  TH1F *PtAxEH = (TH1F*)PtRecH->Clone("PtAxEH");
  PtAxEH->Divide(PtRecH,PtGenH);
  new TCanvas;
  PtAxEH->Draw("PE");

  TFile *SavingFile = new TFile("SavingFile.root", "RECREATE");
  PtAxEH->Write();
  SavingFile->Close();

}
