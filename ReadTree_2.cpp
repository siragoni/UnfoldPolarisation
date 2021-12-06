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
  // TFile* file   = new TFile("AnalysisResultsNewMC_LHC16r_GammaLow.root", "READ");
  // TFile* file   = new TFile("AnalysisResultsMC_gamma_16s.root", "READ");
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
  Double_t Mrec;
  Double_t Mgen;
  ReconTree->SetBranchAddress("fTrkTrkM",  &Mrec);
  GenerTree->SetBranchAddress("fMCTrkTrkM",&Mgen);
  Double_t Yrec;
  Double_t Ygen;
  ReconTree->SetBranchAddress("fTrkTrkY",  &Yrec);
  GenerTree->SetBranchAddress("fMCTrkTrkY",&Ygen);


  Int_t nentriesRec = (Int_t) ReconTree->GetEntries();
  Int_t nentriesGen = (Int_t) GenerTree->GetEntries();


  TH1F *PtRecH = new TH1F("PtRecH", "PtRecH", 4000, 0, 20);
  TH1F *PtGenH = new TH1F("PtGenH", "PtGenH", 4000, 0, 20);


  TH1F *PtRecH_0 = new TH1F("PtRecH_0", "PtRecH_0", 4000, 0, 20);
  TH1F *PtRecH_1 = new TH1F("PtRecH_1", "PtRecH_1", 4000, 0, 20);
  TH1F *PtRecH_2 = new TH1F("PtRecH_2", "PtRecH_2", 4000, 0, 20);
  TH1F *PtRecH_3 = new TH1F("PtRecH_3", "PtRecH_3", 4000, 0, 20);
  TH1F *PtRecH_4 = new TH1F("PtRecH_4", "PtRecH_4", 4000, 0, 20);

  TH1F *MassRapRec[20];
  for( Int_t iRapidity = 0; iRapidity < 20; iRapidity++ ) {
    MassRapRec[iRapidity] = new TH1F( Form("MassRapRec_%d", iRapidity), Form("MassRapRec_%d", iRapidity), 100, 0, 5);
  }
  TH1F *MassRapGen[20];
  for( Int_t iRapidity = 0; iRapidity < 20; iRapidity++ ) {
    MassRapGen[iRapidity] = new TH1F( Form("MassRapGen_%d", iRapidity), Form("MassRapGen_%d", iRapidity), 100, 0, 5);
  }


  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    PtRecH   ->Fill(pTrec);
    if(Mrec > 0.8 && Mrec < 1.3) PtRecH_0->Fill(pTrec);
    if(Mrec > 1.3 && Mrec < 1.8) PtRecH_1->Fill(pTrec);
    if(Mrec > 1.8 && Mrec < 2.3) PtRecH_2->Fill(pTrec);
    if(Mrec > 2.3 && Mrec < 2.8) PtRecH_3->Fill(pTrec);
    if(Mrec > 1.0 && Mrec < 2.5) PtRecH_4->Fill(pTrec);
    // if ( Mrec > 1.7 && Mrec < 2.2 ) {
    if ( Mrec > 0. && Mrec < 6. ) {
      if ( Yrec > -4.0 && Yrec < -2.5 ) {
        MassRapRec[0]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.625 ) {
        MassRapRec[1]->Fill(Mrec);
      }
      if ( Yrec > -3.625 && Yrec < -3.25 ) {
        MassRapRec[2]->Fill(Mrec);
      }
      if ( Yrec > -3.25  && Yrec < -2.875 ) {
        MassRapRec[3]->Fill(Mrec);
      }
      if ( Yrec > -2.875 && Yrec < -2.5 ) {
        MassRapRec[4]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.5 ) {
        MassRapRec[5]->Fill(Mrec);
      }
      if ( Yrec > -3.5   && Yrec < -3.0 ) {
        MassRapRec[6]->Fill(Mrec);
      }
      if ( Yrec > -3.0   && Yrec < -2.5 ) {
        MassRapRec[7]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.25 ) {
        MassRapRec[8]->Fill(Mrec);
      }
      if ( Yrec > -3.25   && Yrec < -2.5 ) {
        MassRapRec[9]->Fill(Mrec);
      }
    }
    if ( Mrec > 2.2 && Mrec < 2.7 ) {
      if ( Yrec > -4.0 && Yrec < -2.5 ) {
        MassRapRec[10]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.625 ) {
        MassRapRec[11]->Fill(Mrec);
      }
      if ( Yrec > -3.625 && Yrec < -3.25 ) {
        MassRapRec[12]->Fill(Mrec);
      }
      if ( Yrec > -3.25  && Yrec < -2.875 ) {
        MassRapRec[13]->Fill(Mrec);
      }
      if ( Yrec > -2.875 && Yrec < -2.5 ) {
        MassRapRec[14]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.5 ) {
        MassRapRec[15]->Fill(Mrec);
      }
      if ( Yrec > -3.5   && Yrec < -3.0 ) {
        MassRapRec[16]->Fill(Mrec);
      }
      if ( Yrec > -3.0   && Yrec < -2.5 ) {
        MassRapRec[17]->Fill(Mrec);
      }
      if ( Yrec > -4.0   && Yrec < -3.25 ) {
        MassRapRec[18]->Fill(Mrec);
      }
      if ( Yrec > -3.25   && Yrec < -2.5 ) {
        MassRapRec[19]->Fill(Mrec);
      }
    }

  }


  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesGen; i++) {
    GenerTree->GetEntry(i);
    PtGenH   ->Fill(pTgen);
    if ( Mgen > 1.7 && Mgen < 2.2 ) {
      if ( Ygen > -4.0 && Ygen < -2.5 ) {
        MassRapGen[0]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.625 ) {
        MassRapGen[1]->Fill(Mgen);
      }
      if ( Ygen > -3.625 && Ygen < -3.25 ) {
        MassRapGen[2]->Fill(Mgen);
      }
      if ( Ygen > -3.25  && Ygen < -2.875 ) {
        MassRapGen[3]->Fill(Mgen);
      }
      if ( Ygen > -2.875 && Ygen < -2.5 ) {
        MassRapGen[4]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.5 ) {
        MassRapGen[5]->Fill(Mgen);
      }
      if ( Ygen > -3.5   && Ygen < -3.0 ) {
        MassRapGen[6]->Fill(Mgen);
      }
      if ( Ygen > -3.0   && Ygen < -2.5 ) {
        MassRapGen[7]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.25 ) {
        MassRapGen[8]->Fill(Mgen);
      }
      if ( Ygen > -3.25   && Ygen < -2.5 ) {
        MassRapGen[9]->Fill(Mgen);
      }
    }
    if ( Mgen > 2.2 && Mgen < 2.7 ) {
      if ( Ygen > -4.0 && Ygen < -2.5 ) {
        MassRapGen[10]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.625 ) {
        MassRapGen[11]->Fill(Mgen);
      }
      if ( Ygen > -3.625 && Ygen < -3.25 ) {
        MassRapGen[12]->Fill(Mgen);
      }
      if ( Ygen > -3.25  && Ygen < -2.875 ) {
        MassRapGen[13]->Fill(Mgen);
      }
      if ( Ygen > -2.875 && Ygen < -2.5 ) {
        MassRapGen[14]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.5 ) {
        MassRapGen[15]->Fill(Mgen);
      }
      if ( Ygen > -3.5   && Ygen < -3.0 ) {
        MassRapGen[16]->Fill(Mgen);
      }
      if ( Ygen > -3.0   && Ygen < -2.5 ) {
        MassRapGen[17]->Fill(Mgen);
      }
      if ( Ygen > -4.0   && Ygen < -3.25 ) {
        MassRapGen[18]->Fill(Mgen);
      }
      if ( Ygen > -3.25   && Ygen < -2.5 ) {
        MassRapGen[19]->Fill(Mgen);
      }
    }

  }


  for (size_t i = 0; i < 20; i++) {
    cout << "AxE M and Y [" << i <<"]  =  " << MassRapRec[i]->GetEntries()/MassRapGen[i]->GetEntries() << " +/- " << (MassRapRec[i]->GetEntries()/MassRapGen[i]->GetEntries())*TMath::Sqrt(1./((Double_t)MassRapRec[i]->GetEntries())) /*<< "% "*/ << endl;
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
  PtRecH_0->Write();
  PtRecH_1->Write();
  PtRecH_2->Write();
  PtRecH_3->Write();
  PtRecH_4->Write();
  SavingFile->Close();

}
