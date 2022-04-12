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

//______________________________________________
Double_t helicity2Dv3(Double_t *x, Double_t *par) {

  Double_t CosSquaredTheta      = x[1] * x[1];
  Double_t partialTheta         = 1. + par[0] * CosSquaredTheta;

  Double_t CosOfTwoPhi          = TMath::Cos( 2. * x[0] );
  Double_t SinSquaredTheta      = 1. - CosSquaredTheta;
  Double_t partialPhi           = par[1] * SinSquaredTheta * CosOfTwoPhi;

  Double_t CosPhi               = TMath::Cos( x[0] );
  // Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredOfTwoTheta = 4. * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t partialMix           = par[2] * SinOfTwoTheta * CosPhi;


  Double_t sumOfTheSubFits      = partialTheta + partialPhi + partialMix;
  Double_t FinalResult          = par[3] * sumOfTheSubFits / ( 3. + par[0] );
  return   FinalResult;

}
//______________________________________________
TF2* Model5 = new TF2("Model5", helicity2Dv3,0., 2.*TMath::Pi(), -1. ,1., 4 );
// TF2* Model5 = new TF2("Model5", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
//______________________________________________
Double_t IntegralFunction(Double_t *x, Double_t *par) {

  Model5->SetParameter(0, par[0]);
  Model5->SetParameter(1, par[1]);
  Model5->SetParameter(2, par[2]);
  Model5->SetParameter(3, par[3]);


  Double_t iPhiBins      = x[0];
  Double_t iCosThetaBins = x[1];
  Double_t M = 6.;
  Double_t IntegralFunction = Model5->Integral(
                                      (((Double_t) iPhiBins     )*2.*TMath::Pi()/M),
                                      (((Double_t) iPhiBins + 1.)*2.*TMath::Pi()/M),
                                      (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)),
                                      (-1.+((Double_t) iCosThetaBins + 1.)*(0.08+0.01/3.))
                                     );
  return   IntegralFunction;

}

//_____________________________________________________________________________
void ParseMC_modulationphi(Double_t LPhi = 0., Double_t Threshold = 10000., Double_t LTheta = 0.){
// ParseMC_modulationphi(1., 100, 0.)
  TFile* file   = new TFile("AnalysisResultsCohWithSinglePt.root", "READ");
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
  TH1F *PhiRecH      = new TH1F("PhiRecH",      "PhiRecH",      250, 0., 2.*TMath::Pi());
  TH1F *PhiGenH      = new TH1F("PhiGenH",      "PhiGenH",      250, 0., 2.*TMath::Pi());
  TH1F *TildePhiRecH = new TH1F("TildePhiRecH", "TildePhiRecH", 25, 0., 2.*TMath::Pi());
  TH1F *TildePhiGenH = new TH1F("TildePhiGenH", "TildePhiGenH", 25, 0., 2.*TMath::Pi());


  // TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -1. ,1., 4 );
  Model->FixParameter(0, LTheta);
  Model->FixParameter(1, LPhi);
  Model->FixParameter(2, 0.);
  Model->FixParameter(3, Threshold);


  TH2F* generatedlimit = new TH2F("generatedlimit", "generatedlimit", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  TH2F* generatedlimit2 = new TH2F("generatedlimit2", "generatedlimit2", 450, 0., 2.*TMath::Pi(), 450, -1., 1.);
  // 0.99 LPhi with 400, 400, thre /= 100
  Double_t PhiCenters2[24];
  Double_t Spacing = TMath::Pi()/24.;
  for (size_t i = 0; i < 24; i++) {
    PhiCenters2[i] = (2.*((Double_t) i) + 1.) *Spacing;
    cout << "Phi centres = " << PhiCenters2[i] << endl;

  }
  Double_t CosThetaCenters2[24];
  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters2[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
    cout << "CosTheta centres = " << CosThetaCenters2[i] << endl;

  }




  TH1F*** InvMassH;
  InvMassH = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassH[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j] = new TH1F(Form("InvMassH_%d_%d", i, j),Form("InvMassH_%d_%d", i, j),2000, 0, 20);
    }
  }


  TH1F*** InvMassHg;
  InvMassHg = new TH1F**[24];
  for (Int_t i = 0; i < 24; i++) {
    InvMassHg[i] = new TH1F*[24];
    for (Int_t j = 0; j < 24; j++) {
      InvMassHg[i][j] = new TH1F(Form("InvMassHg_%d_%d", i, j),Form("InvMassHg_%d_%d", i, j),2000, 0, 20);
    }
  }







  Double_t controlFlag  = 0;
  Double_t controlFlag2 = 0;
  Double_t controlFlag3 = 0;

  Int_t binx = -1;
  Double_t Flag = -999.;


  //read all entries and fill the histograms
  for (Int_t i=0; i<nentriesRec; i++) {
    ReconTree->GetEntry(i);
    GenerTree->GetEntry(i);










      if( (pTgen > 0. && pTgen < 0.25) && (Ygen > -4. && Ygen < -2.5) ) {

        CosThetaHEgen2 = -1.*CosThetaHEgen;
        if((CosThetaHEgen2 < 1.0) && (CosThetaHEgen2 > -1.0)){
        PhiHEgen2 = PhiHEgen+TMath::Pi();









        Int_t biny2 = generatedlimit2->GetYaxis()->FindBin(CosThetaHEgen2);
        Int_t binx2 = generatedlimit2->GetXaxis()->FindBin(PhiHEgen2);
        Flag = -999.;
        if ( generatedlimit2->GetBinContent(binx2,biny2) < ((Threshold/90.) * Model->Eval(PhiHEgen2, CosThetaHEgen2))  ){
          generatedlimit->Fill( PhiHEgen2, CosThetaHEgen2);
          generatedlimit2->Fill( PhiHEgen2, CosThetaHEgen2);
          Flag = 1.;
        }

        if( Flag > 0.5){

          //================================
          // Translating Phi and CosTheta
          // into bin numbers
          //--------------------------------
          Double_t TraslatedCosThetaGeng = 0.5*(CosThetaHEgen2 + 1.)*24.;
          Double_t iCosThetaBins2g = -1;
          Double_t RemainderCosThetag = 0.;
          RemainderCosThetag = modf(TraslatedCosThetaGeng, &iCosThetaBins2g);
          Int_t iCosThetaBinsg = (Int_t)  iCosThetaBins2g;
          //--------------------------------
          // Binning in phi depends
          // on CosTheta
          //--------------------------------
          Double_t Mg = 6.;
          // Double_t TraslatedPhiGen = (PhiHEgen2+TMath::Pi())*M/(2.*TMath::Pi());
          Double_t TraslatedPhiGeng = (PhiHEgen2)*Mg/(2.*TMath::Pi());
          Double_t iPhiBins2g = -1;
          Double_t RemainderPhig = 0.;
          RemainderPhig = modf(TraslatedPhiGeng, &iPhiBins2g);
          Int_t iPhiBinsg = (Int_t)  iPhiBins2g;
          //--------------------------------------
          InvMassHg[iCosThetaBinsg][iPhiBinsg]->Fill(Mrec);


        }







          if( (pTrec > 0. && pTrec < 0.25) && (Yrec > -4. && Yrec < -2.5) && (Flag > 0.5) ) {

                  if(PhiHErec < 0.) PhiHErec = PhiHErec + 2.*TMath::Pi();


                  //================================
                  // Translating Phi and CosTheta
                  // into bin numbers
                  //--------------------------------
                  Double_t TraslatedCosThetaGen = 0.5*(CosThetaHEgen2 + 1.)*24.;
                  Double_t iCosThetaBins2 = -1;
                  Double_t RemainderCosTheta = 0.;
                  RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
                  Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;
                  //--------------------------------
                  // Binning in phi depends
                  // on CosTheta
                  //--------------------------------
                  Double_t M = 6.;
                  // Double_t TraslatedPhiGen = (PhiHEgen2+TMath::Pi())*M/(2.*TMath::Pi());
                  Double_t TraslatedPhiGen = (PhiHErec)*M/(2.*TMath::Pi());
                  Double_t iPhiBins2 = -1;
                  Double_t RemainderPhi = 0.;
                  RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
                  Int_t iPhiBins = (Int_t)  iPhiBins2;
                  //--------------------------------------
                  InvMassH[iCosThetaBins][iPhiBins]->Fill(Mrec);








          }
      }
    }
  }













  TFile *SavingFile = new TFile("ScanModulation/LPhi005/ModulationYields/CohHe_phimodulation.root", "RECREATE");
  for (Int_t i = 0; i < 24; i++) {
    for (Int_t j = 0; j < 24; j++) {
      InvMassH[i][j]->Write();
    }
  }
  for (Int_t i = 0; i < 24; i++) {
    for (Int_t j = 0; j < 24; j++) {
      InvMassHg[i][j]->Write();
    }
  }
  SavingFile       ->Close();














  Int_t Counter[24] = {6,6,6,6,6, 6,6, 12,12, 24,24,24,24,24,24, 12,12, 6,6, 6,6,6,6,6 };

  Double_t CosThetaCenters[24];
  Double_t PhiCenters[24];

  PhiCenters[0] = PhiCenters[1] = PhiCenters[2]  = PhiCenters[3]  = PhiCenters[4] = PhiCenters[23] = PhiCenters[22] = PhiCenters[21] = PhiCenters[20] = PhiCenters[19] = TMath::Pi();
  Double_t Spacing1 = TMath::Pi()/12.;
  PhiCenters[5] = PhiCenters[6] = PhiCenters[18] = PhiCenters[17] = Spacing1;
  Double_t Spacing2 = TMath::Pi()/24.;
  PhiCenters[7] = PhiCenters[8] = PhiCenters[16] = PhiCenters[15] = Spacing2;
  Double_t Spacing3 = TMath::Pi()/48.;
  PhiCenters[9] = PhiCenters[10] = PhiCenters[11] = PhiCenters[12] = PhiCenters[13] = PhiCenters[14] = Spacing3;


  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
  }

  TH1F* RawYields[24];
  for (Int_t i = 0; i < 24; i++) {
    RawYields[i] = new TH1F(Form("h_%d", i), Form("h_%d", i), 100, -0.5, 99.5);
  }
  TH1F* RawYieldsg[24];
  for (Int_t i = 0; i < 24; i++) {
    RawYieldsg[i] = new TH1F(Form("hg_%d", i), Form("hg_%d", i), 100, -0.5, 99.5);
  }
  TFile* parsedMC = new TFile("ScanModulation/LPhi005/ModulationYields/CohHe_phimodulation.root");

  Double_t JPsiPeakValue    = 0;
  Double_t JPsiPeakValueErr = 0;
  Double_t BkgValue         = 0;
  Double_t BkgValueError    = 0;
  Double_t JPsiPeakValue2    = 0;
  Double_t JPsiPeakValueErr2 = 0;
  Double_t BkgValue2         = 0;
  Double_t BkgValueError2    = 0;
  TH1F* fCohJpsiToMu        = 0x0;
  TH1F* fCohJpsiToMu2        = 0x0;
  for (Int_t iCosThetaBins = 4; iCosThetaBins < 20; iCosThetaBins++) {

    for (Int_t iPhiBins = 0; iPhiBins < Counter[iCosThetaBins]; iPhiBins++) {

      JPsiPeakValue    = 0;
      JPsiPeakValueErr = 0;
      BkgValue         = 0;
      BkgValueError    = 0;
      JPsiPeakValue2    = 0;
      JPsiPeakValueErr2 = 0;
      BkgValue2         = 0;
      BkgValueError2    = 0;

      fCohJpsiToMu     = (TH1F*)parsedMC->Get( Form( "InvMassH_%d_%d", iCosThetaBins, iPhiBins ) );
      JPsiPeakValue    = (Double_t) fCohJpsiToMu->GetEntries();
      JPsiPeakValueErr = TMath::Sqrt((Double_t) fCohJpsiToMu->GetEntries());
      fCohJpsiToMu2     = (TH1F*)parsedMC->Get( Form( "InvMassHg_%d_%d", iCosThetaBins, iPhiBins ) );
      JPsiPeakValue2    = (Double_t) fCohJpsiToMu2->GetEntries();
      JPsiPeakValueErr2 = TMath::Sqrt((Double_t) fCohJpsiToMu2->GetEntries());

      RawYields[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue   );
      RawYields[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr);
      RawYieldsg[iCosThetaBins]->Fill(       iPhiBins,   JPsiPeakValue2   );
      RawYieldsg[iCosThetaBins]->SetBinError(iPhiBins+1, JPsiPeakValueErr2);

    }
    TFile f(Form("ScanModulation/LPhi005/ModulationYields/MonteCarloYieldsHe_phimodulation_%d.root", iCosThetaBins),   "recreate");
    RawYields[iCosThetaBins]->Write();
    RawYieldsg[iCosThetaBins]->Write();
    generatedlimit->Write();
    generatedlimit2->Write();
    f.Close();

  }





}
