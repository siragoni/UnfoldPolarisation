#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;
#include <math.h>
#include "TH2D.h"
#include "TF2.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TVirtualFitter.h"
#include "TList.h"

#include <vector>
#include <map>


//_____________________________________________________________________________
Double_t calcChi2(TH1* h1, TH1* h2, Int_t nbins) {
    // Int_t nb = h1->GetNbinsX();
    Int_t nb = nbins;
    std::cout << "chi2 Nbins:" << nb << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nb+1; i++) {
        Double_t b1  = h1->GetBinContent(i);
        Double_t Eb1 = h1->GetBinError(i);
        Double_t Eb2 = h2->GetBinError(i);
        Double_t b2  = h2->GetBinContent(i);
        cout << " partial = " << (b1-b2)*(b1-b2)/Eb1 << endl;
        chi2 += (b1-b2)*(b1-b2)/Eb1;
        // chi2 += (b1-b1)*(b2-b1)/Eb2;
    }
    return chi2;
}
Double_t calcChi2(TH1* h1, TH1* h2, Int_t nbins, Double_t* container) {
  // Int_t nb = h1->GetNbinsX();
  Int_t nb = nbins;
  std::cout << "chi2 Nbins:" << nb << std::endl;
  Double_t chi2 = 0;
  for(int i = 1; i < nb+1; i++) {
      Double_t b1  = h1->GetBinContent(i);
      Double_t Eb1 = h1->GetBinError(i);
      Double_t Eb2 = h2->GetBinError(i);
      Double_t b2  = h2->GetBinContent(i);
      cout << " partial = " << (b1-b2)*(b1-b2)/Eb1 << endl;
      chi2 += (b1-b2)*(b1-b2)/Eb1;
      container[i-1] = (b1-b2)*(b1-b2)/Eb1;
      // chi2 += (b1-b1)*(b2-b1)/Eb2;
  }
  return chi2;
}
//_____________________________________________________________________________
/* - Toy MC
   - Generate flat phi
   - Use STARlight full
   - response matrix
   -
 */
void ToyData(Int_t Stripe = 14){

  // ===================
  // Retrieve REAL data
  // -------------------
  TFile* fData[24];
  TH1F* hdata[24];
  for (Int_t i = 5; i < 20; i++) {
    fData[i] = new TFile(Form("SignalExtraction/RawYieldsHeV3_%d.root", i));
    hdata[i] = (TH1F*)fData[i]->Get(Form("h_%d", i));
  }
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight
  // -------------------
  TFile* f[24];
  RooUnfoldResponse*  h[24];
  for (Int_t i = 5; i < 20; i++) {
    f[i]    = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f[i]->Get("response;1");
  }
  Int_t M = 1;
  if( Stripe == 0  || Stripe == 1  || Stripe == 2  || Stripe == 3  ||
      Stripe == 4  || Stripe == 23 || Stripe == 22 || Stripe == 21 ||
      Stripe == 20 || Stripe == 19 )
  {
    M = 1;
  } else if ( Stripe == 5  || Stripe == 6  || Stripe == 18  || Stripe == 17 ) {
    M = 6;
  } else if ( Stripe == 7  || Stripe == 8  || Stripe == 16  || Stripe == 15 ) {
    M = 12;
  } else if ( Stripe == 9  || Stripe == 10  || Stripe == 11  ||
              Stripe == 12 || Stripe == 13  || Stripe == 14 ) {
    M = 24;
  }
  TH1F* hdata_formatted = new TH1F("hdata_formatted","hdata_formatted", M, 0, 2.*TMath::Pi());
  for (Int_t i = 1; i < M+1; i++) {
    hdata_formatted->SetBinContent(i, hdata[Stripe]->GetBinContent(i));
    hdata_formatted->SetBinError(  i, hdata[Stripe]->GetBinError(i));
  }

  // ===================
  // Folding the model with
  // the response matrix
  // from STARlight and
  // unfolding it
  // -------------------
  RooUnfoldBayes unfold[50];
  TH1D* hunfold_default[50];
  TH1D* hunfold_kerrors[50];
  TH1D* hunfold_kcovariance[50];
  TH1D* hunfold_kcovtoy[50];
  TH1D* hunfold_manual[50];
  TH1D* runfold_default[50];
  TH1D* runfold_kerrors[50];
  TH1D* runfold_kcovariance[50];
  TH1D* runfold_kcovtoy[50];
  TH1D* runfold_manual[50];
  Double_t chi2_default[50];
  Double_t chi2_kerrors[50];
  Double_t chi2_kcovariance[50];
  Double_t chi2_kcovtoy[50];
  Double_t chi2_manual[50];
  Double_t** partials_chi2_default;
  partials_chi2_default = new Double_t*[50];
  for (Int_t i = 0; i < 50; i++) {
    unfold[i]              = RooUnfoldBayes(h[Stripe], hdata_formatted,i+1);
    hunfold_default[i]     = (TH1D*) unfold[i].Hreco();
    hunfold_kerrors[i]     = (TH1D*) unfold[i].Hreco(RooUnfold::kErrors);
    hunfold_kcovariance[i] = (TH1D*) unfold[i].Hreco(RooUnfold::kCovariance);
    hunfold_kcovtoy[i]     = (TH1D*) unfold[i].Hreco(RooUnfold::kCovToy);
    hunfold_manual[i]      = (TH1D*) unfold[i].Hreco();
    runfold_default[i]     = (TH1D*) h[Stripe]->ApplyToTruth(hunfold_default[i]);
    runfold_kerrors[i]     = (TH1D*) h[Stripe]->ApplyToTruth(hunfold_kerrors[i]);
    runfold_kcovariance[i] = (TH1D*) h[Stripe]->ApplyToTruth(hunfold_kcovariance[i]);
    runfold_kcovtoy[i]     = (TH1D*) h[Stripe]->ApplyToTruth(hunfold_kcovtoy[i]);
    runfold_manual[i]      = (TH1D*) h[Stripe]->ApplyToTruth(hunfold_manual[i]);
    partials_chi2_default[i] = new Double_t[M];
    for (Int_t j = 1; j < M+1; j++) {
      Double_t value1  = runfold_manual[i]->GetBinContent(j);
      Double_t value2  = hunfold_manual[i]->GetBinContent(j);
      Double_t Evalue2 = hunfold_manual[i]->GetBinError(j);
      Double_t Evalue1 = value1 * Evalue2 / value2;
      // cout << "value1  = " << value1  << endl;
      // cout << "value2  = " << value2  << endl;
      // cout << "Evalue2 = " << Evalue2 << endl;
      // cout << "Evalue1 = " << Evalue1 << endl;
      // runfold_manual[i]->SetBinError(j, runfold_manual[i]->GetBinContent(j)*hunfold_manual[i]->GetBinError(j)/hunfold_manual[i]->GetBinContent(j));
      runfold_manual[i]->SetBinError(j, Evalue1);
    }
    std::cout << "chi2  " << i << std::endl;
    // chi2_default[i]     = calcChi2(runfold_default[i],     hdata_formatted, M    );
    chi2_default[i]     = calcChi2(runfold_default[i],     hdata_formatted, M, partials_chi2_default[i]    );
    chi2_kerrors[i]     = calcChi2(runfold_kerrors[i],     hdata_formatted, M    );
    chi2_kcovariance[i] = calcChi2(runfold_kcovariance[i], hdata_formatted, M    );
    chi2_kcovtoy[i]     = calcChi2(runfold_kcovtoy[i],     hdata_formatted, M    );
    chi2_manual[i]      = calcChi2(runfold_manual[i],      hdata_formatted, M    );
    std::cout << "finished chi2  " << i << std::endl;
  }
  for (Int_t i = 0; i < 50; i++) {
    std::cout << "============================================================="     << std::endl;
    std::cout << "chi2 ref-rec, default     [" << i << "] = " << chi2_default[i]     << std::endl;
    std::cout << "chi2 ref-rec, kErrors     [" << i << "] = " << chi2_kerrors[i]     << std::endl;
    std::cout << "chi2 ref-rec, kCovariance [" << i << "] = " << chi2_kcovariance[i] << std::endl;
    std::cout << "chi2 ref-rec, kCovToy     [" << i << "] = " << chi2_kcovtoy[i]     << std::endl;
    std::cout << "chi2 ref-rec, Manual      [" << i << "] = " << chi2_manual[i]      << std::endl;
  }

  TFile* file = new TFile(Form("Errors_%d.root", Stripe), "recreate");
  for (Int_t i = 0; i < 50; i++) {
    hunfold_default[i]->Write();
  }
  file->Close();







  std::cout << "============================================================="     << std::endl;
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial0     [" << i << "] = " << partials_chi2_default[i][0]     << std::endl;
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial1     [" << i << "] = " << partials_chi2_default[i][1]     << std::endl;
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial2     [" << i << "] = " << partials_chi2_default[i][2]     << std::endl;
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial3     [" << i << "] = " << partials_chi2_default[i][3]     << std::endl;
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial4     [" << i << "] = " << partials_chi2_default[i][4]     << std::endl;
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "partial5     [" << i << "] = " << partials_chi2_default[i][5]     << std::endl;
  }
  // for (Int_t i = 0; i < 20; i++) {
  //   std::cout << "partial6     [" << i << "] = " << partials_chi2_default[i][6]     << std::endl;
  // }









}
