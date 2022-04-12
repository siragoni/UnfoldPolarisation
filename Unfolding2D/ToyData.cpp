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
Double_t calcChi2(TH1* h1, TH1* h2, int nbins) {
    Int_t nb = h1->GetNbinsX();
    std::cout << "chi2 Nbins:" << nb << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nbins+1; i++) {
        Double_t b1 =h1->GetBinContent(i);
        Double_t b2 =h2->GetBinContent(i);
        //std::cout << "b1:" << b1 << " b2:" << b2 << std::endl;
        chi2 += (b1-b2) *(b1-b2)/b1;
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
  // Generate flat phi
  // -------------------
  TH1F* generatedH_few = new TH1F("generatedH_few", "generatedH_few", 24, 0., 2.*TMath::Pi());
  TH1F* generatedH_lot = new TH1F("generatedH_lot", "generatedH_lot", 24, 0., 2.*TMath::Pi());
  // TF1 *fun = new TF1("fun","sin(x)",0.,2.*TMath::Pi() );
  TF1 *fun = new TF1("fun","1/(2.*TMath::Pi())",0.,2.*TMath::Pi() );
  for(int i =0; i<(  400.*24.); i++) generatedH_few->Fill(fun->GetRandom());
  for(int i =0; i<(40000.*24.); i++) generatedH_lot->Fill(fun->GetRandom());
  // ===================
  // Retrieve the response
  // matrices made with
  // STARlight
  // -------------------
  TFile* f[24];
  TFile* f2[24];
  TFile* fData[24];
  RooUnfoldResponse*  h[24];
  TH1F* hdata[24];
  TFile* cheating = new TFile(Form("RomanExercise10_%d.root", Stripe));
  TH1F* prior = (TH1F*) cheating->Get("response;1");
  new TCanvas;
  prior -> Draw();
  Double_t Determinants[24];
  for (Int_t i = 5; i < 20; i++) {
    f[i]     = new TFile(Form("Unfolding2D/UnfoldedClosureHe_%d.root", i));
    f2[i]    = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    h[i]     = (RooUnfoldResponse*) f2[i]->Get("response;1");
    // h[i] ->SetPriors(prior);
    fData[i] = new TFile(Form("SignalExtraction/RawYieldsHeV3_%d.root", i));
    hdata[i] = (TH1F*)fData[i]->Get(Form("h_%d", i));
    // unfold[i]  = RooUnfoldBayes(h[i], histo[i], 1);
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
  TH1F* gener_real   = (TH1F*)f[14]->Get("PhiGenBinsH_14");
  TH1F* recon_real   = (TH1F*)f[14]->Get("PhiRecBinsH_14");
  TH1D* hfolded_real = (TH1D*)h[14]->ApplyToTruth(gener_real);
  TH1D* hfolded_few  = (TH1D*)h[14]->ApplyToTruth(generatedH_few);
  TH1D* hfolded_lot  = (TH1D*)h[14]->ApplyToTruth(generatedH_lot);
  RooUnfoldBayes unfold_real[20];
  RooUnfoldBayes unfold_data[50];
  RooUnfoldBayes unfold_real2[20];
  RooUnfoldBayes unfold_few[20];
  RooUnfoldBayes unfold_lot[20];
  TH1D* hunfold_data[20];
  TH1D* hunfold_real[50];
  TH1D* hunfold_real2[20];
  TH1D* hunfold_few[20];
  TH1D* hunfold_lot[20];
  Double_t chi2_real[20];
  Double_t chi2_real2[20];
  Double_t chi2_few[20];
  Double_t chi2_lot[20];
  for (Int_t i = 0; i < 50; i++) {
    // unfold_data[i]  = RooUnfoldBayes(h[Stripe], hdata_formatted,i+1);
    if (i == 49) {
      unfold_data[i]  = RooUnfoldBayes(h[Stripe], hdata_formatted,200);
    } else {
      unfold_data[i]  = RooUnfoldBayes(h[Stripe], hdata_formatted,i+1);
    }
    // unfold_data[i].SetPriors(prior);
    // unfold_real[i]  = RooUnfoldBayes(h[14], hfolded_real,i+1);
    // unfold_real2[i] = RooUnfoldBayes(h[14], recon_real,  i+1);
    // unfold_few[i]   = RooUnfoldBayes(h[14], hfolded_few, i+1);
    // unfold_lot[i]   = RooUnfoldBayes(h[14], hfolded_lot, i+1);
    hunfold_data[i] = (TH1D*) unfold_data[i].Hreco();
    // hunfold_real[i] = (TH1D*) unfold_real[i].Hreco();
    // hunfold_real2[i]= (TH1D*) unfold_real2[i].Hreco();
    // hunfold_few[i]  = (TH1D*) unfold_few[i].Hreco();
    // hunfold_lot[i]  = (TH1D*) unfold_few[i].Hreco();
    // chi2_real[i]    = calcChi2(gener_real,     hunfold_real[i], 24);
    // chi2_real2[i]   = calcChi2(gener_real,     hunfold_real2[i],24);
    // chi2_few[i]     = calcChi2(generatedH_few, hunfold_few[i],  24);
    // chi2_lot[i]     = calcChi2(generatedH_few, hunfold_lot[i],  24);
  }
  // for (Int_t i = 0; i < 20; i++) {
  //   std::cout << "============================================================="           << std::endl;
  //   std::cout << "chi2 gen-unfold, real R comp generated [" << i << "] = " << chi2_real[i] << std::endl;
  //   std::cout << "chi2 gen-unfold, real reconstructed    [" << i << "] = " << chi2_real2[i]<< std::endl;
  //   std::cout << "chi2 gen-unfold, few                   [" << i << "] = " << chi2_few[i]  << std::endl;
  //   std::cout << "chi2 gen-unfold, lot                   [" << i << "] = " << chi2_lot[i]  << std::endl;
  // }



  // =============================
  // Refolded exercise
  // =============================
  TH1D* refolded_data[50];
  Double_t chi2_data[50];
  for (size_t i = 0; i < 50; i++) {
    refolded_data[i] = (TH1D*)h[Stripe]->ApplyToTruth(hunfold_data[i]);
    chi2_data[i]     = calcChi2(hdata_formatted,    refolded_data[i], M);
  }
  for (Int_t i = 0; i < 50; i++) {
    std::cout << "============================================================="           << std::endl;
    std::cout << "chi2 recon-refolded[" << i << "] = " << chi2_data[i] << std::endl;
  }








  TFile* file = new TFile(Form("RomanExercise10_%d_3.root", Stripe), "recreate");
  // gener_real     ->Write();
  // generatedH_few ->Write();
  // generatedH_lot ->Write();
  // recon_real     ->Write();
  // hfolded_real   ->Write();
  // hfolded_few    ->Write();
  // hfolded_lot    ->Write();
  // for (Int_t i = 0; i < 20; i++) {
  //   // hunfold_real[i]->Write();
  //   hunfold_real2[i]->Write();
  //   // hunfold_few[i] ->Write();
  //   // hunfold_lot[i] ->Write();
  // }
  hunfold_data[49]->Write();
  hunfold_data[3]->Write();
  file->Close();
















  // //===================
  // // Final check on
  // // chi square
  // // ------------------
  // TH1F* gener_all[24];
  // TH1F* recon_all[24];
  // for (Int_t i = 5; i < 19; i++) {
  //   gener_all[i] = (TH1F*)f[i]->Get(Form("PhiGenBinsH_%d", i));
  //   recon_all[i] = (TH1F*)f[i]->Get(Form("PhiRecBinsH_%d", i));
  // }
  // RooUnfoldBayes unfold_all[24];
  // TH1D* hunfold_all[24];
  // Double_t chi2_all[24];
  // Double_t chi2_total = 0.;
  // for (Int_t i = 5; i < 19; i++) {
  //   unfold_all[i]  = RooUnfoldBayes(h[i], recon_all[i],10);
  //   hunfold_all[i] = (TH1D*) unfold_all[i].Hreco();
  //   Int_t nbins    = recon_all[i]->GetNbinsX();
  //   chi2_all[i]    = calcChi2(gener_all[i], hunfold_all[i], nbins);
  //   chi2_total    += chi2_all[i];
  // }
  // cout << "=========================================" << endl;
  // cout << "++ All Chi2 from closure distributions ++" << endl;
  // cout << "-----------------------------------------" << endl;
  // for (Int_t i = 5; i < 20; i++) {
  //   cout << "Partial Chi2[" << i << "] = " << chi2_all[i] << endl;
  // }
  // cout << "Total Chi2 = " << chi2_total << endl;



  // for (Int_t i = 5; i < 20; i++) {
  //   f[i]     ->Close();
  //   f2[i]    ->Close();
  //   fData[i] ->Close();
  //   delete h[i];
  // }

}
