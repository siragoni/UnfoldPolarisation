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
void ToyMC(){

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
  RooUnfoldResponse*  h[24];
  Double_t Determinants[24];
  for (Int_t i = 5; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_%d.root", i));
    h[i] = (RooUnfoldResponse*) f[i]->Get("response;1");
    // unfold[i]  = RooUnfoldBayes(h[i], histo[i], 1);
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
  RooUnfoldBayes unfold_real2[20];
  RooUnfoldBayes unfold_few[20];
  RooUnfoldBayes unfold_lot[20];
  TH1D* hunfold_real[20];
  TH1D* hunfold_real2[20];
  TH1D* hunfold_few[20];
  TH1D* hunfold_lot[20];
  Double_t chi2_real[20];
  Double_t chi2_real2[20];
  Double_t chi2_few[20];
  Double_t chi2_lot[20];
  for (Int_t i = 0; i < 20; i++) {
    unfold_real[i]  = RooUnfoldBayes(h[14], hfolded_real,i+1);
    unfold_real2[i] = RooUnfoldBayes(h[14], recon_real,  i+1);
    unfold_few[i]   = RooUnfoldBayes(h[14], hfolded_few, i+1);
    unfold_lot[i]   = RooUnfoldBayes(h[14], hfolded_lot, i+1);
    hunfold_real[i] = (TH1D*) unfold_real[i].Hreco();
    hunfold_real2[i]= (TH1D*) unfold_real2[i].Hreco();
    hunfold_few[i]  = (TH1D*) unfold_few[i].Hreco();
    hunfold_lot[i]  = (TH1D*) unfold_few[i].Hreco();
    chi2_real[i]    = calcChi2(gener_real,     hunfold_real[i], 24);
    chi2_real2[i]   = calcChi2(gener_real,     hunfold_real2[i],24);
    chi2_few[i]     = calcChi2(generatedH_few, hunfold_few[i],  24);
    chi2_lot[i]     = calcChi2(generatedH_few, hunfold_lot[i],  24);
  }
  for (Int_t i = 0; i < 20; i++) {
    std::cout << "============================================================="           << std::endl;
    std::cout << "chi2 gen-unfold, real R comp generated [" << i << "] = " << chi2_real[i] << std::endl;
    std::cout << "chi2 gen-unfold, real reconstructed    [" << i << "] = " << chi2_real2[i]<< std::endl;
    std::cout << "chi2 gen-unfold, few                   [" << i << "] = " << chi2_few[i]  << std::endl;
    std::cout << "chi2 gen-unfold, lot                   [" << i << "] = " << chi2_lot[i]  << std::endl;
  }

  TFile* file = new TFile("RomanExercise.root", "recreate");
  gener_real     ->Write();
  generatedH_few ->Write();
  generatedH_lot ->Write();
  recon_real     ->Write();
  hfolded_real   ->Write();
  hfolded_few    ->Write();
  hfolded_lot    ->Write();
  for (Int_t i = 0; i < 20; i++) {
    // hunfold_real[i]->Write();
    hunfold_real2[i]->Write();
    // hunfold_few[i] ->Write();
    // hunfold_lot[i] ->Write();
  }
  file->Close();



  // new TCanvas;
  // gener_real       ->SetLineColor(kRed);
  // hunfold_real2[10]->SetLineColor(kBlue);
  // hunfold_real2[10]->SetLineWidth(5);
  // // hunfold_real2[20]->SetLineColor(kGreen);
  // gener_real       ->Draw();
  // hunfold_real2[10]->Draw("same");
  // hunfold_real2[20]->Draw("same");



  cout << "=====================" << endl;
  cout << "Integral gen real = "  << gener_real    ->Integral() << endl;
  cout << "Integral gen real = "  << generatedH_few->Integral() << endl;
  cout << "Integral gen real = "  << generatedH_lot->Integral() << endl;
  cout << "---------------------" << endl;
  cout << "Integral RECON  real = "  << recon_real  ->Integral() << endl;
  cout << "Integral R comp real = "  << hfolded_real->Integral() << endl;
  cout << "Integral R comp few  = "  << hfolded_few ->Integral() << endl;
  cout << "Integral R comp lot  = "  << hfolded_lot ->Integral() << endl;







  // cout << "---------------------" << endl;
  // TH1F* generated_parser = (TH1F*)f[14]->Get("histo3_14");
  // TH1F* unfolded_closure = (TH1F*)f[14]->Get("histo4_14");
  // cout << "Integral generated saved  = "  << generated_parser ->Integral() << endl;
  // cout << "Integral unfolded saved   = "  << unfolded_closure ->Integral() << endl;



  //===================
  // Final check on
  // chi square
  // ------------------
  TH1F* gener_all[24];
  TH1F* recon_all[24];
  for (Int_t i = 5; i < 19; i++) {
    gener_all[i] = (TH1F*)f[i]->Get(Form("PhiGenBinsH_%d", i));
    recon_all[i] = (TH1F*)f[i]->Get(Form("PhiRecBinsH_%d", i));
  }
  RooUnfoldBayes unfold_all[24];
  TH1D* hunfold_all[24];
  Double_t chi2_all[24];
  Double_t chi2_total = 0.;
  for (Int_t i = 5; i < 19; i++) {
    unfold_all[i]  = RooUnfoldBayes(h[i], recon_all[i],10);
    hunfold_all[i] = (TH1D*) unfold_all[i].Hreco();
    Int_t nbins    = recon_all[i]->GetNbinsX();
    chi2_all[i]    = calcChi2(gener_all[i], hunfold_all[i], nbins);
    chi2_total    += chi2_all[i];
  }
  cout << "=========================================" << endl;
  cout << "++ All Chi2 from closure distributions ++" << endl;
  cout << "-----------------------------------------" << endl;
  for (Int_t i = 5; i < 20; i++) {
    cout << "Partial Chi2[" << i << "] = " << chi2_all[i] << endl;
  }
  cout << "Total Chi2 = " << chi2_total << endl;




}
