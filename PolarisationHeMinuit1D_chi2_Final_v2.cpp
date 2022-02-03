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

#include "TDatime.h"

#include <vector>
#include <map>

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
void PolarisationHeMinuit1D( Int_t SignalRangeSelectionMode = 0, Int_t FitRangeMode = 0 ){
  #ifdef __CINT__
  // if (!TClass::GetDict("RooUnfold")) gSystem->Load("../RooUnfold/libRooUnfold");
    gSystem->Load("../RooUnfold/libRooUnfold");
  #endif

  SignalRangeModeFromBash = SignalRangeSelectionMode;
  TDatime d;
  // TFile* file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  TFile* file1D = 0x0;
  if        ( SignalRangeSelectionMode == 0 ) {
    file1D = new TFile("SavingFileData.root");
    // file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_long.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 1 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1DresultsV2/PolarisationCorrectedHe1Dv2_long.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 2 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_2.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 3 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_3.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 4 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_4.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else if ( SignalRangeSelectionMode == 5 ) {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D_5.root", d.GetYear(), d.GetMonth(), d.GetDay() ) );
  } else {
    file1D = new TFile(Form("pngResults/%d-%2.2d-%2.2d/1Dresults/PolarisationCorrectedHe1D.root",   d.GetYear(), d.GetMonth(), d.GetDay() ) );
  }
  TH1D* CorrectedCosTheta = (TH1D*) file1D->Get("response;5");
  TH1D* CorrectedPhi      = (TH1D*) file1D->Get("response;4");
  TH1D* CorrectedTildePhi = (TH1D*) file1D->Get("response;6");


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
  // gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.2e+04, 3.2e+04);
  gMinuit->DefineParameter(1, "NormalTheta", 0, 0,  -1., 3.2e+04);
  // gMinuit->DefineParameter(2, "NormalisPhi",      8300, 100,  7700, 8700);
  gMinuit->DefineParameter(2, "NormalisPhi",      0, 0,  0, 8700);
  // gMinuit->DefineParameter(3, "LambdaPhi",           0, 0.01,    -2, 2   );
  // gMinuit->DefineParameter(3, "LambdaPhi",           0, 0.01,    -2, 2   );
  gMinuit->DefineParameter(3, "LambdaPhi",           0.01, 0.00,    -2,2   );
  gMinuit->DefineParameter(4, "NormalisTildePhi", 8300, 100,  7700, 8700);
  // gMinuit->DefineParameter(4, "NormalisTildePhi", 3000, 100,  2000, 4000);
  // gMinuit->DefineParameter(5, "LambdaThetaPhi",      0, 0.01,    -2, 2   );
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
  CorrectedCosTheta->GetYaxis()->SetRangeUser(0., CorrectedCosTheta->GetMaximum()*2);
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
  CorrectedPhi->GetYaxis()->SetRangeUser(0., CorrectedPhi->GetMaximum()*2);
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
  CorrectedTildePhi->GetYaxis()->SetRangeUser(0., CorrectedPhi->GetMaximum()*2);
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

  // TFile SavingFile( Form("pngResults/%d-%2.2d-%2.2d/1Dresults/Parameters_SigEx_%d_FitRange_%d_HE.root", d.GetYear(), d.GetMonth(), d.GetDay(), SignalRangeSelectionMode, FitRangeMode), "recreate" );
  // TH1F* SavingParamH = new TH1F( "SavingParamH", "SavingParamH", 10, 0, 10 );
  // SavingParamH->SetBinContent( 1, LambdaTheta );
  // SavingParamH->SetBinContent( 2, LambdaPhi );
  // SavingParamH->SetBinContent( 3, LambdaThetaPhi );
  // SavingParamH->SetBinContent( 6, LambdaThetaErr );
  // SavingParamH->SetBinContent( 7, LambdaPhiErr );
  // SavingParamH->SetBinContent( 8, LambdaThetaPhiErr );
  // SavingParamH->SetBinError( 1, 0 );
  // SavingParamH->SetBinError( 2, 0 );
  // SavingParamH->SetBinError( 3, 0 );
  // SavingParamH->SetBinError( 6, 0 );
  // SavingParamH->SetBinError( 7, 0 );
  // SavingParamH->SetBinError( 8, 0 );
  // SavingParamH     ->Write();
  // CorrectedCosTheta->Write();
  // CorrectedPhi     ->Write();
  // CorrectedTildePhi->Write();
  // SavingFile.Close();












  // CONTOUR PLOTS
  //
  //
  // // gMinuit->GetParameter(0, LambdaTheta,      LambdaThetaErr     );
  // // gMinuit->GetParameter(3, LambdaPhi,        LambdaPhiErr       );
  // // gMinuit->GetParameter(5, LambdaThetaPhi,   LambdaThetaPhiErr  );
  // new TCanvas;
  // TFile ContourFile( "pngResults/ContourFileHeV2.root", "recreate" );
  // gMinuit->SetErrorDef(9); //note 4 and not 2!
  // TGraph *LThetaVsLPthree = (TGraph*)gMinuit->Contour(80,0,3);
  // LThetaVsLPthree->SetName("LThetaVsLPthree");
  // LThetaVsLPthree->SetTitle("LThetaVsLPthree");
  // LThetaVsLPthree->SetFillColor(42);
  // LThetaVsLPthree->Write();
  // gMinuit->SetErrorDef(4); //note 4 and not 2!
  // TGraph *LThetaVsLPtwo = (TGraph*)gMinuit->Contour(80,0,3);
  // LThetaVsLPtwo->SetName("LThetaVsLPtwo");
  // LThetaVsLPtwo->SetTitle("LThetaVsLPtwo");
  // LThetaVsLPtwo->SetFillColor(42);
  // LThetaVsLPtwo->Write();
  // /*Get contour for parameter 0 versus parameter 2 for ERRDEF=1*/
  // gMinuit->SetErrorDef(1);
  // TGraph *LThetaVsLPone = (TGraph*)gMinuit->Contour(80,0,3);
  // LThetaVsLPone->SetName("LThetaVsLPone");
  // LThetaVsLPone->SetTitle("LThetaVsLPone");
  // LThetaVsLPone->SetFillColor(38);
  // LThetaVsLPone->Write();
  //
  //
  //
  // gMinuit->SetErrorDef(9); //note 4 and not 2!
  // TGraph *LThetaVsLTPthree = (TGraph*)gMinuit->Contour(80,0,5);
  // LThetaVsLTPthree->SetName("LThetaVsLTPthree");
  // LThetaVsLTPthree->SetTitle("LThetaVsLTPthree");
  // LThetaVsLTPthree->SetFillColor(42);
  // LThetaVsLTPthree->Write();
  // gMinuit->SetErrorDef(4); //note 4 and not 2!
  // TGraph *LThetaVsLTPtwo = (TGraph*)gMinuit->Contour(80,0,5);
  // LThetaVsLTPtwo->SetName("LThetaVsLTPtwo");
  // LThetaVsLTPtwo->SetTitle("LThetaVsLTPtwo");
  // LThetaVsLTPtwo->SetFillColor(42);
  // LThetaVsLTPtwo->Write();
  // /*Get contour for parameter 0 versus parameter 2 for ERRDEF=1*/
  // gMinuit->SetErrorDef(1);
  // TGraph *LThetaVsLTPone = (TGraph*)gMinuit->Contour(80,0,5);
  // LThetaVsLTPone->SetName("LThetaVsLTPone");
  // LThetaVsLTPone->SetTitle("LThetaVsLTPone");
  // LThetaVsLTPone->SetFillColor(38);
  // LThetaVsLTPone->Write();
  //
  //
  //
  //
  //
  // gMinuit->SetErrorDef(9); //note 4 and not 2!
  // TGraph *LPhiVsLTPthree = (TGraph*)gMinuit->Contour(80,3,5);
  // LPhiVsLTPthree->SetName("LPhiVsLTPthree");
  // LPhiVsLTPthree->SetTitle("LPhiVsLTPthree");
  // LPhiVsLTPthree->SetFillColor(42);
  // LPhiVsLTPthree->Write();
  // gMinuit->SetErrorDef(4); //note 4 and not 2!
  // TGraph *LPhiVsLTPtwo = (TGraph*)gMinuit->Contour(80,3,5);
  // LPhiVsLTPtwo->SetName("LPhiVsLTPtwo");
  // LPhiVsLTPtwo->SetTitle("LPhiVsLTPtwo");
  // LPhiVsLTPtwo->SetFillColor(42);
  // LPhiVsLTPtwo->Write();
  // /*Get contour for parameter 0 versus parameter 2 for ERRDEF=1*/
  // gMinuit->SetErrorDef(1);
  // TGraph *LPhiVsLTPone = (TGraph*)gMinuit->Contour(80,3,5);
  // LPhiVsLTPone->SetName("LPhiVsLTPone");
  // LPhiVsLTPone->SetTitle("LPhiVsLTPone");
  // LPhiVsLTPone->SetFillColor(38);
  // LPhiVsLTPone->Write();
  // ContourFile.Close();

}
//_____________________________________________________________________________
/* - Fit function for the helicity case.
 * - It is basically a parabolic fit...
 * - Only CosTheta Fit.
 */
void fitOnlyCosTheta(){

  TFile* fileNew           = new TFile("pngResults/TH1corr25bins.root");
  TH1F*  CorrectedCosTheta = (TH1F*) fileNew->Get("RawCosTheta2H");

  Double_t CosThetaLowLimit   = -1;
  Double_t CosThetaUpperLimit = +1;

  Int_t nBinsCosTheta = CorrectedCosTheta->GetNbinsX();

  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure
  for (Int_t ix = 7; ix <= nBinsCosTheta-6; ++ix) {
    // coords.push_back( std::make_pair(xaxis1->GetBinCenter(ix), yaxis1->GetBinCenter(iy) ) );
    coords.push_back( CorrectedCosTheta->GetXaxis()->GetBinCenter(ix) );
    values.push_back( CorrectedCosTheta->GetBinContent(ix)            );
    errors.push_back( CorrectedCosTheta->GetBinError(ix)              );
  }

  for(int i=0; i < 18; i++){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }

  TMinuit *gMinuit = new TMinuit(2);
  gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.1, -2, 2);
  gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.58e+04, 2.62e+04);
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  // gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,LambdaPhi,NormalTheta,NormalisPhi;
  Double_t LambdaThetaErr,LambdaPhiErr,NormalThetaErr,NormalisPhiErr;
  gMinuit->GetParameter(0,LambdaTheta,LambdaThetaErr);
  gMinuit->GetParameter(1,NormalTheta,NormalThetaErr);
  printf("LambdaTheta: %+.7f +- %.7f\n",LambdaTheta,LambdaThetaErr);
  printf("NormalTheta: %+.7f +- %.7f\n",NormalTheta,NormalThetaErr);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  TF1* Model = new TF1("Model", "[1]*(1+[0]*x*x)/(3+[0])", -1 ,1 );
  new TCanvas;
  CorrectedCosTheta->Draw();
  Model->SetParameter( 0, LambdaTheta );
  Model->SetParameter( 1, NormalTheta );
  Model->SetNpx(500);
  Model->Draw("same");

}
//_____________________________________________________________________________
/* - Fit function for the helicity case.
 * - Flat Phi ONLY...
 * -
 */
void fitOnlyPhi(){

  TFile* file1D = new TFile("pngResults/TH1corr.root");
  TH1F* CorrectedPhi      = (TH1F*) file1D->Get("RawPhiH");

  Double_t PhiLowLimit        = -3.14;
  Double_t PhiUpperLimit      = +3.14;

  Int_t nBinsPhi      = CorrectedPhi     ->GetNbinsX();


  /// reset data structure
  // coords = std::vector<std::pair<double,double> >();
  coords = std::vector<Double_t>();
  values = std::vector<Double_t>();
  errors = std::vector<Double_t>();
  /// fill data structure
  for (Int_t iy = 1; iy <= nBinsPhi; ++iy) {
    coords.push_back( CorrectedPhi     ->GetXaxis()->GetBinCenter(iy) );
    values.push_back( CorrectedPhi     ->GetBinContent(iy)            );
    errors.push_back( CorrectedPhi     ->GetBinError(iy)              );
  }

  for(int i=0; i < 50; i++){
    cout << i << "  " << coords[i] << "  " << values[i] << endl;
  }


  TMinuit *gMinuit = new TMinuit(2);
  gMinuit->SetFCN(FcnForMinimisation);
  // gMinuit->DefineParameter(0, "LambdaTheta", 1., 0.1, -2, 2);
  // gMinuit->DefineParameter(1, "NormalTheta", 2.60e+04, 100,  2.58e+04, 2.62e+04);
  // gMinuit->DefineParameter(1, "NormalTheta", 2.52497e+04, 100,  2.50e+04, 2.55e+04);
  gMinuit->DefineParameter(0, "NormalisPhi", 4000, 100,  3500, 5000);
  gMinuit->DefineParameter(1, "LambdaPhi",   0, 0.1, -2, 2);
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  // gMinuit->Command("MIGRAD");
  // gMinuit->Command("MINOS");
  Double_t LambdaTheta,LambdaPhi,NormalTheta,NormalisPhi;
  Double_t LambdaThetaErr,LambdaPhiErr,NormalThetaErr,NormalisPhiErr;
  // gMinuit->GetParameter(0,LambdaTheta,LambdaThetaErr);
  // gMinuit->GetParameter(1,NormalTheta,NormalThetaErr);
  gMinuit->GetParameter(0,NormalisPhi,NormalisPhiErr);
  gMinuit->GetParameter(1,LambdaPhi,  LambdaPhiErr);
  printf("LambdaTheta: %+.7f +- %.7f\n",LambdaTheta,LambdaThetaErr);
  printf("NormalTheta: %+.7f +- %.7f\n",NormalTheta,NormalThetaErr);
  printf("NormalisPhi: %+.7f +- %.7f\n",NormalisPhi,NormalisPhiErr);
  printf("LambdaPhi  : %+.7f +- %.7f\n",LambdaPhi,  LambdaPhiErr);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);



  TF1* Model2 = new TF1("Model2", "[1]*(1+2*[2]*cos(2*x)/(3+[0]))", -3.1 ,3.1 );
  new TCanvas;
  CorrectedPhi->Draw();
  Model2->SetParameter( 0, 1.13220 );
  Model2->SetParameter( 1, NormalisPhi );
  Model2->SetParameter( 2, LambdaPhi );
  Model2->SetNpx(500);
  Model2->Draw("same");


}
