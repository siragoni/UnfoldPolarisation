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

// ==========================
// Available methods
// --------------------------
Double_t helicity2Dv3              (Double_t *x, Double_t *par);
void     FcnForMinimisation        (Int_t &, Double_t *, Double_t &fval, Double_t *p, Int_t );
void     FcnForMinimisationWithBins(Int_t &, Double_t *, Double_t &fval, Double_t *p, Int_t );
Double_t IntegralFunction          (Double_t *x, Double_t *par);
void     PolarisationHeMinuit2D    (Int_t FitRangeMode = 0, Int_t Iterations = 1 );
Double_t calcChi2                  (TH1* h1, TH1* h2, int nbins);
void     DrawResidualV3            (Int_t Iterations = 0);
// ==========================
// Global data
// --------------------------
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsY;
std::vector< Double_t > values;
std::vector< Double_t > errors;
Double_t GlobalChi = 0;
TF2* Model5 = new TF2("Model5", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
//______________________________________________
void BeautifyPad(){
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);
}
//______________________________________________
void BeautifyHisto(TH1* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");
}
//______________________________________________
void BeautifyHisto2D(TH2* histogram){
  histogram->SetTitle("");
  histogram->GetXaxis()->SetTitleOffset(1.15);
  histogram->GetYaxis()->SetTitleOffset(1.45);
  histogram->GetXaxis()->SetTitleSize(0.045);
  histogram->GetYaxis()->SetTitleSize(0.045);
  histogram->GetXaxis()->SetLabelSize(0.045);
  histogram->GetYaxis()->SetLabelSize(0.045);
  histogram->GetXaxis()->SetTitleFont(42);
  histogram->GetYaxis()->SetTitleFont(42);
  histogram->GetXaxis()->SetLabelFont(42);
  histogram->GetYaxis()->SetLabelFont(42);
}
//______________________________________________
Double_t helicity2Dv3(Double_t *x, Double_t *par) {

  Double_t CosSquaredTheta      = x[1] * x[1];
  Double_t partialTheta         = 1 + par[0] * CosSquaredTheta;

  Double_t CosOfTwoPhi          = TMath::Cos( 2 * x[0] );
  Double_t SinSquaredTheta      = 1 - CosSquaredTheta;
  Double_t partialPhi           = par[1] * SinSquaredTheta * CosOfTwoPhi;

  Double_t CosPhi               = TMath::Cos( x[0] );
  // Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredOfTwoTheta = 4. * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t partialMix           = par[2] * SinOfTwoTheta * CosPhi;


  Double_t sumOfTheSubFits      = partialTheta + partialPhi + partialMix;
  Double_t FinalResult          = par[3] * sumOfTheSubFits / ( 3 + par[0] );
  return   FinalResult;

}
//_____________________________________________________________________________
/* - Data needed for the global fit.
 * - They have to be global to be visible.
 */
void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  Int_t n = coordsX.size();
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    x[0] = coordsX[i];
    x[1] = coordsY[i];
    if ( values[i] != 0 ) {
      tmp = ( values[i] - helicity2Dv3( x, p ) ) / (50.);
      // tmp = ( values[i] - helicity2Dv3( x, p ) ) / (0.00001*values[i]);
      // tmp = ( values[i] - helicity2Dv3( x, p ) ) / errors[i];
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  fval      = chi2;
  GlobalChi = chi2;
}
//______________________________________________
Double_t IntegralFunction(Double_t *x, Double_t *par) {

  Model5->SetParameter(0, par[0]);
  Model5->SetParameter(1, par[1]);
  Model5->SetParameter(2, par[2]);
  Model5->SetParameter(3, par[3]);


  Double_t iPhiBins      = x[0];
  Double_t iCosThetaBins = x[1];
  Double_t M = 1.;
  if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
      iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
      iCosThetaBins == 20 || iCosThetaBins == 19 )
  {
    M = 1.;
  } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
    M = 6.;
  } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
    M = 12.;
  } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
              iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
    M = 24.;
  }

  Double_t IntegralFunction = Model5->Integral(
                                      (((Double_t) iPhiBins     )*2.*TMath::Pi()/M),
                                      (((Double_t) iPhiBins + 1.)*2.*TMath::Pi()/M),
                                      (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)),
                                      (-1.+((Double_t) iCosThetaBins + 1.)*(0.08+0.01/3.))
                                     );

  IntegralFunction *= M;
  return   IntegralFunction;

}
//______________________________________________
void FcnForMinimisationWithBins(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  Int_t n = coordsX.size();
  Double_t chi2 = 0;
  Double_t tmp,x[2],y[2];
  for ( Int_t i = 0; i < n; ++i ) {
    x[0] = coordsX[i];
    x[1] = coordsY[i];

    // cout << "====================================" << endl;
    // cout << "CosTheta Orig = " << x[1] << endl;
    // cout << "Phi Orig = " << x[0] << endl;
    Double_t TraslatedCosThetaGen = 0.5*(coordsY[i] + 1.)*24.;
    Double_t iCosThetaBins2 = -1;
    Double_t RemainderCosTheta = 0.;
    RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
    Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;

    Double_t M = 1.;
    if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
        iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
        iCosThetaBins == 20 || iCosThetaBins == 19 )
    {
      M = 1.;
    } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
      M = 6.;
    } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
      M = 12.;
    } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
      M = 24.;
    }
    Double_t TraslatedPhiGen = coordsX[i]*M/(2.*TMath::Pi());
    Double_t iPhiBins2 = -1;
    Double_t RemainderPhi = 0.;
    RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
    Int_t iPhiBins = (Int_t)  iPhiBins2;


    if ( values[i] != 0 ) {
      y[0] = iPhiBins;
      y[1] = iCosThetaBins;
      // tmp = ( values[i] - IntegralFunction(y,p)) / errors[i];
      tmp = ( values[i] - IntegralFunction(y,p)) / (50.);
      // tmp = ( values[i] - IntegralFunction(y,p)) / (0.00001*values[i]);
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  fval      = chi2;
  GlobalChi = chi2;
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void PolarisationHeMinuit2D(Int_t FitRangeMode = 0, Int_t Iterations = 1 ){
  TFile* f[24];
  TH1F*  h[24];
  TH1F*  generated[24];
  Double_t chi2_all[24];
  Double_t chi2_total = 0.;
  for (Int_t i = 4; i < 20; i++) {
    f[i]         = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    // h[i]         = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i]         = (TH1F*) f[i]->Get("response;2"); // actual
    // h[i]         = (TH1F*) f[i]->Get(Form("PhiGenBinsH_%d", i));
    // generated[i] = (TH1F*) f[i]->Get(Form("PhiGenBinsH_%d", i));
    // Int_t nbins  = generated[i]->GetNbinsX();
    // chi2_all[i]  = calcChi2(generated[i], h[i], nbins);
    // chi2_total  += chi2_all[i];
    // h[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i)); // already multiplied by M
  }
  // cout << "=========================================" << endl;
  // cout << "++ All Chi2 from closure distributions ++" << endl;
  // cout << "-----------------------------------------" << endl;
  // for (Int_t i = 5; i < 20; i++) {
  //   cout << "Partial Chi2[" << i << "] = " << chi2_all[i] << endl;
  // }
  // cout << "Total Chi2 = " << chi2_total << endl;
  // cout << "=========================================" << endl;



  // ===============================================================
  // Draw Residuals V3
  // ---------------------------------------------------------------
  // Int_t Counter = 1;
  // Int_t M = 1;
  // TH1F* DifferencesH = new TH1F("DifferencesH","DifferencesH", 300, -0.5, 299.5);
  // for (Int_t i = 4; i < 20; i++) {
  //   if( i == 0  || i == 1  || i == 2  || i == 3  ||
  //       i == 4  || i == 23 || i == 22 || i == 21 ||
  //       i == 20 || i == 19 )
  //   {
  //     M = 1;
  //   } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
  //     M = 6;
  //   } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
  //     M = 12;
  //   } else if ( i == 9  || i == 10  || i == 11  ||
  //               i == 12 || i == 13  || i == 14 ) {
  //     M = 24;
  //
  //   }
  //   for (Int_t j = 1; j < M+1; j++) {
  //     DifferencesH->SetBinContent(Counter, (h[i]->GetBinContent(j)-generated[i]->GetBinContent(j))/generated[i]->GetBinContent(j));
  //     DifferencesH->SetBinError(  Counter, 0.000000000000000000);
  //     Counter++;
  //   }
  // }
  // TH1F* histogram = (TH1F*) DifferencesH->Clone("histogram");
  // new TCanvas;
  // BeautifyPad();
  // BeautifyHisto(histogram);
  // histogram->GetXaxis()->SetTitle("Bin number");
  // histogram->GetYaxis()->SetTitle("(Unfolded REC - GEN)/GEN [a.u.]");
  // // histogram->GetYaxis()->SetRangeUser(-0.10,0.10);
  // histogram->GetYaxis()->SetRangeUser(histogram->GetMaximum()*(-2),histogram->GetMaximum()*(10.));
  // histogram->GetXaxis()->SetRangeUser(-0.5, 217.5);
  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  // latex5->DrawLatex(0.31,0.94,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  // latex5->DrawLatex(0.31,0.74,"This thesis");
  // if (Iterations == 1){latex5->DrawLatex(0.31,0.80,"1 iteration");}
  // else {latex5->DrawLatex(0.31,0.80,Form("%d iterations", Iterations));}
  // latex5->DrawLatex(0.2,0.64,Form("#chi^{2} = #sum (Unfolded REC - GEN)^{2}/GEN = %0.02f", chi2_total));
  // gPad->SaveAs(Form("Unfolding2D/plotsThesis/residuals-closure-%d.pdf", Iterations), "recreate");



  // =========================================
  // Drawing Refolded Distribution
  // -----------------------------------------
  new TCanvas;
  BeautifyPad();
  TH1F* UnlikeSignDimuon2    = (TH1F*)f[10]->Get("histo5");
  Double_t ChiSquareTestRef = 0;
  TH1F* UnlikeSignDimuon3    = (TH1F*)f[10]->Get("histo6");
  TH1F* UnlikeSignDimuon    = new TH1F("UnlikeSignDimuon","UnlikeSignDimuon", 300, -0.5, 299.5);
  BeautifyHisto(UnlikeSignDimuon);
  for (size_t i = 0; i < 230; i++) {
    if (UnlikeSignDimuon3->GetBinContent(i+1) > 0.0000000001){
      ChiSquareTestRef += UnlikeSignDimuon3->GetBinContent(i+1);
    }
  }
  Int_t j = 1;
  for (Int_t i = 0; i < 600; i++) {
    if (UnlikeSignDimuon2->GetBinContent(i+1) > 0.0000000001){
      UnlikeSignDimuon->SetBinError(j, 0.);
      UnlikeSignDimuon->SetBinContent(j, UnlikeSignDimuon2->GetBinContent(i+1));
      j++;
    }
  }
  UnlikeSignDimuon->GetXaxis()->SetTitle("Bin number");
  UnlikeSignDimuon->GetYaxis()->SetTitle("Ratio of refolded to REC modulation");
  // UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.99,1.01);
  // UnlikeSignDimuon->GetXaxis()->SetRangeUser(-0.5 , 229.5);
  UnlikeSignDimuon->GetYaxis()->SetRangeUser(UnlikeSignDimuon->GetMaximum()*(0.999),UnlikeSignDimuon->GetMaximum()*(1.001));
  UnlikeSignDimuon->GetXaxis()->SetRangeUser(-0.5, 215.5);
  UnlikeSignDimuon->Draw("");
  TLatex* latex10 = new TLatex();
  latex10->SetTextSize(0.045);
  latex10->SetTextFont(42);
  latex10->SetTextAlign(11);
  latex10->SetNDC();
  latex10->DrawLatex(0.31,0.94,"This thesis, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.72,"1 iteration");}
  else {latex10->DrawLatex(0.31,0.72,Form("%d iterations", Iterations));}
  latex10->DrawLatex(0.2,0.80,Form("#chi^{2} = #sum (Refolded - REC)^{2}/REC = %0.02f", ChiSquareTestRef));

  gPad->SaveAs(Form("Unfolding2D/plotsThesis/residual-refolded-data-%d.pdf", Iterations), "recreate");






  // =========================================
  // Filling the 2D histogram
  // -----------------------------------------
  TH2F* histo2 = new TH2F("histo2", "histo2", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  Double_t PhiCenters[24];
  Double_t Spacing = TMath::Pi()/24.;
  for (size_t i = 0; i < 24; i++) {
    PhiCenters[i] = (2.*((Double_t) i) + 1.) *Spacing;
    cout << "Phi centres = " << PhiCenters[i] << endl;

  }
  Double_t CosThetaCenters[24];
  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
    cout << "CosTheta centres = " << CosThetaCenters[i] << endl;

  }
  /// reset data structure
  coordsX = std::vector<Double_t>();
  coordsY = std::vector<Double_t>();
  values  = std::vector<Double_t>();
  errors  = std::vector<Double_t>();
  /// fill data structure
  Double_t Mn = 1.;
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 8; iCosThetaBins < 17; iCosThetaBins++) {
    if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
        iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
        iCosThetaBins == 20 || iCosThetaBins == 19 )
    {
      Mn = 1.;
    } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
      Mn = 6.;
    } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
      Mn = 12.;
    } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
      Mn = 24.;

    }

    if        ( FitRangeMode == 1 ) {
      if (iCosThetaBins      == 5)   continue;
    } else if ( FitRangeMode == 2 ) {
      if (iCosThetaBins      == 18)  continue;
    } else if ( FitRangeMode == 3 ) {
      if ( (iCosThetaBins    == 5) || (iCosThetaBins == 19) )  continue;
    } else {
    }
    for (size_t iPhiBins = 0; iPhiBins < 24; iPhiBins++) {



      Int_t binx = -1;

      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        binx = 1;
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        if(iPhiBins < 4){
          binx = 1;
        } else if (iPhiBins < 8){
          binx = 2;
        } else if (iPhiBins < 12){
          binx = 3;
        } else if (iPhiBins < 16){
          binx = 4;
        } else if (iPhiBins < 20){
          binx = 5;
        } else if (iPhiBins < 24){
          binx = 6;
        }
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        if(iPhiBins < 2){
          binx = 1;
        } else if (iPhiBins < 4){
          binx = 2;
        } else if (iPhiBins < 6){
          binx = 3;
        } else if (iPhiBins < 8){
          binx = 4;
        } else if (iPhiBins < 10){
          binx = 5;
        } else if (iPhiBins < 12){
          binx = 6;
        } else if (iPhiBins < 14){
          binx = 7;
        } else if (iPhiBins < 16){
          binx = 8;
        } else if (iPhiBins < 18){
          binx = 9;
        } else if (iPhiBins < 20){
          binx = 10;
        } else if (iPhiBins < 22){
          binx = 11;
        } else if (iPhiBins < 24){
          binx = 12;
        }
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
        binx = iPhiBins+1;
      }

      cout << "binx  = " << binx << endl;
      Int_t biny2 = histo2->GetYaxis()->FindBin(CosThetaCenters[iCosThetaBins]);
      Int_t binx2 = histo2->GetXaxis()->FindBin(PhiCenters[iPhiBins]);

      // histo->Fill( PhiCenters[iPhiBins], CosThetaCenters[iCosThetaBins], h[iCosThetaBins]->GetBinContent(binx) );
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx)*Mn );
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx)*Mn );

      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        if (iPhiBins == 0) {
          coordsX.push_back( TMath::Pi() );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        }
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        cout << "iPhibins" << iPhiBins << endl;
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/12. + 4.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/12. + 8.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/12. + 12.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/12. + 16.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/12. + 20.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        }
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 2){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*4.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 6){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*6.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*8.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 10){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*10.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*12.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 14){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*14.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*16.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 18){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*18.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*20.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 22){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*22.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        } else if (iPhiBins == 24){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*24.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
        }
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 )
      {
        coordsX.push_back( PhiCenters[iPhiBins] );
        coordsY.push_back( CosThetaCenters[iCosThetaBins] );
        // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
        // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        values.push_back(  h[iCosThetaBins]->GetBinContent(binx)*Mn        );
        errors.push_back(  h[iCosThetaBins]->GetBinError(binx)*Mn        );
      }


    }
  }










  /* - NB:
   * - COMPUTING THE
   * - CHI SQUARE
   */
  int ndf  = coordsX.size()-4;
  int ndf2 = coordsY.size()-4;
  int ndf3 = values.size()-4;
  int ndf4 = errors.size()-4;

  cout << "ndf  = " << ndf  << endl;
  cout << "ndf2 = " << ndf2 << endl;
  cout << "ndf3 = " << ndf3 << endl;
  cout << "ndf4 = " << ndf4 << endl;






  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationWithBins);

  // gMinuit->SetFCN(FcnForMinimisationWithBinsCovariance);
  // gMinuit->SetFCN(FcnForMinimisationCovarianceSimple);
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.,    -2, 2        );
  gMinuit->DefineParameter(0, "LambdaTheta",        0.6, 0.01,    -2, 2        );
  gMinuit->DefineParameter(1, "LambdaPhi",           0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(3, "Normalisation",   2420000, 100, 1, 1000000000    );
  // gMinuit->DefineParameter(3, "Normalisation",   (4.00207e+07*3./(4.*TMath::Pi())), 0, 1, 10000000    ); // from projection
  // gMinuit->DefineParameter(3, "Normalisation",   (394082.*3./(4.*TMath::Pi())), 100, ((394082.-2.*19165.5)*3./(4.*TMath::Pi())), ((394082.+2.*19165.5)*3./(4.*TMath::Pi()))    ); // from projection
  gMinuit->Command("SIMPLEX");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MIGRAD");
  gMinuit->Command("MINOS");
  Double_t LambdaTheta,    LambdaPhi,        LambdaThetaPhi,    Normalisation;
  Double_t LambdaThetaErr, LambdaPhiErr,     LambdaThetaPhiErr, NormalisationErr;
  gMinuit->GetParameter(0, LambdaTheta,      LambdaThetaErr     );
  gMinuit->GetParameter(1, LambdaPhi,        LambdaPhiErr       );
  gMinuit->GetParameter(2, LambdaThetaPhi,   LambdaThetaPhiErr  );
  gMinuit->GetParameter(3, Normalisation,    NormalisationErr   );
  printf("LambdaTheta     : %+.7f +- %.7f\n",LambdaTheta,     LambdaThetaErr     );
  printf("LambdaPhi       : %+.7f +- %.7f\n",LambdaPhi,       LambdaPhiErr       );
  printf("LambdaThetaPhi  : %+.7f +- %.7f\n",LambdaThetaPhi,  LambdaThetaPhiErr  );
  printf("Normalisation   : %+.7f +- %.7f\n",Normalisation,   NormalisationErr   );
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(3,amin);

  gStyle->SetOptStat(0);

  cout << "OK1" << endl << flush;

















  TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  Model->FixParameter(0, LambdaTheta);
  Model->FixParameter(1, LambdaPhi);
  Model->FixParameter(2, LambdaThetaPhi);
  Model->SetParameter(3, Normalisation);

  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  BeautifyHisto2D(histo2);
  histo2->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );









    const Int_t NRGBs = 6;
    const Int_t NCont = 999;

    Double_t stops[NRGBs] = { 0.00, 0.1, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.99, 0.0, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.0, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.99, 0.0, 1.00, 0.12, 0.00, 0.00 };


    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    gStyle->SetOptStat(0);

    //here the actually interesting code starts
    // const Double_t min = 8000;
    // const Double_t min = 1500000;
    // const Double_t max = 5000000;
    const Double_t min = histo2->GetMaximum()*0.3;
    const Double_t max = histo2->GetMaximum()*1.3;

    const Int_t nLevels = 999;
    Double_t levels[nLevels];


    for(int i = 1; i < nLevels; i++) {
      levels[i] = min + (max - min) / (nLevels - 1) * (i);
    }
    levels[0] = 0.01;

    histo2->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    histo2->DrawClone("col");// draw "axes", "contents", "statistics box"
    histo2->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
    histo2->Draw("colz same"); // draw the "color palette"
    // Model->Draw("same");



  TLatex* latex = new TLatex();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.12,0.94,"This thesis, Data, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->DrawLatex(0.7,0.94,"Helicity");
  latex->SetTextSize(0.04);
  // LambdaTheta = 0.600;
  // LambdaPhi = 0.000;
  // LambdaThetaPhi = -0.004;
  // LambdaThetaErr = 0.006;
  // LambdaPhiErr = 0.001;
  latex->DrawLatex(0.6,0.84,Form("#lambda_{#theta} = %0.3f #pm %0.3f", LambdaTheta,      LambdaThetaErr ));
  latex->DrawLatex(0.6,0.79,Form("#lambda_{#varphi} = %0.3f #pm %0.3f", LambdaPhi,        LambdaPhiErr ));
  latex->DrawLatex(0.6,0.74,Form("#lambda_{#theta#varphi} = %0.3f #pm %0.3f", LambdaThetaPhi,   LambdaThetaPhiErr ));
  // latex->DrawLatex(0.6,0.15,Form("Normalisation = %0.3f #pm %0.3f ", Normalisation,    NormalisationErr ));
  latex->DrawLatex(0.15,0.15,Form("#tilde{#chi^{2}} = %0.3f / %d = %0.3f ", GlobalChi,    ndf2, GlobalChi/((Double_t)ndf2) ));
  // gPad->SaveAs("Unfolding2D/plotsThesis/Inverted2Dmaps.pdf",  "RECREATE");
  gPad->SaveAs(Form("Unfolding2D/plotsThesis/Inverted2Dmaps-data-%d.pdf", Iterations), "recreate");


  new TCanvas;
  histo2->Draw("surf same");
  cout << "Fitting model for Normalisation " << endl;
  histo2->Fit("Model", "R");
  // gPad->SaveAs("Unfolding2D/plotsThesis/3D.pdf",  "RECREATE");
  gPad->SaveAs(Form("Unfolding2D/plotsThesis/3D-data-%d.pdf", Iterations), "recreate");




  // TF2* Model2 = new TF2("Model2", IntegralFunction,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  // Model2->FixParameter(0, LambdaTheta);
  // Model2->FixParameter(1, LambdaPhi);
  // Model2->FixParameter(2, LambdaThetaPhi);
  // Model2->FixParameter(3, Normalisation);
  // Model5->FixParameter(0, LambdaTheta);
  // Model5->FixParameter(1, LambdaPhi);
  // Model5->FixParameter(2, LambdaThetaPhi);
  // Model5->FixParameter(3, Normalisation);

  TH1F* pulls2     = new TH1F("pulls2",     "pulls2",     500, -0.5, 499.5);
  TH1F* valuehis2  = new TH1F("valuehis2",  "valuehis2",  500, -0.5, 499.5);
  TH1F* modelhis2  = new TH1F("modelhis2",  "modelhis2",  500, -0.5, 499.5);
  TH1F* residuals2 = new TH1F("residuals2", "residuals2", 500, -0.5, 499.5);
  for (size_t i = 0; i < ndf; i++) {
    // pulls2    ->SetBinContent(i+1,  values[i] - Model2->Eval(coordsX[i], coordsY[i]));
    // valuehis2 ->SetBinContent(i+1,  values[i] );
    // modelhis2 ->SetBinContent(i+1,  Model2->Eval(coordsX[i], coordsY[i]) );
    // residuals2->SetBinContent(i+1, (values[i] - Model2->Eval(coordsX[i], coordsY[i]))/errors[i]);
    pulls2    ->SetBinContent(i+1,  values[i] - Model5->Eval(coordsX[i], coordsY[i]));
    valuehis2 ->SetBinContent(i+1,  values[i] );
    modelhis2 ->SetBinContent(i+1,  Model5->Eval(coordsX[i], coordsY[i]) );
    residuals2->SetBinContent(i+1, (values[i] - Model5->Eval(coordsX[i], coordsY[i]))/errors[i]);
    // pulls2    ->SetBinContent(i+1,  values[i] - Model->Eval(coordsX[i], coordsY[i]));
    // valuehis2 ->SetBinContent(i+1,  values[i] );
    // modelhis2 ->SetBinContent(i+1,  Model->Eval(coordsX[i], coordsY[i]) );
    // residuals2->SetBinContent(i+1, (values[i] - Model->Eval(coordsX[i], coordsY[i]))/errors[i]);

  }
  new TCanvas;
  pulls2->Draw();
  new TCanvas;
  residuals2->Draw();
  new TCanvas;
  valuehis2->Draw();
  new TCanvas;
  modelhis2->Draw();

}
//_____________________________________________________________________________
Double_t calcChi2(TH1* h1, TH1* h2, int nbins) {
    // Int_t nb = h1->GetNbinsX();
    // std::cout << "chi2 Nbins:" << nb << std::endl;
    std::cout << "chi2 Nbins:" << nbins << std::endl;
    Double_t chi2 = 0;
    for(int i = 1; i < nbins+1; i++) {
        Double_t b1 =h1->GetBinContent(i);
        Double_t b2 =h2->GetBinContent(i);
        //std::cout << "b1:" << b1 << " b2:" << b2 << std::endl;
        chi2 += (b1-b2) *(b1-b2)/b1;
    }
    return chi2;
}
