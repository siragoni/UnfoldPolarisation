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

Double_t GlobalChi = 0;


//_____________________________________________________________________________
/* - Codign in the fit functions.
   - The fit is modelled as the sum of 3 different functions:
   - 1) CosTheta only
   - 2) Phi      only
   - 3) Mix of the two
   -
 */
Double_t CosTheta(Double_t *x, Double_t *par) {
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t returnValue     = 1 + par[0] * CosSquaredTheta;
  return   returnValue;
}
//______________________________________________
Double_t Phi(Double_t *x, Double_t *par) {
  Double_t CosOfTwoPhi     = TMath::Cos( 2 * x[1] );
  Double_t CosSquaredTheta = x[0] * x[0];
  Double_t SinSquaredTheta = 1 - CosSquaredTheta;

  Double_t returnValue = par[0] * SinSquaredTheta * CosOfTwoPhi;
  return   returnValue;
}
//______________________________________________
Double_t Mix(Double_t *x, Double_t *par) {
  Double_t CosPhi               = TMath::Cos( x[1] );
  Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredTheta      = 1 - CosSquaredTheta;
  Double_t SinSquaredOfTwoTheta = 4 * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t returnValue          = par[0] * SinOfTwoTheta * CosPhi;
  return   returnValue;
}
//______________________________________________
Double_t helicity2D(Double_t *x, Double_t *par) {
   Double_t *lambdaTheta    = &par[0];
   Double_t *lambdaPhi      = &par[1];
   Double_t *lambdaThetaPhi = &par[2];
   Double_t sumOfTheSubFits = CosTheta( x, lambdaTheta) + Phi( x, lambdaPhi ) + Mix( x, lambdaThetaPhi );
   Double_t FinalResult     = par[3] * sumOfTheSubFits / ( 3 + par[0] );
   return   FinalResult;
}
//______________________________________________
Double_t helicity2Dv2(Double_t *x, Double_t *par) {

  Double_t CosSquaredTheta      = x[0] * x[0];
  Double_t partialTheta         = 1 + par[0] * CosSquaredTheta;

  Double_t CosOfTwoPhi          = TMath::Cos( 2 * x[1] );
  Double_t SinSquaredTheta      = 1 - CosSquaredTheta;
  Double_t partialPhi           = par[1] * SinSquaredTheta * CosOfTwoPhi;

  Double_t CosPhi               = TMath::Cos( x[1] );
  // Double_t CosSquaredTheta      = x[0]*x[0];
  Double_t SinSquaredOfTwoTheta = 4 * CosSquaredTheta * SinSquaredTheta;
  Double_t SinOfTwoTheta        = TMath::Sqrt(SinSquaredOfTwoTheta);
  Double_t partialMix           = par[2] * SinOfTwoTheta * CosPhi;


  Double_t sumOfTheSubFits      = partialTheta + partialPhi + partialMix;
  Double_t FinalResult          = par[3] * sumOfTheSubFits / ( 3 + par[0] );
  return   FinalResult;

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
  Double_t SinSquaredOfTwoTheta = 4 * CosSquaredTheta * SinSquaredTheta;
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
// std::vector< std::pair< Double_t, Double_t > > coords;
std::vector< Double_t > coordsX;
std::vector< Double_t > coordsY;
std::vector< Double_t > values;
std::vector< Double_t > errors;
std::vector< Double_t > errorsBin;

void FcnForMinimisation(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  // cout << "HI" << flush << endl;
  Int_t n = coordsX.size();
  // cout << "HI2" << flush << endl;

  // Int_t n = 40;
  Double_t chi2 = 0;
  Double_t tmp,x[2];
  for ( Int_t i = 0; i < n; ++i ) {
    x[1] = coordsX[i];
    x[0] = coordsY[i];
    if ( values[i] != 0 ) {
      tmp = ( values[i] - helicity2Dv2( x, p ) ) / errors[i];
      // tmp = ( values[i] - helicity2Dv2( x, p ) ) / (values[i]*TMath::Sqrt(errors[i]*errors[i]/(values[i]*values[i])+errorsBin[i]*errorsBin[i]));
    } else {
      tmp = 0;
    }
    chi2 += tmp*tmp;
  }
  fval      = chi2;
  GlobalChi = chi2;
}
//_____________________________________________________________________________
// TF2* Model5 = new TF2("Model5", helicity2Dv3,0., 2.*TMath::Pi(), -1. ,1., 4 );
TF2* Model5 = new TF2("Model5", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
//______________________________________________
Double_t IntegralFunction(Double_t *x, Double_t *par) {

  Model5->SetParameter(0, par[0]);
  Model5->SetParameter(1, par[1]);
  Model5->SetParameter(2, par[2]);
  Model5->SetParameter(3, par[3]);
  // cout << "Model5 eval " << Model5->Eval(0, 0) << endl;


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

  // cout << "iPhiBins      = " << iPhiBins << endl;
  // cout << "iCosThetaBins = " << iCosThetaBins << endl;
  // cout << "Phi      = " << ((Double_t) iPhiBins     )*2.*TMath::Pi()/M <<  endl;
  // cout << "CosTheta = " << (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)) << endl;

  // Double_t IntegralFunction = Model5->Integral(
  //                                     (((Double_t) iPhiBins     )*2.*TMath::Pi()/M),
  //                                     (((Double_t) iPhiBins + 1.)*2.*TMath::Pi()/M),
  //                                     (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)),
  //                                     (-1.+((Double_t) iCosThetaBins + 1.)*(0.08+0.01/3.))
  //                                    );
  Double_t IntegralFunction = Model5->Integral(
                                      (((Double_t) iPhiBins     )*2.*TMath::Pi()/M),
                                      (((Double_t) iPhiBins + 1.)*2.*TMath::Pi()/M),
                                      (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)),
                                      (-1.+((Double_t) iCosThetaBins + 1.)*(0.08+0.01/3.))
                                     );

  IntegralFunction *= M;
  // IntegralFunction /= 2.*TMath::Pi()/M;
  // IntegralFunction /= 2.*TMath::Pi()/(M*(0.08+0.01/3.));
  return   IntegralFunction;

}

void FcnForMinimisationWithBins(Int_t & /*nPar*/, Double_t * /*grad*/ , Double_t &fval, Double_t *p, Int_t /*iflag */  )
{
  // cout << "HI" << flush << endl;
  Int_t n = coordsX.size();
  // cout << "HI2" << flush << endl;
  // Int_t n = 40;
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
      // // tmp = ( values[i] - helicity2Dv2( x, p ) ) / errors[i];
      // Model5->SetParameter(0, 1);
      // Model5->SetParameter(1, 0);
      // Model5->SetParameter(2, 0);
      // Model5->SetParameter(3, 251000);
      // // Model5->SetParameter(0, p[0]);
      // // Model5->SetParameter(1, p[1]);
      // // Model5->SetParameter(2, p[2]);
      // // Model5->SetParameter(3, p[3]);
      // // Double_t IntegralFunction = 0;
      // Double_t IntegralFunction = Model5->Integral(
      //                                     (((Double_t) iPhiBins     )*2.*TMath::Pi()/M),
      //                                     (((Double_t) iPhiBins + 1.)*2.*TMath::Pi()/M),
      //                                     (-1.+((Double_t) iCosThetaBins     )*(0.08+0.01/3.)),
      //                                     (-1.+((Double_t) iCosThetaBins + 1.)*(0.08+0.01/3.))
      //                                    );
      y[0] = iPhiBins;
      y[1] = iCosThetaBins;
      tmp = ( values[i] - IntegralFunction(y,p)) / errors[i];
      // tmp = ( values[i] - IntegralFunction(y,p)) / values[i];

      // tmp = ( values[i] - IntegralFunction) / errors[i];
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
void PolarisationHeMinuit2D(Int_t FitRangeMode = 0 ){
// void PolarisationHeMinuit2D( Int_t SignalRangeSelectionMode = 0, Int_t FitRangeMode = 0 ){


  TFile* f[24];
  TH1F*  h[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo3_%d", i));
  }

  Double_t integrals[24];
  Double_t uncertainties[24];
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
    for (size_t j = 0; j < 50; j++) {


      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        if((j+1)>1) continue;
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        if((j+1)>6) continue;

      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        if((j+1)>12) continue;

      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 )
      {
        if((j+1)>24) continue;

      }



      integrals[iCosThetaBins] += h[iCosThetaBins]->GetBinContent(j+1);
      uncertainties[iCosThetaBins] += h[iCosThetaBins]->GetBinError(j+1);

      cout << "(iCosThetaBins,iPhibins) = " << "(" << iCosThetaBins << "," << j << "), " << h[iCosThetaBins]->GetBinContent(j+1) << "+/-" << h[iCosThetaBins]->GetBinError(j+1) << endl;
      cout << "(iCosThetaBins,iPhibins) = " << "(" << iCosThetaBins << "," << j << "), " << integrals[iCosThetaBins] << "(+/-)" << uncertainties[iCosThetaBins] << endl;

    }
  }
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {


      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        integrals[iCosThetaBins]      /= 1.;
        uncertainties[iCosThetaBins]  /= 1.;
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        integrals[iCosThetaBins]      /= 6.;
        uncertainties[iCosThetaBins]  /= 6.;

      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        integrals[iCosThetaBins]      /= 12.;
        uncertainties[iCosThetaBins]  /= 12.;

      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 )
      {
        integrals[iCosThetaBins]      /= 24.;
        uncertainties[iCosThetaBins]  /= 24.;

      }




  }

  TH1F* projection = new TH1F("proj", "proj", 24, -1, 1);
  for (size_t i = 5; i < 19; i++) {
    projection->SetBinContent(i+1, integrals[i]);
    projection->SetBinError(i+1, uncertainties[i]);

  }
  new TCanvas;
  projection->Draw();

  Double_t FullIntegral    = 0;
  Double_t FullIntegralErr = 0;
  for (size_t i = 5; i < 19; i++) {
    FullIntegral    += integrals[i];
    FullIntegralErr += uncertainties[i];

  }



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
  // for (size_t iCosThetaBins = 9; iCosThetaBins < 18; iCosThetaBins++) {
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 17; iCosThetaBins++) {
    // if        ( FitRangeMode == 1 ) {
    //   if (iCosThetaBins      == 5)   continue;
    // } else if ( FitRangeMode == 2 ) {
    //   if (iCosThetaBins      == 18)  continue;
    // } else if ( FitRangeMode == 3 ) {
    //   if ( (iCosThetaBins    == 5) || (iCosThetaBins == 19) )  continue;
    // } else {
    // }
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
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) );
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        errorsBin.push_back(  TMath::Pi()/TMath::Sqrt(12.)     );
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        errorsBin.push_back(  TMath::Pi()/(6.*TMath::Sqrt(12.))    );
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        errorsBin.push_back(  TMath::Pi()/(12.*TMath::Sqrt(12.))    );
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
        errorsBin.push_back(  TMath::Pi()/(24.*TMath::Sqrt(12.))    );
      }
      // errorsBin.push_back( 0.    );


      // coordsX.push_back( ((TAxis*) Distr2D->GetXaxis())->GetBinCenter( ix ) );
      // coordsY.push_back( CosThetaCenters[iCosThetaBins] );
      // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
      // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        if (iPhiBins == 0) {
          coordsX.push_back( TMath::Pi() );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        }
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        cout << "iPhibins" << iPhiBins << endl;
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/12. + 4.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/12. + 8.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/12. + 12.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/12. + 16.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/12. + 20.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        }
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 2){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*4.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 6){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*6.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*8.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 10){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*10.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*12.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 14){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*14.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*16.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 18){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*18.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*20.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 22){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*22.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        } else if (iPhiBins == 24){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*24.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
        }
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 )
      {
        coordsX.push_back( PhiCenters[iPhiBins] );
        coordsY.push_back( CosThetaCenters[iCosThetaBins] );
        values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
        errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
      }


    }
  }

  //
  // // Filling the data
  // if        ( FitRangeMode != 1 ||  FitRangeMode != 3) {
  //       coordsX.push_back( TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 2.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 4.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 6.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 8.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 10.*TMath::Pi()/12. );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       coordsY.push_back( CosThetaCenters[5] );
  //       for (size_t i = 0; i < 6; i++) {
  //         values.push_back(  h[5]->GetBinContent(i+1) );
  //         errors.push_back(  h[5]->GetBinError(i+1)   );
  //       }
  // }
  // coordsX.push_back( TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 2.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 4.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 6.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 8.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 10.*TMath::Pi()/12. );
  // coordsY.push_back( CosThetaCenters[6] );
  // coordsY.push_back( CosThetaCenters[6] );
  // coordsY.push_back( CosThetaCenters[6] );
  // coordsY.push_back( CosThetaCenters[6] );
  // coordsY.push_back( CosThetaCenters[6] );
  // coordsY.push_back( CosThetaCenters[6] );
  // for (size_t i = 0; i < 6; i++) {
  //   values.push_back(  h[6]->GetBinContent(i+1) );
  //   errors.push_back(  h[6]->GetBinError(i+1)   );
  // }
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  // coordsX.push_back( TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 2.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 4.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 6.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 8.*TMath::Pi()/12. );
  // coordsX.push_back( TMath::Pi()/12. + 10.*TMath::Pi()/12. );
  // coordsY.push_back( CosThetaCenters[17] );
  // coordsY.push_back( CosThetaCenters[17] );
  // coordsY.push_back( CosThetaCenters[17] );
  // coordsY.push_back( CosThetaCenters[17] );
  // coordsY.push_back( CosThetaCenters[17] );
  // coordsY.push_back( CosThetaCenters[17] );
  // for (size_t i = 0; i < 6; i++) {
  //   values.push_back(  h[17]->GetBinContent(i+1) );
  //   errors.push_back(  h[17]->GetBinError(i+1)   );
  // }
  // if        ( FitRangeMode != 2 ||  FitRangeMode != 3) {
  //       coordsX.push_back( TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 2.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 4.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 6.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 8.*TMath::Pi()/12. );
  //       coordsX.push_back( TMath::Pi()/12. + 10.*TMath::Pi()/12. );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       coordsY.push_back( CosThetaCenters[18] );
  //       for (size_t i = 0; i < 6; i++) {
  //         values.push_back(  h[18]->GetBinContent(i+1) );
  //         errors.push_back(  h[18]->GetBinError(i+1)   );
  //       }
  // }









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


  cout << "Coordinates:" << endl;
  for( Int_t i = 0; i < coordsX.size(); i++ ){
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
  }





  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationWithBins);
  gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.,    -2, 2        );
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.1,    -2, 2        );
  gMinuit->DefineParameter(1, "LambdaPhi",           0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.1,    -2, 2        );
  // gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.,    -2, 2        );
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
  Model->FixParameter(3, Normalisation);
  TF2* Model2 = new TF2("Model2", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  Model2->FixParameter(0, 1.);
  Model2->FixParameter(1, LambdaPhi);
  Model2->FixParameter(2, LambdaThetaPhi);
  Model2->FixParameter(3, Normalisation);

  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  // Distr2D->GetXaxis()->SetTitleOffset(1.15);
  // // Distr2D->GetYaxis()->SetTitleOffset(1.25);
  // Distr2D->GetYaxis()->SetTitleOffset(1.);
  // Distr2D->GetXaxis()->SetTitleSize(0.045);
  // Distr2D->GetYaxis()->SetTitleSize(0.045);
  // Distr2D->GetXaxis()->SetLabelSize(0.045);
  // Distr2D->GetYaxis()->SetLabelSize(0.045);
  // Distr2D->GetXaxis()->SetTitleFont(42);
  // Distr2D->GetYaxis()->SetTitleFont(42);
  // Distr2D->GetXaxis()->SetLabelFont(42);
  // Distr2D->GetYaxis()->SetLabelFont(42);
  // Distr2D->GetXaxis()->SetNdivisions(408);

  histo2->GetXaxis()->SetTitleOffset(1.15);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.25);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.);
  histo2->GetYaxis()->SetTitleOffset(0.9);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelSize(0.045);
  histo2->GetXaxis()->SetTitleSize(0.055);
  histo2->GetYaxis()->SetTitleSize(0.055);
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleFont(42);
  histo2->GetYaxis()->SetTitleFont(42);
  histo2->GetXaxis()->SetLabelFont(42);
  histo2->GetYaxis()->SetLabelFont(42);
  histo2->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
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
    const Double_t min = histo2->GetMaximum()*0.001;
    const Double_t max = histo2->GetMaximum()*1.1;
    // const Double_t min = 600000;
    // const Double_t max = 900000;

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
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");
  latex->DrawLatex(0.3,0.8,"1 iteration, this thesis");

  gPad->SaveAs("Unfolding2D/Inverted2Dmaps_phimodulation.pdf",  "RECREATE");


  new TCanvas;
  histo2->Draw("surf same");
  Model->Draw("surf same");




  for( Int_t i = 0; i < coordsX.size(); i++ ){
    Double_t TraslatedCosThetaGen = 0.5*(coordsY[i] + 1.)*24.;
    Double_t iCosThetaBins2 = -1;
    Double_t RemainderCosTheta = 0.;
    RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
    Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;

    Double_t Mn = 1.;
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
    Double_t TraslatedPhiGen = coordsX[i]*Mn/(2.*TMath::Pi());
    Double_t iPhiBins2 = -1;
    Double_t RemainderPhi = 0.;
    RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
    Int_t iPhiBins = (Int_t)  iPhiBins2;


    Double_t xn[2] = {(Double_t)iPhiBins, (Double_t)iCosThetaBins};
    Double_t yn[4] = {LambdaTheta,LambdaPhi,LambdaThetaPhi,Normalisation};

    cout << "==============================================" << endl;
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
    cout << "As computed in the function? " << endl;
    cout << "CosThetaF = " << (-1.+((Double_t) iCosThetaBins  +0.5   )*(0.08+0.01/3.)) << ", PhiF =  " << (((Double_t) iPhiBins +0.5    )*2.*TMath::Pi()/Mn) <<  endl;
    cout << "iCosThetaBins = " << iCosThetaBins << endl << ",   iPhiBins = " << iPhiBins << endl;
    cout << "TraslatedPhiGen = " << TraslatedPhiGen << endl << ",   Mn = " << Mn << endl;
    cout << "IntegralFunction = " << IntegralFunction(xn,yn) << endl;
    cout << "Model            = " << Model->Eval(coordsX[i], coordsY[i]) << endl;
  }
  printf("LambdaTheta     : %+.7f +- %.7f\n",LambdaTheta,     LambdaThetaErr     );
  printf("LambdaPhi       : %+.7f +- %.7f\n",LambdaPhi,       LambdaPhiErr       );
  printf("LambdaThetaPhi  : %+.7f +- %.7f\n",LambdaThetaPhi,  LambdaThetaPhiErr  );
  printf("Normalisation   : %+.7f +- %.7f\n",Normalisation,   NormalisationErr   );


  TH1F* pulls     = new TH1F("pulls",     "pulls",     500, -0.5, 499.5);
  TH1F* valuehis  = new TH1F("valuehis",  "valuehis",  500, -0.5, 499.5);
  TH1F* modelhis  = new TH1F("modelhis",  "modelhis",  500, -0.5, 499.5);
  TH1F* residuals = new TH1F("residuals", "residuals", 500, -0.5, 499.5);
  for (size_t i = 0; i < ndf; i++) {
    pulls    ->SetBinContent(i+1,  values[i] - Model->Eval(coordsX[i], coordsY[i]));
    valuehis ->SetBinContent(i+1,  values[i] );
    modelhis ->SetBinContent(i+1,  Model->Eval(coordsX[i], coordsY[i]) );
    residuals->SetBinContent(i+1, (values[i] - Model->Eval(coordsX[i], coordsY[i]))/errors[i]);

  }
  new TCanvas;
  pulls->Draw();
  new TCanvas;
  residuals->Draw();
  new TCanvas;
  valuehis->Draw();
  new TCanvas;
  modelhis->Draw();







  TH1F* pulls2     = new TH1F("pulls2",     "pulls2",     500, -0.5, 499.5);
  TH1F* valuehis2  = new TH1F("valuehis2",  "valuehis2",  500, -0.5, 499.5);
  TH1F* modelhis2  = new TH1F("modelhis2",  "modelhis2",  500, -0.5, 499.5);
  TH1F* residuals2 = new TH1F("residuals2", "residuals2", 500, -0.5, 499.5);
  for (size_t i = 0; i < ndf; i++) {
    pulls2    ->SetBinContent(i+1,  values[i] - Model2->Eval(coordsX[i], coordsY[i]));
    valuehis2 ->SetBinContent(i+1,  values[i] );
    modelhis2 ->SetBinContent(i+1,  Model2->Eval(coordsX[i], coordsY[i]) );
    residuals2->SetBinContent(i+1, (values[i] - Model2->Eval(coordsX[i], coordsY[i]))/errors[i]);

  }
  new TCanvas;
  pulls2->Draw();
  new TCanvas;
  residuals2->Draw();
  new TCanvas;
  valuehis2->Draw();
  new TCanvas;
  modelhis2->Draw();


  std::cout << "Full integral      = " << FullIntegral << '\n';
  std::cout << "Full integral Err  = " << FullIntegralErr << '\n';



  std::cout << "Norm * 4pi/3      = " << Normalisation*4.*TMath::Pi()/3. << endl;
  std::cout << "Err(Norm * 4pi/3) = " << NormalisationErr*4.*TMath::Pi()/3. << endl;
}
//______________________________________________
void Draw2D(){

  TFile* f[24];
  TH1F*  h[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    h[i] = (TH1F*) f[i]->Get("histo");
  }



  TH2F* histo = new TH2F("histo2", "histo2", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  Double_t PhiCenters[24];
  Double_t Spacing = TMath::Pi()/24.;
  for (size_t i = 0; i < 24; i++) {
    PhiCenters[i] = (2.*((Double_t) i) + 1.) *Spacing;
    // cout << "Phi centres = " << PhiCenters[i] << endl;

  }
  Double_t CosThetaCenters[24];
  for (size_t i = 0; i < 24; i++) {
    CosThetaCenters[i] = -1. + (2.*(Double_t) i + 1.)*(0.08+0.01/3.)*0.5;
    // cout << "CosTheta centres = " << CosThetaCenters[i] << endl;

  }



  for (size_t iCosThetaBins = 4; iCosThetaBins < 20; iCosThetaBins++) {
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

      // cout << "binx  = " << binx << endl;
      Int_t biny2 = histo->GetYaxis()->FindBin(CosThetaCenters[iCosThetaBins]);
      Int_t binx2 = histo->GetXaxis()->FindBin(PhiCenters[iPhiBins]);

      // histo->Fill( PhiCenters[iPhiBins], CosThetaCenters[iCosThetaBins], h[iCosThetaBins]->GetBinContent(binx) );
      histo->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) );
      histo->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );

    }
  }











  TF2* Model = new TF2("Model", helicity2D, -0.6 ,0.6, -3.14, 3.14, 4 );
  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  // Distr2D->GetXaxis()->SetTitleOffset(1.15);
  // // Distr2D->GetYaxis()->SetTitleOffset(1.25);
  // Distr2D->GetYaxis()->SetTitleOffset(1.);
  // Distr2D->GetXaxis()->SetTitleSize(0.045);
  // Distr2D->GetYaxis()->SetTitleSize(0.045);
  // Distr2D->GetXaxis()->SetLabelSize(0.045);
  // Distr2D->GetYaxis()->SetLabelSize(0.045);
  // Distr2D->GetXaxis()->SetTitleFont(42);
  // Distr2D->GetYaxis()->SetTitleFont(42);
  // Distr2D->GetXaxis()->SetLabelFont(42);
  // Distr2D->GetYaxis()->SetLabelFont(42);
  // Distr2D->GetXaxis()->SetNdivisions(408);

  histo->GetXaxis()->SetTitleOffset(1.15);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.25);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.);
  histo->GetYaxis()->SetTitleOffset(0.9);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelSize(0.045);
  histo->GetXaxis()->SetTitleSize(0.055);
  histo->GetYaxis()->SetTitleSize(0.055);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleFont(42);
  histo->GetYaxis()->SetTitleFont(42);
  histo->GetXaxis()->SetLabelFont(42);
  histo->GetYaxis()->SetLabelFont(42);
  histo->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
  histo->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );









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
    const Double_t min = 100000;
    const Double_t max = 850000;

    const Int_t nLevels = 999;
    Double_t levels[nLevels];


    for(int i = 1; i < nLevels; i++) {
      levels[i] = min + (max - min) / (nLevels - 1) * (i);
    }
    levels[0] = 0.01;

    histo->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    histo->DrawClone("col");// draw "axes", "contents", "statistics box"
    histo->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
    histo->Draw("colz same"); // draw the "color palette"




  TLatex* latex = new TLatex();
  latex->SetTextSize(0.055);
  latex->SetTextFont(42);
  latex->SetTextAlign(11);
  latex->SetNDC();
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");

  gPad->SaveAs("Unfolding2D/Inverted2Dmaps_phimodulation.pdf",  "RECREATE");




}



//______________________________________________
void DrawResidual(){

  TFile* f[24];
  TH1F*  h[24];
  TH1F*  h2[24];
  TH1F*  h3[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    h[i]  = (TH1F*) f[i]->Get(Form("histo3_%d", i));
    h2[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    h3[i] = (TH1F*) h2[i]->Clone(Form("histo4_%d", i));
    // h3[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    h3[i]->Sumw2();
    h2[i]->Sumw2();
    h[i]->Sumw2();
    h3[i] -> Add(h[i], -1.);
    // for (Int_t j = 0; j < 30; j++) {
    //   if
    // }
    h3[i] ->Divide(h2[i]);
    // h3[i] ->Divide(h[i]);
  }

  Int_t M = 1;
  TH1F* histogram = new TH1F("histogram","histogram", 1000, -0.5,999.5);
  for (Int_t iC = 4; iC < 20; iC++) {
    if( iC == 0  || iC == 1  || iC == 2  || iC == 3  ||
        iC == 4  || iC == 23 || iC == 22 || iC == 21 ||
        iC == 20 || iC == 19 )
    {
      M = 1;
    } else if ( iC == 5  || iC == 6  || iC == 18  || iC == 17 ) {
      M = 6;
    } else if ( iC == 7  || iC == 8  || iC == 16  || iC == 15 ) {
      M = 12;
    } else if ( iC == 9  || iC == 10  || iC == 11  ||
                iC == 12 || iC == 13  || iC == 14 ) {
      M = 24;
    }
    for (Int_t j = 0; j < M; j++) {
      // if( iC == 0  || iC == 1  || iC == 2  || iC == 3  ||
      //     iC == 4  || iC == 23 || iC == 22 || iC == 21 ||
      //     iC == 20 || iC == 19 )
      // {
      //   M = 1;
      // } else if ( iC == 5  || iC == 6  || iC == 18  || iC == 17 ) {
      //   M = 6;
      // } else if ( iC == 7  || iC == 8  || iC == 16  || iC == 15 ) {
      //   M = 12;
      // } else if ( iC == 9  || iC == 10  || iC == 11  ||
      //             iC == 12 || iC == 13  || iC == 14 ) {
      //   M = 24;
      // }
      histogram->SetBinContent( iC*30+j, h3[iC]->GetBinContent( j+1 ) );
      histogram->SetBinError(   iC*30+j, 0. );
      // histogram->SetBinError(   iC*30+j, h3[iC]->GetBinError( j+1 ) );
    }
  }
  new TCanvas;
  histogram->GetYaxis()->SetRangeUser(0.95,1.05);
  histogram->Draw("*H");


  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);

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

  histogram->GetXaxis()->SetTitle("Bin number");
  histogram->GetYaxis()->SetTitle("(Unfolded REC - GEN)/GEN [a.u.]");
  histogram->GetYaxis()->SetRangeUser(-0.006,0.004);
  histogram->GetXaxis()->SetRangeUser(100, 600);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");

  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis");

  gPad->SaveAs("residuals-closure_phimodulation.pdf", "recreate");


}
//_____________________________________________________________________________
/* -
   -
 */
void DrawDifferences2D(){


  TFile* f[24];
  TH1F*  h[24];
  TH1F*  h2[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i]  = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    h2[i] = (TH1F*) f[i]->Get(Form("histo3_%d", i));
  }





  TH2F* histo2 = new TH2F("histo2", "histo2", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  TH2F* histo3 = new TH2F("histo3", "histo3", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  TH2F* histo4 = new TH2F("histo4", "histo4", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
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


  // for (size_t iCosThetaBins = 7; iCosThetaBins < 18; iCosThetaBins++) {
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 17; iCosThetaBins++) {
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
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) );
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );
      histo3->SetBinContent( binx2, biny2, h2[iCosThetaBins]->GetBinContent(binx) );
      histo3->SetBinError( binx2, biny2, h2[iCosThetaBins]->GetBinError(binx) );
      histo4->SetBinContent( binx2, biny2, h2[iCosThetaBins]->GetBinContent(binx) );
      histo4->SetBinError( binx2, biny2, h2[iCosThetaBins]->GetBinError(binx) );

    }
  }










  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);

  histo2->GetXaxis()->SetTitleOffset(1.15);
  histo2->GetYaxis()->SetTitleOffset(0.9);
  histo2->GetXaxis()->SetTitleSize(0.055);
  histo2->GetYaxis()->SetTitleSize(0.055);
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleFont(42);
  histo2->GetYaxis()->SetTitleFont(42);
  histo2->GetXaxis()->SetLabelFont(42);
  histo2->GetYaxis()->SetLabelFont(42);
  histo2->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
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
    // const Double_t min = 1500000;
    // const Double_t max = 5000000;
    const Double_t min = histo2->GetMaximum()*0.6;
    const Double_t max = histo2->GetMaximum()*1.1;

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
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity, this thesis");

  gPad->SaveAs("Unfolding2D/unfolded-reconstructed_phimodulation.pdf",  "RECREATE");





















  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);

  histo3->GetXaxis()->SetTitleOffset(1.15);
  histo3->GetYaxis()->SetTitleOffset(0.9);
  histo3->GetXaxis()->SetTitleSize(0.055);
  histo3->GetYaxis()->SetTitleSize(0.055);
  histo3->GetXaxis()->SetLabelSize(0.05);
  histo3->GetYaxis()->SetLabelSize(0.05);
  histo3->GetXaxis()->SetTitleFont(42);
  histo3->GetYaxis()->SetTitleFont(42);
  histo3->GetXaxis()->SetLabelFont(42);
  histo3->GetYaxis()->SetLabelFont(42);
  histo3->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
  histo3->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptStat(0);
  histo3->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
  histo3->DrawClone("col");// draw "axes", "contents", "statistics box"
  histo3->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
  histo3->Draw("colz same"); // draw the "color palette"
  // Model->Draw("same");



  TLatex* latex3 = new TLatex();
  latex3->SetTextSize(0.055);
  latex3->SetTextFont(42);
  latex3->SetTextAlign(11);
  latex3->SetNDC();
  latex3->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex3->SetTextSize(0.055);
  latex3->DrawLatex(0.7,0.94,"Helicity generated, this thesis");

  gPad->SaveAs("Unfolding2D/generated_phimodulation.pdf",  "RECREATE");






  histo4->Add(histo2, -1.);
  histo4->Divide(histo2);
  histo4->Scale(-1.);
  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  histo4->GetXaxis()->SetTitleOffset(1.15);
  histo4->GetYaxis()->SetTitleOffset(0.9);
  histo4->GetXaxis()->SetTitleSize(0.055);
  histo4->GetYaxis()->SetTitleSize(0.055);
  histo4->GetXaxis()->SetLabelSize(0.05);
  histo4->GetYaxis()->SetLabelSize(0.05);
  histo4->GetXaxis()->SetTitleFont(42);
  histo4->GetYaxis()->SetTitleFont(42);
  histo4->GetXaxis()->SetLabelFont(42);
  histo4->GetYaxis()->SetLabelFont(42);
  histo4->GetXaxis()->SetNdivisions(408);
  cout << "OK2" << endl << flush;
  histo4->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );
  gStyle->SetNumberContours(NCont);
  gStyle->SetOptStat(0);
  const Double_t min2 = histo4->GetMaximum()*(-2.);
  const Double_t max2 = histo4->GetMaximum()*1.2;
  for(int i = 1; i < nLevels; i++) {
    levels[i] = min2 + (max2 - min2) / (nLevels - 1) * (i);
  }
  levels[0] = 0.01;
  // histo4->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
  histo4->GetZaxis()->SetRangeUser(min2, max2); // ... set the range ...
  Int_t myPalette[20] ;  // 1=black / 9=blue / 0=white
  // myPalette[0] =  0 ; // -10 ->  0
  // myPalette[1] =  1 ; //   0 -> 10
  // myPalette[2] =  2 ; //  10 -> 20
  // myPalette[3] =  3 ; //  20 -> 30
  // myPalette[4] =  4 ; //  30 -> 40
  // myPalette[5] =  5 ; //  40 -> 50
  for (Int_t i = 0; i < 20; i++) {
    myPalette[i] =  ((Double_t)i+1.+10)*7./30. ;
  }
  gStyle->SetPalette(20,myPalette) ;
  histo4->DrawClone("col");// draw "axes", "contents", "statistics box"
  // histo4->GetZaxis()->SetRangeUser(min2, max2); // ... set the range ...
  // histo4->Draw("colz same"); // draw the "color palette"
  histo4->Draw("COLZ1"); // draw the "color palette"
  // Model->Draw("same");

  TLatex* latex4 = new TLatex();
  latex4->SetTextSize(0.055);
  latex4->SetTextFont(42);
  latex4->SetTextAlign(11);
  latex4->SetNDC();
  latex4->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex4->SetTextSize(0.055);
  latex4->DrawLatex(0.7,0.94,"Helicity residual, this thesis");

  gPad->SaveAs("Unfolding2D/residuals-2D_phimodulation.pdf",  "RECREATE");
  new TCanvas;
  histo4->Draw("surf");
}
//_____________________________________________________________________________
/* -
   -
 */
void DrawUncertainties(){


  TFile* f[24];
  TH1F*  h[24];
  TH1F*  h2[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i]  = (TH1F*) f[i]->Get(Form("histo4_%d", i));
  }





  TH2F* histo2 = new TH2F("histo2", "histo2", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  TH2F* histo3 = new TH2F("histo3", "histo3", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
  TH2F* histo4 = new TH2F("histo4", "histo4", 24, 0., 2.*TMath::Pi(), 24, -1., 1.);
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


  // for (size_t iCosThetaBins = 7; iCosThetaBins < 18; iCosThetaBins++) {
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 17; iCosThetaBins++) {
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
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) );
      histo2->SetBinError( binx2, biny2, 0.00000000000001 );
      histo3->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );
      histo3->SetBinError( binx2, biny2, 0.00000000000001 );
      histo4->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );
      histo4->SetBinError( binx2, biny2, 0.00000000000001 );

    }
  }
















  histo4->Divide(histo2);
  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  histo4->GetXaxis()->SetTitleOffset(1.15);
  histo4->GetYaxis()->SetTitleOffset(0.9);
  histo4->GetXaxis()->SetTitleSize(0.055);
  histo4->GetYaxis()->SetTitleSize(0.055);
  histo4->GetXaxis()->SetLabelSize(0.05);
  histo4->GetYaxis()->SetLabelSize(0.05);
  histo4->GetXaxis()->SetTitleFont(42);
  histo4->GetYaxis()->SetTitleFont(42);
  histo4->GetXaxis()->SetLabelFont(42);
  histo4->GetYaxis()->SetLabelFont(42);
  histo4->GetXaxis()->SetNdivisions(408);
  cout << "OK2" << endl << flush;
  histo4->SetTitle( "; #it{#varphi} ;cos(#it{#theta})" );
  gStyle->SetOptStat(0);
  const Double_t min2 = histo4->GetMaximum()*0.001;
  const Double_t max2 = histo4->GetMaximum()*1.2;
  histo4->GetZaxis()->SetRangeUser(min2, max2); // ... set the range ...
  Int_t myPalette[20] ;  // 1=black / 9=blue / 0=white
  for (Int_t i = 0; i < 20; i++) {
    myPalette[i] =  ((Double_t)i+1.+10)*7./30. ;
  }
  gStyle->SetPalette(20,myPalette) ;
  histo4->DrawClone("col");// draw "axes", "contents", "statistics box"
  // histo4->GetZaxis()->SetRangeUser(min2, max2); // ... set the range ...
  // histo4->Draw("colz same"); // draw the "color palette"
  histo4->Draw("COLZ1"); // draw the "color palette"
  // Model->Draw("same");

  TLatex* latex4 = new TLatex();
  latex4->SetTextSize(0.055);
  latex4->SetTextFont(42);
  latex4->SetTextAlign(11);
  latex4->SetNDC();
  latex4->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex4->SetTextSize(0.055);
  latex4->DrawLatex(0.7,0.94,"Helicity uncertainties, this thesis");

  gPad->SaveAs("Unfolding2D/uncertainties_phimodulation.pdf",  "RECREATE");
  new TCanvas;
  histo4->Draw("surf");
}
//_____________________________________________________________________________
void DrawRefolded(Int_t Iterations = 1)
{
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);


  // Double_t ChiSquareTest = 0;
  // TFile* file = new TFile("Unfolding2D/UnfoldedClosureHe_phimodulation_10.root");
  // TH1F* UnlikeSignDimuon2    = (TH1F*)file->Get("histo5");
  // TH1F* UnlikeSignDimuon3    = (TH1F*)file->Get("histo6");
  // for (size_t i = 0; i < 230; i++) {
  //   if (UnlikeSignDimuon3->GetBinContent(i+1) > 0.0000000001){
  //     ChiSquareTest += UnlikeSignDimuon3->GetBinContent(i+1);
  //   }
  // }
  TFile* file = new TFile("Unfolding2D/UnfoldedClosureHe_phimodulation_10.root");
  TH1F* UnlikeSignDimuon2    = (TH1F*)file->Get("histo5");
  Double_t ChiSquareTest = 0;
  TH1F* UnlikeSignDimuon3    = (TH1F*)file->Get("histo6");
  for (size_t i = 0; i < 230; i++) {
    if (UnlikeSignDimuon3->GetBinContent(i+1) > 0.0000000001){
      ChiSquareTest += UnlikeSignDimuon3->GetBinContent(i+1);
    }
  }
  TH1F* UnlikeSignDimuon    = new TH1F("UnlikeSignDimuon","UnlikeSignDimuon", 300, -0.5, 299.5);
  Int_t j = 1;
  for (Int_t i = 0; i < 600; i++) {
    if (UnlikeSignDimuon2->GetBinContent(i+1) > 0.0000000001){
      UnlikeSignDimuon->SetBinError(j, 0.);
      UnlikeSignDimuon->SetBinContent(j, UnlikeSignDimuon2->GetBinContent(i+1));
      j++;
    }
  }
  UnlikeSignDimuon->SetTitle("");
  UnlikeSignDimuon->GetXaxis()->SetTitleOffset(1.15);
  UnlikeSignDimuon->GetYaxis()->SetTitleOffset(1.25);
  UnlikeSignDimuon->GetXaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon->GetYaxis()->SetTitleSize(0.045);
  UnlikeSignDimuon->GetXaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon->GetYaxis()->SetLabelSize(0.045);
  UnlikeSignDimuon->GetXaxis()->SetTitleFont(42);
  UnlikeSignDimuon->GetYaxis()->SetTitleFont(42);
  UnlikeSignDimuon->GetXaxis()->SetLabelFont(42);
  UnlikeSignDimuon->GetYaxis()->SetLabelFont(42);

  UnlikeSignDimuon->GetXaxis()->SetTitle("Bin number");
  UnlikeSignDimuon->GetYaxis()->SetTitle("Ratio of refolded to REC modulation");
  UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.9,1.1);
  UnlikeSignDimuon->GetXaxis()->SetRangeUser(-0.5 , 229.5);
  UnlikeSignDimuon->SetLineWidth(5);
  UnlikeSignDimuon->SetLineColor(2);
  // UnlikeSignDimuon->SetFillColor(kRed-3);
  // UnlikeSignDimuon->SetFillStyle(1001);
  // UnlikeSignDimuon->Draw("ep");
  UnlikeSignDimuon->Draw("");




  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"LHC18l7, PbPb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.80,"1 iteration");}
  else {latex5->DrawLatex(0.31,0.70,Form("This thesis, %d iterations", Iterations));}
  latex5->DrawLatex(0.2,0.80,Form("#chi^{2} = #sum (Refolded - REC)^{2}/REC = %0.02f", ChiSquareTest));

  gPad->SaveAs(Form("Unfolding2D/residuals-closure-modulation-refolded-%d.pdf", Iterations), "recreate");
}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void FitGeneratedModulation(){


  TFile* f = new TFile("SignalExtraction/SimpleClosure_phimodulation.root");
  TH2F*  histo2 = (TH2F*) f->Get("generatedlimit");


  /// reset data structure
  coordsX = std::vector<Double_t>();
  coordsY = std::vector<Double_t>();
  values  = std::vector<Double_t>();
  errors  = std::vector<Double_t>();
  /// fill data structure
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 18; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 0; iCosThetaBins < 23; iCosThetaBins++) {
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
    for (size_t iPhiBins = 0; iPhiBins < 24; iPhiBins++) {



      Int_t binx = -1;
      cout << "binx  = " << binx << endl;
      Int_t biny2 = histo2->GetYaxis()->FindBin(-1.+((Double_t)iCosThetaBins + 0.5) * (0.08+0.01/3.));
      Int_t binx2 = histo2->GetXaxis()->FindBin(((Double_t)iPhiBins + 0.5) * (2.*TMath::Pi())/24.);

      coordsX.push_back( ((TAxis*) histo2->GetXaxis())->GetBinCenter( binx2 ) );
      coordsY.push_back( ((TAxis*) histo2->GetYaxis())->GetBinCenter( biny2 ) );
      values.push_back(  histo2->GetBinContent( binx2, biny2)        );
      errors.push_back(  histo2->GetBinError( binx2, biny2)       );
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        errorsBin.push_back(  TMath::Pi()/TMath::Sqrt(12.)     );
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        errorsBin.push_back(  TMath::Pi()/(6.*TMath::Sqrt(12.))    );
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        errorsBin.push_back(  TMath::Pi()/(12.*TMath::Sqrt(12.))    );
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
        errorsBin.push_back(  TMath::Pi()/(24.*TMath::Sqrt(12.))    );
      }
      // errorsBin.push_back(0);



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


  cout << "Coordinates:" << endl;
  for( Int_t i = 0; i < coordsX.size(); i++ ){
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
  }





  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationWithBins);
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.,    -2, 2        );
  // gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.1,    -2, 2        );
  gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.,    -2, 2        );
  gMinuit->DefineParameter(1, "LambdaPhi",           0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.1,    -2, 2        );
  // gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.,    -2, 2        );
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

















  TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -1. ,1., 4 );
  Model->FixParameter(0, LambdaTheta);
  Model->FixParameter(1, LambdaPhi);
  Model->FixParameter(2, LambdaThetaPhi);
  Model->FixParameter(3, Normalisation);
  TF2* Model2 = new TF2("Model2", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  Model2->FixParameter(0, 1.);
  Model2->FixParameter(1, LambdaPhi);
  Model2->FixParameter(2, LambdaThetaPhi);
  Model2->FixParameter(3, Normalisation);

  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  histo2->GetXaxis()->SetTitleOffset(1.15);
  histo2->GetYaxis()->SetTitleOffset(0.9);
  histo2->GetXaxis()->SetTitleSize(0.055);
  histo2->GetYaxis()->SetTitleSize(0.055);
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleFont(42);
  histo2->GetYaxis()->SetTitleFont(42);
  histo2->GetXaxis()->SetLabelFont(42);
  histo2->GetYaxis()->SetLabelFont(42);
  histo2->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
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
    const Double_t min = histo2->GetMaximum()*0.001;
    const Double_t max = histo2->GetMaximum()*1.1;
    // const Double_t min = 600000;
    // const Double_t max = 900000;

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
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");

  gPad->SaveAs("Unfolding2D/generated_phimodulation.pdf",  "RECREATE");


  new TCanvas;
  histo2->Draw("surf same");
  Model->Draw("surf same");

  new TCanvas;
  Model->Draw("surf same");

}
//_____________________________________________________________________________
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void FitGeneratedModulationWithSameBinning(Int_t FitRangeMode = 0 ){
// void PolarisationHeMinuit2D( Int_t SignalRangeSelectionMode = 0, Int_t FitRangeMode = 0 ){


  TFile* f[24];
  TH1F*  h[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_phimodulation_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i] = (TH1F*) f[i]->Get(Form("hg_%d", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo3_%d", i));
  }




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
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 18; iCosThetaBins++) {
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 9; iCosThetaBins < 15; iCosThetaBins++) {
    // if        ( FitRangeMode == 1 ) {
    //   if (iCosThetaBins      == 5)   continue;
    // } else if ( FitRangeMode == 2 ) {
    //   if (iCosThetaBins      == 18)  continue;
    // } else if ( FitRangeMode == 3 ) {
    //   if ( (iCosThetaBins    == 5) || (iCosThetaBins == 19) )  continue;
    // } else {
    // }
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
      Double_t M = 1.;
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
      //   M = 1.;
      // } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
      //   M = 6.;
      // } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
      //   M = 12.;
      // } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
      //             iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
      //   M = 24.;
        M = 1.;
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        M = 6.;
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        M = 12.;
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
        M = 24.;

      }

      // histo->Fill( PhiCenters[iPhiBins], CosThetaCenters[iCosThetaBins], h[iCosThetaBins]->GetBinContent(binx) );
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) * M );
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) * M );


      // coordsX.push_back( ((TAxis*) Distr2D->GetXaxis())->GetBinCenter( ix ) );
      // coordsY.push_back( CosThetaCenters[iCosThetaBins] );
      // values.push_back(  h[iCosThetaBins]->GetBinContent(binx)        );
      // errors.push_back(  h[iCosThetaBins]->GetBinError(binx)        );
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        if (iPhiBins == 0) {
          coordsX.push_back( TMath::Pi() );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx) *M       );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)  *M      );
        }
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        cout << "iPhibins" << iPhiBins << endl;
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx) *M       );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)  *M      );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/12. + 4.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx) *M       );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)  *M      );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/12. + 8.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)  *M      );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/12. + 12.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)  *M      );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/12. + 16.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)   *M     );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/12. + 20.*TMath::Pi()/12. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        }
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        if(iPhiBins == 0){
          coordsX.push_back( 2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 2){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*2.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)    *M    );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 4){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*4.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 6){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*6.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)   *M     );
        } else if (iPhiBins == 8){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*8.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)    *M    );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)     *M   );
        } else if (iPhiBins == 10){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*10.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 12){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*12.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 14){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*14.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)   *M     );
        } else if (iPhiBins == 16){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*16.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 18){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*18.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)     *M   );
        } else if (iPhiBins == 20){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*20.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)   *M     );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)   *M     );
        } else if (iPhiBins == 22){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*22.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)    *M    );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)    *M    );
        } else if (iPhiBins == 24){
          coordsX.push_back( 2.*TMath::Pi()/24. + 2.*24.*TMath::Pi()/24. );
          coordsY.push_back( CosThetaCenters[iCosThetaBins] );
          values.push_back(  h[iCosThetaBins]->GetBinContent(binx)    *M    );
          errors.push_back(  h[iCosThetaBins]->GetBinError(binx)     *M   );
        }
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 )
      {
        coordsX.push_back( PhiCenters[iPhiBins] );
        coordsY.push_back( CosThetaCenters[iCosThetaBins] );
        values.push_back(  h[iCosThetaBins]->GetBinContent(binx)  *M      );
        errors.push_back(  h[iCosThetaBins]->GetBinError(binx)     *M   );
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


  cout << "Coordinates:" << endl;
  for( Int_t i = 0; i < coordsX.size(); i++ ){
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
  }





  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationWithBins);
  gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.,    -2, 2        );
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.1,    -2, 2        );
  // gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.,    -2, 2        );
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

















  TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -0.5 ,0.5, 4 );
  Model->FixParameter(0, LambdaTheta);
  Model->FixParameter(1, LambdaPhi);
  Model->FixParameter(2, LambdaThetaPhi);
  Model->FixParameter(3, Normalisation);
  TF2* Model2 = new TF2("Model2", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  Model2->FixParameter(0, 1.);
  Model2->FixParameter(1, LambdaPhi);
  Model2->FixParameter(2, LambdaThetaPhi);
  Model2->FixParameter(3, Normalisation);

  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  // Distr2D->GetXaxis()->SetTitleOffset(1.15);
  // // Distr2D->GetYaxis()->SetTitleOffset(1.25);
  // Distr2D->GetYaxis()->SetTitleOffset(1.);
  // Distr2D->GetXaxis()->SetTitleSize(0.045);
  // Distr2D->GetYaxis()->SetTitleSize(0.045);
  // Distr2D->GetXaxis()->SetLabelSize(0.045);
  // Distr2D->GetYaxis()->SetLabelSize(0.045);
  // Distr2D->GetXaxis()->SetTitleFont(42);
  // Distr2D->GetYaxis()->SetTitleFont(42);
  // Distr2D->GetXaxis()->SetLabelFont(42);
  // Distr2D->GetYaxis()->SetLabelFont(42);
  // Distr2D->GetXaxis()->SetNdivisions(408);

  histo2->GetXaxis()->SetTitleOffset(1.15);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.25);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleOffset(1.);
  histo2->GetYaxis()->SetTitleOffset(0.9);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetTitleSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetXaxis()->SetLabelSize(0.045);
  // helicity2DafterSignalExtractionErrors->GetYaxis()->SetLabelSize(0.045);
  histo2->GetXaxis()->SetTitleSize(0.055);
  histo2->GetYaxis()->SetTitleSize(0.055);
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleFont(42);
  histo2->GetYaxis()->SetTitleFont(42);
  histo2->GetXaxis()->SetLabelFont(42);
  histo2->GetYaxis()->SetLabelFont(42);
  histo2->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
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
    const Double_t min = histo2->GetMaximum()*0.001;
    const Double_t max = histo2->GetMaximum()*1.1;
    // const Double_t min = 600000;
    // const Double_t max = 900000;

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
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");

  gPad->SaveAs("Unfolding2D/Inverted2Dmaps_phimodulation.pdf",  "RECREATE");


  new TCanvas;
  histo2->Draw("surf same");
  Model->Draw("surf same");


  for( Int_t i = 0; i < coordsX.size(); i++ ){
    Double_t TraslatedCosThetaGen = 0.5*(coordsY[i] + 1.)*24.;
    Double_t iCosThetaBins2 = -1;
    Double_t RemainderCosTheta = 0.;
    RemainderCosTheta = modf(TraslatedCosThetaGen, &iCosThetaBins2);
    Int_t iCosThetaBins = (Int_t)  iCosThetaBins2;

    Double_t Mn = 1.;
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
    Double_t TraslatedPhiGen = coordsX[i]*Mn/(2.*TMath::Pi());
    Double_t iPhiBins2 = -1;
    Double_t RemainderPhi = 0.;
    RemainderPhi = modf(TraslatedPhiGen, &iPhiBins2);
    Int_t iPhiBins = (Int_t)  iPhiBins2;


    Double_t xn[2] = {(Double_t)iPhiBins, (Double_t)iCosThetaBins};
    Double_t yn[4] = {LambdaTheta,LambdaPhi,LambdaThetaPhi,Normalisation};

    cout << "==============================================" << endl;
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
    cout << "As computed in the function? " << endl;
    cout << "CosThetaF = " << (-1.+((Double_t) iCosThetaBins  +0.5   )*(0.08+0.01/3.)) << ", PhiF =  " << (((Double_t) iPhiBins +0.5    )*2.*TMath::Pi()/Mn) <<  endl;
    cout << "iCosThetaBins = " << iCosThetaBins << endl << ",   iPhiBins = " << iPhiBins << endl;
    cout << "TraslatedPhiGen = " << TraslatedPhiGen << endl << ",   Mn = " << Mn << endl;
    cout << "IntegralFunction = " << IntegralFunction(xn,yn) << endl;
    cout << "Model            = " << Model->Eval(coordsX[i], coordsY[i]) << endl;
  }
  printf("LambdaTheta     : %+.7f +- %.7f\n",LambdaTheta,     LambdaThetaErr     );
  printf("LambdaPhi       : %+.7f +- %.7f\n",LambdaPhi,       LambdaPhiErr       );
  printf("LambdaThetaPhi  : %+.7f +- %.7f\n",LambdaThetaPhi,  LambdaThetaPhiErr  );
  printf("Normalisation   : %+.7f +- %.7f\n",Normalisation,   NormalisationErr   );





  TH1F* pulls     = new TH1F("pulls",     "pulls",     500, -0.5, 499.5);
  TH1F* valuehis  = new TH1F("valuehis",  "valuehis",  500, -0.5, 499.5);
  TH1F* modelhis  = new TH1F("modelhis",  "modelhis",  500, -0.5, 499.5);
  TH1F* residuals = new TH1F("residuals", "residuals", 500, -0.5, 499.5);
  for (size_t i = 0; i < ndf; i++) {
    pulls    ->SetBinContent(i+1,  values[i] - Model->Eval(coordsX[i], coordsY[i]));
    valuehis ->SetBinContent(i+1,  values[i] );
    modelhis ->SetBinContent(i+1,  Model->Eval(coordsX[i], coordsY[i]) );
    residuals->SetBinContent(i+1, (values[i] - Model->Eval(coordsX[i], coordsY[i]))/errors[i]);

  }
  new TCanvas;
  pulls->Draw();
  new TCanvas;
  residuals->Draw();
  new TCanvas;
  valuehis->Draw();
  new TCanvas;
  modelhis->Draw();







  TH1F* pulls2     = new TH1F("pulls2",     "pulls2",     500, -0.5, 499.5);
  TH1F* valuehis2  = new TH1F("valuehis2",  "valuehis2",  500, -0.5, 499.5);
  TH1F* modelhis2  = new TH1F("modelhis2",  "modelhis2",  500, -0.5, 499.5);
  TH1F* residuals2 = new TH1F("residuals2", "residuals2", 500, -0.5, 499.5);
  for (size_t i = 0; i < ndf; i++) {
    pulls2    ->SetBinContent(i+1,  values[i] - Model2->Eval(coordsX[i], coordsY[i]));
    valuehis2 ->SetBinContent(i+1,  values[i] );
    modelhis2 ->SetBinContent(i+1,  Model2->Eval(coordsX[i], coordsY[i]) );
    residuals2->SetBinContent(i+1, (values[i] - Model2->Eval(coordsX[i], coordsY[i]))/errors[i]);

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
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void FitGeneratedModulationFinerBinning(){


  TFile* f = new TFile("SignalExtraction/SimpleClosure_phimodulation.root");
  TH2F*  histo2 = (TH2F*) f->Get("generatedlimit2");


  /// reset data structure
  coordsX = std::vector<Double_t>();
  coordsY = std::vector<Double_t>();
  values  = std::vector<Double_t>();
  errors  = std::vector<Double_t>();
  /// fill data structure
  // for (size_t iCosThetaBins = 7; iCosThetaBins < 18; iCosThetaBins++) {
  for (size_t iCosThetaBins = 20; iCosThetaBins < 80; iCosThetaBins++) {
    for (size_t iPhiBins = 0; iPhiBins < 100; iPhiBins++) {



      Int_t binx = -1;
      cout << "binx  = " << binx << endl;
      Int_t biny2 = histo2->GetYaxis()->FindBin(-1.+((Double_t)iCosThetaBins + 0.5) * (0.02));
      Int_t binx2 = histo2->GetXaxis()->FindBin(((Double_t)iPhiBins + 0.5) * (2.*TMath::Pi())/100.);

      coordsX.push_back( ((TAxis*) histo2->GetXaxis())->GetBinCenter( binx2 ) );
      coordsY.push_back( ((TAxis*) histo2->GetYaxis())->GetBinCenter( biny2 ) );
      values.push_back(  histo2->GetBinContent( binx2, biny2)        );
      errors.push_back(  histo2->GetBinError( binx2, biny2)       );


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


  cout << "Coordinates:" << endl;
  for( Int_t i = 0; i < coordsX.size(); i++ ){
    cout << "CosTheta = " << coordsY[i] << ", Phi =  " << coordsX[i] << ", Value =  " << values[i] << "+/-" << errors[i] << endl;
  }





  TMinuit *gMinuit = new TMinuit(4);
  // gMinuit->SetFCN(FcnForMinimisation);
  gMinuit->SetFCN(FcnForMinimisationWithBins);
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.,    -2, 2        );
  gMinuit->DefineParameter(0, "LambdaTheta",        0., 0.,    -2, 2        );
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.1,    -2, 2        );
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

















  TF2* Model = new TF2("Model", helicity2Dv3,0., 2.*TMath::Pi(), -1. ,1., 4 );
  Model->FixParameter(0, LambdaTheta);
  Model->FixParameter(1, LambdaPhi);
  Model->FixParameter(2, LambdaThetaPhi);
  Model->FixParameter(3, Normalisation);
  TF2* Model2 = new TF2("Model2", helicity2Dv3,0., 2.*TMath::Pi(), -0.6 ,0.6, 4 );
  Model2->FixParameter(0, 1.);
  Model2->FixParameter(1, LambdaPhi);
  Model2->FixParameter(2, LambdaThetaPhi);
  Model2->FixParameter(3, Normalisation);

  new TCanvas;
  gPad->SetMargin(0.13,0.13,0.12,0.12);
  histo2->GetXaxis()->SetTitleOffset(1.15);
  histo2->GetYaxis()->SetTitleOffset(0.9);
  histo2->GetXaxis()->SetTitleSize(0.055);
  histo2->GetYaxis()->SetTitleSize(0.055);
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleFont(42);
  histo2->GetYaxis()->SetTitleFont(42);
  histo2->GetXaxis()->SetLabelFont(42);
  histo2->GetYaxis()->SetLabelFont(42);
  histo2->GetXaxis()->SetNdivisions(408);

  cout << "OK2" << endl << flush;
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
    const Double_t min = histo2->GetMaximum()*0.001;
    const Double_t max = histo2->GetMaximum()*1.1;
    // const Double_t min = 600000;
    // const Double_t max = 900000;

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
  latex->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex->SetTextSize(0.055);
  latex->DrawLatex(0.7,0.94,"Helicity");

  gPad->SaveAs("Unfolding2D/generated_phimodulation.pdf",  "RECREATE");


  new TCanvas;
  histo2->Draw("surf same");
Model->Draw("surf same");

}
//______________________________________________
void DrawResidualV2(){

  TFile* f[24];
  TFile* f2[24];
  TH1F*  h[24];
  TH1F*  h2[24];
  TH1F*  h3[24];
  TH1F*  h8[24];
  Int_t Counter = 1;
  Double_t Values[250];
  Double_t ChiSquareTest = 0;
  Double_t M = 1.;
  TH1F* DifferencesH = new TH1F("DifferencesH","DifferencesH", 300, -0.5, 299.5);
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    f2[i] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_phimodulation_%d.root", i));
    // h[i]  = (TH1F*) f[i]->Get(Form("histo3_%d", i));
    h[i]  = (TH1F*) f2[i]->Get(Form("hg_%d", i));
    h[i]->Sumw2();
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      M = 1.;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      M = 6.;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      M = 12.;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      M = 24.;

    }
    h[i]->Scale(M);
    h2[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    h3[i] = (TH1F*) h[i]->Clone(Form("histo7_%d", i));
    h3[i]->Sumw2();
    h2[i]->Sumw2();
    h[i]->Sumw2();
    h3[i] -> Add(h[i], -1.);
    h3[i]->Scale(-1.);
    h8[i] = (TH1F*) h3[i]->Clone(Form("histo8_%d", i));
    h8[i]->Sumw2();
    // h8[i] ->Divide(h3[i],h2[i], 1,1, "B");
    h8[i] ->Divide(h3[i],h[i], 1,1, "B");
    // for (Int_t j = 0; j < 30; j++) {
    //   if
    // }
    // if (i == 4 || i == 5) {
      for (Int_t j = 1; j < 31; j++) {
        if ((h2[i]->GetBinContent(j) > 0.00000000000001) && ((h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))/h[i]->GetBinContent(j) > 0.5)) {
        // if (h2[i]->GetBinContent(j) > 0.00000000000001) {
          cout << "=========================================================== " << endl ;
          cout << "i = " << i << ", j = " << j << endl;
          cout << "Bin content unfolded rec = " <<  h2[i]->GetBinContent(j) << endl ;
          cout << "Bin content generated    = " <<  h[i]->GetBinContent(j) << endl ;
          cout << "Bin content diff         = " <<  (h2[i]->GetBinContent(j)-h[i]->GetBinContent(j)) << endl ;
          cout << "Bin content diff/gen     = " <<  (h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))/h[i]->GetBinContent(j) << endl ;
          cout << "=========================================================== " << endl ;
          cout << "h2-h"      << h3[i]->GetBinContent(j) << endl;
          cout << "(h2-h)/h2" << h8[i]->GetBinContent(j) << endl;
          Counter++;
        }
      // }
    }
    // h3[i] ->Divide(h2[i]);
    h3[i] ->Divide(h[i]);
    for (Int_t j = 1; j < 31; j++) {
      if (h[i]->GetBinContent(j) > 0.00000000000001) {
        DifferencesH->SetBinContent(Counter, (h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))/h[i]->GetBinContent(j));
        // DifferencesH->SetBinError(Counter, h8[i]->GetBinError(j));
        DifferencesH->SetBinError(Counter, 0.000000000000000000);
        Values[Counter] = (h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))*(h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))/h[i]->GetBinContent(j);
        cout << "Values[" << Counter << "] = " << Values[Counter] << endl;
        ChiSquareTest += Values[Counter];
        Counter++;
      }
    }
  }
  TH1F* histogram = (TH1F*) DifferencesH->Clone("histogram");

  // Int_t M = 1;
  // TH1F* histogram = new TH1F("histogram","histogram", 1000, -0.5,999.5);
  // for (Int_t iC = 4; iC < 20; iC++) {
  //   if( iC == 0  || iC == 1  || iC == 2  || iC == 3  ||
  //       iC == 4  || iC == 23 || iC == 22 || iC == 21 ||
  //       iC == 20 || iC == 19 )
  //   {
  //     M = 1;
  //   } else if ( iC == 5  || iC == 6  || iC == 18  || iC == 17 ) {
  //     M = 6;
  //   } else if ( iC == 7  || iC == 8  || iC == 16  || iC == 15 ) {
  //     M = 12;
  //   } else if ( iC == 9  || iC == 10  || iC == 11  ||
  //               iC == 12 || iC == 13  || iC == 14 ) {
  //     M = 24;
  //   }
  //   for (Int_t j = 0; j < M; j++) {
  //     // if( iC == 0  || iC == 1  || iC == 2  || iC == 3  ||
  //     //     iC == 4  || iC == 23 || iC == 22 || iC == 21 ||
  //     //     iC == 20 || iC == 19 )
  //     // {
  //     //   M = 1;
  //     // } else if ( iC == 5  || iC == 6  || iC == 18  || iC == 17 ) {
  //     //   M = 6;
  //     // } else if ( iC == 7  || iC == 8  || iC == 16  || iC == 15 ) {
  //     //   M = 12;
  //     // } else if ( iC == 9  || iC == 10  || iC == 11  ||
  //     //             iC == 12 || iC == 13  || iC == 14 ) {
  //     //   M = 24;
  //     // }
  //     histogram->SetBinContent( iC*30+j, h3[iC]->GetBinContent( j+1 ) );
  //     histogram->SetBinError(   iC*30+j, 0. );
  //     // histogram->SetBinError(   iC*30+j, h3[iC]->GetBinError( j+1 ) );
  //   }
  // }
  new TCanvas;
  histogram->GetYaxis()->SetRangeUser(0.95,1.05);
  histogram->Draw("*H");


  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);

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

  histogram->GetXaxis()->SetTitle("Bin number");
  histogram->GetYaxis()->SetTitle("(Unfolded REC - GEN)/GEN [a.u.]");
  // histogram->GetYaxis()->SetRangeUser(-0.006,0.004);
  histogram->GetYaxis()->SetRangeUser(-0.010,0.010);
  histogram->GetXaxis()->SetRangeUser(-0.5, 229.5);
  // histogram->GetXaxis()->SetRangeUser(100, 600);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");

  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis");
  latex5->DrawLatex(0.31,0.80,"1 iteration");
  latex5->DrawLatex(0.2,0.64,Form("#chi^{2} = #sum (Unfolded REC - GEN)^{2}/GEN = %0.02f", ChiSquareTest));

  gPad->SaveAs("residuals-closure-v2-modulation-1.pdf", "recreate");


}
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
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void DrawDifferencesGenUnfolded(Int_t FitRangeMode = 0 ){


  TFile* f[24];
  TFile* f2[24];
  TH1F*  h[24];
  TH1F*  h2[24];
  Int_t Counter = 1;
  Double_t Values[250];
  Double_t ChiSquareTest = 0;
  TH1F* DifferencesH = new TH1F("DifferencesH","DifferencesH", 300, -0.5, 299.5);
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_phimodulation_%d.root", i));
    h[i] = (TH1F*) f[i]->Get(Form("hg_%d", i));
    f2[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    // h2[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h2[i] = (TH1F*) f2[i]->Get(Form("histo4_%d", i));
    for (Int_t j = 1; j < 30; j++) {
      if (  h2[i]->GetBinContent(j) > 0.0001 ) {
        DifferencesH->SetBinContent(Counter, (h2[i]->GetBinContent(j) - h[i]->GetBinContent(j))/h2[i]->GetBinContent(j) );
        Values[Counter] = (h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))*(h2[i]->GetBinContent(j)-h[i]->GetBinContent(j))/h[i]->GetBinContent(j);
        cout << "Values[" << Counter << "] = " << Values[Counter] << endl;
        ChiSquareTest += Values[Counter];
        Counter++;
      }
    }
  }


  TH1F* histogram = (TH1F*) DifferencesH->Clone("histogram");
  new TCanvas;
  histogram->GetYaxis()->SetRangeUser(0.95,1.05);
  histogram->Draw("*H");


  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);

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

  histogram->GetXaxis()->SetTitle("Bin number");
  histogram->GetYaxis()->SetTitle("(Unfolded REC - GEN)/GEN [a.u.]");
  // histogram->GetYaxis()->SetRangeUser(-0.006,0.004);
  histogram->GetYaxis()->SetRangeUser(0.50,1.30);
  histogram->GetXaxis()->SetRangeUser(-0.5, 229.5);
  // histogram->GetXaxis()->SetRangeUser(100, 600);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");

  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis");
  latex5->DrawLatex(0.31,0.80,"1 iteration");
  latex5->DrawLatex(0.2,0.64,Form("#chi^{2} = #sum (Unfolded REC - GEN)^{2}/GEN = %0.f", ChiSquareTest));


  gPad->SaveAs("Unfolding2D/differences-modulation-2d-1.pdf", "recreate");
  histogram->GetXaxis()->SetRangeUser(29.5, 180.5);
  histogram->GetYaxis()->SetRangeUser(0.85,1.20);
  gPad->SaveAs("Unfolding2D/differences-modulation-2d-1-zoom.pdf", "recreate");



}
//______________________________________________
void DrawResidualV3(Int_t Iterations){

  TFile* f[24];
  TFile* f2[24];
  TH1F*  GenH[24];
  TH1F*  UnfoldedH[24];
  TH1F*  GenCopyH[24];
  Int_t Counter = 1;
  Double_t Values[250];
  Double_t ChiSquareTest = 0;
  Double_t ChiSquareTest2 = 0;
  Double_t M = 1.;
  TH1F* DifferencesH = new TH1F("DifferencesH","DifferencesH", 300, -0.5, 299.5);
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_phimodulation_%d.root", i));
    f2[i] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_phimodulation_%d.root", i));
    GenH[i]  = (TH1F*) f2[i]->Get(Form("hg_%d", i));
    GenH[i]->Sumw2();
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      M = 1.;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      M = 6.;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      M = 12.;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      M = 24.;

    }
    GenH[i]->Scale(M);
    UnfoldedH[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
    UnfoldedH[i]->Sumw2();
    GenCopyH[i] = (TH1F*) GenH[i]->Clone(Form("histo7_%d", i));
    GenCopyH[i]->Sumw2();
    UnfoldedH[i]->Sumw2();
    // ChiSquareTest2 += calcChi2(GenCopyH[i], UnfoldedH[i], N);
    // GenCopyH[i] -> Add(GenH[i], -1.);
    // GenCopyH[i]->Scale(-1.);
    // h8[i] ->Divide(GenCopyH[i],GenH[i], 1,1, "B");
    // for (Int_t j = 1; j < 31; j++) {
    //     if ((UnfoldedH[i]->GetBinContent(j) > 0.00000000000001) && ((UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j) > 0.5)) {
    //       cout << "=========================================================== " << endl ;
    //       cout << "i = " << i << ", j = " << j << endl;
    //       cout << "Bin content unfolded rec = " <<  UnfoldedH[i]->GetBinContent(j) << endl ;
    //       cout << "Bin content generated    = " <<  GenH[i]->GetBinContent(j) << endl ;
    //       cout << "Bin content diff         = " <<  (UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j)) << endl ;
    //       cout << "Bin content diff/gen     = " <<  (UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j) << endl ;
    //       cout << "=========================================================== " << endl ;
    //       cout << "h2-h"      << GenCopyH[i]->GetBinContent(j) << endl;
    //       cout << "(h2-h)/h2" << (UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j) << endl;
    //       Counter++;
    //     }
    // }
    for (Int_t j = 1; j < 31; j++) {
      if ((UnfoldedH[i]->GetBinContent(j) > 0.00000000000001) && GenH[i]->GetBinContent(j) > 0.00000000000001) {
        if( TMath::Abs((UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j)) > 0.5  )  continue;
        DifferencesH->SetBinContent(Counter, (UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j));
        DifferencesH->SetBinError(Counter, 0.000000000000000000);
        Values[Counter] = (UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))*(UnfoldedH[i]->GetBinContent(j)-GenH[i]->GetBinContent(j))/GenH[i]->GetBinContent(j);
        cout << "Values[" << Counter << "] = " << Values[Counter] << endl;
        ChiSquareTest += Values[Counter];
        Counter++;
      }
    }
  }
  TH1F* histogram = (TH1F*) DifferencesH->Clone("histogram");





  new TCanvas;
  // histogram->GetYaxis()->SetRangeUser(0.95,1.05);
  histogram->Draw("*H");


  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);

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

  histogram->GetXaxis()->SetTitle("Bin number");
  histogram->GetYaxis()->SetTitle("(Unfolded REC - GEN)/GEN [a.u.]");
  // histogram->GetYaxis()->SetRangeUser(-0.006,0.004);
  // histogram->GetYaxis()->SetRangeUser(-0.010,0.010);
  histogram->GetXaxis()->SetRangeUser(-0.5, 229.5);
  // histogram->GetXaxis()->SetRangeUser(100, 600);
  histogram->SetLineWidth(5);
  histogram->SetLineColor(2);
  histogram->Draw("");

  TLatex* latex5 = new TLatex();
  latex5->SetTextSize(0.045);
  latex5->SetTextFont(42);
  latex5->SetTextAlign(11);
  latex5->SetNDC();
  latex5->DrawLatex(0.31,0.94,"LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis");
  if (Iterations == 1){latex5->DrawLatex(0.31,0.80,"1 iteration");}
  else {latex5->DrawLatex(0.31,0.80,Form("%d iterations", Iterations));}
  latex5->DrawLatex(0.2,0.64,Form("#chi^{2} = #sum (Unfolded REC - GEN)^{2}/GEN = %0.02f", ChiSquareTest));


  gPad->SaveAs(Form("Unfolding2D/residuals-closure-v3-modulation-%d.pdf", Iterations), "recreate");


}
