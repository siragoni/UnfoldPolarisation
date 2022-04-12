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
    f[i] = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    h[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
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
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 6; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 6; iCosThetaBins < 18; iCosThetaBins++) {
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
  // gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.,    -2, 2        );
  gMinuit->DefineParameter(0, "LambdaTheta",        1., 0.1,    -2, 2        );
  gMinuit->DefineParameter(1, "LambdaPhi",           0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(2, "LambdaThetaPhi",      0, 0.1,    -2, 2        );
  gMinuit->DefineParameter(3, "Normalisation",   2420000, 100, 1, 10000000    );
  // gMinuit->DefineParameter(3, "Normalisation",   (394082.*3./(4.*TMath::Pi())), 0, 1, 10000000    ); // from projection
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
    // const Double_t min = 10000;
    // const Double_t max = 50000;
    const Double_t min = histo2->GetMaximum()*0.2;
    const Double_t max = histo2->GetMaximum()*1.2;

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

  gPad->SaveAs("Unfolding2D/Inverted2Dmaps.pdf",  "RECREATE");


  new TCanvas;
  histo2->Draw("surf same");
  histo2->Fit("Model", "R");




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
    f[i] = new TFile(Form("Unfolding2D/UnfoldHe_%d.root", i));
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

  gPad->SaveAs("Unfolding2D/Inverted2Dmaps.pdf",  "RECREATE");




}



//______________________________________________
void DrawRaw(){


  TFile* f[24];
  TH1F*  h[24];
  for (Int_t i = 4; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldHeV2_%d.root", i));
    h[i] = (TH1F*) f[i]->Get(Form("histo_%d", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo2_%d", i));
    // h[i] = (TH1F*) f[i]->Get(Form("histo4_%d", i));
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

  Double_t M = 1.;
  for (size_t iCosThetaBins = 5; iCosThetaBins < 19; iCosThetaBins++) {
  // for (size_t iCosThetaBins = 6; iCosThetaBins < 18; iCosThetaBins++) {
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
      if( iCosThetaBins == 0  || iCosThetaBins == 1  || iCosThetaBins == 2  || iCosThetaBins == 3  ||
          iCosThetaBins == 4  || iCosThetaBins == 23 || iCosThetaBins == 22 || iCosThetaBins == 21 ||
          iCosThetaBins == 20 || iCosThetaBins == 19 )
      {
        M = 1;
      } else if ( iCosThetaBins == 5  || iCosThetaBins == 6  || iCosThetaBins == 18  || iCosThetaBins == 17 ) {
        M = 6;
      } else if ( iCosThetaBins == 7  || iCosThetaBins == 8  || iCosThetaBins == 16  || iCosThetaBins == 15 ) {
        M = 12;
      } else if ( iCosThetaBins == 9  || iCosThetaBins == 10  || iCosThetaBins == 11  ||
                  iCosThetaBins == 12 || iCosThetaBins == 13  || iCosThetaBins == 14 ) {
        M = 24;
      }


      // histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx) );
      // histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx) );
      histo2->SetBinContent( binx2, biny2, h[iCosThetaBins]->GetBinContent(binx)*((Double_t) M));
      histo2->SetBinError( binx2, biny2, h[iCosThetaBins]->GetBinError(binx)*((Double_t) M));
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
    // const Double_t min = 8000;
    const Double_t min = histo2->GetMinimum()*0.9;
    const Double_t max = histo2->GetMaximum()*1.1;
    // const Double_t min = 10000;
    // const Double_t max = 50000;
    const Int_t nLevels = 999;
    Double_t levels[nLevels];


    for(int i = 1; i < nLevels; i++) {
      levels[i] = min + (max - min) / (nLevels - 1) * (i);
    }
    levels[0] = 0.01;

    histo2->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    histo2->DrawClone("col");// draw "axes", "contents", "statistics box"
    histo2->GetZaxis()->SetRangeUser(min, max); // ... set the range ...
    // histo2->GetZaxis()->SetRangeUser(histo2->GetMinimum()*0.9, histo2->GetMaximum()*1.1); // ... set the range ...
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

  gPad->SaveAs("Unfolding2D/HeRaw2D-colz.pdf",  "RECREATE");


  new TCanvas;
  histo2->GetXaxis()->SetTitleOffset(1.3);
  histo2->GetYaxis()->SetTitleOffset(1.3);
  histo2->Draw("surf same");
  TLatex* latex2 = new TLatex();
  latex2->SetTextSize(0.055);
  latex2->SetTextFont(42);
  latex2->SetTextAlign(11);
  latex2->SetNDC();
  latex2->DrawLatex(0.12,0.94,"ALICE, Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latex2->SetTextSize(0.055);
  latex2->DrawLatex(0.7,0.94,"Helicity, this thesis");
  gPad->SaveAs("Unfolding2D/HeRaw2D-surf.pdf",  "RECREATE");


}
//_____________________________________________________________________________
void DrawRefolded()
{
  new TCanvas;
  gPad->SetMargin(0.13,0.10,0.12,0.10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetOptStat(0);



  TFile* file = new TFile("Unfolding2D/UnfoldHeV2_10.root");
  TH1F* UnlikeSignDimuon2    = (TH1F*)file->Get("histo5");
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
  UnlikeSignDimuon->GetYaxis()->SetTitle("Ratio of refolded data to REC data");
  UnlikeSignDimuon->GetYaxis()->SetRangeUser(0.5,2.4);
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
  latex5->DrawLatex(0.31,0.94,"ALICE LHC18l7, PbPb #sqrt{s_{NN}} = 5.02 TeV");
  latex5->DrawLatex(0.31,0.74,"This thesis, 12 iterations");

  gPad->SaveAs("Unfolding2D/residual-refolded-data-he-12.pdf", "recreate");
}
