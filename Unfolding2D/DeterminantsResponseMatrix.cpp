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
/* - Fit function for the helicity case. It is basically a parabolic fit...
   -
 */
void Determinants(){


  TFile* f[24];
  RooUnfoldResponse*  h[24];
  Double_t Determinants[24];
  // TMatrixD ValuesResponseMatrix[24];
  TMatrixD ValuesResponseMatrix(6,6);
  // TMatrixD ValuesResponseMatrix9(24,24);
  TMatrixD ValuesResponseMatrix9(6,6);
  RooUnfoldBayes    unfold[24];
  TMatrixD CovarianceUncert[24];
  TMatrixD CovarianceUncert2(6,6);
  TVectorD *VectorUncert[24];
  TMatrixD CovarianceUncertT[24];
  TFile* fileDataRaw[24];
  TH1F *Data[24];
  for (Int_t iC = 4; iC < 20; iC++) {
    fileDataRaw[iC] = new TFile(Form("SignalExtraction/MonteCarloYieldsHe_%d.root", iC));
    // fileDataRaw[iC] = new TFile(Form("SignalExtraction/RawYieldsHe_%d.root", iC));
    Data[iC]        = (TH1F*)fileDataRaw[iC]->Get(Form("h_%d", iC));
  }
  TH1F* histo[24];
  Int_t N = 1;
  for (Int_t i = 4; i < 20; i++) {
    if( i == 0  || i == 1  || i == 2  || i == 3  ||
        i == 4  || i == 23 || i == 22 || i == 21 ||
        i == 20 || i == 19 )
    {
      N = 1;
    } else if ( i == 5  || i == 6  || i == 18  || i == 17 ) {
      N = 6;
    } else if ( i == 7  || i == 8  || i == 16  || i == 15 ) {
      N = 12;
    } else if ( i == 9  || i == 10  || i == 11  ||
                i == 12 || i == 13  || i == 14 ) {
      N = 24;
    }

    histo[i] = new TH1F(Form("histo_%d", i), Form("histo_%d", i), N, 0., 2.*TMath::Pi());
    for(Int_t j = 1; j < 30; j++){
      histo[i]->SetBinContent(j, Data[i]->GetBinContent(j));
      histo[i]->SetBinError(  j, Data[i]->GetBinError(j)  );
    }
  }


  TArrayD Elements[36];
  for (Int_t i = 5; i < 20; i++) {
    f[i] = new TFile(Form("Unfolding2D/UnfoldedClosureHe_%d.root", i));
    h[i] = (RooUnfoldResponse*) f[i]->Get("response;1");
    // (h[i]->Mresponse()).InvertFast(&Determinants[i]);
    Determinants[i] = (h[i]->Mresponse()).Determinant();
    // ValuesResponseMatrix[i] = TMatrixD((h[i]->Mresponse()));
    if(i == 5)ValuesResponseMatrix = ((h[i]->Mresponse()));
    if(i == 6)ValuesResponseMatrix9 = ((h[i]->Mresponse()));
    // cout << "Determinants[" << i << "] = " << Determinants[i] << endl;
    unfold[i]  = RooUnfoldBayes(h[i], histo[i], 1);
    cout << "CovarianceUncert[" << i << "] "  << endl;
    // CovarianceUncert[i] = unfold[i].Ereco(RooUnfold::kCovariance);
    if(i == 5)CovarianceUncert2 = unfold[i].Ereco(RooUnfold::kCovariance);
    // (unfold[i].Ereco(RooUnfold::kCovariance)).Print();
    // if(i == 5){
    //   for (Int_t iC = 0; iC < 36; iC++) {
    //     Double_t row = -1;
    //     Double_t columns = 0.;
    //     columns = modf(iC, &row);
    //     Int_t columnsI = (Int_t)  columns;
    //     Elements[iC] = (unfold[i].Ereco(RooUnfold::kCovariance))(row, columnsI);
    //   }
    // }
    (unfold[i].ErecoV(RooUnfold::kCovariance)).Print();
    // CovarianceUncert[i].Transpose(CovarianceUncertT[i]);


    // h[i]->Print();
    // ValuesResponseMatrix[i] = (TMatrixD) h[i]->Mresponse();
    // ValuesResponseMatrix[i].InvertFast(&Determinants[i]);
    // cout << "Determinants[" << i << "] = " << Determinants[i] << endl;

  }

  // for (Int_t i = 5; i < 20; i++) {
  //   VectorUncert[i].Print();
  //
  // }

  TFile* f3 = new TFile("saving.root", "RECREATE");
  h[12]->Write();
  f3->Close();
  // for (size_t i = 0; i < 36; i++) {
  // cout << Elements[i];
  // }


  cout << "================" << endl;
  cout << "For Roman " << endl;
  (h[5]->Mresponse()).Print();

  cout << "================" << endl;
  cout << "For Roman 2" << endl;
  ValuesResponseMatrix.Print();





  cout << "================" << endl;
  cout << "For Roman eigenvalue" << endl;

  TMatrixDEigen eigen(ValuesResponseMatrix);
  auto eigenvalues = eigen.GetEigenValues();
  eigenvalues.Print();
  auto eigenvectors = eigen.GetEigenVectors();



  TMatrixDEigen eigen2(ValuesResponseMatrix9);
  auto eigenvalues2 = eigen2.GetEigenValues();
  eigenvalues2.Print();
  auto eigenvectors2 = eigen2.GetEigenVectors();


  cout << "================" << endl;
  cout << "For Roman covariance" << endl;
  CovarianceUncert2.Print();
}
