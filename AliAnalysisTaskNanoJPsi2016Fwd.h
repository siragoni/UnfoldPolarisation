/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskNanoJPsi2016Fwd_H
#define AliAnalysisTaskNanoJPsi2016Fwd_H

#include "AliAnalysisTaskSE.h"

class AliMuonTrackCuts; 	// Include class for standard muon tack cuts
class AliAODTrack;


class AliAnalysisTaskNanoJPsi2016Fwd : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskNanoJPsi2016Fwd();
                                AliAnalysisTaskNanoJPsi2016Fwd(const char *name);
        virtual                 ~AliAnalysisTaskNanoJPsi2016Fwd();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);
        virtual void   			NotifyRun();// Implement the Notify run to search for the new parameters at each new runs



        Double_t                CosThetaCollinsSoper( TLorentzVector muonPositive,
                                              TLorentzVector muonNegative,
                                              TLorentzVector possibleJPsi );
Double_t                CosThetaHelicityFrame( TLorentzVector muonPositive,
                                               TLorentzVector muonNegative,
                                               TLorentzVector possibleJPsi );
Double_t                CosPhiCollinsSoper( TLorentzVector muonPositive,
                                            TLorentzVector muonNegative,
                                            TLorentzVector possibleJPsi );
Double_t                CosPhiHelicityFrame( TLorentzVector muonPositive,
                                             TLorentzVector muonNegative,
                                             TLorentzVector possibleJPsi );

	void TrkTrkKine(AliAODTrack *Track1, AliAODTrack *Track2, Double_t mass);
	void CheckTrigger(Bool_t *trig);
	void SetPeriod(Int_t period);
	void SetMC(Bool_t flag);
	Bool_t GoodMUONTrack(Int_t iTrack);
	Int_t GetMassHypothesis(Int_t *idxPosTrk, Int_t *idxNegTrk);
	void GetMCInfo();
        AliMuonTrackCuts* 		fMuonTrackCuts; 					// Use the class as a data member

    private:
	Double_t gMuonMass;
	Int_t fPeriod;
	Bool_t fIsMC;

        AliAODEvent*            fAOD;       //! input event
        TList*                  fOutputList; //! output list
        TH1F*                   fCounterH; //! counter for events passing each cut
        TH1F*                   fMuonTrackCounterH; //! counter for tracks passing each cut
	TH1F*                   fTriggerCounterFwdH; //! counter for triggersper run
	TH2F*                   fGoodTrkFwdH; //! number of good tracks found in event

	TTree *fAnaTree; //! analysis tree
	TTree *fAnaTreeMC; //! analysis tree for MC
	Int_t fRunNum;
	UInt_t fL0inputs;
	Int_t fTracklets;
	Int_t fAnaType;
	Int_t fGoodPosTrk;
	Int_t fGoodNegTrk;
	Bool_t fZNAfired;
	Bool_t fZNCfired;
	Double_t fZNCEnergy;
	Double_t fZNAEnergy;
	Double_t fZPCEnergy;
	Double_t fZPAEnergy;
	Double_t fZNATDC[4];
	Double_t fZNCTDC[4];
	Double_t fZPATDC[4];
	Double_t fZPCTDC[4];
	Int_t fV0ADecision;
	Int_t fV0CDecision;
	Double_t fV0ATime;
	Double_t fV0CTime;
	Float_t fV0AMultiplicity[32];
	Float_t fV0CMultiplicity[32];
	Int_t fADADecision;
	Int_t fADCDecision;
	Double_t fADATime;
	Double_t fADCTime;
	Float_t fADAMultiplicity[8];
	Float_t fADCMultiplicity[8];
  	TBits fIR1Map;
  	TBits fIR2Map;
	Double_t fTrkTrkPt;
	Double_t fTrkTrkPhi;
	Double_t fTrkTrkY;
	Double_t fTrkTrkM;
	Double_t fTrkPt1;
	Double_t fTrkPt2;
	Double_t fTrkEta1;
	Double_t fTrkEta2;
	Double_t fTrkPhi1;
	Double_t fTrkPhi2;
	Double_t fTrkQ1;
	Double_t fTrkQ2;
	Double_t fTrkRabs1;
	Double_t fTrkRabs2;
	Double_t fMCTrkTrkPt;
	Double_t fMCTrkTrkPhi;
	Double_t fMCTrkTrkY;
	Double_t fMCTrkTrkM;
	Double_t fMCTrkPt1;
	Double_t fMCTrkPt2;
	Double_t fMCTrkEta1;
	Double_t fMCTrkEta2;
	Double_t fMCTrkPhi1;
	Double_t fMCTrkPhi2;

  Double_t fMCCosThetaHE;
  Double_t fMCCosThetaCS;
  Double_t fMCPhiHE;
  Double_t fMCPhiCS;
  Double_t fMCTildePhiHEpos;
  Double_t fMCTildePhiHEneg;
  Double_t fMCTildePhiCSpos;
  Double_t fMCTildePhiCSneg;
  Double_t fCosThetaHE;
  Double_t fCosThetaCS;
  Double_t fPhiHE;
  Double_t fPhiCS;
  Double_t fTildePhiHEpos;
  Double_t fTildePhiHEneg;
  Double_t fTildePhiCSpos;
  Double_t fTildePhiCSneg;

	Int_t fMCTrkQ1;
	Int_t fMCTrkQ2;
  Bool_t                  fV0Hits[64];
  Int_t                   fV0TotalNCells;


        AliAnalysisTaskNanoJPsi2016Fwd(const AliAnalysisTaskNanoJPsi2016Fwd&); // not implemented
        AliAnalysisTaskNanoJPsi2016Fwd& operator=(const AliAnalysisTaskNanoJPsi2016Fwd&); // not implemented

        ClassDef(AliAnalysisTaskNanoJPsi2016Fwd, 1);
};

#endif
