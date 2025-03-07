// ePIC DVCS analysis class

// include guard
//#ifndef ePIC_DVCS_TASK_H
#pragma once
#define ePIC_DVCS_TASK_H

// Include detector resolution scripts (from A. Jentsch)
#include "detectorResolution.h"
// Include ePIC event kinematic utilities
#include "ePIC_IncKinUtils.h"
#include "ePIC_ExcKinUtils.h"

// ROOT::Math aliases
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;

class ePIC_DVCS_TASK{
 private:
  // TStrings of stored campaign attributes (date, beam energy and beam setting)
  TString sCamp{};
  TString sEnergy{};
  TString sSett{};

  // Input and output files (file list string for input, single ROOT file for output
  TString sInList{};
  TFile*  fOutFile{};

  // Boolean behaviour selectors
  Bool_t kUsePID{kFALSE}; // Default to using assumptions about DVCS event
  Bool_t kUseExplicitMatch{kFALSE}; // Default to not using explicit matching between associated and reco tracks
  Bool_t kUseEventBeams{kFALSE}; // Default to using beams averaged over file

  // Global variables
  Float_t fMass_proton{0.938272};
  Float_t fMass_electron{0.000511};

  // Cut variables
  Float_t fPMax_p{};     // Maximum scattered proton momentum
  Float_t fPMax_e{};     // Maximum scattered electron momentum
  Float_t fMinQ2{};      // Miniumum Q2
  Float_t fMaxt_RP{};    // Maximum t expected from Roman Pots
  Float_t fMax_Emiss{};  // Maximum missing energy from 3-particle final state
  Float_t fMax_M2miss{}; // Maximum missing mass squared from 3-particle state
  Float_t fxB_Tail{};

  // Kinematic variables used in cut functions
  Float_t fQ2{0.};
  Float_t fxB{1e-10};
  Float_t ft{0};
  Float_t fEmiss{0};
  Float_t fM2miss{0};

  // Objects for undoing afterburn boost
  Float_t fXAngle{-0.025}; // Crossing angle in radians
  Float_t fRotX{};
  RotationX rotAboutX;
  Float_t fRotY{};
  RotationY rotAboutY;
  MomVector vBoostToCoM;
  MomVector vBoostToHoF;

 public:
  // Default constructor
  ePIC_DVCS_TASK();
  // Specific constructor
  ePIC_DVCS_TASK(TString camp, TString energy, TString sett);

  // Set stored campaign attributes
  void setDate(TString string)   { sCamp = string; }
  void setEnergy(TString string) { sEnergy = string; }
  void setSetting(TString string){ sSett = string; }
  // Set input file list and output ROOT file
  void setInFileList(TString name);
  void setOutFile(TString name);

  void setUsePID(Bool_t usePID){ kUsePID = usePID; }
  void setUseExplicitMatch(Bool_t useExplicitMatch){ kUseExplicitMatch = useExplicitMatch; }
  void setUseEventBeams (Bool_t useEventBeams){ kUseEventBeams = useEventBeams; }

  // Set momenta cuts automatically
  void setMomCuts(Float_t factor = 1.);
  // Set other cuts - used in run macro
  void setMin_Q2(Float_t cut)    { fMinQ2 = cut; }
  void setMax_tRP(Float_t cut)   { fMaxt_RP = cut; }
  void setMax_Emiss(Float_t cut) { fMax_Emiss = cut; }
  void setMax_M2miss(Float_t cut){ fMax_M2miss = cut; }

  // Apply cuts
  Bool_t applyCuts_Electron(std::vector<P3EVector> scate);
  Bool_t applyCuts_Photon(std::vector<P3EVector> scatg);
  Bool_t applyCuts_Proton(std::vector<P3EVector> scatp, TString sProtonDet="all");
  Bool_t applyCuts_DVCS(TString sProtonDet="all");
  Bool_t applyCuts_All(P3EVector beame, P3EVector beamp, std::vector<P3EVector> scate, std::vector<P3EVector> scatp, std::vector<P3EVector> scatg, TString sProtonDet="all");

  // Undo afterburner boost
  void undoAfterburnAndCalc(P3EVector& p, P3EVector& k); // Undo procedure AND calculate boost vectors
  void undoAfterburn(P3EVector& p); // Undo procedure ONLY
  
  // Calculation of kinematic quantities
  // Event kinematics: t, Q2, xB, Phi
  Double_t calcTrentoPhi_qp(P3EVector k, P3EVector kprime, P3EVector pprime);
  Double_t calcTrentoPhi_pg(P3EVector k, P3EVector kprime, P3EVector pprime, P3EVector gprime);
  Double_t calcTrentoPhi_qg(P3EVector k, P3EVector kprime, P3EVector gprime);
  Double_t calcTrentoPhi_4Vec(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime);
  Double_t calcgT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j);
  Double_t calcepsT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j);
  Int_t LeviCivita(int i, int j, int k, int l);
  Double_t calcPhiQPQG(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime);
  Double_t calcConeAngle(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime, P3EVector gprime);

  void doAnalysis();
};

