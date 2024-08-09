// ePIC DVCS analysis class

// include guard
//#ifndef ePIC_DVCS_TASK_H
#pragma once
#define ePIC_DVCS_TASK_H

// Include detector resolution scripts (from A. Jentsch)
#include "detectorResolution.h"

// ROOT::Math aliases
using ROOT::Math::VectorUtil::boost;
using ROOT::Math::VectorUtil::Angle;
using ROOT::Math::RotationX;
using ROOT::Math::RotationY;
using P3MVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzMVector>;
using P3EVector=ROOT::Math::LorentzVector<ROOT::Math::PxPyPzEVector>;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

class ePIC_DVCS_TASK{
 private:
  // TStrings of stored campaign attributes (date, beam energy and beam setting)
  TString sCamp{};
  TString sEnergy{};
  TString sSett{};

  // Input and output files (file list string for input, single ROOT file for output
  TString sInList{};
  TFile*  fOutFile{};

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

  // Set momenta cuts automatically
  void setMomCuts(Float_t factor = 1.);
  // Set other cuts - used in run macro
  void setMin_Q2(Float_t cut)    { fMinQ2 = cut; }
  void setMax_tRP(Float_t cut)   { fMaxt_RP = cut; }
  void setMax_Emiss(Float_t cut) { fMax_Emiss = cut; }
  void setMax_M2miss(Float_t cut){ fMax_M2miss = cut; }

  // Apply cuts
  Bool_t applyCuts_Electron(std::vector<P3MVector> scate);
  Bool_t applyCuts_Photon(std::vector<P3MVector> scatg);
  Bool_t applyCuts_Proton(std::vector<P3MVector> scatp, TString sProtonDet="all");
  Bool_t applyCuts_DVCS(TString sProtonDet="all");
  Bool_t applyCuts_All(P3MVector beame, P3MVector beamp, std::vector<P3MVector> scate, std::vector<P3MVector> scatp, std::vector<P3MVector> scatg, TString sProtonDet="all");

  // Undo afterburner boost
  void undoAfterburnAndCalc(P3MVector& p, P3MVector& k); // Undo procedure AND calculate boost vectors
  void undoAfterburn(P3MVector& p); // Undo procedure ONLY

  // Calculation of kinematic quantities
  // Event kinematics: t, Q2, xB, Phi
  Double_t calcT(P3MVector p, P3MVector pprime);
  Double_t calcTNoP(P3MVector k, P3MVector kprime, P3MVector g);
  Double_t calcQ2(P3MVector k, P3MVector kprime);
  Double_t calcBjorkenX(P3MVector k, P3MVector kprime, P3MVector p);
  Double_t calcBjorkenY(P3MVector k, P3MVector kprime, P3MVector p);
  Double_t calcTrentoPhi(P3MVector k, P3MVector kprime, P3MVector pprime);
  Double_t calcPhiH(P3MVector k, P3MVector kprime, P3MVector p, P3MVector pprime);
  Double_t calcPhiQPQG(P3MVector k, P3MVector p, P3MVector kprime, P3MVector gprime);
  Double_t calcConeAngle(P3MVector k, P3MVector p, P3MVector kprime, P3MVector pprime, P3MVector gprime);
  // Missing mass/energy/momentum: ab->cdf
  Double_t calcPMiss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);
  Double_t calcPtMiss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);
  Double_t calcEMiss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);
  Double_t calcM2Miss_3Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d, P3MVector f);
  // ab->cd
  Double_t calcPMiss_2Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d);
  Double_t calcEMiss_2Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d);
  Double_t calcM2Miss_2Body(P3MVector a, P3MVector b, P3MVector c, P3MVector d);
  

  void doAnalysis();
};

