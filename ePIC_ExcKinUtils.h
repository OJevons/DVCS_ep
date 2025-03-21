//---------------------------------------------------------------------------------------
//
// Utility functions; calculation of kinematic quantities for exclusive ePIC analyses
//
// Author: O. Jevons, 27/02/25
//
//---------------------------------------------------------------------------------------

#include "TMath.h"

// Aliases for common 3/4-vector types
using P3EVector=ROOT::Math::PxPyPzEVector;
using P3MVector=ROOT::Math::PxPyPzMVector;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//
// NOTE: Templating applied for brevity
//       4-vector functions are valid for any type which contains the operators E(), P(), Pt() and M2()
//       This includes TLorentzVector (legacy) and ALL variants of ROOT::Math::LorentzVector class
//
//-----------------------------------------------------------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy from momentum and mass
// 1. Using vector structures for momentum
// Works for ANY structure which contains Mag2() operator
template<typename P>
Double_t calcE(const P& mom, const Float_t& M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 2. Using separate floats for momentum components
Double_t calcE(const Float_t& px, const Float_t& py, const Float_t& pz, const Float_t& M){ 
  return TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(M,2)); 
}


// Calculate Mandelstam t - BABE method using tRECO conventions
// Uses incoming proton BEam and scattered BAryon 4-vectors
// Another way of saying t = -(p' - p)^2
template <typename V>
Double_t calcT_BABE(const V& be, const V& ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eX method using tRECO conventions
// Uses difference between the beam and scattered electron and all of the (non-scattered) final state
// e.g. for DVCS, X is a photon; for electroproduction, it is the species of interest; etc...
// 1. Providing separate vectors for beam and scattered electrons
template <typename V>
Double_t calcT_eX(const V& e, const V& ep, const V& X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
template <typename V>
Double_t calcT_eX(const V& q, const V& X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
template <typename V>
Double_t calcPMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){ 
  return (a+b-c-d-f).P(); 
}
// Missing transverse momentum
template <typename V>
Double_t calcPtMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  return (a+b-c-d-f).Pt(); 
}
// Missing energy
template <typename V>
Double_t calcEMiss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  return (a+b-c-d-f).E(); 
}
// Missing mass (squared)
template <typename V>
Double_t calcM2Miss_3Body(const V& a, const V& b, const V& c, const V& d, const V& f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}

// 2-body final state: ab->cd
// Missing momentum
template <typename V>
Double_t calcPMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).P();
}
// Missing transverse momentum
template <typename V>
Double_t calcPtMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).Pt();
}
// Missing energy
template <typename V>
Double_t calcEMiss_2Body(const V& a, const V& b, const V& c, const V& d){
  return (a+b-c-d).E();
}
// Missing mass (squared)
template <typename V>
Double_t calcM2Miss_2Body(const V& a, const V& b, const V& c, const V& d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
