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
// FUNCTION DECLARATIONS
//
// NOTE: 4-vector functions defined in triplicate
//       1) Using TLorentzVector objects (legacy)
//       2) Using ROOT::Math::LorentzVector with PxPyPzE4D coordinate system (updated alternative to TLorentzVector, using 3 momentum components and energy)
//       3) Using ROOT::Math::LorentzVector with PxPyPzM4D coordinate system (similar to (2), but using mass instead of energy)
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy
Double_t calcE(TVector3 mom, Float_t M);
Double_t calcE(Float_t px, Float_t py, Float_t pz, Float_t M);
Double_t calcE(ROOT::Math::XYZVector mom, Float_t M);

// Calculate Mandelstam t - BABE method by tRECO convention
Double_t calcT_BABE(TLorentzVector be, TLorentzVector ba);
Double_t calcT_BABE(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> be, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ba);
Double_t calcT_BABE(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> be, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ba);

// Calculate Mandelstam t - eX method by tRECO convention
// Using separate vectors for beam and scattered electron
Double_t calcT_eX(TLorentzVector e, TLorentzVector ep, TLorentzVector X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
// Giving virtual photon vector directly 
Double_t calcT_eX(TLorentzVector q, TLorentzVector X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);

/*
// Calculate Phi_h
Double_t calcTrentoPhi_qp(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp);
Double_t calcTrentoPhi_pg(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcTrentoPhi_qg(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcTrentoPhi_4Vec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp);
Double_t calcgT_ij(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, Int_t i, Int_t j);
Double_t calcepsT_ij(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, Int_t i, Int_t j);
Int_t LeviCivita(int i, int j, int k, int l);

// Other event angles
Double_t calcPhiQPQG(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcConeAngle(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
*/

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
Double_t calcPMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f);
Double_t calcPMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcPMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f);
// Missing transverse momentum
Double_t calcPtMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f);
Double_t calcPtMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcPtMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f);
// Missing energy
Double_t calcEMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f);
Double_t calcEMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcEMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f);
// Missing mass (squared)
Double_t calcM2Miss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f);
Double_t calcM2Miss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcM2Miss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f);
// 2-body final state: ab->cd
// Missing momentum
Double_t calcPMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d);
Double_t calcPMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcPMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d);
// Missing transverse momentum
Double_t calcPtMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d);
Double_t calcPtMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcPtMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d);
// Missing energy
Double_t calcEMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d);
Double_t calcEMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcEMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d);
// Missing mass (squared)
Double_t calcM2Miss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d);
Double_t calcM2Miss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcM2Miss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d);


//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy from momentum and mass
// 1. Using TV3ctor3 (legacy) for momentum vector
Double_t calcE(TVector3 mom, Float_t M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 2. Using ROOT::Math::XYZVector (updated alternative to TVector3) for momentum vector
Double_t calcE(ROOT::Math::XYZVector mom, Float_t M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 3. Using separate floats for momentum components
Double_t calcE(Float_t px, Float_t py, Float_t pz, Float_t M){ 
  return TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(M,2)); 
}


// Calculate Mandelstam t - BABE method using tRECO conventions
// Uses incoming proton BEam and scattered BAryon 4-vectors
// Another way of saying t = -(p' - p)^2
Double_t calcT_BABE(TLorentzVector be, TLorentzVector ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_BABE(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> be,
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_BABE(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> be,
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eX method using tRECO conventions
// Uses difference between the beam and scattered electron and all of the (non-scattered) final state
// e.g. for DVCS, X is a photon; for electroproduction, it is the species of interest; etc...
// 1. Providing separate vectors for beam and scattered electrons
Double_t calcT_eX(TLorentzVector e, TLorentzVector ep, TLorentzVector X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
Double_t calcT_eX(TLorentzVector q, TLorentzVector X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
Double_t calcPMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f){ 
  return (a+b-c-d-f).P(); 
}
Double_t calcPMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f){
  return (a+b-c-d-f).P(); 
}
Double_t calcPMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f){
  return (a+b-c-d-f).P(); 
}
// Missing transverse momentum
Double_t calcPtMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f){
  return (a+b-c-d-f).Perp(); 
}
Double_t calcPtMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f){
  return (a+b-c-d-f).Pt(); 
}
Double_t calcPtMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f){
  return (a+b-c-d-f).Pt(); 
}
// Missing energy
Double_t calcEMiss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f){
  return (a+b-c-d-f).E(); 
}
Double_t calcEMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f){
  return (a+b-c-d-f).E(); 
}
Double_t calcEMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f){
  return (a+b-c-d-f).E(); 
}
// Missing mass (squared)
Double_t calcM2Miss_3Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d, TLorentzVector f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
// 2-body final state: ab->cd
// Missing momentum
Double_t calcPMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d){
  return (a+b-c-d).P();
}
Double_t calcPMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d){
  return (a+b-c-d).P();
}
Double_t calcPMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d){
  return (a+b-c-d).P();
}
// Missing transverse momentum
Double_t calcPtMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d){
  return (a+b-c-d).Perp();
}
Double_t calcPtMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d){
  return (a+b-c-d).Pt();
}
Double_t calcPtMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d){
  return (a+b-c-d).Pt();
}
// Missing energy
Double_t calcEMiss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d){
  return (a+b-c-d).E();
}
Double_t calcEMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d){
  return (a+b-c-d).E();
}
Double_t calcEMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			 ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d){
  return (a+b-c-d).E();
}
// Missing mass (squared)
Double_t calcM2Miss_2Body(TLorentzVector a, TLorentzVector b, TLorentzVector c, TLorentzVector d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> a, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> b, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> c, 
			  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
