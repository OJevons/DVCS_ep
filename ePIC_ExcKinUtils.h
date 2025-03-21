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
Double_t calcE(const TVector3& mom, const Float_t& M);
Double_t calcE(const Float_t& px, const Float_t& py, const Float_t& pz, const Float_t& M);
Double_t calcE(const ROOT::Math::XYZVector& mom, const Float_t& M);

// Calculate Mandelstam t - BABE method by tRECO convention
Double_t calcT_BABE(const TLorentzVector& be, const TLorentzVector& ba);
Double_t calcT_BABE(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& be, 
		    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ba);
Double_t calcT_BABE(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& be, 
		    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& ba);

// Calculate Mandelstam t - eX method by tRECO convention
// Using separate vectors for beam and scattered electron
Double_t calcT_eX(const TLorentzVector& e, const TLorentzVector& ep, const TLorentzVector& X);
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& e, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& ep, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& X);
// Giving virtual photon vector directly 
Double_t calcT_eX(const TLorentzVector& q, const TLorentzVector& X);
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& q, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& q, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& X);

/*
// Calculate Phi_h
Double_t calcTrentoPhi_qp(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& pp);
Double_t calcTrentoPhi_pg(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& pp, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
Double_t calcTrentoPhi_qg(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
Double_t calcTrentoPhi_4Vec(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& p, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& pp);
Double_t calcgT_ij(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& q, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& p, Int_t i, Int_t j);
Double_t calcepsT_ij(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& q, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& p, Int_t i, Int_t j);
Int_t LeviCivita(int i, int j, int k, int l);

// Other event angles
Double_t calcPhiQPQG(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& p, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
Double_t calcConeAngle(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& p, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& pp, const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X);
*/

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
Double_t calcPMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f);
Double_t calcPMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f);
Double_t calcPMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f);
// Missing transverse momentum
Double_t calcPtMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f);
Double_t calcPtMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f);
Double_t calcPtMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f);
// Missing energy
Double_t calcEMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f);
Double_t calcEMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f);
Double_t calcEMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f);
// Missing mass (squared)
Double_t calcM2Miss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f);
Double_t calcM2Miss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f);
Double_t calcM2Miss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f);
// 2-body final state: ab->cd
// Missing momentum
Double_t calcPMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d);
Double_t calcPMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d);
Double_t calcPMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d);
// Missing transverse momentum
Double_t calcPtMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d);
Double_t calcPtMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d);
Double_t calcPtMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d);
// Missing energy
Double_t calcEMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d);
Double_t calcEMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d);
Double_t calcEMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d);
// Missing mass (squared)
Double_t calcM2Miss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d);
Double_t calcM2Miss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d);
Double_t calcM2Miss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d);


//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy from momentum and mass
// 1. Using TV3ctor3 (legacy) for momentum vector
Double_t calcE(const TVector3& mom, const Float_t& M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 2. Using ROOT::Math::XYZVector (updated alternative to TVector3) for momentum vector
Double_t calcE(const ROOT::Math::XYZVector& mom, const Float_t& M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 3. Using separate floats for momentum components
Double_t calcE(const Float_t& px, const Float_t& py, const Float_t& pz, const Float_t& M){ 
  return TMath::Sqrt(TMath::Power(px,2) + TMath::Power(py,2) + TMath::Power(pz,2) + TMath::Power(M,2)); 
}


// Calculate Mandelstam t - BABE method using tRECO conventions
// Uses incoming proton BEam and scattered BAryon 4-vectors
// Another way of saying t = -(p' - p)^2
Double_t calcT_BABE(const TLorentzVector& be, const TLorentzVector& ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_BABE(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& be,
		    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_BABE(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& be,
		    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& ba){
  double t = (ba - be).M2();
  
  return TMath::Abs(t);
}

// Calculate Mandelstam t - eX method using tRECO conventions
// Uses difference between the beam and scattered electron and all of the (non-scattered) final state
// e.g. for DVCS, X is a photon; for electroproduction, it is the species of interest; etc...
// 1. Providing separate vectors for beam and scattered electrons
Double_t calcT_eX(const TLorentzVector& e, const TLorentzVector& ep, const TLorentzVector& X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& e, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& ep, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& e, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& ep, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// 2. Giving virtual photon vector directly
Double_t calcT_eX(const TLorentzVector& q, const TLorentzVector& X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& q, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
Double_t calcT_eX(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& q, 
		  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}

// Calculate missing kinematics (mass/energy/momentum)
// 3-body final state: ab->cdf
// Missing momentum
Double_t calcPMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f){ 
  return (a+b-c-d-f).P(); 
}
Double_t calcPMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f){
  return (a+b-c-d-f).P(); 
}
Double_t calcPMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f){
  return (a+b-c-d-f).P(); 
}
// Missing transverse momentum
Double_t calcPtMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f){
  return (a+b-c-d-f).Perp(); 
}
Double_t calcPtMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f){
  return (a+b-c-d-f).Pt(); 
}
Double_t calcPtMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f){
  return (a+b-c-d-f).Pt(); 
}
// Missing energy
Double_t calcEMiss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f){
  return (a+b-c-d-f).E(); 
}
Double_t calcEMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f){
  return (a+b-c-d-f).E(); 
}
Double_t calcEMiss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f){
  return (a+b-c-d-f).E(); 
}
// Missing mass (squared)
Double_t calcM2Miss_3Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d, const TLorentzVector& f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_3Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
// 2-body final state: ab->cd
// Missing momentum
Double_t calcPMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d){
  return (a+b-c-d).P();
}
Double_t calcPMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d){
  return (a+b-c-d).P();
}
Double_t calcPMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d){
  return (a+b-c-d).P();
}
// Missing transverse momentum
Double_t calcPtMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d){
  return (a+b-c-d).Perp();
}
Double_t calcPtMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d){
  return (a+b-c-d).Pt();
}
Double_t calcPtMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d){
  return (a+b-c-d).Pt();
}
// Missing energy
Double_t calcEMiss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d){
  return (a+b-c-d).E();
}
Double_t calcEMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d){
  return (a+b-c-d).E();
}
Double_t calcEMiss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			 const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d){
  return (a+b-c-d).E();
}
// Missing mass (squared)
Double_t calcM2Miss_2Body(const TLorentzVector& a, const TLorentzVector& b, const TLorentzVector& c, const TLorentzVector& d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>& d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
Double_t calcM2Miss_2Body(const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& a, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& b, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& c, 
			  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>& d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}
