//------------------------------------------------------------------------------------
//
// Utility functions; calculation of kinematic quantities for ePIC analyses
//
// Author: O. Jevons, 27/02/25
//
//------------------------------------------------------------------------------------

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
Double_t calcT_eX(TLorentzVector e, TLorentzVector ep, TLorentzVector X);
Double_t calcT_eX(TLorentzVector q, TLorentzVector X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);

// Calculate DIS kinematics - electron method
// e + p -> e' + X
Double_t calcQ2_Elec(TLorentzVector e,TLorentzVector  ep);
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep);
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep);
Double_t calcX_Elec(TLorentzVector e,TLorentzVector  ep,TLorentzVector  p);
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p);
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p);
Double_t calcY_Elec(TLorentzVector e, TLorentzVector ep, TLorentzVector p);
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p);
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p);
void calcKin_Elec(TLorentzVector e, TLorentzVector p, TLorentzVector ep, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  Float_t &Q2, Float_t &x, Float_t &y);

// Calculate DIS kinematics - JB method
// SIDIS/Exclusive case:  e + p -> e' + p' + X
// Separating the scattered proton from the rest of the hadronic final state
Double_t calcQ2_JB(TLorentzVector e, TLorentzVector pp, TLorentzVector X);
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcX_JB(TLorentzVector e, TLorentzVector p, TLorentzVector pp, TLorentzVector X);
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcY_JB(TLorentzVector e, TLorentzVector pp, TLorentzVector X);
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
void calcKin_JB(TLorentzVector e, TLorentzVector p, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y);
// DIS case: e + p -> e' + HFS
// Combining the whole HFS into 1 vector
Double_t calcQ2_JB(TLorentzVector e, TLorentzVector HFS);
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcX_JB(TLorentzVector e, TLorentzVector p, TLorentzVector HFS);
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcY_JB(TLorentzVector e, TLorentzVector HFS);
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e,  
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e,  
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
void calcKin_JB(TLorentzVector e, TLorentzVector p, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y);


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


// Missing mass/energy/momentum: ab->cdf
Double_t calcPMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcPtMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcEMiss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
Double_t calcM2Miss_3Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> f);
// ab->cd
Double_t calcPMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcEMiss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
Double_t calcM2Miss_2Body(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> a, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> b, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> c, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> d);
*/

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate energy from momentum and mass
// 1. Using TVector3 (legacy) for the momentum vector
Double_t calcE(TVector3 mom, Float_t M){ 
  return TMath::Sqrt(mom.Mag2() + TMath::Power(M,2)); 
}
// 2. Using ROOT::Math::XYZVector (updated alternative to TVector3) for the momentum vector
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
// Using TLorentzVectors (legacy) for e, e', X
Double_t calcT_eX(TLorentzVector e, TLorentzVector ep, TLorentzVector X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// Using LorentzVectors (legacy) for q, x
Double_t calcT_eX(TLorentzVector q, TLorentzVector X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
// Using ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D> for e, e', X
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// Using ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D> for q, X
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}
// Using ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D> for e, e', X
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  double t = (e - ep - X).M2();
  
  return TMath::Abs(t);
}
// Using ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D> for q, X
Double_t calcT_eX(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> q, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  double t = (q - X).M2();
  
  return TMath::Abs(t);
}

// Calculate inclusive kinematics (Q2, x, y) with the electron method
// Q2
Double_t calcQ2_Elec(TLorentzVector e, TLorentzVector ep){
  double q2 = (e - ep).M2();
  double Q2 = -q2;

  return Q2;
}
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep){
  double q2 = (e - ep).M2();
  double Q2 = -q2;

  return Q2;
}
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep){
  double q2 = (e - ep).M2();
  double Q2 = -q2;
  
  return Q2;
}

// x
Double_t calcX_Elec(TLorentzVector e,TLorentzVector  ep,TLorentzVector  p){
  TLorentzVector q = e-ep;
  double q2 = q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;
  return xB;
}
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p){
  P3EVector q = e-ep;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;
  return xB;
}
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p){
  P3MVector q = e-ep;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double x = q2/denom;
  return x;
}
Double_t calcY_Elec(TLorentzVector e, TLorentzVector ep, TLorentzVector p){
  TLorentzVector q = e - ep;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  double y = num/den;
  return y;
}
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p){
  P3EVector q = e - ep;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  double y = num/den;
  return y;
}
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p){
  P3MVector q = e - ep;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  double y = num/den;
  return y;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables
void calcKin_Elec(TLorentzVector e, TLorentzVector p, TLorentzVector ep, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;

  // Calculate kinematics
  TLorentzVector q = e-ep;
  Float_t Q2_e = -q.Dot(q);
  
  Float_t q_dot_p = q.Dot(p);
  Float_t x_e = Q2_e/(2*q_dot_p);

  Float_t e_dot_p = e.Dot(p);
  Float_t y_e = q_dot_p/e_dot_p;

  // Export variables
  Q2 = Q2_e;
  x = x_e;
  y = y_e;
}
void calcKin_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Calculate kinematics
  P3EVector q = e-ep;
  Float_t Q2_e = -q.Dot(q);
  
  Float_t q_dot_p = q.Dot(p);
  Float_t x_e = Q2_e/(2*q_dot_p);
  
  Float_t e_dot_p = e.Dot(p);
  Float_t y_e = q_dot_p/e_dot_p;
  
  // Export variables
  Q2 = Q2_e;
  x = x_e;
  y = y_e;
}
void calcKin_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Calculate kinematics
  P3MVector q = e-ep;
  Float_t Q2_e = -q.Dot(q);
  
  Float_t q_dot_p = q.Dot(p);
  Float_t x_e = Q2_e/(2*q_dot_p);

  Float_t e_dot_p = e.Dot(p);
  Float_t y_e = q_dot_p/e_dot_p;
  
  // Export variables
  Q2 = Q2_e;
  x = x_e;
  y = y_e;
}

// Calculate inclusive kinematics (Q2, x, y) with the electron method
// SIDIS case ep -> e'p'X
// Q2
Double_t calcQ2_JB(TLorentzVector e, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
// x
Double_t calcX_JB(TLorentzVector e, TLorentzVector p, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
// y
Double_t calcY_JB(TLorentzVector e, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_JB(TLorentzVector e, TLorentzVector p, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
// DIS case ep -> e'X
// Q2
Double_t calcQ2_JB(TLorentzVector e, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
Double_t calcQ2_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
// x
Double_t calcX_JB(TLorentzVector e, TLorentzVector p, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
Double_t calcX_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
// y
Double_t calcY_JB(TLorentzVector e, TLorentzVector HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e,  
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e,  
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_JB(TLorentzVector e, TLorentzVector p, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
void calcKin_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}
