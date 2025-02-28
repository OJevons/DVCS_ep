//---------------------------------------------------------------------------------------
//
// Utility functions; calculation of kinematic quantities for inclusive ePIC analyses
//
// Q2, x and y for all 5 methods used for the InclusiveKinematics branch of EICrecon
//     - Electron, JB, DA, Sigma and eSigma
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
//
//
// NOTE 2: Order of vectors (whichever are required) is ALWAYS as follows (for process ep -> e'p'X - change scattered baryon as necessary)
//
//         e -> p -> e' -> p' -> X
//
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate DIS kinematics - electron method
// e + p -> e' + X
Double_t calcQ2_Elec(TLorentzVector e,TLorentzVector  ep);
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep);
Double_t calcQ2_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep);
Double_t calcX_Elec(TLorentzVector e, TLorentzVector  p,TLorentzVector  ep);
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep);
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep);
Double_t calcY_Elec(TLorentzVector e, TLorentzVector p, TLorentzVector ep);
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep);
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep);
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

// Calculate DIS kinematics - DA method
// SIDIS/Exclusive case:  e + p -> e' + p' + X
// Separating the scattered proton from the rest of the hadronic final state
Double_t calcQ2_DA(TLorentzVector e, TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcX_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcY_DA(TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
void calcKin_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y);
// DIS case:  e + p -> e' + HFS
// Combining whole HFS into one vector
Double_t calcQ2_DA(TLorentzVector e, TLorentzVector ep, TLorentzVector HFS);
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep,  
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcX_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS);
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcY_DA(TLorentzVector ep, TLorentzVector HFS);
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
void calcKin_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y);

// Calculate DIS kinematics - Sigma method
// SIDIS/Exclusive case:  e + p -> e' + p' + X
// Separating the scattered proton from the rest of the hadronic final state
Double_t calcQ2_Sigma(TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcX_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcY_Sigma(TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
void calcKin_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		   Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		   Float_t &Q2, Float_t &x, Float_t &y);
// DIS case:  e + p -> e' + HFS
// Combining whole HFS into one vector
Double_t calcQ2_Sigma(TLorentzVector ep, TLorentzVector HFS);
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcX_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS);
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcY_Sigma(TLorentzVector ep, TLorentzVector HFS);
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep,
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
void calcKin_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		   Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		   Float_t &Q2, Float_t &x, Float_t &y);

// Calculate DIS kinematics - eSigma method
// SIDIS/Exclusive case:  e + p -> e' + p' + X
// Separating the scattered proton from the rest of the hadronic final state
Double_t calcQ2_ESigma(TLorentzVector e, TLorentzVector ep);
Double_t calcQ2_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		       ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep);
Double_t calcQ2_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		       ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep);
Double_t calcX_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
Double_t calcY_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X);
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X);
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X);
void calcKin_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		    Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		    Float_t &Q2, Float_t &x, Float_t &y);
// DIS case:  e + p -> e' + HFS
// Combining whole HFS into one vector
// Don't need to redefine Q2 - Only uses beam and scattered electron
Double_t calcX_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS);
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
Double_t calcY_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS);
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS);
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS);
void calcKin_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		    Float_t &Q2, Float_t &x, Float_t &y);
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		    Float_t &Q2, Float_t &x, Float_t &y);

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//-----------------------------------------------------------------------------------------------------------------------------

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
Double_t calcX_Elec(TLorentzVector e, TLorentzVector  p, TLorentzVector  ep){
  TLorentzVector q = e-ep;
  double q2 = q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;
  return xB;
}
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep){
  P3EVector q = e-ep;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;
  return xB;
}
Double_t calcX_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep){
  P3MVector q = e-ep;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double x = q2/denom;
  return x;
}
Double_t calcY_Elec(TLorentzVector e, TLorentzVector p, TLorentzVector ep){
  TLorentzVector q = e - ep;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  double y = num/den;
  return y;
}
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep){
  P3EVector q = e - ep;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  double y = num/den;
  return y;
}
Double_t calcY_Elec(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep){
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

// Calculate inclusive kinematics (Q2, x, y) with the JB method
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
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();

  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
Double_t calcY_JB(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
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

// Calculate inclusive kinematics (Q2, x, y) with the DA method
// SIDIS case ep -> e'p'X
// Q2
Double_t calcQ2_DA(TLorentzVector e, TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
// x
Double_t calcX_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
// y
Double_t calcY_DA(TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}
// DIS case ep -> e'X
// Q2
Double_t calcQ2_DA(TLorentzVector e, TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
Double_t calcQ2_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep,  
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
// x
Double_t calcX_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
Double_t calcX_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
// y
Double_t calcY_DA(TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
Double_t calcY_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_DA(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}
void calcKin_DA(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}

// Calculate inclusive kinematics (Q2, x, y) with the Sigma method
// SIDIS case ep -> e'p'X
// Q2
Double_t calcQ2_Sigma(TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
// x
Double_t calcX_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
// y
Double_t calcY_Sigma(TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
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
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
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
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}
// DIS case ep -> e'X
// Q2
Double_t calcQ2_Sigma(TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
Double_t calcQ2_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
// x
Double_t calcX_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
Double_t calcX_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
// y
Double_t calcY_Sigma(TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
Double_t calcY_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		     ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_Sigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		   Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}
void calcKin_Sigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		   ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		   Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}

// Calculate inclusive kinematics (Q2, x, y) with the e-Sigma method
// SIDIS case ep -> e'p'X
// Q2
Double_t calcQ2_ESigma(TLorentzVector e, TLorentzVector ep){
  // Intermediate variables
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  return Q2_e;
}
Double_t calcQ2_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		       ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep){
  // Intermediate variables
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  return Q2_e;
}
Double_t calcQ2_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		       ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep){
  // Intermediate variables
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  return Q2_e;
}
// x
Double_t calcX_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
// y
Double_t calcY_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector pp, TLorentzVector X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pp, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> X, 
		    Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> pp, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> X, 
		    Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
// DIS case ep -> e'X
// NO NEED TO REDEFINE Q2 - ONLY DEPENDS ON BEAM AND SCATTERED ELECTRON
// x
Double_t calcX_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
Double_t calcX_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;  
}
// y
Double_t calcY_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
Double_t calcY_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
void calcKin_ESigma(TLorentzVector e, TLorentzVector p, TLorentzVector ep, TLorentzVector HFS, Float_t &Q2, Float_t &x, Float_t &y){
   // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Perp() * ep.Perp();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> HFS, 
		    Float_t &Q2, Float_t &x, Float_t &y){
   // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
void calcKin_ESigma(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> e, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> p, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> ep, 
		    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>> HFS, 
		    Float_t &Q2, Float_t &x, Float_t &y){
   // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
