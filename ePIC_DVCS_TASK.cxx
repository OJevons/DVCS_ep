// ePIC DVCS analysis class definition

// ROOT Includes
#include <TSystem.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <Math/Vector3D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <Math/LorentzVector.h>

// Class header include
#include "ePIC_DVCS_TASK.h"

//----------------------------------------------------
//----------------------------------------------------
//                    CONSTRUCTORS
//----------------------------------------------------
//----------------------------------------------------

// Default constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(){
}

// Specific constructor
ePIC_DVCS_TASK::ePIC_DVCS_TASK(TString camp, TString energy, TString sett){
  // Set stored campaign attributes
  setDate(camp);
  setEnergy(energy);
  setSetting(sett);

  // Set maximum particle momenta from beam momentum setting
  //setMomCuts(1.1);
  setMomCuts(2);
}


//----------------------------------------------------
//----------------------------------------------------
//                    SETTERS
//----------------------------------------------------
//----------------------------------------------------

// Set input file list
void ePIC_DVCS_TASK::setInFileList(TString name){
  sInList = name;
  std::cout<<"Input file list used: "<<name<<std::endl;
}

// Set output file name and create new
void ePIC_DVCS_TASK::setOutFile(TString name){
  std::cout<<"Output ROOT file: "<<name<<std::endl;
  fOutFile = new TFile(name,"RECREATE");
}

// Automatically set momentum cuts from energy string
// Apply multiplicative factor to cuts (default 1)
void ePIC_DVCS_TASK::setMomCuts(Float_t factor = 1.){
  // Use beam energy string to set maximum momenta
  if(sEnergy == "5x41"){
    fPMax_p=41.0*factor;
    fPMax_e=5.0*factor;
  }
  else if(sEnergy == "10x100"){
    fPMax_p=100.0*factor;
    fPMax_e=10.0*factor;
  }
  else if(sEnergy == "18x275"){
    fPMax_p=275.0*factor;
    fPMax_e=18.0*factor;
  }
  else{
    fPMax_p=100.0*factor;
    fPMax_e=10.0*factor;
  }

}


//----------------------------------------------------
//----------------------------------------------------
//                    APPLY CUTS
//----------------------------------------------------
//----------------------------------------------------

// Single particle cuts - electron
Bool_t ePIC_DVCS_TASK::applyCuts_Electron(std::vector<P3EVector> scate){
   Bool_t passCuts{kTRUE};
   
   // EVENT CUTS
   // Require single particle in final state
   if(scate.size() != 1) passCuts = kFALSE;
   // Return out of function if array is not filled
   if(!passCuts) return passCuts;

   // KINEMATIC CUTS
   // 1. Momentum
   if(scate[0].P() > fPMax_e) passCuts = kFALSE;
   // 2. Q2
   if(fQ2 < fMinQ2) passCuts = kFALSE;

   return passCuts;
}

// Single particle cuts - photon
Bool_t ePIC_DVCS_TASK::applyCuts_Photon(std::vector<P3EVector> scatg){
   Bool_t passCuts{kTRUE};
   
   // EVENT CUTS
   // Require single particle in final state
   if(scatg.size() != 1) passCuts = kFALSE;
   // Return out of function if array is not filled
   if(!passCuts) return passCuts;

   //----------------------------------
   // INSERT ANY OTHER PHOTON CUTS HERE
   //----------------------------------

   return passCuts;
}

// Single particle cuts - proton
Bool_t ePIC_DVCS_TASK::applyCuts_Proton(std::vector<P3EVector> scatp, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // EVENT CUTS
  // Require single particle in final state
  if(scatp.size() != 1) passCuts = kFALSE;
  // Return out of function if array is not filled
  if(!passCuts) return passCuts;
  
  // KINEMATIC CUTS
  // 1. Momentum
  if(scatp[0].P() > fPMax_p) passCuts = kFALSE;
  
  // 2. Scattered proton theta (ensure within B0, Roman Pots or 'all')
  // If invalid detector name used, consider all
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  Float_t fMinPTheta{0.};
  Float_t fMaxPTheta{0.};
  // Need to know beam proton energy for minimum theta in RP
  Float_t beamP{0};
  if(sEnergy == "5x41"){
    beamP=41.0;
  }
  else if(sEnergy == "10x100"){
    beamP=100.0;
  }
  else if(sEnergy == "18x275"){
    beamP=275.0;
  }
  else{
    beamP=100.0;
  }
  // RP momentum acceptance < 200 MeV
  
  // B0 angular acceptance: 5.5 mrad - 20 mrad
  if(sProtonDet == "B0"){
    fMinPTheta = 0.0055;
    fMaxPTheta = 0.02;
  }
  // RP angluar acceptance: < 5.0 mrad
  else if(sProtonDet == "RP"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.005;
  }
  // Full FF proton acceptance: < 20 mrad
  else if(sProtonDet == "all"){
    fMinPTheta = 0.;
    fMaxPTheta = 0.02;
  }
  if(scatp[0].Theta()<fMinPTheta || scatp[0].Theta()>fMaxPTheta) passCuts = kFALSE;

  return passCuts;
}

// Event-level cuts (DVCS kinematics)
Bool_t ePIC_DVCS_TASK::applyCuts_DVCS(TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};

  // 1. MAXIMUM T CUT FOR ROMAN POTS
  if(sProtonDet != "B0" && sProtonDet != "RP" && sProtonDet != "all") sProtonDet="all";
  if(sProtonDet == "RP" && ft > fMaxt_RP) passCuts = kFALSE;

  // 2. BJORKEN X CUT (removing tail from reconstructed histogram)
  if(TMath::Log10(fxB) < fxB_Tail) passCuts = kFALSE;
  
  // 3. MAXIMUM MISSING MASS^2
  if(TMath::Abs(fM2miss) > fMax_M2miss) passCuts = kFALSE;

  return passCuts;
}

// Combination of all cuts
Bool_t ePIC_DVCS_TASK::applyCuts_All(P3EVector beame, P3EVector beamp, vector<P3EVector> scate, vector<P3EVector> scatp, vector<P3EVector> scatg, TString sProtonDet="all"){
  Bool_t passCuts{kTRUE};
  
  // 1. Electron cuts
  // Need to calculate Q2 first - set to zero if no detected electron
  if(scate.size() == 0) fQ2 = 0;
  else fQ2 = calcQ2_Elec(beame, scate[0]);
  passCuts = applyCuts_Electron(scate);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 2. Photon cuts
  passCuts = applyCuts_Photon(scatg);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 3. Proton cuts
  passCuts = applyCuts_Proton(scatp, sProtonDet);
  // Exit from function if failure
  if(!passCuts) return passCuts;

  // 4. Event cuts
  // Need to calculate t, xB and MM2 first (e'p'y final state guaranteed by this point)
  fxB = calcX_Elec(beame, scate[0], beamp);
  fM2miss = calcM2Miss_3Body(beame, beamp, scate[0], scatp[0], scatg[0]);
  ft = calcT_BABE(beamp, scatp[0]);
  passCuts = applyCuts_DVCS(sProtonDet);

  return passCuts;
}

//----------------------------------------------------
//----------------------------------------------------
//            UNDO AFTERBURNER PROCEDURE
//----------------------------------------------------
//----------------------------------------------------

// Undo AB and calculate boost vectors - DO THIS FIRST FOR EACH EVENT
// USE BEAM VECTORS
void ePIC_DVCS_TASK::undoAfterburnAndCalc(P3EVector& p, P3EVector& k){
  // Holding vectors for beam - undoing crossing angle ONLY
  P3EVector p_beam(fXAngle*p.E(), 0., p.E(), p.E());
  P3EVector e_beam(0., 0., -k.E(), k.E());
  
  // Define boost vector to CoM frame
  P3EVector CoM_boost = p_beam+e_beam;
  vBoostToCoM.SetXYZ(-CoM_boost.X()/CoM_boost.E(), -CoM_boost.Y()/CoM_boost.E(), -CoM_boost.Z()/CoM_boost.E());
  
  // Apply boost to beam vectors
  p_beam = boost(p_beam, vBoostToCoM);
  e_beam = boost(e_beam, vBoostToCoM);
  
  // Calculate rotation angles and create rotation objects
  fRotY = -1.0*TMath::ATan2(p_beam.X(), p_beam.Z());
  fRotX = 1.0*TMath::ATan2(p_beam.Y(), p_beam.Z());

  rotAboutY = RotationY(fRotY);
  rotAboutX = RotationX(fRotX);

  // Apply rotation to beam vectors
  p_beam = rotAboutY(p_beam);
  p_beam = rotAboutX(p_beam);
  e_beam = rotAboutY(e_beam);
  e_beam = rotAboutX(e_beam);

  // Define boost vector back to head-on frame
  P3EVector HoF_boost(0., 0., CoM_boost.Z(), CoM_boost.E());
  vBoostToHoF.SetXYZ(HoF_boost.X()/HoF_boost.E(), HoF_boost.Y()/HoF_boost.E(), HoF_boost.Z()/HoF_boost.E());

  // Apply boost back to head on frame to beam vectors
  p_beam = boost(p_beam, vBoostToHoF);
  e_beam = boost(e_beam, vBoostToHoF);

  // Make changes to input vectors
  p.SetPxPyPzE(p_beam.X(), p_beam.Y(), p_beam.Z(), calcE(p_beam.Vect(),fMass_proton));
  k.SetPxPyPzE(e_beam.X(), e_beam.Y(), e_beam.Z(), calcE(e_beam.Vect(),fMass_electron));
}

// Undo afterburn procedure only
void ePIC_DVCS_TASK::undoAfterburn(P3EVector& a){
  Float_t mass = a.M();
  
  // Undo AB procedure for single vector, a^{mu}
  a = boost(a, vBoostToCoM); // BOOST TO COM FRAME
  a = rotAboutY(a);          // ROTATE TO Z-AXIS
  a = rotAboutX(a);          // ROTATE TO Z-AXIS
  a = boost(a, vBoostToHoF); // BOOST BACK TO HEAD ON FRAME

  a.SetPxPyPzE(a.X(), a.Y(), a.Z(), calcE(a.Vect(),mass));
}

//----------------------------------------------------
//----------------------------------------------------
//              KINEMATIC CALCULATIONS
//----------------------------------------------------
//----------------------------------------------------

// Calculate t from beam and scattered proton vectors
// t = (p' - p)^2
/*Double_t ePIC_DVCS_TASK::calcT(P3EVector p, P3EVector pprime){
  double t = (pprime - p).M2();
  
  return TMath::Abs(t);
}

// Calculate t without using proton information
// t = (k - k' - g)^2
Double_t ePIC_DVCS_TASK::calcTNoP(P3EVector k, P3EVector kprime, P3EVector g){
  double t = (k - kprime - g).M2();
  
  return TMath::Abs(t);
}

// Calculate Q2 from beam and scattered electron vectors
// Q2 = (k - k')^2
Double_t ePIC_DVCS_TASK::calcQ2(P3EVector k, P3EVector kprime){
  double q2 = (k - kprime).M2();
  double Q2 = -q2;

  return Q2;
}

// Calculate Bjorken x from both beam vectors and scattered electron vector
// xB = Q2 / 2(q.p)
Double_t ePIC_DVCS_TASK::calcBjorkenX(P3EVector k, P3EVector kprime, P3EVector p){
  P3EVector q = k - kprime;
  double q2 = -q.M2();
  double denom = 2*q.Dot(p);

  double xB = q2/denom;



  return xB;
}

// Calculate Bjorken y from both beam vectors and scattered electron vector
// y = (p.q)/(p.k)
Double_t ePIC_DVCS_TASK::calcBjorkenY(P3EVector k, P3EVector kprime, P3EVector p){
  P3EVector q = k - kprime;
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(k);

  double yB = num/den;

  return yB;
  }*/

// Calculate angle between hadronic and leptonic planes (Trento phi)
// Using planes defined by [k, q] and [q, p']
// Source: Bachetta, A. et al; Phys. Rev. D (2004); eq. 16
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qp(P3EVector k, P3EVector kprime, P3EVector pprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector pp3 = pprime.Vect();
  MomVector qhat3 = (k3-kp3).Unit();

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = qhat3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and scattered proton
  MomVector hNorm = qhat3.Cross(pp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If p' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,pp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

// Alternate Trento Phi calculation (using final state photon and proton)
// Angle between planes defined by [k, k'] and [g, p']
Double_t ePIC_DVCS_TASK::calcTrentoPhi_pg(P3EVector k, P3EVector kprime, P3EVector pprime, P3EVector gprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector pp3 = pprime.Vect();
  MomVector gp3 = gprime.Vect();

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = kp3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using scattered proton and final state photon
  MomVector hNorm = gp3.Cross(pp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If p' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,pp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

// Alternate Trento Phi calculation (using final state photon)
// Angle between planes defined by [k, k'] and [q, g]
Double_t ePIC_DVCS_TASK::calcTrentoPhi_qg(P3EVector k, P3EVector kprime, P3EVector gprime){  
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector q3 = k3-kp3;

  // Define leptonic plane using beam and scattered electron
  MomVector lNorm = kp3.Cross(k3);
  lNorm /= lNorm.R();
  // Define hadronic plane using q vector and final state photon
  MomVector hNorm = q3.Cross(gp3);
  hNorm /= hNorm.R();

  // Angle() function just returns magnitude of angle
  // If g' vector has a component parallel to the leptonic normal, should be positive. If opposite, negative.
  return TMath::Sign(1.,gp3.Dot(lNorm))*Angle(lNorm,hNorm);
}

// Alternate Trento Phi calculation using particle 4-vectors
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.3
Double_t ePIC_DVCS_TASK::calcTrentoPhi_4Vec(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime){
  // First, extract 4-vector components into arrays
  Float_t ppcomp[4]{0.,0.,0.,0.};
  Float_t kcomp[4]{0.,0.,0.,0.};
  pprime.GetCoordinates(ppcomp);
  k.GetCoordinates(kcomp);
  
  // Need q vector for gT calculations
  P3EVector q = k-kprime;
  
  // Need kT^2 and p'T^2 for denominator
  Float_t kTcomp[4]{0.,0.,0.,0.};
  Float_t ppTcomp[4]{0.,0.,0.,0.};
  
  // Loop over indices of 4-vectors
  Float_t numCos{0};
  Float_t numSin{0};
  
  for(int i{0}; i<4; i++){
    for(int j{0}; j<4; j++){
      kTcomp[i] += calcgT_ij(q,p,i,j)*kcomp[j];
      ppTcomp[i] += calcgT_ij(q,p,i,j)*ppcomp[j];

      numCos += kcomp[i]*ppcomp[j]*calcgT_ij(q,p,i,j);
      numSin += kcomp[i]*ppcomp[j]*calcepsT_ij(q,p,i,j);
    } // Loop over index j
  } // Loop over index i

  P3EVector kT(kTcomp[0],kTcomp[1],kTcomp[2],kTcomp[3]);
  P3EVector ppT(ppTcomp[0],ppTcomp[1],ppTcomp[2],ppTcomp[3]);

  Float_t kTSq = (kTcomp[3]*kTcomp[3])-(kTcomp[0]*kTcomp[0])-(kTcomp[1]*kTcomp[1])-(kTcomp[2]*kTcomp[2]);
  Float_t ppTSq = (ppTcomp[3]*ppTcomp[3])-(ppTcomp[0]*ppTcomp[0])-(ppTcomp[1]*ppTcomp[1])-(ppTcomp[2]*ppTcomp[2]);

  Float_t den = TMath::Sqrt(kTSq*ppTSq);
  Float_t cosphi = -numCos/den;
  Float_t sinphi = -numSin/den;

  return TMath::ACos(cosphi)*TMath::Sign(1.,TMath::ASin(sinphi));
}

// Calculate value of transverse metric tensor for 4-vector phi calculation
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.4
Double_t ePIC_DVCS_TASK::calcgT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j){
  // Extract 4-vector components into arrays
  Float_t pcomp[4]{0.,0.,0.,0.};
  Float_t qcomp[4]{0.,0.,0.,0.};
  p.GetCoordinates(pcomp);
  q.GetCoordinates(qcomp);

  // First term comes from standard metric tensor
  Int_t gij{0};
  
  // Metric tensor only contains diagonal components
  if(i==j && i==3) gij = 1;
  else if(i==j) gij = -1;
  
  // Second term comes from cross-product of 4-vector terms
  Float_t fQ2 = -q.M2();
  Float_t fpq = p.Dot(q);
  Float_t fGamma2 = (TMath::Power(fMass_proton,2)*fQ2)/TMath::Power(fpq,2);
  Float_t fPQTerm = qcomp[i]*pcomp[j] + pcomp[i]*qcomp[j];
  fPQTerm /= (1+fGamma2)*fpq;
  
  // Final term comes from products of terms in the same vector
  Float_t fQQTerm = (qcomp[i]*qcomp[j])/fQ2;
  fQQTerm -= (pcomp[i]*pcomp[j])/p.M2();
  fQQTerm *= fGamma2/(1+fGamma2);

  Float_t gT = gij - fPQTerm + fQQTerm;
  return gT;
}

// Calculate value of transverse Levi-Civita tensor for 4-vector phi calculation
// Source: Bachetta, A. et al; JHEP 02 (2007); eq. 2.4
Double_t ePIC_DVCS_TASK::calcepsT_ij(P3EVector q, P3EVector p, Int_t i, Int_t j){
  // Extract 4-vector components into arrays
  Float_t pcomp[4]{0.,0.,0.,0.};
  Float_t qcomp[4]{0.,0.,0.,0.};
  p.GetCoordinates(pcomp);
  q.GetCoordinates(qcomp);

  // Need gamma2 and p.q for calculation
  Float_t fQ2 = -q.M2();
  Float_t fpq = p.Dot(q);
  Float_t fGamma2 = (TMath::Power(fMass_proton,2)*fQ2)/TMath::Power(fpq,2);

  Float_t epsT{0};
  
  // Sum over all combinations of k and l
  for(int k{0}; k<4; k++){
    for(int l{0}; l<4; l++){
      epsT += LeviCivita(i,j,k,l)*(1./fpq)*(1./TMath::Sqrt(1+fGamma2))*pcomp[k]*qcomp[l];
    }
  }

  return epsT;
}

// Need Levi-Civita symbol for epsT calculation
// Using form of epsilon using Pi notation
Int_t ePIC_DVCS_TASK::LeviCivita(int i, int j, int k, int l){
  int ai[4]{i+1,j+1,k+1,l+1};

  int eps{1};

  for(int ind1{3}; ind1>0; ind1--){
    for(int ind2{0}; ind2<ind1; ind2++){
      // epsilon = 0 if any 2 indices are the same
      if(ai[ind1] == ai[ind2]) eps*=0;
      else eps*=TMath::Sign(1.,ai[ind1]-ai[ind2]);
    }
  }
  return eps;
}

// Calculate angle between planes of qp and qg
Double_t ePIC_DVCS_TASK::calcPhiQPQG(P3EVector k, P3EVector p, P3EVector kprime, P3EVector gprime){
  MomVector p3 = p.Vect();
  MomVector k3 = k.Vect();
  MomVector kp3 = kprime.Vect();
  MomVector gp3 = gprime.Vect();
  MomVector q3 = k3-kp3;

  MomVector qp = q3.Cross(p3);
  MomVector qg = q3.Cross(gp3);

  return Angle(qp,qg);
}

// Calculation of cone angle (angle between measured photon and expected photon)
Double_t ePIC_DVCS_TASK::calcConeAngle(P3EVector k, P3EVector p, P3EVector kprime, P3EVector pprime, P3EVector gprime){
  // Initial state vectors
  MomVector p3 = p.Vect(); // Proton beam
  MomVector k3 = k.Vect(); // Electron beam
  MomVector pi = p3+k3;    // Total initial momenta
  // Final state vectors
  MomVector pp3 = pprime.Vect(); // Scattered proton
  MomVector kp3 = kprime.Vect(); // Scattered electron
  MomVector gp3 = gprime.Vect(); // Real photon
  MomVector pfkp = pp3+kp3;      // Sum of scattered proton and scattered electron momenta

  // Return angle between real photon and expected photon (calc. from missing momentum of ep system)
  MomVector gExpected = pi-pfkp;

  return Angle(gp3,gExpected);
}

// Calculation of missing momentum for process ab->cdf
/*Double_t ePIC_DVCS_TASK::calcPMiss_3Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d, P3EVector f){
  Float_t fPMiss = (a+b-c-d-f).P();
  return fPMiss;
}

Double_t ePIC_DVCS_TASK::calcPtMiss_3Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d, P3EVector f){
  Float_t fPtMiss = (a+b-c-d-f).Pt();
  return fPtMiss;
}

// Calculation of missing energy for process ab->cdf
Double_t ePIC_DVCS_TASK::calcEMiss_3Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d, P3EVector f){
  Float_t fEMiss = (a+b-c-d-f).E();
  return fEMiss;
}

// Calculation of missing mass (squared) for process ab->cdf
Double_t ePIC_DVCS_TASK::calcM2Miss_3Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d, P3EVector f){
  Float_t fEMiss = (a+b-c-d-f).E();
  Float_t fPMiss = (a+b-c-d-f).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
}

// Calculation of missing momentum for process ab->cd
Double_t ePIC_DVCS_TASK::calcPMiss_2Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d){
  Float_t fPMiss = (a+b-c-d).P();
  return fPMiss;
}

// Calculation of missing energy for process ab->cd
Double_t ePIC_DVCS_TASK::calcEMiss_2Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d){
  Float_t fEMiss = (a+b-c-d).E();
  return fEMiss;
}

// Calculation of missing mass (squared) for process ab->cd
Double_t ePIC_DVCS_TASK::calcM2Miss_2Body(P3EVector a, P3EVector b, P3EVector c, P3EVector d){
  Float_t fEMiss = (a+b-c-d).E();
  Float_t fPMiss = (a+b-c-d).P();

  Float_t fM2Miss = TMath::Power(fEMiss,2) - TMath::Power(fPMiss,2);
  return fM2Miss;
  }*/

//----------------------------------------------------
//----------------------------------------------------
//                     DO ANALYSIS
//----------------------------------------------------
//----------------------------------------------------

void ePIC_DVCS_TASK::doAnalysis(){

  //---------------------------------------------------------
  // Setup: Load input file list
  //---------------------------------------------------------
  // IF TESTING, LOAD TEST FILE LIST
  if(sSett != "hiAcc" && sSett != "hiDiv") sInList="./filelists/inputFileList_test.list";
  //if(sSett != "hiAcc" && sSett != "hiDiv") sInList="./filelists/inputFileList_single.list";

  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;

  //---------------------------------------------------------
  // Setup: Declare histograms
  //---------------------------------------------------------
  //Particle kinematics - MC generated protons
  TH1D* h_eta_MCp   = new TH1D("eta_MCp",";#eta_{p'}(MC)", 275, -11.0, 11.0);
  TH1D* h_pt_MCp    = new TH1D("pt_MCp", ";p_{T, p'}(MC) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCp     = new TH1D("p_MCp", ";p_{p'}(MC) [GeV/c]", (Int_t)2*fPMax_p, 0.0, 2*fPMax_p);
  TH1D* h_theta_MCp = new TH1D("theta_MCp", ";#theta_{p'}(MC) [mrad]", 100, 0.0, 25.0);
  TH1D* h_phi_MCp   = new TH1D("phi_MCp", ";#phi_{p'}(MC) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCp     = new TH1D("E_MCp", ";E_{p'}(MC) [GeV]", (Int_t)1.5*fPMax_p, 0.0, 1.5*fPMax_p);

  //Particle kinematics - MC generated electrons
  TH1D* h_eta_MCe   = new TH1D("eta_MCe",";#eta_{e'}(MC)", 275, -11.0, 11.0);
  TH1D* h_pt_MCe    = new TH1D("pt_MCe", ";p_{T, e'}(MC) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCe     = new TH1D("p_MCe", ";p_{e'}(MC) [GeV/c]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_theta_MCe = new TH1D("theta_MCe", ";#theta_{e'}(MC) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_MCe   = new TH1D("phi_MCe", ";#phi_{e'}(MC) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCe     = new TH1D("E_MCe", ";E_{e'}(MC) [GeV]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);

  //Particle kinematics - MC generated photons
  TH1D* h_eta_MCg   = new TH1D("eta_MCg",";#eta_{#gamma}(MC)", 275, -11.0, 11.0);
  TH1D* h_pt_MCg    = new TH1D("pt_MCg", ";p_{T, #gamma}(MC) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCg     = new TH1D("p_MCg", ";p_{#gamma}(MC) [GeV/c]", 150, 0.0, 150.);
  TH1D* h_theta_MCg = new TH1D("theta_MCg", ";#theta_{#gamma}(MC) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_MCg   = new TH1D("phi_MCg", ";#phi_{#gamma}(MC) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCg     = new TH1D("E_MCg", ";E_{#gamma}(MC) [GeV]", 50, 0.0, 50.);

  //Particle kinematics - Associated MC scattered protons
  TH1D* h_eta_MCAp   = new TH1D("eta_MCAp",";#eta_{p'}(MC|Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_MCAp    = new TH1D("pt_MCAp", ";p_{T, p'}(MC|Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCAp     = new TH1D("p_MCAp", ";p_{p'}(MC|Reco) [GeV/c]", (Int_t)2*fPMax_p, 0.0, 2*fPMax_p);
  TH1D* h_theta_MCAp = new TH1D("theta_MCAp", ";#theta_{p'}(MC|Reco) [mrad]", 100, 0.0, 25.0);
  TH1D* h_phi_MCAp   = new TH1D("phi_MCAp", ";#phi_{p'}(MC|Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCAp     = new TH1D("E_MCAp", ";E_{p'}(MC|Reco) [GeV]", (Int_t)1.5*fPMax_p, 0.0, 1.5*fPMax_p);

  //Particle kinematics - Associated MC scattered electrons
  TH1D* h_eta_MCAe   = new TH1D("eta_MCAe",";#eta_{e'}(MC|Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_MCAe    = new TH1D("pt_MCAe", ";p_{T, e'}(MC|Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCAe     = new TH1D("p_MCAe", ";p_{e'}(MC|Reco) [GeV/c]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_theta_MCAe = new TH1D("theta_MCAe", ";#theta_{e'}(MC|Reco) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_MCAe   = new TH1D("phi_MCAe", ";#phi_{e'}(MC|Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCAe     = new TH1D("E_MCAe", ";E_{e'}(MC|Reco) [GeV]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);

  //Particle kinematics - Associated MC photons
  TH1D* h_eta_MCAg   = new TH1D("eta_MCAg",";#eta_{#gamma}(MC|Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_MCAg    = new TH1D("pt_MCAg", ";p_{T, #gamma}(MC|Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_MCAg     = new TH1D("p_MCAg", ";p_{#gamma}(MC|Reco) [GeV/c]", 150, 0.0, 150.);
  TH1D* h_theta_MCAg = new TH1D("theta_MCAg", ";#theta_{#gamma}(MC|Reco) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_MCAg   = new TH1D("phi_MCAg", ";#phi_{#gamma}(MC|Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_MCAg     = new TH1D("E_MCAg", ";E_{#gamma}(MC|Reco) [GeV]", 50, 0.0, 50.);

  //Particle kinematics - Reconstructed scattered protons - Barrel and B0
  TH1D* h_eta_RPp   = new TH1D("eta_RPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_RPp    = new TH1D("pt_RPp", ";p_{T, p'}(Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_RPp     = new TH1D("p_RPp", ";p_{p'}(Reco) [GeV/c]", (Int_t)2*fPMax_p, 0.0, 2*fPMax_p);
  TH1D* h_theta_RPp = new TH1D("theta_RPp", ";#theta_{p'}(Reco) [mrad]", 100, 0.0, 25.0);
  TH1D* h_phi_RPp   = new TH1D("phi_RPp", ";#phi_{p'}(Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_RPp     = new TH1D("E_RPp", ";E_{p'}(Reco) [GeV]", (Int_t)1.5*fPMax_p, 0.0, 1.5*fPMax_p);

  //Particle kinematics - Reconstructed scattered protons - Roman Pots
  TH1D* h_eta_RPPp   = new TH1D("eta_RPPp",";#eta_{p'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_RPPp    = new TH1D("pt_RPPp", ";p_{T, p'}(Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_RPPp     = new TH1D("p_RPPp", ";p_{p'}(Reco) [GeV/c]", (Int_t)2*fPMax_p, 0.0, 2*fPMax_p);
  TH1D* h_theta_RPPp = new TH1D("theta_RPPp", ";#theta_{p'}(Reco) [mrad]", 100, 0.0, 25.0);
  TH1D* h_phi_RPPp   = new TH1D("phi_RPPp", ";#phi_{p'}(Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_RPPp     = new TH1D("E_RPPp", ";E_{p'}(Reco) [GeV]", (Int_t)1.5*fPMax_p, 0.0, 1.5*fPMax_p);

  //Particle kinematics - Reconstructed scattered electrons - Barrel
  TH1D* h_eta_RPe   = new TH1D("eta_RPe",";#eta_{e'}(Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_RPe    = new TH1D("pt_RPe", ";p_{T, e'}(Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_RPe     = new TH1D("p_RPe", ";p_{e'}(Reco) [GeV/c]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_theta_RPe = new TH1D("theta_RPe", ";#theta_{e'}(Reco) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_RPe   = new TH1D("phi_RPe", ";#phi_{e'}(Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_RPe     = new TH1D("E_RPe", ";E_{e'}(Reco) [GeV]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);

  //Particle kinematics - Reconstructed photons - Barrel
  TH1D* h_eta_RPg   = new TH1D("eta_RPg",";#eta_{#gamma}(Reco)", 275, -11.0, 11.0);
  TH1D* h_pt_RPg    = new TH1D("pt_RPg", ";p_{T, #gamma}(Reco) [GeV/c]", 200, 0.0, 4.0);
  TH1D* h_p_RPg     = new TH1D("p_RPg", ";p_{#gamma}(Reco) [GeV/c]", 150, 0.0, 150.0);
  TH1D* h_theta_RPg = new TH1D("theta_RPg", ";#theta_{#gamma}(Reco) [rad]", 100, 0.0, 3.2);
  TH1D* h_phi_RPg   = new TH1D("phi_RPg", ";#phi_{#gamma}(Reco) [rad]", 100, -3.2, 3.2);
  TH1D* h_E_RPg     = new TH1D("E_RPg", ";E_{#gamma}(Reco) [GeV]", 50, 0.0, 50.);

  // DVCS event kinematics - Generated particles
  TH1D* h_t_MC    = new TH1D("t_MC"   , ";|t|(MC) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_Q2_MC   = new TH1D("Q2_MC"  , ";Q^{2}(MC) [(GeV/c^{2})^{2}]", 500, 0.0, 10.0);
  TH1D* h_xB_MC   = new TH1D("xB_MC"  , ";log_{10}(x_{B})(MC)"    , 100, -5, 0);
  TH1D* h_y_MC    = new TH1D("y_MC"   , ";y(MC)"                  , 100, 0, 1);
  TH1D* h_TPhi_MC = new TH1D("tphi_MC", ";#phi_{h}(MC) [rad]"     , 175, -3.5, 3.5);

  // DVCS event kinematics - Associated MC particles
  TH1D* h_t_MCA    = new TH1D("t_MCA"   , ";|t|(MC|Reco) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_Q2_MCA   = new TH1D("Q2_MCA"  , ";Q^{2}(MC|Reco) [(GeV/c^{2})^{2}]", 500, 0.0, 10.0);
  TH1D* h_xB_MCA   = new TH1D("xB_MCA"  , ";log_{10}(x_{B})(MC|Reco)"    , 100, -5, 0);
  TH1D* h_y_MCA    = new TH1D("y_MCA"   , ";y(MCA)"                      , 100, 0, 1);
  TH1D* h_TPhi_MCA = new TH1D("tphi_MCA", ";#phi_{h}(MC|Reco) [rad]"     , 175, -3.5, 3.5);
 
  // DVCS event kinematics - Reconstructed particles
  TH1D* h_t_RP     = new TH1D("t_RP"    , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_t_RPP    = new TH1D("t_RPP"   , ";|t|(Reco) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0); // Roman Pots
  TH1D* h_Q2_RP    = new TH1D("Q2_RP"   , ";Q^{2}(Reco) [(GeV/c^{2})^{2}]", 500, 0.0, 10.0);
  TH1D* h_xB_RP    = new TH1D("xB_RP"   , ";log_{10}(x_{B})(Reco)"    , 100, -5, 0);
  TH1D* h_y_RP     = new TH1D("y_RP"    , ";y(Reco)"                  , 100, 0, 1);
  TH1D* h_TPhi_RP  = new TH1D("tphi_RP" , ";#phi_{h}(Reco) [rad]"     , 175, -3.5, 3.5);
  TH1D* h_TPhi_RPP = new TH1D("tphi_RPP", ";#phi_{h}(Reco) [rad]"     , 175, -3.5, 3.5); // Roman Pots

  // 2D plots, associated MC vs reconstructed
  TH2D* h_p_2d = new TH2D("h_p_2d",";p_{p'}(MC|Reco) [GeV/#it{c}];p_{p'}(Reco) [GeV/#it{c}]", 150, 0., 150., 110, 0., 110.);
  TH2D* h_pt_2d = new TH2D("h_pt_2d",";p_{T, p'}(MC|Reco) [GeV/#it{c}];p_{T, p'}(Reco) [GeV/#it{c}]", 100, 0., 2., 100, 0., 2.);
  TH2D* h_t_2d = new TH2D("h_t_2d",";|t|(MC|Reco) [(GeV/#it{c})^{2}];|t|(Reco) [(GeV/#it{c})^{2}]", 100, 0., 2., 100, 0., 2.);
  TH2D* h_tRes_2d = new TH2D("h_tres_2d",";|t|_{MC} [(GeV/#it{c})^{2}];#delta t [(GeV/#it{c})^{2}]", 100, 0., 2., 100, -5., 5.);
  TH1D* h_extracted_t_resolution;
  
  // Angular distances between final state particles
  TH1D* h_dphi_MCeg = new TH1D("dphi_MCeg",";#delta#phi_{e'#gamma} [rad]"      ,100,-3.2,3.2);
  TH1D* h_dtheta_MCeg = new TH1D("dtheta_MCeg",";#delta#theta_{e'#gamma} [rad]",100,-3.2,3.2);
  TH1D* h_dphi_MCAeg = new TH1D("dphi_MCAeg",";#delta#phi_{e'#gamma} [rad]"    ,100,-3.2,3.2);
  TH1D* h_dtheta_MCAeg = new TH1D("dtheta_MCAeg",";#delta#theta_{e'#gamma} [rad]",100,-3.2,3.2);
  TH1D* h_dphi_RPeg = new TH1D("dphi_RPeg",";#delta#phi_{e'#gamma} [rad]"        ,100,-3.2,3.2);
  TH1D* h_dtheta_RPeg = new TH1D("dtheta_RPeg",";#delta#theta_{e'#gamma} [rad]"  ,100,-3.2,3.2);

  // Occupancies - B0 tracker
  TH2D* h_B0_occupancy_map_layer_0 = new TH2D("B0_occupancy_map_0", "B0_occupancy_map_0;x [mm];y [mm]", 100, -400, 0, 100, -170, 170);
  TH2D* h_B0_occupancy_map_layer_1 = new TH2D("B0_occupancy_map_1", "B0_occupancy_map_1;x [mm];y [mm]", 100, -400, 0, 100, -170, 170);
  TH2D* h_B0_occupancy_map_layer_2 = new TH2D("B0_occupancy_map_2", "B0_occupancy_map_2;x [mm];y [mm]", 100, -400, 0, 100, -170, 170);
  TH2D* h_B0_occupancy_map_layer_3 = new TH2D("B0_occupancy_map_3", "B0_occupancy_map_3;x [mm];y [mm]", 100, -400, 0, 100, -170, 170);
  TH1D* h_B0_hit_energy_deposit = new TH1D("B0_tracker_hit_energy_deposit", ";Deposited Energy [keV]", 100, 0.0, 500.0);
  
  // Occupancies - B0 ECAL
  TH2D* h_B0_emcal_occupancy_map = new TH2D("B0_emcal_occupancy_map", ";hit x [mm];hit y [mm]", 100, -400, 0, 100, -170, 170);
  TH1D* h_B0_emcal_cluster_energy = new TH1D("B0_emcal_cluster_energy", ";Cluster Energy [GeV]", 100, 0.0, 100.0);

  // Occupancies - Roman Pots
  //TH2D* h_rp_occupancy_map = new TH2D("Roman_pots_occupancy_map", ";hit x [mm];hit y [mm]", 100, -150, 150, 100, -80, 80);
  TH2D* h_rp_occupancy_map = new TH2D("Roman_pots_occupancy_map", ";hit x [mm];hit y [mm]", 100, -1300, -1000, 100, -80, 80);
  TH1D* h_rp_z = new TH1D("Roman_pots_z", ";hit z [mm]", 50, 25990, 26005);

  // B0 momentum resolution
  TH1D* h_b0_pt_resolution = new TH1D("b0_pt_resolution", ";#Delta p_{T} [GeV/c]", 100, -2.0, 2.0);
  TH2D* h_b0_pt_resolution_percent = new TH2D("b0_deltaPt_over_pt_vs_pt", ";P_{T, MC} [GeV/c]; #Delta p_{T}/p_{T, MC} [percent/100]", 100, 0.0, 2.0, 100, -1.0, 1.0);	
  TH2D* h_b0_p_resolution_percent = new TH2D("b0_deltaP_over_p_vs_p", ";Three-Momentum,  p_{MC} [GeV/c]; #Delta p/p_{MC} [percent/100]", 550, 0., 110., 100, -1.0, 1.0);	
  TH1D* h_b0_extracted_pt_resolution;
  TH1D* h_b0_extracted_p_resolution;

  // Kinematic conservation - full final state
  TH1D* h_Emiss3_MC = new TH1D("Emiss3_MC",";E_{miss}(MC) [GeV]",110,-5,50);
  TH1D* h_Pmiss3_MC = new TH1D("Pmiss3_MC",";p_{miss}(MC) [GeV/#it{c}]",110,-5,50);
  TH1D* h_Ptmiss3_MC = new TH1D("Ptmiss3_MC",";p_{T,miss}(MC) [GeV/#it{c}]",60,-2,10);
  TH1D* h_M2miss3_MC = new TH1D("M2miss3_MC",";M^{2}_{miss}(MC) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  TH1D* h_Emiss3_MCA = new TH1D("Emiss3_MCA",";E_{miss}(MC|Reco) [GeV]",110,-5,50);
  TH1D* h_Pmiss3_MCA = new TH1D("Pmiss3_MCA",";p_{miss}(MC|Reco) [GeV/#it{c}]",110,-5,50);
  TH1D* h_Ptmiss3_MCA = new TH1D("Ptmiss3_MCA",";p_{T,miss}(MC|Reco) [GeV/#it{c}]",60,-2,10);
  TH1D* h_M2miss3_MCA = new TH1D("M2miss3_MCA",";M^{2}_{miss}(MC|Reco) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  TH1D* h_Emiss3_RP = new TH1D("Emiss3_RP",";E_{miss}(Reco) [GeV]",110,-5,50);
  TH1D* h_Pmiss3_RP = new TH1D("Pmiss3_RP",";p_{miss}(Reco) [GeV/#it{c}]",110,-5,50);
  TH1D* h_Ptmiss3_RP = new TH1D("Ptmiss3_RP",";p_{T,miss}(Reco) [GeV/#it{c}]",60,-2,10);
  TH1D* h_M2miss3_RP = new TH1D("M2miss3_RP",";M^{2}_{miss}(Reco) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  // Kinematic conservation - partial final state
  TH1D* h_M2miss2ep_MC = new TH1D("M2miss2ep_MC",";M^{2}_{miss, e'p'}(MC) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  TH1D* h_M2miss2ep_MCA = new TH1D("M2miss2ep_MCA",";M^{2}_{miss, e'p'}(MCA) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  TH1D* h_M2miss2ep_RP = new TH1D("M2miss2ep_RP",";M^{2}_{miss, e'p'}(Reco) [(GeV/#it{c}^{2})^{2}]",200,-50,50);
  TH1D* h_M2miss2eg_MC = new TH1D("M2miss2eg_MC",";M^{2}_{miss, e'#gamma}(MC) [(GeV/#it{c}^{2})^{2}]",300,-500,1000);
  TH1D* h_M2miss2eg_MCA = new TH1D("M2miss2eg_MCA",";M^{2}_{miss, e'#gamma}(MCA) [(GeV/#it{c}^{2})^{2}]",300,-500,1000);
  TH1D* h_M2miss2eg_RP = new TH1D("M2miss2eg_RP",";M^{2}_{miss, e'#gamma}(Reco) [(GeV/#it{c}^{2})^{2}]",300,-500,1000);
  
  // Exclusivity variables - qpqg angle and cone angle
  TH1D* h_QPQG_MC  = new TH1D("qpqg_MC",";#theta_{[#gamma*p'][#gamma*#gamma]} [rad]",200,0,3.2);
  TH1D* h_QPQG_MCA = new TH1D("qpqg_MCA",";#theta_{[#gamma*p'][#gamma*#gamma]}[rad]",200,0,3.2);
  TH1D* h_QPQG_RP  = new TH1D("qpqg_RP",";#theta_{[#gamma*p'][#gamma*#gamma]} [rad]",200,0,3.2);
  TH1D* h_Cone_MC  = new TH1D("cone_MC",";#theta_{#gamma#gamma} [rad]",200,0,3.2);
  TH1D* h_Cone_MCA = new TH1D("cone_MCA",";#theta_{#gamma#gamma} [rad]",200,0,3.2);
  TH1D* h_Cone_RP  = new TH1D("cone_RP",";#theta_{#gamma#gamma} [rad]",200,0,3.2);
  TH1D* h_Cone_RPP = new TH1D("cone_RPP",";#theta_{#gamma#gamma} [rad]",200,0,3.2);

  // Photon resolution plots - energy and angle (reco - gen)
  TH1D* h_PhotRes_E = new TH1D("photres_E",";E_{#gamma}(Reco)-E_{#gamma}(MC) [GeV]",50,0,10);
  TH1D* h_PhotRes_theta = new TH1D("photres_theta",";#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [rad]",600,-1.5,1.5);
  TH2D* h_PhotRes2D_theta = new TH2D("photres2d_theta",";#theta_{#gamma, MC} [rad]; #delta#theta_{#gamma}",320,0,3.2,600,-1.5,1.5);
  TH2D* h_PhotRes2D_thetaReco = new TH2D("photres2d_thetareco",";#theta_{#gamma, Reco} [rad]; #delta#theta_{#gamma}",320,0,3.2,600,-1.5,1.5);
  
  // 2D kinematic distributions
  // x - always on x-axis
  // Q2 - always on y-axis
  // t - fill as appropriate
  TH2D* h_2D_xVt_MC = new TH2D("2d_xvt_mc",";log_{10}(x_{B,MC});|t|_{MC} [GeV^{2}]",100,-5.,0.,100,0.,2.);
  TH2D* h_2D_xVt_MCA = new TH2D("2d_xvt_mca",";log_{10}(x_{B,MC|Reco});|t|_{MC|Reco} [GeV^{2}]",100,-5.,0.,100,0.,2.);
  TH2D* h_2D_xVt_RP = new TH2D("2d_xvt_rp",";log_{10}(x_{B,Reco});|t|_{Reco} [GeV^{2}]",100,-5.,0.,100,0.,2.);
  TH2D* h_2D_xVQ2_MC = new TH2D("2d_xvq2_mc",";log_{10}(x_{B,MC});Q^{2}_{MC} [GeV^{2}]",100,-5.,0.,450,1.,10.);
  TH2D* h_2D_xVQ2_MCA = new TH2D("2d_xvq2_mca",";log_{10}(x_{B,MC|Reco});Q^{2}_{MC|Reco} [GeV^{2}]",100,-5.,0.,450,1.,10.);
  TH2D* h_2D_xVQ2_RP = new TH2D("2d_xvq2_rp",";log_{10}(x_{B,Reco});Q^{2}_{Reco} [GeV^{2}]",100,-5.,0.,450,1.,10.);
  TH2D* h_2D_tVQ2_MC = new TH2D("2d_tvq2_mc",";t_{MC} [GeV^{2}];Q^{2}_{MC} [GeV^2]",100,0.,2.,450,1.,10.);
  TH2D* h_2D_tVQ2_MCA = new TH2D("2d_tvq2_mca",";t_{MC|Reco} [GeV^{2}];Q^{2}_{MC|Reco} [GeV^2]",100,0.,2.,450,1.,10.);
  TH2D* h_2D_tVQ2_RP = new TH2D("2d_tvq2_rp",";t_{Reco} [GeV^{2}];Q^{2}_{Reco} [GeV^2]",100,0.,2.,450,1.,10.);

  // 2D coverage distributions
  TH2D* h_2D_EvEta_g = new TH2D("2d_eveta_g",";#eta_{#gamma};E_{#gamma} [GeV]",200,-4.,4.,100,0.,50.);
  TH2D* h_2D_EvEta_e = new TH2D("2d_eveta_e",";#eta_{e'};E_{e'} [GeV]",200,-4.,4.,(Int_t)1.1*fPMax_e, 0., 1.1*fPMax_e);
  TH2D* h_2D_EvEta_p = new TH2D("2d_eveta_p",";#eta_{p'};E_{p'} [GeV]",150,4.,10.,(Int_t)1.1*fPMax_p, 0., 1.1*fPMax_p);
  
  // Inclusive t-distributions (only proton detected)
  TH1D* h_tInc_MC = new TH1D("tInc_MC"   , ";|t|(MC) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_tInc_RP = new TH1D("tInc_RP"   , ";|t|(Reco. - B0) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);
  TH1D* h_tInc_RPP = new TH1D("tInc_RPP"   , ";|t|(Reco. - RP) [(GeV/c^{2})^{2}]"  , 100, 0.0, 2.0);

  // Other calculations for trento phi
  TH1D* h_TPhi_pg_MC  = new TH1D("tphi_pg_MC",  ";#phi_{h}(MC) [rad]"   , 175, -3.5, 3.5);
  TH1D* h_TPhi_pg_RP  = new TH1D("tphi_pg_RP" , ";#phi_{h}(Reco) [rad]" , 175, -3.5, 3.5);
  TH1D* h_TPhi_pg_RPP = new TH1D("tphi_pg_RPP", ";#phi_{h}(Reco) [rad]" , 175, -3.5, 3.5);
  TH1D* h_TPhi_qg_MC  = new TH1D("tphi_qg_MC",  ";#phi_{h}(MC) [rad]"   , 175, -3.5, 3.5);
  TH1D* h_TPhi_qg_RP  = new TH1D("tphi_qg_RP" , ";#phi_{h}(Reco) [rad]" , 175, -3.5, 3.5);
  TH1D* h_TPhi_qg_RPP = new TH1D("tphi_qg_RPP", ";#phi_{h}(Reco) [rad]" , 175, -3.5, 3.5);

  // Bjorken-x - LINEAR
  TH1D* h_xBLin_MC = new TH1D("xblin_mc",";x_{B}",1000000,0.,1.);
  
  // Momenta plots for RP tracks
  TH1D* h_px_RPP = new TH1D("px_RPP",";p_{x} [GeV/#it{c}]",200,-10.,10.);
  TH1D* h_py_RPP = new TH1D("py_RPP",";p_{y} [GeV/#it{c}]",200,-10.,10.);
  TH1D* h_pz_RPP = new TH1D("pz_RPP",";p_{z} [GeV/#it{c}]",(Int_t)3*fPMax_p, -1.5*fPMax_p, 1.5*fPMax_p);

  // 2D eta vs p distributions - species inclusive

  //---------------------------------------------------------
  // Loop over files in list
  //---------------------------------------------------------
  int fileCounter{0};

  // 4-vectors for beam particles - need these defined outside of file loop
  P3EVector beame4(0,0,0,-1);     // Beam electron (generated)
  P3EVector beamp4(0,0,0,-1);     // Beam proton (generated)
  
  // Start file loop
  while(getline(fileListStream,fileName)){
    //---------------------------------------------------------
    // Open 1 file from list at a time, then get event TTree
    //---------------------------------------------------------
    // Get file
    TString tmp{fileName};
    std::cout<<"Input file "<<fileCounter<<" : "<<tmp<<std::endl;
    auto inputRootFile = TFile::Open(tmp);
    if(!inputRootFile){ std::cout<<"MISSING_ROOT_FILE"<<tmp<<endl; continue;}
    fileCounter++;

    // Get TTree
    TTree * evtTree = (TTree*)inputRootFile->Get("events");
    if (!(inputRootFile->GetListOfKeys()->Contains("events"))) continue;
    int numEvents = evtTree->GetEntries();
    std::cout<<"File has "<<numEvents<<" events..."<<std::endl;

    //---------------------------------------------------------
    // Declare TTreeReader, and choose appropriate branches
    //---------------------------------------------------------
    TTreeReader tree_reader(evtTree);
    // MC particles
    TTreeReaderArray<float>  mc_px_array        = {tree_reader, "MCParticles.momentum.x"};
    TTreeReaderArray<float>  mc_py_array        = {tree_reader, "MCParticles.momentum.y"};
    TTreeReaderArray<float>  mc_pz_array        = {tree_reader, "MCParticles.momentum.z"};
    TTreeReaderArray<double> mc_mass_array      = {tree_reader, "MCParticles.mass"};
    TTreeReaderArray<int>    mc_genStatus_array = {tree_reader, "MCParticles.generatorStatus"};
    TTreeReaderArray<int>    mc_pdg_array       = {tree_reader, "MCParticles.PDG"};    
    // Reconstructed/MC particle associations - BARREL (using ReconstructedParticles branch)
    TTreeReaderArray<unsigned int> assoc_rec_id = {tree_reader, "ReconstructedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> assoc_sim_id = {tree_reader, "ReconstructedParticleAssociations.simID"};
    // Reconstructed particles - BARREL (using ReconstructedParticles branch)
    TTreeReaderArray<float>  re_px_array     = {tree_reader, "ReconstructedParticles.momentum.x"};
    TTreeReaderArray<float>  re_py_array     = {tree_reader, "ReconstructedParticles.momentum.y"};
    TTreeReaderArray<float>  re_pz_array     = {tree_reader, "ReconstructedParticles.momentum.z"};
    TTreeReaderArray<float>  re_e_array      = {tree_reader, "ReconstructedParticles.energy"};
    TTreeReaderArray<float>  re_charge_array = {tree_reader, "ReconstructedParticles.charge"};
    TTreeReaderArray<float>  re_mass_array   = {tree_reader, "ReconstructedParticles.mass"};
    TTreeReaderArray<int>    re_pdg_array    = {tree_reader, "ReconstructedParticles.PDG"};
    // Reconstructed/MC particle associations - B0 (using ReconstructedTruthSeededChargedParticles branch)
    TTreeReaderArray<unsigned int> tsassoc_rec_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.recID"};
    TTreeReaderArray<unsigned int> tsassoc_sim_id = {tree_reader, "ReconstructedTruthSeededChargedParticleAssociations.simID"};
    // Reconstructed particles - B0 (using ReconstructedTruthSeededChargedParticles branch)
    TTreeReaderArray<float>  tsre_px_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.x"};
    TTreeReaderArray<float>  tsre_py_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.y"};
    TTreeReaderArray<float>  tsre_pz_array     = {tree_reader, "ReconstructedTruthSeededChargedParticles.momentum.z"};
    TTreeReaderArray<float>  tsre_e_array      = {tree_reader, "ReconstructedTruthSeededChargedParticles.energy"};
    TTreeReaderArray<float>  tsre_charge_array = {tree_reader, "ReconstructedTruthSeededChargedParticles.charge"};
    TTreeReaderArray<float>  tsre_mass_array   = {tree_reader, "ReconstructedTruthSeededChargedParticles.mass"};
    // B0 tracker hits
    TTreeReaderArray<float> b0_hits_x    = {tree_reader, "B0TrackerRecHits.position.x"};
    TTreeReaderArray<float> b0_hits_y    = {tree_reader, "B0TrackerRecHits.position.y"};
    TTreeReaderArray<float> b0_hits_z    = {tree_reader, "B0TrackerRecHits.position.z"};
    TTreeReaderArray<float> b0_hits_eDep = {tree_reader, "B0TrackerRecHits.edep"}; //deposited energy per hit
    // B0 ECAL hits
    TTreeReaderArray<float> b0_cluster_x      = {tree_reader, "B0ECalClusters.position.x"};
    TTreeReaderArray<float> b0_cluster_y      = {tree_reader, "B0ECalClusters.position.y"};
    TTreeReaderArray<float> b0_cluster_z      = {tree_reader, "B0ECalClusters.position.z"};
    TTreeReaderArray<float> b0_cluster_energy = {tree_reader, "B0ECalClusters.energy"}; //deposited energy in cluster
    // RP hits
    TTreeReaderArray<float> global_hit_RP_x = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.x"};
    TTreeReaderArray<float> global_hit_RP_y = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.y"};
    TTreeReaderArray<float> global_hit_RP_z = {tree_reader, "ForwardRomanPotRecParticles.referencePoint.z"};
    TTreeReaderArray<float> rp_px_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.x"};
    TTreeReaderArray<float> rp_py_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.y"};
    TTreeReaderArray<float> rp_pz_array     = {tree_reader, "ForwardRomanPotRecParticles.momentum.z"};
    TTreeReaderArray<float> rp_mass_array   = {tree_reader, "ForwardRomanPotRecParticles.mass"};
    TTreeReaderArray<int>   rp_pdg_array    = {tree_reader, "ForwardRomanPotRecParticles.PDG"};

    tree_reader.SetEntriesRange(0, evtTree->GetEntries());
    
    // If averaging beams from file
    if(!kUseEventBeams){
      // MUST DO THIS FIRST
      // Full run over tree in first file before anything else
      // Calculate beams from average of individual event beam particles
      if(fileCounter == 1){
	// Accumulator variables
	P3EVector beame4_acc(0,0,0,-1);
	P3EVector beamp4_acc(0,0,0,-1);
	
	while (tree_reader.Next()){
	  // Beams for each event
	  P3EVector beame4_evt(0,0,0,-1);
	  P3EVector beamp4_evt(0,0,0,-1);
	  TVector3 mctrk;
	  for(int imc=0;imc<mc_px_array.GetSize();imc++){
	    mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);

	    // Beam particles ==> Generator Status 4
	    if(mc_genStatus_array[imc] == 4){
	      // Proton
	      if(mc_pdg_array[imc] == 2212){ 
		beamp4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	      }
	      // Electron
	      else if(mc_pdg_array[imc] == 11){ 
		beame4_evt.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	      }		
	    } // Found beam particles for event
	  } // End MCParticles loop
	  
	  // Add found beams to accumulator
	  beame4_acc += beame4_evt;
	  beamp4_acc += beamp4_evt;
	} // Tree complete
	
	// Divide by number of events in file
	beame4.SetCoordinates(beame4_acc.X()/evtTree->GetEntries(), beame4_acc.Y()/evtTree->GetEntries(), beame4_acc.Z()/evtTree->GetEntries(), beame4_acc.E()/evtTree->GetEntries());
	beamp4.SetCoordinates(beamp4_acc.X()/evtTree->GetEntries(), beamp4_acc.Y()/evtTree->GetEntries(), beamp4_acc.Z()/evtTree->GetEntries(), beamp4_acc.E()/evtTree->GetEntries());
      } // File complete
      
      // Undo afterburn on beam particles and calculate "postburn" variables
      undoAfterburnAndCalc(beamp4,beame4);
    }

    //---------------------------------------------------------
    // Loop over all entries in event TTree
    //---------------------------------------------------------
    while (tree_reader.Next()){ 
      // 4-vectors for MC raw particles
      vector<P3EVector> scate4_gen;   // Scattered electron (generated)
      vector<P3EVector> scatp4_gen;   // Scattered proton (generated)
      vector<P3EVector> scatg4_gen;   // Scattered photon (generated)
      // 4-vectors for associated MC particles (ONLY SCATTERED)
      vector<P3EVector> scate4_aso;   // Scattered electron (associated MC)
      vector<P3EVector> scatp4_aso;   // Scattered proton (associated MC)
      vector<P3EVector> scatg4_aso;   // Scattered photon (associated MC)
      // 4-vectors for reconstructed particles (SEPARATE PROTONS FOR B0 AND ROMAN POTS)
      vector<P3EVector> scate4_rec;   // Scattered electron (reconstructed)
      vector<P3EVector> scatp4_rec;   // Scattered proton (B0 reconstructed)
      vector<P3EVector> scatp4_rom;   // Scattered proton (Roman Pots reconstructed)
      vector<P3EVector> scatg4_rec;   // Scattered photon (reconstructed)      

      // Holding 3-vectors
      TVector3 mctrk, assoctrk, recotrk;
      
      // Get beams for each event
      if(kUseEventBeams){
	for(int imc=0;imc<mc_px_array.GetSize();imc++){
	  mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
	  
	  // Beam particles ==> Generator Status 4
	  if(mc_genStatus_array[imc] == 4){
	    // Proton
	    if(mc_pdg_array[imc] == 2212){ 
	      beamp4.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	    }
	    // Electron
	    else if(mc_pdg_array[imc] == 11){ 
	      beame4.SetCoordinates(mctrk.X(), mctrk.Y(), mctrk.Z(), calcE(mctrk,mc_mass_array[imc]));
	    }		
	  } // Found beam particles
	}
	// Undo afterburn on beam particles and calculate "postburn" variables
	undoAfterburnAndCalc(beamp4,beame4);
      }

      //---------------------------------------------------------
      // Fill particle holding arrays
      //---------------------------------------------------------
      // 1. MC generated
      for(int imc=0;imc<mc_px_array.GetSize();imc++){
	mctrk.SetXYZ(mc_px_array[imc], mc_py_array[imc], mc_pz_array[imc]);
	P3EVector q_scat(mctrk.X(),mctrk.Y(),mctrk.Z(),calcE(mctrk,mc_mass_array[imc]));
	
	// Undo afterburner
	undoAfterburn(q_scat);

	// Look for scattered particles ==> Generator status 1
	if(mc_genStatus_array[imc] == 1){
	  if(mc_pdg_array[imc] == 2212){
	    scatp4_gen.push_back(q_scat);
	  }
	  if(mc_pdg_array[imc] == 11){
	    scate4_gen.push_back(q_scat);
	  }
	  if(mc_pdg_array[imc] == 22){
	    scatg4_gen.push_back(q_scat);
	  }
	} // Found scattered particles
      }// End of generated particles loop
      
      // 2. and 3. Associated MC tracks and their matched reconstucted tracks
      // ONLY IF USING EXPLICIT TRACK MATCHING
      if(kUseExplicitMatch){
	unsigned int mc_assoc_index = -1;
	// LOOK FOR ELECTRONS AND PHOTONS (using ReconstructedParticleAssociations)
	for(unsigned int iAssoc{0};iAssoc<assoc_rec_id.GetSize();iAssoc++){
	  mc_assoc_index = assoc_sim_id[iAssoc];

	  // If reco track isn't associated to an MC track, then skip
	  if(mc_assoc_index == -1) continue;
	
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]); 
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_assoc);
	  recotrk.SetXYZ(re_px_array[iAssoc], re_py_array[iAssoc], re_pz_array[iAssoc]);
	  P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_reco);

	  // Fill track vectors based on associated PID
	  // Electrons
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 11){ 
	    scate4_aso.push_back(q_assoc); 
	    scate4_rec.push_back(q_reco); 
	  }
	  // Photons
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 22){ 
	    scatg4_aso.push_back(q_assoc); 
	    scatg4_rec.push_back(q_reco); 
	  }
	} // End of associations loop

	mc_assoc_index=-1; // Reset association index
	// THEN LOOK FOR PROTONS (using ReconstructedTruthSeededChargedParticleAssociations)
	for(unsigned int iTSAssoc{0};iTSAssoc<tsassoc_rec_id.GetSize();iTSAssoc++){
	  mc_assoc_index = tsassoc_sim_id[iTSAssoc];

	  // Only care about protons here (PID 2212)
	  if(mc_assoc_index != -1 && mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 2212){
	    assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]);
	    P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	    undoAfterburn(q_assoc);
	    recotrk.SetXYZ(tsre_px_array[iTSAssoc], tsre_py_array[iTSAssoc], tsre_pz_array[iTSAssoc]);
	    P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,mc_mass_array[mc_assoc_index]));
	    undoAfterburn(q_reco);

	    scatp4_aso.push_back(q_assoc); 
	    scatp4_rec.push_back(q_reco); 
	  }
	} // End of truth-seeded association loop
      }

      // 2. and 3. Associated MC and reconstructed tracks
      // IF NOT USING EXPLICIT MATCHING (default behaviour)
      else if(!kUseExplicitMatch){
	// 2. Associated MC (scattered)
	for(unsigned int iAssoc{0};iAssoc<assoc_rec_id.GetSize();iAssoc++){
	  unsigned int mc_assoc_index = assoc_sim_id[iAssoc];
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]);
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  // Undo afterburner
	  undoAfterburn(q_assoc);
	  // Look for scattered particles ==> Generator status 1
	  if(mc_genStatus_array[mc_assoc_index] == 1){
	    if(mc_pdg_array[mc_assoc_index] == 11){ scate4_aso.push_back(q_assoc); }
	    if(mc_pdg_array[mc_assoc_index] == 22){ scatg4_aso.push_back(q_assoc); }
	  } // Found associated particles
	}// End of associated particles loop
	
	// 2a. Associated MC protons (found using different association branch)
	for(unsigned int iTSAssoc{0};iTSAssoc<tsassoc_rec_id.GetSize();iTSAssoc++){
	  unsigned int mc_assoc_index = tsassoc_sim_id[iTSAssoc];
	  assoctrk.SetXYZ(mc_px_array[mc_assoc_index], mc_py_array[mc_assoc_index], mc_pz_array[mc_assoc_index]); 
	  P3EVector q_assoc(assoctrk.X(),assoctrk.Y(),assoctrk.Z(),calcE(assoctrk,mc_mass_array[mc_assoc_index]));
	  undoAfterburn(q_assoc);
	  
	  if(mc_genStatus_array[mc_assoc_index] == 1 && mc_pdg_array[mc_assoc_index] == 2212){ scatp4_aso.push_back(q_assoc); }
	  
	} // End of truth-seeded association loop
      
	// 3. Reconstructed particles
	// Start with ACTS reconstructed particles (barrel and B0)
	for(int ireco{0}; ireco<re_px_array.GetSize(); ireco++){
	  recotrk.SetXYZ(re_px_array[ireco], re_py_array[ireco], re_pz_array[ireco]);
	  // If not using ePIC PID, assume particles based on charge of track
	  // Look for electrons and photons from ReconstructedParticles branch
	  if(!kUsePID){
	    // Negative => ELECTRON
	    if(re_charge_array[ireco] == -1){
	      P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,fMass_electron));
	      undoAfterburn(q_reco);
	      scate4_rec.push_back(q_reco);
	    }
	    // Neutral => REAL PHOTON
	    if(re_charge_array[ireco] == 0){
	      P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),recotrk.Mag());
	      undoAfterburn(q_reco);
	      scatg4_rec.push_back(q_reco);
	    }
	  }
	  // Using ePIC PID
	  else{
	    P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,re_mass_array[ireco]));
	    // Undo afterburner
	    undoAfterburn(q_reco);
	    if(re_pdg_array[ireco] == 11){ scate4_rec.push_back(q_reco); }
	    if(re_pdg_array[ireco] == 22){ scatg4_rec.push_back(q_reco); }
	  }
	}// End of ACTS reconstructed particles loop

	// 3a. ReconstructedB0 protons (ReconstructedTruthSeededChargedParticles branch)
	for(int ireco{0}; ireco<tsre_px_array.GetSize(); ireco++){
	  recotrk.SetXYZ(tsre_px_array[ireco], tsre_py_array[ireco], tsre_pz_array[ireco]);
	  
	  // Using track charge
	  if(!kUsePID){
	    // Positive => PROTON
	    if(re_charge_array[ireco] == 1){
	      P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,fMass_proton));
	      undoAfterburn(q_reco);
	      scatp4_rec.push_back(q_reco);
	    }
	  }
	  
	  // Not using track charge
	  else{
	    P3EVector q_reco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,tsre_mass_array[ireco]));
	    undoAfterburn(q_reco);
	    // Select on protons with eta cut
	    if(q_reco.Eta() > 4.2){ scatp4_rec.push_back(q_reco); }
	  }
	}// End of truth seeded charged particles

      }

      // Add in RP hits - only looking at protons
      // NO NEED TO UNDO AFTERBURNER FOR FF DETECTORS - NOT APPLIED IN FIRST PLACE
      for(int irpreco{0}; irpreco<rp_px_array.GetSize(); irpreco++){
 	recotrk.SetXYZ(rp_px_array[irpreco], rp_py_array[irpreco], rp_pz_array[irpreco]);	
	P3EVector q_rpreco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,rp_mass_array[irpreco]));
	if(rp_pdg_array[irpreco] == 2212){
	  scatp4_rom.push_back(q_rpreco);
	}
      }// End of RP reconstructed particles loop
      // NEED TO COMBINE B0 and RP tracks to avoid double counting
      vector<P3EVector> scatp4_all;   
      for(auto i:scatp4_rom) scatp4_all.push_back(i);
      for(auto j:scatp4_rec) scatp4_all.push_back(j);
            
      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Generated particles
      // Need Q2 for electron cuts
      if(scate4_gen.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_gen[0]);
      // Start with setting xB tail far away from expected point
      fxB_Tail = -10;
      // Fill histograms
      // Electron kinematics
      if(applyCuts_Electron(scate4_gen)){
	h_eta_MCe->Fill(scate4_gen[0].Eta());
	h_pt_MCe->Fill(scate4_gen[0].Pt());
	h_p_MCe->Fill(scate4_gen[0].P());
	h_theta_MCe->Fill(scate4_gen[0].Theta());
	h_phi_MCe->Fill(scate4_gen[0].Phi());
	h_E_MCe->Fill(scate4_gen[0].E());
      }
      // Photon kinematics
      if(applyCuts_Photon(scatg4_gen)){
	h_eta_MCg->Fill(scatg4_gen[0].Eta());
	h_pt_MCg->Fill(scatg4_gen[0].Pt());
	h_p_MCg->Fill(scatg4_gen[0].P());
	h_theta_MCg->Fill(scatg4_gen[0].Theta());
	h_phi_MCg->Fill(scatg4_gen[0].Phi());
	h_E_MCg->Fill(scatg4_gen[0].E());
      }
      // Proton kinematics
      if(applyCuts_Proton(scatp4_gen, "all")){
	h_eta_MCp->Fill(scatp4_gen[0].Eta());
	h_pt_MCp->Fill(scatp4_gen[0].Pt());
	h_p_MCp->Fill(scatp4_gen[0].P());
	h_theta_MCp->Fill(scatp4_gen[0].Theta()*1000);
	h_phi_MCp->Fill(scatp4_gen[0].Phi());
	h_E_MCp->Fill(scatp4_gen[0].E());
      }
      if(applyCuts_Electron(scate4_gen) && applyCuts_Photon(scatg4_gen)){
	// Electron-photon angles
	h_dphi_MCeg->Fill(scate4_gen[0].Phi()-scatg4_gen[0].Phi());
	h_dtheta_MCeg->Fill(scate4_gen[0].Theta()-scatg4_gen[0].Theta());
	
	// Electron-based event kinematics
	h_Q2_MC->Fill(fQ2);
	h_y_MC->Fill(calcY_Elec(beame4, scate4_gen[0], beamp4));
	fxB = calcX_Elec(beame4, scate4_gen[0], beamp4);
	h_xB_MC->Fill(TMath::Log10(fxB));
	h_QPQG_MC->Fill(calcPhiQPQG(beame4, beamp4, scate4_gen[0], scatg4_gen[0]));
	
	h_xBLin_MC->Fill(fxB);

	// Now add proton information
	if(applyCuts_Proton(scatp4_gen, "all")){
	  // And event DVCS cuts
	  // Need to calculate missing mass squared and t for DVCS cuts
	  fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]);
	  ft = calcT_BABE(beamp4, scatp4_gen[0]);
	  if(applyCuts_DVCS("all")){
	    h_t_MC->Fill(calcT_BABE(beamp4, scatp4_gen[0]));
	    // Need to calculate phi_h in the CoM frame
	    MomVector vCoMi = (beame4+beamp4).BoostToCM();
	    MomVector vCoMf = (scate4_gen[0]+scatp4_gen[0]+scatg4_gen[0]).BoostToCM();
	    // Change CoM definition: use qp/p-gamma system
	    // Try using virtual photon (q = k-k')
	    //MomVector vCoMgp = (beamp4+beame4-scate4_gen[0]).BoostToCM();
	    MomVector vCoMgp = (scatg4_gen[0]+scatp4_gen[0]).BoostToCM();
	    MomVector vTargetRest = beamp4.BoostToCM();
	    // Trento Phi - Target rest frame
	    //h_TPhi_MC->Fill(calcTrentoPhi_qp(boost(beame4, vTargetRest), boost(scate4_gen[0], vTargetRest), boost(scatp4_gen[0], vTargetRest)));
	    //h_TPhi_pg_MC->Fill(calcTrentoPhi_pg(boost(beame4, vTargetRest), boost(scate4_gen[0], vTargetRest), boost(scatp4_gen[0], vTargetRest), boost(scatg4_gen[0], vTargetRest)));
	    //h_TPhi_qg_MC->Fill(calcTrentoPhi_qg(boost(beame4, vTargetRest), boost(scate4_gen[0], vTargetRest), boost(scatg4_gen[0], vTargetRest)));
	    // Trento Phi - Centre-of-mass frame
	    h_TPhi_MC->Fill(calcTrentoPhi_qp(boost(beame4, vCoMgp), boost(scate4_gen[0], vCoMgp), boost(scatp4_gen[0], vCoMgp)));
	    h_TPhi_pg_MC->Fill(calcTrentoPhi_pg(boost(beame4, vCoMgp), boost(scate4_gen[0], vCoMgp), boost(scatp4_gen[0], vCoMgp), boost(scatg4_gen[0], vCoMgp)));
	    h_TPhi_qg_MC->Fill(calcTrentoPhi_qg(boost(beame4, vCoMgp), boost(scate4_gen[0], vCoMgp), boost(scatg4_gen[0], vCoMgp)));
	    h_Cone_MC->Fill(calcConeAngle(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
	  }
	}
      } // FILLED MC GENERATED HISTOS

      // Decide minimum xB value from MC generated histogram
      fxB_Tail = h_xB_MC->GetBinLowEdge( h_xB_MC->FindFirstBinAbove(0) );
      // Associated MC
      // Need Q2 for electron cuts
      if(scate4_aso.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_aso[0]);
      // Fill histograms
      // Electron kinematics
      if(applyCuts_Electron(scate4_aso)){
	h_eta_MCAe->Fill(scate4_aso[0].Eta());
	h_pt_MCAe->Fill(scate4_aso[0].Pt());
	h_p_MCAe->Fill(scate4_aso[0].P());
	h_theta_MCAe->Fill(scate4_aso[0].Theta());
	h_phi_MCAe->Fill(scate4_aso[0].Phi());
	h_E_MCAe->Fill(scate4_aso[0].E());
      }
      // Photon kinematics
      if(applyCuts_Photon(scatg4_aso)){
	h_eta_MCAg->Fill(scatg4_aso[0].Eta());
	h_pt_MCAg->Fill(scatg4_aso[0].Pt());
	h_p_MCAg->Fill(scatg4_aso[0].P());
	h_theta_MCAg->Fill(scatg4_aso[0].Theta());
	h_phi_MCAg->Fill(scatg4_aso[0].Phi());	
	h_E_MCAg->Fill(scatg4_aso[0].E());
      }
      // Proton kinematics
      if(applyCuts_Proton(scatp4_aso, "all")){
	h_eta_MCAp->Fill(scatp4_aso[0].Eta());
	h_pt_MCAp->Fill(scatp4_aso[0].Pt());
	h_p_MCAp->Fill(scatp4_aso[0].P());
	h_theta_MCAp->Fill(scatp4_aso[0].Theta()*1000);
	h_phi_MCAp->Fill(scatp4_aso[0].Phi());
	h_E_MCAp->Fill(scatp4_aso[0].E());
      }
      if(applyCuts_Electron(scate4_aso) && applyCuts_Photon(scatg4_aso)){
	// Electron-photon angles
	h_dphi_MCAeg->Fill(scate4_aso[0].Phi()-scatg4_aso[0].Phi());
	h_dtheta_MCAeg->Fill(scate4_aso[0].Theta()-scatg4_aso[0].Theta());

	// Electron-based event kinematics
	h_Q2_MCA->Fill(fQ2);
	h_y_MCA->Fill(calcY_Elec(beame4, scate4_aso[0], beamp4));
	fxB = calcX_Elec(beame4, scate4_aso[0], beamp4);
	h_xB_MCA->Fill(TMath::Log10(fxB));
	h_QPQG_MCA->Fill(calcPhiQPQG(beame4, beamp4, scate4_aso[0], scatg4_aso[0]));

	// Now add proton information
	if(applyCuts_Proton(scatp4_aso, "all")){
	  // And event DVCS cuts
	  // Need to calculate missing mass squared and t for DVCS cuts
	  fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]);
	  ft = calcT_BABE(beamp4, scatp4_aso[0]);
	  if(applyCuts_DVCS("all")){
	    h_t_MCA->Fill(calcT_BABE(beamp4, scatp4_aso[0]));
	    MomVector vCoMi = (beame4+beamp4).BoostToCM();
	    MomVector vCoMf = (scate4_aso[0]+scatp4_aso[0]+scatg4_aso[0]).BoostToCM();
	    MomVector vTargetRest = beamp4.BoostToCM();
	    // Change CoM definition: use qp/p-gamma system
	    MomVector vCoMgp = (beamp4+beame4-scate4_gen[0]).BoostToCM();
	    //h_TPhi_MCA->Fill(calcTrentoPhi_qp(boost(beame4, vCoMgp), boost(scate4_aso[0], vCoMgp), boost(scatp4_aso[0], vCoMgp)));
	    //h_TPhi_MCA->Fill(calcTrentoPhi(beame4, scate4_aso[0], scatp4_aso[0]));
	    h_TPhi_MCA->Fill(calcTrentoPhi_qp(boost(beame4, vTargetRest), boost(scate4_aso[0], vTargetRest), boost(scatp4_aso[0], vTargetRest)));
	    h_Cone_MCA->Fill(calcConeAngle(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	  }
	}
      } // FILLED ASSOCIATED MC HISTOS

      // Reconstructed particles - Barrel/B0
      // Need Q2 for electron cuts
      if(scate4_rec.size() == 0) fQ2 = 0;
      else fQ2 = calcQ2_Elec(beame4, scate4_rec[0]);
      // Fill histograms
      // Electron kinematics
      if(applyCuts_Electron(scate4_rec) && applyCuts_Photon(scatg4_rec)){
	h_eta_RPe->Fill(scate4_rec[0].Eta());
	h_pt_RPe->Fill(scate4_rec[0].Pt());
	h_p_RPe->Fill(scate4_rec[0].P());
	h_theta_RPe->Fill(scate4_rec[0].Theta());
	h_phi_RPe->Fill(scate4_rec[0].Phi());
	h_E_RPe->Fill(scate4_rec[0].E());

	
	h_2D_EvEta_e->Fill(scate4_rec[0].Eta(), scate4_rec[0].E());
      }
      // Photon kinematics
      if(applyCuts_Electron(scate4_rec) && applyCuts_Photon(scatg4_rec)){
	h_eta_RPg->Fill(scatg4_rec[0].Eta());
	h_pt_RPg->Fill(scatg4_rec[0].Pt());
	h_p_RPg->Fill(scatg4_rec[0].P());
	h_theta_RPg->Fill(scatg4_rec[0].Theta());
	h_phi_RPg->Fill(scatg4_rec[0].Phi());
	h_E_RPg->Fill(scatg4_rec[0].E());

	h_2D_EvEta_g->Fill(scatg4_rec[0].Eta(), scatg4_rec[0].E());
      }
      // Proton kinematics - B0
      if(applyCuts_Proton(scatp4_rec, "B0")){
	h_eta_RPp->Fill(scatp4_rec[0].Eta());
	h_pt_RPp->Fill(scatp4_rec[0].Pt());
	h_p_RPp->Fill(scatp4_rec[0].P());
	h_theta_RPp->Fill(scatp4_rec[0].Theta()*1000);
	h_phi_RPp->Fill(scatp4_rec[0].Phi());
	h_E_RPp->Fill(scatp4_rec[0].E());

	h_2D_EvEta_p->Fill(scatp4_rec[0].Eta(), scatp4_rec[0].E());
      }
      // Proton kinematics - RP
      if(applyCuts_Proton(scatp4_rom, "RP")){
	h_eta_RPPp->Fill(scatp4_rom[0].Eta());
	h_pt_RPPp->Fill(scatp4_rom[0].Pt());
	h_p_RPPp->Fill(scatp4_rom[0].P());
	h_theta_RPPp->Fill(scatp4_rom[0].Theta()*1000);
	h_phi_RPPp->Fill(scatp4_rom[0].Phi());
	h_E_RPPp->Fill(scatp4_rom[0].E());

	h_px_RPP->Fill(scatp4_rom[0].Px());
	h_py_RPP->Fill(scatp4_rom[0].Py());
	h_pz_RPP->Fill(scatp4_rom[0].Pz());

	h_2D_EvEta_p->Fill(scatp4_rom[0].Eta(), scatp4_rom[0].E());
      }
      if(applyCuts_Electron(scate4_rec) && applyCuts_Photon(scatg4_rec)){
	// Electron-photon angles
	h_dphi_RPeg->Fill(scate4_rec[0].Phi()-scatg4_rec[0].Phi());
	h_dtheta_RPeg->Fill(scate4_rec[0].Theta()-scatg4_rec[0].Theta());
	
	// Electron-based event kinematics
	h_Q2_RP->Fill(fQ2);
	h_y_RP->Fill(calcY_Elec(beame4, scate4_rec[0], beamp4));
	fxB = calcX_Elec(beame4, scate4_rec[0], beamp4);
	h_xB_RP->Fill(TMath::Log10(fxB));
	h_QPQG_RP->Fill(calcPhiQPQG(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
	
	// Now add proton information - B0
	if(applyCuts_Proton(scatp4_all, "B0")){
	  // And event DVCS cuts
	  // Need to calculate missing mass squared and t for DVCS cuts
	  fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_all[0], scatg4_rec[0]);
	  ft = calcT_BABE(beamp4, scatp4_all[0]);
	  
	  // Add proton information - RP
	  if(applyCuts_Proton(scatp4_all, "RP")){
	    // And event DVCS cuts
	    // Need to calculate missing mass squared and t for DVCS cuts
	    fM2miss = calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_all[0], scatg4_rec[0]);
	    ft = calcT_BABE(beamp4, scatp4_all[0]);
	    if(applyCuts_DVCS("RP")){
	      h_t_RPP->Fill(calcT_BABE(beamp4, scatp4_all[0]));
	      MomVector vCoMi = (beame4+beamp4).BoostToCM();
	      MomVector vCoMf = (scate4_rec[0]+scatp4_all[0]+scatg4_rec[0]).BoostToCM();
	      MomVector vTargetRest = beamp4.BoostToCM();
	      // Change CoM definition: use qp/p-gamma system
	      MomVector vCoMgp = (beamp4+beame4-scate4_gen[0]).BoostToCM();
	      // Trento Phi - target rest frame
	      /*h_TPhi_RPP->Fill(calcTrentoPhi_qp(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatp4_all[0], vTargetRest)));
		h_TPhi_pg_RPP->Fill(calcTrentoPhi_pg(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatp4_all[0], vTargetRest), boost(scatg4_rec[0], vTargetRest)));	   
		h_TPhi_qg_RPP->Fill(calcTrentoPhi_qg(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatg4_rec[0], vTargetRest)));*/
	      h_TPhi_RPP->Fill(calcTrentoPhi_qp(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatp4_all[0], vCoMgp)));
	      h_TPhi_pg_RPP->Fill(calcTrentoPhi_pg(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatp4_all[0], vCoMgp), boost(scatg4_rec[0], vCoMgp))); 
	      h_TPhi_qg_RPP->Fill(calcTrentoPhi_qg(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatg4_rec[0], vCoMgp)));
	      h_Cone_RPP->Fill(calcConeAngle(beame4, beamp4, scate4_rec[0], scatp4_all[0], scatg4_rec[0]));
	    }
	  }
	  
	  // Only move onto B0 if RP tracks are empty
	  if(applyCuts_DVCS("B0")){
	    h_t_RP->Fill(calcT_BABE(beamp4, scatp4_all[0]));
	    MomVector vCoMi = (beame4+beamp4).BoostToCM();
	    MomVector vCoMf = (scate4_rec[0]+scatp4_all[0]+scatg4_rec[0]).BoostToCM();
	    MomVector vTargetRest = beamp4.BoostToCM();
	    // Change CoM definition: use qp/p-gamma system
	    MomVector vCoMgp = (beamp4+beame4-scate4_gen[0]).BoostToCM();
	    //h_TPhi_RP->Fill(calcTrentoPhi_qp(boost(beame4, vCoMi), boost(scate4_rec[0], vCoMf), boost(scatp4_rec[0], vCoMf)));
	    // Trento Phi - target rest frame
	    /*h_TPhi_RP->Fill(calcTrentoPhi_qp(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatp4_rec[0], vTargetRest)));
	    h_TPhi_pg_RP->Fill(calcTrentoPhi_pg(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatp4_rec[0], vTargetRest), boost(scatg4_rec[0], vTargetRest)));	   
	    h_TPhi_qg_RP->Fill(calcTrentoPhi_qg(boost(beame4, vTargetRest), boost(scate4_rec[0], vTargetRest), boost(scatg4_rec[0], vTargetRest)));*/
	    // Trento Phi - centre-of-mass frame
	    h_TPhi_RP->Fill(calcTrentoPhi_qp(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatp4_all[0], vCoMgp)));
	    h_TPhi_pg_RP->Fill(calcTrentoPhi_pg(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatp4_all[0], vCoMgp), boost(scatg4_rec[0], vCoMgp)));	   
	    h_TPhi_qg_RP->Fill(calcTrentoPhi_qg(boost(beame4, vCoMgp), boost(scate4_rec[0], vCoMgp), boost(scatg4_rec[0], vCoMgp)));
	    h_Cone_RP->Fill(calcConeAngle(beame4, beamp4, scate4_rec[0], scatp4_all[0], scatg4_rec[0]));
	  }
	}
      } // FILLED RECONSTRUCTED HISTOS

      // 2D histograms for proton momenta
      if(scatp4_aso.size()==1 && scatp4_rec.size()==1){
	h_p_2d->Fill(scatp4_aso[0].P(), scatp4_rec[0].P());
	h_pt_2d->Fill(scatp4_aso[0].Pt(), scatp4_rec[0].Pt());
	h_t_2d->Fill(calcT_BABE(beamp4, scatp4_aso[0]), calcT_BABE(beamp4, scatp4_rec[0]));
      }

      // Fill 2D t-resolution
      // ASSUME THAT GENERATED POSITIVE TRACK MATCHES RECONSTRUCTED POSITIVE TRACK
      Float_t t_rec{0}, t_gen{0};
      if(applyCuts_Proton(scatp4_all, "B0") && applyCuts_Proton(scatp4_gen, "all")){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_all[0]);
	h_tRes_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
      }
      if(applyCuts_Proton(scatp4_all, "RP") && applyCuts_Proton(scatp4_gen, "all")){
	t_gen = calcT_BABE(beamp4, scatp4_gen[0]);
	t_rec = calcT_BABE(beamp4, scatp4_all[0]);
	h_tRes_2d->Fill(t_gen, TMath::Abs(t_rec-t_gen));
      }

      // Standalone histograms - B0 occupancies
      double hit_x = -9999.;
      double hit_y = -9999.;
      double hit_z = -9999.;
      double hit_deposited_energy = -9999.;
      // B0 ECAL
      for(int b0cluster = 0; b0cluster < b0_cluster_x.GetSize(); b0cluster++){	
	hit_x = b0_cluster_x[b0cluster];
	hit_y = b0_cluster_y[b0cluster];
	hit_z = b0_cluster_z[b0cluster];
	hit_deposited_energy = b0_cluster_energy[b0cluster]*1.246; //poor man's calibration constant, for now
							
	h_B0_emcal_occupancy_map->Fill(hit_x, hit_y);
	h_B0_emcal_cluster_energy->Fill(hit_deposited_energy);
      }
      // B0 tracker layers
      for(int b0hit = 0; b0hit < b0_hits_x.GetSize(); b0hit++){
	hit_x = b0_hits_x[b0hit];
	hit_y = b0_hits_y[b0hit];
	hit_z = b0_hits_z[b0hit];
	hit_deposited_energy = b0_hits_eDep[b0hit]*1e6; //convert GeV --> keV
	h_B0_hit_energy_deposit->Fill(hit_deposited_energy);
				
	if(hit_deposited_energy < 10.0){ continue; } //threshold value -- 10 keV, arbitrary for now
			
	//ACLGAD layout
	if(hit_z > 5700 && hit_z < 5990){ h_B0_occupancy_map_layer_0->Fill(hit_x, hit_y); }
	if(hit_z > 6100 && hit_z < 6200){ h_B0_occupancy_map_layer_1->Fill(hit_x, hit_y); }
	if(hit_z > 6400 && hit_z < 6500){ h_B0_occupancy_map_layer_2->Fill(hit_x, hit_y); }
	if(hit_z > 6700 && hit_z < 6750){ h_B0_occupancy_map_layer_3->Fill(hit_x, hit_y); }
      }

      // Standalone histograms - p/pT resolutions (ONLY FOR B0)
      // NOTE: Can ONLY match tracks to associated MC -> IF MCA ARRAY IS SMALLER THAN RECO ARRAY, CANNOT USE
      // Force at least 1 proton in MCA and reco arrays
      if( scatp4_aso.size()!=0 && scatp4_rec.size()!=0 && scatp4_aso.size()==scatp4_rec.size() ){
	for(int i{0}; i<scatp4_rec.size(); i++){
	  double delPt = scatp4_rec[i].Pt() - scatp4_aso[i].Pt();
	  double delP = scatp4_rec[i].P() - scatp4_aso[i].P();
	  
	  //if(scatp4_aso[i].P() >= 80 && scatp4_aso[i].P() <= 100){
	  h_b0_pt_resolution->Fill(delPt/scatp4_aso[i].Pt());
	  h_b0_pt_resolution_percent->Fill(scatp4_aso[i].Pt(), delPt/scatp4_aso[i].Pt());
	  h_b0_p_resolution_percent->Fill(scatp4_aso[i].P(), delP/scatp4_aso[i].P());
	  //}
	}
      }

      // Standalone histograms - Roman Pots occupancies
      for(int iRPPart = 0; iRPPart < global_hit_RP_x.GetSize(); iRPPart++){
    	h_rp_occupancy_map->Fill(global_hit_RP_x[iRPPart], global_hit_RP_y[iRPPart]);
	h_rp_z->Fill(global_hit_RP_z[iRPPart]);
      }

      //------------------------------------------------------------
      // Standalone histograms -  Exclusivity variables
      //------------------------------------------------------------
      // 1. FULL FINAL STATE
      // Generated particles
      if( !applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all") ){ }
      else{
	h_Emiss3_MC->Fill(calcEMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
	h_Pmiss3_MC->Fill(calcPMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
	h_Ptmiss3_MC->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
	h_M2miss3_MC->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0], scatg4_gen[0]));
      }
      //  Associated MC
      if( !applyCuts_All(beame4, beamp4, scate4_aso, scatp4_aso, scatg4_aso, "all") ){ }
      else{
	h_Emiss3_MCA->Fill(calcEMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	h_Pmiss3_MCA->Fill(calcPMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	h_Ptmiss3_MCA->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
	h_M2miss3_MCA->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0], scatg4_aso[0]));
      }
      //  Reconstructed - B0
      if( !applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0") ){ }
      else{
	h_Emiss3_RP->Fill(calcEMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
	h_Pmiss3_RP->Fill(calcPMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
	h_Ptmiss3_RP->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
	h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0], scatg4_rec[0]));
      }
      //  Reconstructed - Roman Pots
      if( !applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP") ){ }
      else{
	h_Emiss3_RP->Fill(calcEMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
	h_Pmiss3_RP->Fill(calcPMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
	h_Ptmiss3_RP->Fill(calcPtMiss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
	h_M2miss3_RP->Fill(calcM2Miss_3Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0], scatg4_rec[0]));
      }

      // 2. PARTIAL FINAL STATE: e'p'
      // Generated MC
      if(scate4_gen.size()==1 && scatp4_gen.size()==1 && scate4_gen[0].P()<=fPMax_e && scatp4_gen[0].P()<=fPMax_p && scatp4_gen[0].Theta()<0.02){
	h_M2miss2ep_MC->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_gen[0], scatp4_gen[0]));
      }
      // Associated MC
      if(scate4_aso.size()==1 && scatp4_aso.size()==1 && scate4_aso[0].P()<=fPMax_e && scatp4_aso[0].P()<=fPMax_p && scatp4_aso[0].Theta()<0.02){
	h_M2miss2ep_MCA->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_aso[0], scatp4_aso[0]));
      }
      // Reconstructed (B0 protons)
      if(scate4_rec.size()==1 && scatp4_rec.size()==1 && scate4_rec[0].P()<=fPMax_e && scatp4_rec[0].P()<=fPMax_p && scatp4_rec[0].Theta()<0.02 && scatp4_gen[0].Theta()>0.0055){
	h_M2miss2ep_RP->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_rec[0], scatp4_rec[0]));
      }
      // Reconstructed (Roman Pot protons)
      if(scate4_rec.size()==1 && scatp4_rom.size()==1 && scate4_rec[0].P()<=fPMax_e && scatp4_rom[0].P()<=fPMax_p && scatp4_rom[0].Theta()<0.0055){
	h_M2miss2ep_RP->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_rec[0], scatp4_rom[0]));
      }
      // 3. PARTIAL FINAL STATE: e'gamma
      // Generated MC
      if(scate4_gen.size()==1 && scatg4_gen.size()==1 && scate4_gen[0].P()<=fPMax_e){
	h_M2miss2eg_MC->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_gen[0], scatg4_gen[0]));
      }
      // Associated MC
      if(scate4_aso.size()==1 && scatg4_aso.size()==1 && scate4_aso[0].P()<=fPMax_e){
	h_M2miss2eg_MCA->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_aso[0], scatg4_aso[0]));
      }
      // Reconstructed
      if(scate4_rec.size()==1 && scatg4_rec.size()==1 && scate4_rec[0].P()<=fPMax_e){
	h_M2miss2eg_RP->Fill(calcM2Miss_2Body(beame4, beamp4, scate4_rec[0], scatg4_rec[0]));
      }

      //------------------------------------------------------------
      // Standalone histograms -  Photon resolutions
      //------------------------------------------------------------
      if(applyCuts_Photon(scatg4_rec) && applyCuts_Photon(scatg4_aso)){
	h_PhotRes_E->Fill(scatg4_rec[0].E() - scatg4_aso[0].E());
	h_PhotRes_theta->Fill(scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
	h_PhotRes2D_theta->Fill(scatg4_aso[0].Theta(), scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
	h_PhotRes2D_thetaReco->Fill(scatg4_rec[0].Theta(), scatg4_rec[0].Theta() - scatg4_aso[0].Theta());
      }
      
      //------------------------------------------------------------
      // Standalone histograms -  2D distributions (t vs xB)
      //------------------------------------------------------------
      if(applyCuts_All(beame4, beamp4, scate4_gen, scatp4_gen, scatg4_gen, "all")){
	h_2D_xVt_MC->Fill(TMath::Log10(calcX_Elec(beame4, scate4_gen[0], beamp4)), calcT_BABE(beamp4, scatp4_gen[0]));
	h_2D_tVQ2_MC->Fill(calcT_BABE(beamp4, scatp4_gen[0]), calcQ2_Elec(beame4, scate4_gen[0]));
      }
      if(applyCuts_All(beame4, beamp4, scate4_aso, scatp4_aso, scatg4_aso, "all")){
	h_2D_xVt_MCA->Fill(TMath::Log10(calcX_Elec(beame4, scate4_aso[0], beamp4)), calcT_BABE(beamp4, scatp4_aso[0]));
	h_2D_tVQ2_MCA->Fill(calcT_BABE(beamp4, scatp4_aso[0]), calcQ2_Elec(beame4, scate4_aso[0]));
      }
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rec, scatg4_rec, "B0")){
	h_2D_xVt_RP->Fill(TMath::Log10(calcX_Elec(beame4, scate4_rec[0], beamp4)), calcT_BABE(beamp4, scatp4_rec[0]));
	h_2D_tVQ2_RP->Fill(calcT_BABE(beamp4, scatp4_rec[0]), calcQ2_Elec(beame4, scate4_rec[0]));
      }
      if(applyCuts_All(beame4, beamp4, scate4_rec, scatp4_rom, scatg4_rec, "RP")){
	h_2D_xVt_RP->Fill(TMath::Log10(calcX_Elec(beame4, scate4_rec[0], beamp4)), calcT_BABE(beamp4, scatp4_rom[0]));
	h_2D_tVQ2_RP->Fill(calcT_BABE(beamp4, scatp4_rom[0]), calcQ2_Elec(beame4, scate4_rec[0]));
      }

      //------------------------------------------------------------
      // Standalone histograms -  2D distributions (Q2 vs xB)
      //------------------------------------------------------------
      if(applyCuts_Electron(scate4_gen)){
	h_2D_xVQ2_MC->Fill(TMath::Log10(calcX_Elec(beame4, scate4_gen[0], beamp4)), calcQ2_Elec(beame4, scate4_gen[0]));
      }
      if(applyCuts_Electron(scate4_aso)){
	h_2D_xVQ2_MCA->Fill(TMath::Log10(calcX_Elec(beame4, scate4_aso[0], beamp4)), calcQ2_Elec(beame4, scate4_aso[0]));
      }
      if(applyCuts_Electron(scate4_rec)){
	h_2D_xVQ2_RP->Fill(TMath::Log10(calcX_Elec(beame4, scate4_rec[0], beamp4)), calcQ2_Elec(beame4, scate4_rec[0]));
      }

      //------------------------------------------------------------
      // Standalone histograms -  Inclusive t-distribution
      //------------------------------------------------------------
      if(applyCuts_Proton(scatp4_gen, "all")) h_tInc_MC->Fill(calcT_BABE(beamp4,scatp4_gen[0]));
      if(applyCuts_Proton(scatp4_rec,"B0"))   h_tInc_RP->Fill(calcT_BABE(beamp4,scatp4_rec[0]));
      if(applyCuts_Proton(scatp4_rom,"RP") && calcT_BABE(beamp4,scatp4_rom[0])<fMaxt_RP)  h_tInc_RPP->Fill(calcT_BABE(beamp4,scatp4_rom[0]));

    } // END EVENT/TTREEREADER LOOP

    inputRootFile->Close();
  } // END OF FILE LOOP

  //------------------------------------------------------------
  // Calculate B0 resolutions
  //------------------------------------------------------------
  // pT resolution
  h_b0_extracted_pt_resolution = extractResolution("b0_extracted_pt_resolution", h_b0_pt_resolution_percent);
  h_b0_extracted_pt_resolution->GetXaxis()->SetTitle("p_{T, MC} [GeV/c]");
  h_b0_extracted_pt_resolution->GetYaxis()->SetTitle("#Delta p_{T}/p_{T, MC}");  

  // p resolution
  h_b0_extracted_p_resolution = extractResolution("b0_extracted_p_resolution", h_b0_p_resolution_percent);
  h_b0_extracted_p_resolution->GetXaxis()->SetTitle("p_{MC} [GeV/c]");
  h_b0_extracted_p_resolution->GetYaxis()->SetTitle("#Delta p/p_{MC}");

  // t resolution (also includes RP)
  h_extracted_t_resolution = extractResolution("extracted_t_resolution",h_tRes_2d);
  h_extracted_t_resolution->GetXaxis()->SetTitle("|t|_{MC} [GeV^{2}]");
  h_extracted_t_resolution->GetYaxis()->SetTitle("#delta t [GeV^{2}]");
  
  // Check that errors aren't too large - if are, set bin content to -1
  for(int binp{1}; binp<h_b0_extracted_p_resolution->GetNbinsX(); binp++){
    if(h_b0_extracted_p_resolution->GetBinError(binp) > 0.02){
      h_b0_extracted_p_resolution->SetBinContent(binp,-1);
      h_b0_extracted_p_resolution->SetBinError(binp,0);
    }
  }
  for(int binpt{1}; binpt<h_b0_extracted_pt_resolution->GetNbinsX(); binpt++){
    if(h_b0_extracted_pt_resolution->GetBinError(binpt) > 0.02){
      h_b0_extracted_pt_resolution->SetBinContent(binpt,-1);
      h_b0_extracted_pt_resolution->SetBinError(binpt,0);
    }
  }


  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  fOutFile->cd();
  
  // MC generated
  h_eta_MCp->Write();
  h_pt_MCp->Write();
  h_p_MCp->Write();
  h_theta_MCp->Write();
  h_phi_MCp->Write();
  h_E_MCp->Write();
  h_eta_MCe->Write();
  h_pt_MCe->Write();
  h_p_MCe->Write();
  h_theta_MCe->Write();
  h_phi_MCe->Write();
  h_E_MCe->Write();
  h_eta_MCg->Write();
  h_pt_MCg->Write();
  h_p_MCg->Write();
  h_theta_MCg->Write();
  h_phi_MCg->Write();
  h_E_MCg->Write();
  h_t_MC->Write();
  h_Q2_MC->Write();
  h_xB_MC->Write();
  h_y_MC->Write();
  h_TPhi_MC->Write();
  h_QPQG_MC->Write();
  h_Cone_MC->Write();
  h_dphi_MCeg->Write();
  h_dtheta_MCeg->Write();
  // MC associated
  h_eta_MCAp->Write();
  h_pt_MCAp->Write();
  h_p_MCAp->Write();
  h_theta_MCAp->Write();
  h_phi_MCAp->Write();
  h_E_MCAp->Write();
  h_eta_MCAe->Write();
  h_pt_MCAe->Write();
  h_p_MCAe->Write();
  h_theta_MCAe->Write();
  h_phi_MCAe->Write();
  h_E_MCAe->Write();
  h_eta_MCAg->Write();
  h_pt_MCAg->Write();
  h_p_MCAg->Write();
  h_theta_MCAg->Write();
  h_phi_MCAg->Write();
  h_E_MCAg->Write();
  h_t_MCA->Write();
  h_Q2_MCA->Write();
  h_xB_MCA->Write();
  h_y_MCA->Write();
  h_TPhi_MCA->Write();
  h_QPQG_MCA->Write();
  h_Cone_MCA->Write();
  h_dphi_MCAeg->Write();
  h_dtheta_MCAeg->Write();
  // Reconstructed (Barrel/B0)
  h_eta_RPp->Write();
  h_pt_RPp->Write();
  h_p_RPp->Write();
  h_theta_RPp->Write();
  h_phi_RPp->Write();
  h_E_RPp->Write();
  h_eta_RPe->Write();
  h_pt_RPe->Write();
  h_p_RPe->Write();
  h_theta_RPe->Write();
  h_phi_RPe->Write();
  h_E_RPe->Write();
  h_eta_RPg->Write();
  h_pt_RPg->Write();
  h_p_RPg->Write();
  h_theta_RPg->Write();
  h_phi_RPg->Write();
  h_E_RPg->Write();
  h_t_RP->Write();
  h_Q2_RP->Write();
  h_xB_RP->Write();
  h_y_RP->Write();
  h_TPhi_RP->Write();
  h_QPQG_RP->Write();
  h_Cone_RP->Write();
  h_dphi_RPeg->Write();
  h_dtheta_RPeg->Write();
  // Reconstructed (Roman Pots)
  h_eta_RPPp->Write();
  h_pt_RPPp->Write();
  h_p_RPPp->Write();
  h_theta_RPPp->Write();
  h_phi_RPPp->Write();
  h_E_RPPp->Write();
  h_t_RPP->Write();
  h_TPhi_RPP->Write();
  h_Cone_RPP->Write();
  // 2D distributions
  h_p_2d->Write();
  h_pt_2d->Write();
  h_t_2d->Write();
  h_tRes_2d->Write();
  h_extracted_t_resolution->Write();
  // B0 occupancies (tracker and ECAL)
  h_B0_occupancy_map_layer_0->Write();
  h_B0_occupancy_map_layer_1->Write();
  h_B0_occupancy_map_layer_2->Write();
  h_B0_occupancy_map_layer_3->Write();
  h_B0_hit_energy_deposit->Write();
  h_B0_emcal_occupancy_map->Write();
  h_B0_emcal_cluster_energy->Write();
  // Roman Pots occupancies
  h_rp_occupancy_map->Write();
  h_rp_z->Write();
  // B0 momentum resolution
  h_b0_pt_resolution->Write();
  h_b0_pt_resolution_percent->Write();
  h_b0_p_resolution_percent->Write();
  h_b0_extracted_pt_resolution->Write();
  h_b0_extracted_p_resolution->Write();
  // Exclusivity variables
  h_Emiss3_MC->Write();
  h_Pmiss3_MC->Write();
  h_Ptmiss3_MC->Write();
  h_M2miss3_MC->Write();
  h_Emiss3_MCA->Write();
  h_Pmiss3_MCA->Write();
  h_Ptmiss3_MCA->Write();
  h_M2miss3_MCA->Write();
  h_Emiss3_RP->Write();
  h_Pmiss3_RP->Write();
  h_Ptmiss3_RP->Write();
  h_M2miss3_RP->Write();
  h_M2miss2ep_MC->Write();
  h_M2miss2ep_MCA->Write();
  h_M2miss2ep_RP->Write();
  h_M2miss2eg_MC->Write();
  h_M2miss2eg_MCA->Write();
  h_M2miss2eg_RP->Write();
  // Photon resolutions
  h_PhotRes_E->Write();
  h_PhotRes_theta->Write();
  h_PhotRes2D_theta->Write();
  h_PhotRes2D_thetaReco->Write();
  // 2D distributions
  h_2D_xVt_MC->Write();
  h_2D_xVt_MCA->Write();
  h_2D_xVt_RP->Write();
  h_2D_xVQ2_MC->Write();
  h_2D_xVQ2_MCA->Write();
  h_2D_xVQ2_RP->Write();
  h_2D_tVQ2_MC->Write();
  h_2D_tVQ2_MCA->Write();
  h_2D_tVQ2_RP->Write();
  h_2D_EvEta_e->Write();
  h_2D_EvEta_p->Write();
  h_2D_EvEta_g->Write();
  // Inclusive t-dist
  h_tInc_MC->Write();
  h_tInc_RP->Write();
  h_tInc_RPP->Write();
  // Other Trento phi calculations
  h_TPhi_pg_MC->Write();
  h_TPhi_qg_MC->Write();
  h_TPhi_pg_RP->Write();
  h_TPhi_qg_RP->Write();
  h_TPhi_pg_RPP->Write();
  h_TPhi_qg_RPP->Write();
  // OTHER MISC.
  h_xBLin_MC->Write();
  h_px_RPP->Write();
  h_py_RPP->Write();
  h_pz_RPP->Write();

  fOutFile->Close();

  return;
}
