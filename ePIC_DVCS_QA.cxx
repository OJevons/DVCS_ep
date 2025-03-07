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
//                     DO ANALYSIS
//----------------------------------------------------
//----------------------------------------------------

void ePIC_DVCS_TASK::doAnalysis(){

  //---------------------------------------------------------
  // Setup: Load input file list
  //---------------------------------------------------------
  ifstream fileListStream;
  fileListStream.open(sInList);
  string fileName;
  TFile* inputRootFile;

  //---------------------------------------------------------
  // Setup: Declare QA histograms
  //---------------------------------------------------------
  // Inclusive kinematic quantities
  // Q2
  TH1D* h_Q2_Truth = new TH1D("q2_truth",";Q^{2}(MC) [(GeV/c^{2})^{2}]"      , 500, 0., 10.);
  TH1D* h_Q2_Acc   = new TH1D("q2_acc"  ,";Q^{2}(MC|Reco.) [(GeV/c^{2})^{2}]", 500, 0., 10.);
  TH1D* h_Q2_Reco  = new TH1D("q2_reco" ,";Q^{2}(Reco.) [(GeV/c^{2})^{2}]"   , 500, 0., 10.);
  TH2D* h_Q2_Resp  = new TH2D("q2_resp" ,";Q^{2}(MC|Reco.) [(GeV/c^{2})^{2}];Q^{2}(Reco.) [(GeV/c^{2})^{2}]", 500, 0., 10., 500, 0., 10.);
  TH1D* h_Q2_Pur   = new TH1D("q2_pur"  ,";Q^{2}(Reco.) [(GeV/c^{2})^{2}];#frac{MC|Reco.}{Reco.}", 500, 0., 10.);
  TH2D* h_dQ2vQ2   = new TH2D("dq2vq2"  ,";Q^{2}(MC|Reco.) [(GeV/c^{2})^{2}];#Delta Q^{2} [(GeV/c^{2})^{2}]", 500, 0., 10., 200, -10., 10.);
  // x
  TH1D* h_xB_Truth = new TH1D("xb_truth",";x_{B}(MC)"      , 1e4, 0., 1.);
  TH1D* h_xB_Acc   = new TH1D("xb_acc"  ,";x_{B}(MC|Reco.)", 1e4, 0., 1.);
  TH1D* h_xB_Reco  = new TH1D("xb_reco" ,";x_{B}(Reco.)"   , 1e4, 0., 1.);
  TH2D* h_xB_Resp  = new TH2D("xb_resp" ,";x_{B}(MC|Reco.);x_{B}(Reco.)", 1e4, 0., 1., 1e4, 0., 1.);
  TH1D* h_xB_Pur   = new TH1D("xb_pur"  ,";x_{B}(Reco.);#frac{MC|Reco.}{Reco.}", 1e4, 0., 1.);
  TH2D* h_dxBvxB   = new TH2D("dxbvxb"  ,";x_{B}(MC|Reco.);#Delta x_{B}",1e4, 0., 1., 200, -1., 1.);
  // y
  TH1D* h_y_Truth = new TH1D("y_truth",";y(MC)"      , 100, 0., 1.);
  TH1D* h_y_Acc   = new TH1D("y_acc"  ,";y(MC|Reco.)", 100, 0., 1.);
  TH1D* h_y_Reco  = new TH1D("y_reco" ,";y(Reco.)"   , 100, 0., 1.);
  TH2D* h_y_Resp  = new TH2D("y_resp" ,";y(MC|Reco.);y(Reco.)", 100, 0., 1., 100, 0., 1.); 
  TH1D* h_y_Pur   = new TH1D("y_pur"  ,";y(Reco.);#frac{MC|Reco.}{Reco.}", 100, 0., 1.);
  TH2D* h_dyvy    = new TH2D("dyvy"   ,";y(MC|Reco.);#Delta y",100, 0., 1., 200, -1., 1.);
  // Exclusive kinematic quantities
  TH1D* h_t_Truth  = new TH1D("t_truth" ,";|t|(MC) [(GeV/c^{2})^{2}]"           , 100, 0., 2.);
  TH1D* h_t_B0Acc  = new TH1D("t_b0acc" ,";|t|(MC|Reco. - B0) [(GeV/c^{2})^{2}]", 100, 0., 2.);
  TH1D* h_t_RPAcc  = new TH1D("t_rpacc" ,";|t|(MC|Reco. - RP) [(GeV/c^{2})^{2}]", 100, 0., 2.);
  TH1D* h_t_B0Reco = new TH1D("t_b0reco",";|t|(Reco. - B0) [(GeV/c^{2})^{2}]"   , 100, 0., 2.);
  TH1D* h_t_RPReco = new TH1D("t_rpreco",";|t|(Reco. - RP) [(GeV/c^{2})^{2}]"   , 100, 0., 2.);
  TH2D* h_t_B0Resp = new TH2D("t_b0resp",";|t|(MC|Reco. - B0) [(GeV/c^{2})^{2}];|t|(Reco. - B0) [(GeV/c^{2})^{2}]", 100, 0., 2., 100, 0., 2.);
  TH2D* h_t_RPResp = new TH2D("t_rpresp",";|t|(MC|Reco. - RP) [(GeV/c^{2})^{2}];|t|(Reco. - RP) [(GeV/c^{2})^{2}]", 100, 0., 2., 100, 0., 2.);
  TH1D* h_t_B0Pur  = new TH1D("t_b0pur" ,";|t|(Reco. - B0) [(GeV/c^{2})^{2}];#frac{MC|Reco.}{Reco.}", 100, 0., 2.);
  TH1D* h_t_RPPur  = new TH1D("t_rppur" ,";|t|(Reco. - RP) [(GeV/c^{2})^{2}];#frac{MC|Reco.}{Reco.}", 100, 0., 2.);
  TH2D* h_B0dtvt   = new TH2D("b0dtvt"  ,";|t|(MC|Reco. - B0) [(GeV/c^{2})^{2}];#Delta t [(GeV/c^{2})^{2}]", 100, 0., 2., 400, -2., 2.);
  TH2D* h_RPdtvt   = new TH2D("rpdtvt"  ,";|t|(MC|Reco. - RP) [(GeV/c^{2})^{2}];#Delta t [(GeV/c^{2})^{2}]", 100, 0., 2., 400, -2., 2.);
  //Single particle kinematics - protons
  TH1D* h_theta_p_Truth  = new TH1D("theta_p_truth" ,";#theta_{p'}(MC) [mrad]"           , 100, 0., 25.);
  TH1D* h_theta_p_B0Acc  = new TH1D("theta_p_b0acc" ,";#theta_{p'}(MC|Reco. - B0) [mrad]", 100, 0., 25.);
  TH1D* h_theta_p_RPAcc  = new TH1D("theta_p_rpacc" ,";#theta_{p'}(MC|Reco. - RP) [mrad]", 100, 0., 25.);
  TH1D* h_theta_p_B0Reco = new TH1D("theta_p_b0reco",";#theta_{p'}(Reco. - B0) [mrad]"   , 100, 0., 25.);
  TH1D* h_theta_p_RPReco = new TH1D("theta_p_rpreco",";#theta_{p'}(Reco. - RP) [mrad]"   , 100, 0., 25.);
  TH2D* h_theta_p_B0Resp = new TH2D("theta_p_b0resp",";#theta_{p'}(MC|Reco. - B0) [mrad];#theta_{p'}(Reco. - B0) [mrad]", 100, 0., 25., 100, 0., 25.);
  TH2D* h_theta_p_RPResp = new TH2D("theta_p_rpresp",";#theta_{p'}(MC|Reco. - RP) [mrad];#theta_{p'}(Reco. - RP) [mrad]", 100, 0., 25., 100, 0., 25.);
  TH1D* h_theta_p_B0Pur  = new TH1D("theta_p_b0pur" ,";#theta_{p'}(Reco. - B0) [mrad];#frac{MC|Reco.}{Reco.}", 100, 0., 25.);
  TH1D* h_theta_p_RPPur  = new TH1D("theta_p_rppur" ,";#theta_{p'}(Reco. - RP) [mrad];#frac{MC|Reco.}{Reco.}", 100, 0., 25.);
  TH1D* h_E_p_Truth  = new TH1D("E_p_truth" ,";E_{p'}(MC) [GeV]"           , (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_B0Acc  = new TH1D("E_p_b0acc" ,";E_{p'}(MC|Reco. - B0) [GeV]", (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_RPAcc  = new TH1D("E_p_rpacc" ,";E_{p'}(MC|Reco. - RP) [GeV]", (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_B0Reco = new TH1D("E_p_b0reco",";E_{p'}(Reco. - B0) [GeV]"   , (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_RPReco = new TH1D("E_p_rpreco",";E_{p'}(Reco. - RP) [GeV]"   , (Int_t)fPMax_p, 0., fPMax_p);
  TH2D* h_E_p_B0Resp = new TH2D("E_p_b0resp",";E_{p'}(MC|Reco. - B0) [GeV];E_{p'}(Reco. - B0) [GeV]", (Int_t)fPMax_p, 0., fPMax_p, (Int_t)fPMax_p, 0., fPMax_p);
  TH2D* h_E_p_RPResp = new TH2D("E_p_rpresp",";E_{p'}(MC|Reco. - RP) [GeV];E_{p'}(Reco. - RP) [GeV]", (Int_t)fPMax_p, 0., fPMax_p, (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_B0Pur  = new TH1D("E_p_b0pur" ,";E_{p'}(Reco. - B0) [GeV];#frac{MC|Reco.}{Reco.}", (Int_t)fPMax_p, 0., fPMax_p);
  TH1D* h_E_p_RPPur  = new TH1D("E_p_rppur" ,";E_{p'}(Reco. - RP) [GeV];#frac{MC|Reco.}{Reco.}", (Int_t)fPMax_p, 0., fPMax_p);
  //Single particle kinematics - electrons
  TH1D* h_theta_e_Truth = new TH1D("theta_e_truth",";#theta_{e'}(MC) [rad]"      , 100, 0., 3.2);
  TH1D* h_theta_e_Acc   = new TH1D("theta_e_acc"  ,";#theta_{e'}(MC|Reco.) [rad]", 100, 0., 3.2);
  TH1D* h_theta_e_Reco  = new TH1D("theta_e_reco" ,";#theta_{e'}(Reco.) [rad]"   , 100, 0., 3.2);
  TH2D* h_theta_e_Resp  = new TH2D("theta_e_resp" ,";#theta_{e'}(MC|Reco.) [rad];#theta_{e'}(Reco.) [rad]", 100, 0., 3.2, 100, 0., 3.2);
  TH1D* h_theta_e_Pur   = new TH1D("theta_e_pur"  ,";#theta_{e'}(Reco.) [rad];#frac{MC|Reco.}{Reco.}", 100, 0., 3.2);
  TH1D* h_E_e_Truth = new TH1D("E_e_truth",";E_{e'}(MC) [GeV]"      , (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_E_e_Acc   = new TH1D("E_e_acc"  ,";E_{e'}(MC|Reco.) [GeV]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_E_e_Reco  = new TH1D("E_e_reco" ,";E_{e'}(Reco.) [GeV]"   , (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH2D* h_E_e_Resp  = new TH2D("E_e_resp" ,";E_{e'}(MC|Reco.) [GeV];E_{e'}(Reco.) [GeV]", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e, (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  TH1D* h_E_e_Pur   = new TH1D("E_e_pur"  ,";E_{e'}(Reco.) [GeV];#frac{MC|Reco.}{Reco.}", (Int_t)2*fPMax_e, 0.0, 2*fPMax_e);
  //Single particle kinematics - photons
  TH1D* h_theta_g_Truth = new TH1D("theta_g_truth",";#theta_{#gamma}(MC) [rad]"     , 100, 0., 3.2);
  TH1D* h_theta_g_Acc   = new TH1D("theta_g_acc"  ,";#theta_{#gamma}(MC|Reco) [rad]", 100, 0., 3.2);
  TH1D* h_theta_g_Reco  = new TH1D("theta_g_reco" ,";#theta_{#gamma}(Reco) [rad]"   , 100, 0., 3.2);
  TH2D* h_theta_g_Resp  = new TH2D("theta_g_resp" ,";#theta_{#gamma}(MC|Reco.) [rad];#theta_{#gamma}(Reco.) [rad]", 100, 0., 3.2, 100, 0., 3.2);
  TH1D* h_theta_g_Pur   = new TH1D("theta_g_pur"  ,";#theta_{#gamma}(Reco.) [rad];#frac{MC|Reco.}{Reco.}", 100, 0., 3.2);
  TH1D* h_E_g_Truth = new TH1D("E_g_truth",";E_{#gamma}(MC) [GeV]", 50, 0., 50.);
  TH1D* h_E_g_Acc   = new TH1D("E_g_acc"  ,";E_{#gamma}(MC|Reco.) [GeV]", 50, 0., 50.);
  TH1D* h_E_g_Reco  = new TH1D("E_g_reco" ,";E_{#gamma}(Reco.) [GeV]", 50, 0., 50.);
  TH2D* h_E_g_Resp  = new TH2D("E_g_resp" ,";E_{#gamma}(MC|Reco.) [GeV];E_{#gamma}(Reco.) [GeV]", 50, 0., 50., 50, 0., 50.);
  TH1D* h_E_g_Pur   = new TH1D("E_g_pur"  ,";E_{#gamma}(Reco.) [GeV];#frac{MC|Reco.}{Reco.}", 50, 0., 50.);
  // HFS quantities
  TH1D* h_HFSSigma_Truth = new TH1D("hfssigma_truth",";#Sigma_{h}(MC) [GeV]"      , 100, 0., 10.);
  TH1D* h_HFSSigma_Acc   = new TH1D("hfssigma_acc"  ,";#Sigma_{h}(MC|Reco.) [GeV]", 100, 0., 10.);
  TH1D* h_HFSSigma_Reco  = new TH1D("hfssigma_reco" ,";#Sigma_{h}(Reco.) [GeV]"   , 100, 0., 10.);
  TH2D* h_HFSSigma_Resp  = new TH2D("hfssigma_resp" ,";#Sigma_{h}(MC|Reco.) [GeV];#Sigma_{h}(Reco.) [GeV]", 100, 0., 10., 100, 0., 10.);
  TH1D* h_HFSSigma_Pur   = new TH1D("hfssigma_pur"  ,";#Sigma_{h}(Reco.) [GeV];#frac{MC|Reco.}{Reco.}", 100, 0., 10.);
  TH1D* h_HFSpT2_Truth = new TH1D("hfspt2_truth",";p_{T,h}^{2}(MC) [(GeV/c^{2})^{2}]"      , 100, 0., 10.);
  TH1D* h_HFSpT2_Acc   = new TH1D("hfspt2_acc"  ,";p_{T,h}^{2}(MC|Reco.) [(GeV/c^{2})^{2}]", 100, 0., 10.);
  TH1D* h_HFSpT2_Reco  = new TH1D("hfspt2_reco" ,";p_{T,h}^{2}(Reco.) [(GeV/c^{2})^{2}]"   , 100, 0., 10.);
  TH2D* h_HFSpT2_Resp  = new TH2D("hfspt2_resp" ,";p_{T,h}^{2}(MC|Reco.) [(GeV/c^{2})^{2}];p_{T,h}^{2}(Reco.) [(GeV/c^{2})^{2}]", 100, 0., 10., 100, 0., 10.);
  TH1D* h_HFSpT2_Pur   = new TH1D("hfspt2_pur"  ,";p_{T,h}^{2}(Reco.) [(GeV/c^{2})^{2}];#frac{MC|Reco}{Reco.}", 100, 0., 10.);
  TH1D* h_FullSigma_Truth = new TH1D("fullsigma_truth",";#Sigma(MC) [GeV]"      , 100, 0., 10.);
  TH1D* h_FullSigma_Acc   = new TH1D("fullsigma_acc"  ,";#Sigma(MC|Reco.) [GeV]", 100, 0., 10.);
  TH1D* h_FullSigma_Reco  = new TH1D("fullsigma_reco" ,";#Sigma(Reco.) [GeV]"   , 100, 0., 10.);
  TH2D* h_FullSigma_Resp  = new TH2D("fullsigma_resp" ,";#Sigma(MC|Reco.) [GeV];#Sigma(Reco.) [GeV]", 100, 0., 10., 100, 0., 10.);
  TH1D* h_FullSigma_Pur   = new TH1D("fullsigma_pur"  ,";#Sigma(Reco.) [GeV];#frac{MC|Reco}{Reco.}", 100, 0., 10.);
  // E/p for electron/photon
  TH1D* h_EOverp_e_Truth = new TH1D("eoverp_e_truth",";E/p(MC)", 200, 0., 2.);
  TH1D* h_EOverp_e_Reco  = new TH1D("eoverp_e_reco" ,";E/p(Reco.)", 200, 0., 2.);
  TH2D* h_EOverp_e_Resp  = new TH2D("eoverp_e_resp" ,";E/p(MC|Reco.);E/p(Reco.)", 200, 0., 2., 200, 0., 2.);
  TH1D* h_EOverp_g_Truth = new TH1D("eoverp_g_truth",";E/p(MC)", 200, 0., 2.);
  TH1D* h_EOverp_g_Reco  = new TH1D("eoverp_g_reco" ,";E/p(Reco.)", 200, 0., 2.);
  TH2D* h_EOverp_g_Resp  = new TH2D("eoverp_g_resp" ,";E/p(MC|Reco.);E/p(Reco.)", 200, 0., 2., 200, 0., 2.);
  
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
    // RP hits
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
      // USING EXPLICIT TRACK MATCHING
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
        
      // Add in RP hits - only looking at protons
      // NO NEED TO UNDO AFTERBURNER FOR FF DETECTORS - NOT APPLIED IN FIRST PLACE
      for(int irpreco{0}; irpreco<rp_px_array.GetSize(); irpreco++){
 	recotrk.SetXYZ(rp_px_array[irpreco], rp_py_array[irpreco], rp_pz_array[irpreco]);	
	// Assume track is a proton
	P3EVector q_rpreco(recotrk.X(),recotrk.Y(),recotrk.Z(),calcE(recotrk,fMass_proton));
	if(rp_pdg_array[irpreco] == 2212){
	  scatp4_rom.push_back(q_rpreco);
	}
      }// End of RP reconstructed particles loop
            
      //---------------------------------------------------------
      // Fill histograms
      //---------------------------------------------------------
      // Fill histograms
      // Inclusive electrons
      // MC - electron energy, theta, E/p, Q2, x and y
      if(scate4_gen.size() == 1){
	h_Q2_Truth->Fill(calcQ2_Elec(beame4,scate4_gen[0]));
	h_xB_Truth->Fill(calcX_Elec(beame4,beamp4,scate4_gen[0]));
	h_y_Truth->Fill(calcY_Elec(beame4,beamp4,scate4_gen[0]));
	h_theta_e_Truth->Fill(scate4_gen[0].Theta());
	h_E_e_Truth->Fill(scate4_gen[0].E());
	h_EOverp_e_Truth->Fill(scate4_gen[0].E()/scate4_gen[0].P());
      }
      // MC accepted
      if(scate4_aso.size() == 1){
	h_Q2_Acc->Fill(calcQ2_Elec(beame4,scate4_aso[0]));
	h_xB_Acc->Fill(calcX_Elec(beame4,beamp4,scate4_aso[0]));
	h_y_Acc->Fill(calcY_Elec(beame4,beamp4,scate4_aso[0]));
	h_theta_e_Acc->Fill(scate4_aso[0].Theta());
	h_E_e_Acc->Fill(scate4_aso[0].E());
      }
      // Reconstructed
      if(scate4_rec.size() == 1){
	h_Q2_Reco->Fill(calcQ2_Elec(beame4,scate4_rec[0]));
	h_xB_Reco->Fill(calcX_Elec(beame4,beamp4,scate4_rec[0]));
	h_y_Reco->Fill(calcY_Elec(beame4,beamp4,scate4_rec[0]));
	h_theta_e_Reco->Fill(scate4_rec[0].Theta());
	h_E_e_Reco->Fill(scate4_rec[0].E());
	h_EOverp_e_Reco->Fill(scate4_rec[0].E()/scate4_rec[0].P());
      }
      // Reco. vs MCA 
      if(scate4_aso.size() == 1 && scate4_rec.size() == 1){
	Float_t Q2_acc = calcQ2_Elec(beame4,scate4_aso[0]);
	Float_t Q2_rec = calcQ2_Elec(beame4,scate4_rec[0]);
	Float_t x_acc = calcX_Elec(beame4,beamp4,scate4_aso[0]);
	Float_t x_rec = calcX_Elec(beame4,beamp4,scate4_rec[0]);
	Float_t y_acc = calcY_Elec(beame4,beamp4,scate4_aso[0]);
	Float_t y_rec = calcY_Elec(beame4,beamp4,scate4_rec[0]);
	// 2D response histograms
	h_Q2_Resp->Fill(Q2_acc, Q2_rec);
	h_xB_Resp->Fill(x_acc, x_rec);
	h_y_Resp->Fill(y_acc, y_rec);
	h_theta_e_Resp->Fill(scate4_aso[0].Theta(),scate4_rec[0].Theta());
	h_E_e_Resp->Fill(scate4_aso[0].E(),scate4_rec[0].E());
	h_EOverp_e_Resp->Fill(scate4_aso[0].E()/scate4_aso[0].P(), scate4_rec[0].E()/scate4_rec[0].P());

	// 2D difference histograms
	h_dQ2vQ2->Fill(Q2_acc, Q2_rec-Q2_acc);
	h_dxBvxB->Fill( x_acc, x_rec-x_acc);
	h_dyvy->Fill(y_acc, y_rec-y_acc);

	// Purities
	// Fill bin if reco matches MC
	if(h_Q2_Acc->FindBin(Q2_acc) == h_Q2_Reco->FindBin(Q2_rec)) h_Q2_Pur->Fill(Q2_rec);
	if(h_xB_Acc->FindBin(x_acc) == h_xB_Reco->FindBin(x_rec))   h_xB_Pur->Fill(x_rec);
	if(h_y_Acc->FindBin(y_acc) == h_y_Reco->FindBin(y_rec))	    h_y_Pur->Fill(y_rec);
	if(h_theta_e_Acc->FindBin(scate4_aso[0].Theta()) == h_theta_e_Reco->FindBin(scate4_rec[0].Theta())) h_theta_e_Pur->Fill(scate4_rec[0].Theta());
	if(h_E_e_Acc->FindBin(scate4_aso[0].E()) == h_E_e_Reco->FindBin(scate4_rec[0].E())) h_E_e_Pur->Fill(scate4_rec[0].E());
      }

      // Inclusive photons
      // MC - photon energy, theta, E/p
      if(scatg4_gen.size() == 1){
	h_theta_g_Truth->Fill(scatg4_gen[0].Theta());
	h_E_g_Truth->Fill(scatg4_gen[0].E());
	Float_t EoverP = scatg4_gen[0].E()/scatg4_gen[0].P();
	h_EOverp_g_Truth->Fill(EoverP);
      }
      // MC accepted
      if(scatg4_aso.size() == 1){
	h_theta_g_Acc->Fill(scatg4_aso[0].Theta());
	h_E_g_Acc->Fill(scatg4_aso[0].E());
      }
      // Reconstructed
      if(scatg4_rec.size() == 1){
	h_theta_g_Reco->Fill(scatg4_rec[0].Theta());
	h_E_g_Reco->Fill(scatg4_rec[0].E());
	Float_t EoverP_rec = scatg4_rec[0].E()/scatg4_rec[0].P();
	h_EOverp_g_Reco->Fill(EoverP_rec);
      }
      // Reco. vs MCA 
      if(scatg4_aso.size() == 1 && scatg4_rec.size() == 1){
	// Responses
	h_theta_g_Resp->Fill(scatg4_aso[0].Theta(),scatg4_rec[0].Theta());
	h_E_g_Resp->Fill(scatg4_aso[0].E(),scatg4_rec[0].E());
	Float_t EoverP_acc = scatg4_aso[0].E()/scatg4_aso[0].P();
	Float_t EoverP_rec = scatg4_rec[0].E()/scatg4_rec[0].P();
	h_EOverp_g_Resp->Fill(EoverP_acc,EoverP_rec);

	// Purities
	// Fill bin if reco matches MC
	if(h_theta_g_Acc->FindBin(scatg4_aso[0].Theta()) == h_theta_g_Reco->FindBin(scatg4_rec[0].Theta())) h_theta_g_Pur->Fill(scatg4_rec[0].Theta());
	if(h_E_g_Acc->FindBin(scatg4_aso[0].E()) == h_E_g_Reco->FindBin(scatg4_rec[0].E())) h_E_g_Pur->Fill(scatg4_rec[0].E());
      }

      // Inclusive protons
      // MC - photon energy, theta
      if(scatp4_gen.size() == 1){
	h_theta_p_Truth->Fill(scatp4_gen[0].Theta()*1000);
	h_E_p_Truth->Fill(scatp4_gen[0].E());
      }
      // MC accepted - only works for B0
      if(scatp4_aso.size() == 1){
	h_theta_p_B0Acc->Fill(scatp4_aso[0].Theta()*1000);
	h_E_p_B0Acc->Fill(scatp4_aso[0].E());
      }
      // Reconstructed - B0
      if(scatp4_rec.size() == 1 && scatp4_rom.size() == 0){
	h_theta_p_B0Reco->Fill(scatp4_rec[0].Theta()*1000);
	h_E_p_B0Reco->Fill(scatp4_rec[0].E());

	// And 2D response/purity
	if(scatp4_aso.size() == 1){
	  h_theta_p_B0Resp->Fill(scatp4_aso[0].Theta()*1000,scatp4_rec[0].Theta()*1000);
	  h_E_p_B0Resp->Fill(scatp4_aso[0].E(),scatp4_rec[0].E());

	  if(h_theta_p_B0Acc->FindBin(scatp4_aso[0].Theta()*1000) == h_theta_p_B0Reco->FindBin(scatp4_rec[0].Theta()*1000)) h_theta_p_B0Pur->Fill(scatp4_rec[0].Theta()*1000);
	  if(h_E_p_B0Acc->FindBin(scatp4_aso[0].E()) == h_E_p_B0Reco->FindBin(scatp4_rec[0].E())) h_E_p_B0Pur->Fill(scatp4_rec[0].E());
	}
      }
      // Reconstructed - RP
      if(scatp4_rom.size() == 1 && scatp4_rec.size() == 0){
	h_theta_p_RPReco->Fill(scatp4_rom[0].Theta()*1000);
	h_E_p_RPReco->Fill(scatp4_rom[0].E());
	
	// And 2D response/purity
	if(scatp4_gen.size() == 1){
	  h_theta_p_RPAcc->Fill(scatp4_gen[0].Theta()*1000);
	  h_E_p_RPAcc->Fill(scatp4_gen[0].E());

	  h_theta_p_RPResp->Fill(scatp4_gen[0].Theta()*1000,scatp4_rom[0].Theta()*1000);
	  h_E_p_RPResp->Fill(scatp4_gen[0].E(),scatp4_rom[0].E());

	  if(h_theta_p_RPAcc->FindBin(scatp4_gen[0].Theta()*1000) == h_theta_p_RPReco->FindBin(scatp4_rom[0].Theta()*1000)) h_theta_p_RPPur->Fill(scatp4_rom[0].Theta()*1000);
	  if(h_E_p_RPAcc->FindBin(scatp4_gen[0].E()) == h_E_p_RPReco->FindBin(scatp4_rom[0].E())) h_E_p_RPPur->Fill(scatp4_rom[0].E());
	}
      }

      // Event variables - t
      // MC
      if(scate4_gen.size() == 1 && scatg4_gen.size() == 1 && scatp4_gen.size() == 1) h_t_Truth->Fill(calcT_BABE(beamp4,scatp4_gen[0]));
      // MC accepted - B0 only
      if(scate4_aso.size() == 1 && scatg4_aso.size() == 1 && scatp4_aso.size() == 1) h_t_B0Acc->Fill(calcT_BABE(beamp4,scatp4_aso[0]));
      // Reconstructed - B0 only
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rec.size() == 1 && scatp4_rom.size() == 0){
	h_t_B0Reco->Fill(calcT_BABE(beamp4,scatp4_rec[0]));

	// And 2D response/purity
	if(scate4_aso.size() == 1 && scatg4_aso.size() == 1 && scatp4_aso.size() == 1){
	  Float_t t_acc = calcT_BABE(beamp4,scatp4_aso[0]);
	  Float_t t_rec = calcT_BABE(beamp4,scatp4_rec[0]);
	  h_t_B0Resp->Fill(t_acc,t_rec);
	  h_B0dtvt->Fill(t_acc,t_rec-t_acc);

	  if(h_t_B0Acc->FindBin(t_acc) == h_t_B0Reco->FindBin(t_rec)) h_t_B0Pur->Fill(t_rec);
	}
      }
      // Reconstructed - RP only
      if(scate4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rom.size() == 1 && scatp4_rec.size() == 0){
	h_t_RPReco->Fill(calcT_BABE(beamp4,scatp4_rom[0]));

	// And 2D response/purity (and MC accepted)
	if(scate4_aso.size() == 1 && scatg4_aso.size() == 1 && scatp4_gen.size() == 1){
	  Float_t t_acc = calcT_BABE(beamp4,scatp4_gen[0]);
	  Float_t t_rec = calcT_BABE(beamp4,scatp4_rom[0]);
	  
	  h_t_RPAcc->Fill(t_acc);
	  h_t_RPReco->Fill(t_rec);
	  h_t_RPResp->Fill(t_acc,t_rec);
	  h_RPdtvt->Fill(t_acc,t_rec-t_acc);
	  
	  if(h_t_RPAcc->FindBin(t_acc) == h_t_RPReco->FindBin(t_rec)) h_t_RPPur->Fill(t_rec);
	}
      }

      // HFS quantities - Sigma, PT^2
      // And full event sigma
      // MC
      if(scatp4_gen.size() == 1 && scatg4_gen.size() == 1 && scatp4_gen.size() == 1){
	Float_t sigma = (scatg4_gen[0]+scatp4_gen[0]).E() - (scatg4_gen[0]+scatp4_gen[0]).Pz();
	Float_t pT2 = TMath::Power((scatg4_gen[0]+scatp4_gen[0]).Px(), 2) + TMath::Power((scatg4_gen[0]+scatp4_gen[0]).Py(), 2);

	h_HFSSigma_Truth->Fill(sigma);
	h_HFSpT2_Truth->Fill(pT2);
	// And full sigma
	if(scate4_gen.size() == 1){
	  Float_t fullsig = (scatg4_gen[0]+scatp4_gen[0]+scate4_gen[0]).E() - (scatg4_gen[0]+scatp4_gen[0]+scate4_gen[0]).Pz();
	  h_FullSigma_Truth->Fill(fullsig);
	}
      }
      // MC accepted - for B0 only
      if(scatp4_aso.size() == 1 && scatg4_aso.size() == 1 && scatp4_aso.size() == 1){
	Float_t sigma = (scatg4_aso[0]+scatp4_aso[0]).E() - (scatg4_aso[0]+scatp4_aso[0]).Pz();
	Float_t pT2 = TMath::Power((scatg4_aso[0]+scatp4_aso[0]).Px(), 2) + TMath::Power((scatg4_aso[0]+scatp4_aso[0]).Py(), 2);

	h_HFSSigma_Acc->Fill(sigma);
	h_HFSpT2_Acc->Fill(pT2);
	// And full sigma
	if(scate4_aso.size() == 1){
	  Float_t fullsig = (scatg4_aso[0]+scatp4_aso[0]+scate4_aso[0]).E() - (scatg4_aso[0]+scatp4_aso[0]+scate4_aso[0]).Pz();
	  h_FullSigma_Acc->Fill(fullsig);
	}
      }
      // Reconstructed - for B0 only
      if(scatp4_rec.size() == 1 && scatg4_rec.size() == 1 && scatp4_rom.size() == 0){
	Float_t sigma_rec = (scatg4_rec[0]+scatp4_rec[0]).E() - (scatg4_rec[0]+scatp4_rec[0]).Pz();
	Float_t pT2_rec = TMath::Power((scatg4_rec[0]+scatp4_rec[0]).Px(), 2) + TMath::Power((scatg4_rec[0]+scatp4_rec[0]).Py(), 2);

	h_HFSSigma_Reco->Fill(sigma_rec);
	h_HFSpT2_Reco->Fill(pT2_rec);

	// 2D response/purity
	if(scatp4_aso.size() == 1 && scatg4_aso.size() == 1 && scatp4_aso.size() == 1){
	  Float_t sigma_acc = (scatg4_aso[0]+scatp4_aso[0]).E() - (scatg4_aso[0]+scatp4_aso[0]).Pz();
	  Float_t pT2_acc = TMath::Power((scatg4_aso[0]+scatp4_aso[0]).Px(), 2) + TMath::Power((scatg4_aso[0]+scatp4_aso[0]).Py(), 2);

	  h_HFSSigma_Resp->Fill(sigma_acc,sigma_rec);
	  h_HFSpT2_Resp->Fill(pT2_acc,pT2_rec);

	  if(h_HFSSigma_Acc->FindBin(sigma_acc) == h_HFSSigma_Reco->FindBin(sigma_rec)) h_HFSSigma_Pur->Fill(sigma_rec);
	  if(h_HFSpT2_Acc->FindBin(pT2_acc) == h_HFSpT2_Reco->FindBin(pT2_rec)) h_HFSpT2_Pur->Fill(pT2_rec);
	}

	// And full sigma
	if(scate4_rec.size() == 1){
	  Float_t fullsig_rec = (scatg4_rec[0]+scatp4_rec[0]+scate4_rec[0]).E() - (scatg4_rec[0]+scatp4_rec[0]+scate4_rec[0]).Pz();
	  h_FullSigma_Reco->Fill(fullsig_rec);

	  // 2D response/purity
	  if(scate4_aso.size() == 1){
	    Float_t fullsig_acc = (scatg4_aso[0]+scatp4_aso[0]+scate4_aso[0]).E() - (scatg4_aso[0]+scatp4_aso[0]+scate4_aso[0]).Pz();
	    
	    h_FullSigma_Resp->Fill(fullsig_acc,fullsig_rec);
	    if(h_FullSigma_Acc->FindBin(fullsig_acc) == h_FullSigma_Reco->FindBin(fullsig_rec)) h_FullSigma_Pur->Fill(fullsig_rec);
	  }
	}
      }
      // Reconstructed - for RP only
      if(scatp4_rom.size() == 1 && scatg4_rec.size() == 1 && scatp4_rec.size() == 0){
	Float_t sigma_rec = (scatg4_rec[0]+scatp4_rom[0]).E() - (scatg4_rec[0]+scatp4_rom[0]).Pz();
	Float_t pT2_rec = TMath::Power((scatg4_rec[0]+scatp4_rom[0]).Px(), 2) + TMath::Power((scatg4_rec[0]+scatp4_rom[0]).Py(), 2);

	h_HFSSigma_Reco->Fill(sigma_rec);
	h_HFSpT2_Reco->Fill(pT2_rec);

	// 2D response/purity
	if(scatp4_gen.size() == 1 && scatg4_aso.size() == 1){
	  Float_t sigma_acc = (scatg4_aso[0]+scatp4_gen[0]).E() - (scatg4_aso[0]+scatp4_gen[0]).Pz();
	  Float_t pT2_acc = TMath::Power((scatg4_aso[0]+scatp4_gen[0]).Px(), 2) + TMath::Power((scatg4_aso[0]+scatp4_gen[0]).Py(), 2);

	  h_HFSSigma_Resp->Fill(sigma_acc,sigma_rec);
	  h_HFSpT2_Resp->Fill(pT2_acc,pT2_rec);

	  if(h_HFSSigma_Acc->FindBin(sigma_acc) == h_HFSSigma_Reco->FindBin(sigma_rec)) h_HFSSigma_Pur->Fill(sigma_rec);
	  if(h_HFSpT2_Acc->FindBin(pT2_acc) == h_HFSpT2_Reco->FindBin(pT2_rec)) h_HFSpT2_Pur->Fill(pT2_rec);
	}

	// And full sigma
	if(scate4_rec.size() == 1){
	  Float_t fullsig_rec = (scatg4_rec[0]+scatp4_rom[0]+scate4_rec[0]).E() - (scatg4_rec[0]+scatp4_rom[0]+scate4_rec[0]).Pz();
	  h_FullSigma_Reco->Fill(fullsig_rec);

	  // 2D response/purity
	  if(scate4_aso.size() == 1){
	    Float_t fullsig_acc = (scatg4_aso[0]+scatp4_gen[0]+scate4_aso[0]).E() - (scatg4_aso[0]+scatp4_gen[0]+scate4_aso[0]).Pz();
	    
	    h_FullSigma_Resp->Fill(fullsig_acc,fullsig_rec);
	    if(h_FullSigma_Acc->FindBin(fullsig_acc) == h_FullSigma_Reco->FindBin(fullsig_rec)) h_FullSigma_Pur->Fill(fullsig_rec);
	  }
	}
      }

    } // END EVENT/TTREEREADER LOOP

    inputRootFile->Close();
  } // END OF FILE LOOP
  
  // For purity/bin migration histograms, divide by reco.
  h_Q2_Pur->Divide(h_Q2_Reco);
  h_xB_Pur->Divide(h_xB_Reco);
  h_y_Pur->Divide(h_y_Reco);
  h_t_B0Pur->Divide(h_t_B0Reco);
  h_t_RPPur->Divide(h_t_RPReco);
  h_theta_p_B0Pur->Divide(h_theta_p_B0Reco);
  h_theta_p_RPPur->Divide(h_theta_p_RPReco);
  h_E_p_B0Pur->Divide(h_E_p_B0Reco);
  h_E_p_RPPur->Divide(h_E_p_RPReco);
  h_theta_e_Pur->Divide(h_theta_e_Reco);
  h_E_e_Pur->Divide(h_E_e_Reco);
  h_theta_g_Pur->Divide(h_theta_g_Reco);
  h_E_g_Pur->Divide(h_E_g_Reco);
  h_HFSSigma_Pur->Divide(h_HFSSigma_Reco);
  h_HFSpT2_Pur->Divide(h_HFSpT2_Reco);
  h_FullSigma_Pur->Divide(h_FullSigma_Reco);
  
  //------------------------------------------------------------
  // Write to output file
  //------------------------------------------------------------
  fOutFile->cd();
  h_Q2_Truth->Write();
  h_Q2_Acc->Write();
  h_Q2_Reco->Write();
  h_Q2_Resp->Write();
  h_Q2_Pur->Write();
  h_dQ2vQ2->Write();
  h_xB_Truth->Write();
  h_xB_Acc->Write();
  h_xB_Reco->Write();
  h_xB_Resp->Write();
  h_xB_Pur->Write();
  h_dxBvxB->Write();
  h_y_Truth->Write();
  h_y_Acc->Write();
  h_y_Reco->Write();
  h_y_Resp->Write();
  h_y_Pur->Write();
  h_dyvy->Write();
  h_t_Truth->Write();
  h_t_B0Acc->Write();
  h_t_RPAcc->Write();
  h_t_B0Reco->Write();
  h_t_RPReco->Write();
  h_t_B0Resp->Write();
  h_t_RPResp->Write();
  h_t_B0Pur->Write();
  h_t_RPPur->Write();
  h_B0dtvt->Write();
  h_RPdtvt->Write();
  h_theta_p_Truth->Write();
  h_theta_p_B0Acc->Write();
  h_theta_p_RPAcc->Write();
  h_theta_p_B0Reco->Write();
  h_theta_p_RPReco->Write();
  h_theta_p_B0Resp->Write();
  h_theta_p_RPResp->Write();
  h_theta_p_B0Pur->Write();
  h_theta_p_RPPur->Write();
  h_E_p_Truth->Write();
  h_E_p_B0Acc->Write();
  h_E_p_RPAcc->Write();
  h_E_p_B0Reco->Write();
  h_E_p_RPReco->Write();
  h_E_p_B0Resp->Write();
  h_E_p_RPResp->Write();
  h_E_p_B0Pur->Write();
  h_E_p_RPPur->Write();
  h_theta_e_Truth->Write();
  h_theta_e_Acc->Write();
  h_theta_e_Reco->Write();
  h_theta_e_Resp->Write();
  h_theta_e_Pur->Write();
  h_E_e_Truth->Write();
  h_E_e_Acc->Write();
  h_E_e_Reco->Write();
  h_E_e_Resp->Write();
  h_E_e_Pur->Write();
  h_theta_g_Truth->Write();
  h_theta_g_Acc->Write();
  h_theta_g_Reco->Write();
  h_theta_g_Resp->Write();
  h_theta_g_Pur->Write();
  h_E_g_Truth->Write();
  h_E_g_Acc->Write();
  h_E_g_Reco->Write();
  h_E_g_Resp->Write();
  h_E_g_Pur->Write();
  h_HFSSigma_Truth->Write();
  h_HFSSigma_Acc->Write();
  h_HFSSigma_Reco->Write();
  h_HFSSigma_Resp->Write();
  h_HFSSigma_Pur->Write();
  h_HFSpT2_Truth->Write();
  h_HFSpT2_Acc->Write();
  h_HFSpT2_Reco->Write();
  h_HFSpT2_Resp->Write();
  h_HFSpT2_Pur->Write();
  h_FullSigma_Truth->Write();
  h_FullSigma_Acc->Write();
  h_FullSigma_Reco->Write();
  h_FullSigma_Resp->Write();
  h_FullSigma_Pur->Write();
  h_EOverp_e_Truth->Write();
  h_EOverp_e_Reco->Write();
  h_EOverp_e_Resp->Write();
  h_EOverp_g_Truth->Write();
  h_EOverp_g_Reco->Write();
  h_EOverp_g_Resp->Write();

  fOutFile->Close();

  return;
}
