using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


//---------------------------------------------------------------------
// MAIN
//---------------------------------------------------------------------
void DVCSPlots(TString campaign, TString energy, TString setting){
  //---------------------------------------------------------------------
  // Get data file from DVCSAnalysis script
  //---------------------------------------------------------------------
  TString inFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_" + campaign + "_" + energy + "_" + setting + ".root";
  //TString inFileName = "/scratch/ePIC_DVCS_" + campaign + "_" + energy + "_" + setting + ".root";
  cout<<"Input file: "<<inFileName<<endl;
  TFile* inFile = new TFile(inFileName);

  //---------------------------------------------------------------------
  // Extract histograms
  //---------------------------------------------------------------------
  // Particle kinematics - MC generated protons, electrons and photons
  TH1D* h_eta_MCp   = (TH1D*)inFile->Get("eta_MCp");
  TH1D* h_pt_MCp    = (TH1D*)inFile->Get("pt_MCp");
  TH1D* h_p_MCp     = (TH1D*)inFile->Get("p_MCp");
  TH1D* h_theta_MCp = (TH1D*)inFile->Get("theta_MCp");
  TH1D* h_phi_MCp   = (TH1D*)inFile->Get("phi_MCp");
  TH1D* h_E_MCp     = (TH1D*)inFile->Get("E_MCp");
  TH1D* h_eta_MCe   = (TH1D*)inFile->Get("eta_MCe");
  TH1D* h_pt_MCe    = (TH1D*)inFile->Get("pt_MCe");
  TH1D* h_p_MCe     = (TH1D*)inFile->Get("p_MCe");
  TH1D* h_theta_MCe = (TH1D*)inFile->Get("theta_MCe");
  TH1D* h_phi_MCe   = (TH1D*)inFile->Get("phi_MCe");
  TH1D* h_E_MCe     = (TH1D*)inFile->Get("E_MCe");
  TH1D* h_eta_MCg   = (TH1D*)inFile->Get("eta_MCg");
  TH1D* h_pt_MCg    = (TH1D*)inFile->Get("pt_MCg");
  TH1D* h_p_MCg     = (TH1D*)inFile->Get("p_MCg");
  TH1D* h_theta_MCg = (TH1D*)inFile->Get("theta_MCg");
  TH1D* h_phi_MCg   = (TH1D*)inFile->Get("phi_MCg");
  TH1D* h_E_MCg     = (TH1D*)inFile->Get("E_MCg");
  // DVCS event kinematics - MC generated particles
  TH1D* h_t_MC    = (TH1D*)inFile->Get("t_MC");
  TH1D* h_Q2_MC   = (TH1D*)inFile->Get("Q2_MC");
  TH1D* h_xB_MC   = (TH1D*)inFile->Get("xB_MC");
  TH1D* h_y_MC    = (TH1D*)inFile->Get("y_MC");
  TH1D* h_TPhi_MC = (TH1D*)inFile->Get("tphi_MC");
  TH1D* h_QPQG_MC = (TH1D*)inFile->Get("qpqg_MC");
  TH1D* h_Cone_MC = (TH1D*)inFile->Get("cone_MC");
  // Particle kinematics - Associated MC protons, electrons and photons
  TH1D* h_eta_MCAp   = (TH1D*)inFile->Get("eta_MCAp");
  TH1D* h_pt_MCAp    = (TH1D*)inFile->Get("pt_MCAp");
  TH1D* h_p_MCAp     = (TH1D*)inFile->Get("p_MCAp");
  TH1D* h_theta_MCAp = (TH1D*)inFile->Get("theta_MCAp");
  TH1D* h_phi_MCAp   = (TH1D*)inFile->Get("phi_MCAp");
  TH1D* h_E_MCAp     = (TH1D*)inFile->Get("E_MCAp");
  TH1D* h_eta_MCAe   = (TH1D*)inFile->Get("eta_MCAe");
  TH1D* h_pt_MCAe    = (TH1D*)inFile->Get("pt_MCAe");
  TH1D* h_p_MCAe     = (TH1D*)inFile->Get("p_MCAe");
  TH1D* h_theta_MCAe = (TH1D*)inFile->Get("theta_MCAe");
  TH1D* h_phi_MCAe   = (TH1D*)inFile->Get("phi_MCAe");
  TH1D* h_E_MCAe     = (TH1D*)inFile->Get("E_MCAe");
  TH1D* h_eta_MCAg   = (TH1D*)inFile->Get("eta_MCAg");
  TH1D* h_pt_MCAg    = (TH1D*)inFile->Get("pt_MCAg");
  TH1D* h_p_MCAg     = (TH1D*)inFile->Get("p_MCAg");
  TH1D* h_theta_MCAg = (TH1D*)inFile->Get("theta_MCAg");
  TH1D* h_phi_MCAg   = (TH1D*)inFile->Get("phi_MCAg");
  TH1D* h_E_MCAg     = (TH1D*)inFile->Get("E_MCAg");
  // DVCS event kinematics - Associated MC particles
  TH1D* h_t_MCA    = (TH1D*)inFile->Get("t_MCA");
  TH1D* h_Q2_MCA   = (TH1D*)inFile->Get("Q2_MCA");
  TH1D* h_xB_MCA   = (TH1D*)inFile->Get("xB_MCA");
  TH1D* h_y_MCA   = (TH1D*)inFile->Get("y_MCA");
  TH1D* h_TPhi_MCA = (TH1D*)inFile->Get("tphi_MCA");
  TH1D* h_QPQG_MCA = (TH1D*)inFile->Get("qpqg_MCA");
  TH1D* h_Cone_MCA = (TH1D*)inFile->Get("cone_MCA");
  // Particle kinematics - Reconstructed protons (B0 ONLY), electrons and photons
  TH1D* h_eta_RPp   = (TH1D*)inFile->Get("eta_RPp");
  TH1D* h_pt_RPp    = (TH1D*)inFile->Get("pt_RPp");
  TH1D* h_p_RPp     = (TH1D*)inFile->Get("p_RPp");
  TH1D* h_theta_RPp = (TH1D*)inFile->Get("theta_RPp");
  TH1D* h_phi_RPp   = (TH1D*)inFile->Get("phi_RPp");
  TH1D* h_E_RPp     = (TH1D*)inFile->Get("E_RPp");
  TH1D* h_eta_RPe   = (TH1D*)inFile->Get("eta_RPe");
  TH1D* h_pt_RPe    = (TH1D*)inFile->Get("pt_RPe");
  TH1D* h_p_RPe     = (TH1D*)inFile->Get("p_RPe");
  TH1D* h_theta_RPe = (TH1D*)inFile->Get("theta_RPe");
  TH1D* h_phi_RPe   = (TH1D*)inFile->Get("phi_RPe");
  TH1D* h_E_RPe     = (TH1D*)inFile->Get("E_RPe");
  TH1D* h_eta_RPg   = (TH1D*)inFile->Get("eta_RPg");
  TH1D* h_pt_RPg    = (TH1D*)inFile->Get("pt_RPg");
  TH1D* h_p_RPg     = (TH1D*)inFile->Get("p_RPg");
  TH1D* h_theta_RPg = (TH1D*)inFile->Get("theta_RPg");
  TH1D* h_phi_RPg   = (TH1D*)inFile->Get("phi_RPg");
  TH1D* h_E_RPg     = (TH1D*)inFile->Get("E_RPg");
  // DVCS event kinematics - Reconstructed particles (BARREL/B0 ONLY)
  TH1D* h_t_RP    = (TH1D*)inFile->Get("t_RP");
  TH1D* h_Q2_RP   = (TH1D*)inFile->Get("Q2_RP");
  TH1D* h_xB_RP   = (TH1D*)inFile->Get("xB_RP");
  TH1D* h_y_RP   = (TH1D*)inFile->Get("y_RP");
  TH1D* h_TPhi_RP = (TH1D*)inFile->Get("tphi_RP");
  TH1D* h_QPQG_RP = (TH1D*)inFile->Get("qpqg_RP");
  TH1D* h_Cone_RP = (TH1D*)inFile->Get("cone_RP");
  // Particle kinematics - Reconstructed protons (ROMAN POTS)
  TH1D* h_eta_RPPp   = (TH1D*)inFile->Get("eta_RPPp");
  TH1D* h_pt_RPPp    = (TH1D*)inFile->Get("pt_RPPp");
  TH1D* h_p_RPPp     = (TH1D*)inFile->Get("p_RPPp");
  TH1D* h_theta_RPPp = (TH1D*)inFile->Get("theta_RPPp");
  TH1D* h_phi_RPPp   = (TH1D*)inFile->Get("phi_RPPp");
  TH1D* h_E_RPPp     = (TH1D*)inFile->Get("E_RPPp");
  // DVCS event kinematics - Reconstructed particles (ROMAN POTS ONLY)
  TH1D* h_t_RPP    = (TH1D*)inFile->Get("t_RPP");
  TH1D* h_TPhi_RPP = (TH1D*)inFile->Get("tphi_RPP");
  TH1D* h_Cone_RPP = (TH1D*)inFile->Get("cone_RPP");
  // 2D distibutions of kinematic variables (MC vs recon.)
  TH2D* h_p_2d = (TH2D*)inFile->Get("h_p_2d");
  TH2D* h_pt_2d = (TH2D*)inFile->Get("h_pt_2d");
  TH2D* h_t_2d = (TH2D*)inFile->Get("h_t_2d");

  // Angular distances between final state particles
  TH1D* h_dphi_MCeg    = (TH1D*)inFile->Get("dphi_MCeg");
  TH1D* h_dtheta_MCeg  = (TH1D*)inFile->Get("dtheta_MCeg");
  TH1D* h_dphi_MCAeg   = (TH1D*)inFile->Get("dphi_MCAeg");
  TH1D* h_dtheta_MCAeg = (TH1D*)inFile->Get("dtheta_MCAeg");
  TH1D* h_dphi_RPeg    = (TH1D*)inFile->Get("dphi_RPeg");
  TH1D* h_dtheta_RPeg  = (TH1D*)inFile->Get("dtheta_RPeg");

  // Occupancies - B0 tracker
  TH2D* h_B0_occupancy_map_layer_0 = (TH2D*)inFile->Get("B0_occupancy_map_0");
  TH2D* h_B0_occupancy_map_layer_1 = (TH2D*)inFile->Get("B0_occupancy_map_1");
  TH2D* h_B0_occupancy_map_layer_2 = (TH2D*)inFile->Get("B0_occupancy_map_2");
  TH2D* h_B0_occupancy_map_layer_3 = (TH2D*)inFile->Get("B0_occupancy_map_3");
  TH1D* h_B0_hit_energy_deposit    = (TH1D*)inFile->Get("B0_tracker_hit_energy_deposit");
  
  // Occupancies - B0 ECAL
  TH2D* h_B0_emcal_occupancy_map  = (TH2D*)inFile->Get("B0_emcal_occupancy_map");
  TH1D* h_B0_emcal_cluster_energy = (TH1D*)inFile->Get("B0_emcal_cluster_energy");

  // Occupancies - Roman Pots
  TH2D* h_rp_occupancy_map = (TH2D*)inFile->Get("Roman_pots_occupancy_map");
  TH1D* h_rp_z             = (TH1D*)inFile->Get("Roman_pots_z");

  // B0 momentum resolution
  TH1D* h_b0_pt_resolution         = (TH1D*)inFile->Get("b0_pt_resolution");
  TH2D* h_b0_pt_resolution_percent = (TH2D*)inFile->Get("b0_deltaPt_over_pt_vs_pt");	
  TH2D* h_b0_p_resolution_percent  = (TH2D*)inFile->Get("b0_deltaP_over_p_vs_p");	
  TH1D* h_b0_extracted_pt_resolution = (TH1D*)inFile->Get("b0_extracted_pt_resolution");
  TH1D* h_b0_extracted_p_resolution  = (TH1D*)inFile->Get("b0_extracted_p_resolution");

  // Exclusivity variables
  TH1D* h_Emiss3_MC = (TH1D*)inFile->Get("Emiss3_MC");
  TH1D* h_Pmiss3_MC = (TH1D*)inFile->Get("Pmiss3_MC");
  TH1D* h_Ptmiss3_MC = (TH1D*)inFile->Get("Ptmiss3_MC");
  TH1D* h_M2miss3_MC = (TH1D*)inFile->Get("M2miss3_MC");
  TH1D* h_Emiss3_MCA = (TH1D*)inFile->Get("Emiss3_MCA");
  TH1D* h_Pmiss3_MCA = (TH1D*)inFile->Get("Pmiss3_MCA");
  TH1D* h_Ptmiss3_MCA = (TH1D*)inFile->Get("Ptmiss3_MCA");
  TH1D* h_M2miss3_MCA = (TH1D*)inFile->Get("M2miss3_MCA");
  TH1D* h_Emiss3_RP = (TH1D*)inFile->Get("Emiss3_RP");
  TH1D* h_Pmiss3_RP = (TH1D*)inFile->Get("Pmiss3_RP");
  TH1D* h_Ptmiss3_RP = (TH1D*)inFile->Get("Ptmiss3_RP");
  TH1D* h_M2miss3_RP = (TH1D*)inFile->Get("M2miss3_RP");
  TH1D* h_M2miss2ep_MC = (TH1D*)inFile->Get("M2miss2ep_MC");
  TH1D* h_M2miss2ep_MCA = (TH1D*)inFile->Get("M2miss2ep_MCA");
  TH1D* h_M2miss2ep_RP = (TH1D*)inFile->Get("M2miss2ep_RP");
  TH1D* h_M2miss2eg_MC = (TH1D*)inFile->Get("M2miss2eg_MC");
  TH1D* h_M2miss2eg_MCA = (TH1D*)inFile->Get("M2miss2eg_MCA");
  TH1D* h_M2miss2eg_RP = (TH1D*)inFile->Get("M2miss2eg_RP");
  
  //--------------------------------------------------------------------------------------
  // Set histogram line colours: MC raw - black; MC associated - red; Reconstructed - blue
  //--------------------------------------------------------------------------------------
  // MC generated
  h_eta_MCp->SetLineColor(kBlack);
  h_pt_MCp->SetLineColor(kBlack);
  h_p_MCp->SetLineColor(kBlack);
  h_theta_MCp->SetLineColor(kBlack);
  h_phi_MCp->SetLineColor(kBlack);
  h_E_MCp->SetLineColor(kBlack);
  h_eta_MCe->SetLineColor(kBlack);
  h_pt_MCe->SetLineColor(kBlack);
  h_p_MCe->SetLineColor(kBlack);
  h_theta_MCe->SetLineColor(kBlack);
  h_phi_MCe->SetLineColor(kBlack);
  h_E_MCe->SetLineColor(kBlack);
  h_eta_MCg->SetLineColor(kBlack);
  h_pt_MCg->SetLineColor(kBlack);
  h_p_MCg->SetLineColor(kBlack);
  h_theta_MCg->SetLineColor(kBlack);
  h_phi_MCg->SetLineColor(kBlack);
  h_E_MCg->SetLineColor(kBlack);
  h_t_MC->SetLineColor(kBlack);
  h_Q2_MC->SetLineColor(kBlack);
  h_xB_MC->SetLineColor(kBlack);
  h_y_MC->SetLineColor(kBlack);
  h_TPhi_MC->SetLineColor(kBlack);
  h_dphi_MCeg->SetLineColor(kBlack);
  h_dtheta_MCeg->SetLineColor(kBlack);
  h_Emiss3_MC->SetLineColor(kBlack);
  h_Pmiss3_MC->SetLineColor(kBlack);
  h_Ptmiss3_MC->SetLineColor(kBlack);
  h_M2miss3_MC->SetLineColor(kBlack);
  h_M2miss2ep_MC->SetLineColor(kBlack);
  h_M2miss2eg_MC->SetLineColor(kBlack);
  h_QPQG_MC->SetLineColor(kBlack);
  h_Cone_MC->SetLineColor(kBlack);
  // Associated MC
  h_eta_MCAp->SetLineColor(kRed);
  h_pt_MCAp->SetLineColor(kRed);
  h_p_MCAp->SetLineColor(kRed);
  h_theta_MCAp->SetLineColor(kRed);
  h_phi_MCAp->SetLineColor(kRed);
  h_E_MCAp->SetLineColor(kRed);
  h_eta_MCAe->SetLineColor(kRed);
  h_pt_MCAe->SetLineColor(kRed);
  h_p_MCAe->SetLineColor(kRed);
  h_theta_MCAe->SetLineColor(kRed);
  h_phi_MCAe->SetLineColor(kRed);
  h_E_MCAe->SetLineColor(kRed);
  h_eta_MCAg->SetLineColor(kRed);
  h_pt_MCAg->SetLineColor(kRed);
  h_p_MCAg->SetLineColor(kRed);
  h_theta_MCAg->SetLineColor(kRed);
  h_phi_MCAg->SetLineColor(kRed);
  h_E_MCAg->SetLineColor(kRed);
  h_t_MCA->SetLineColor(kRed);
  h_Q2_MCA->SetLineColor(kRed);
  h_xB_MCA->SetLineColor(kRed);
  h_y_MCA->SetLineColor(kRed);
  h_TPhi_MCA->SetLineColor(kRed);
  h_dphi_MCAeg->SetLineColor(kRed);
  h_dtheta_MCAeg->SetLineColor(kRed);
  h_Emiss3_MCA->SetLineColor(kRed);
  h_Pmiss3_MCA->SetLineColor(kRed);
  h_Ptmiss3_MCA->SetLineColor(kRed);
  h_M2miss3_MCA->SetLineColor(kRed);
  h_M2miss2ep_MCA->SetLineColor(kRed);
  h_M2miss2eg_MCA->SetLineColor(kRed);
  h_QPQG_MCA->SetLineColor(kRed);
  h_Cone_MCA->SetLineColor(kRed);
  // Reconstructed
  h_eta_RPp->SetLineColor(kBlue);
  h_pt_RPp->SetLineColor(kBlue);
  h_p_RPp->SetLineColor(kBlue);
  h_theta_RPp->SetLineColor(kBlue);
  h_phi_RPp->SetLineColor(kBlue);
  h_E_RPp->SetLineColor(kBlue);
  h_eta_RPPp->SetLineColor(kCyan);
  h_pt_RPPp->SetLineColor(kCyan);
  h_p_RPPp->SetLineColor(kCyan);
  h_theta_RPPp->SetLineColor(kCyan);
  h_phi_RPPp->SetLineColor(kCyan);
  h_E_RPPp->SetLineColor(kCyan);
  h_eta_RPe->SetLineColor(kBlue);
  h_pt_RPe->SetLineColor(kBlue);
  h_p_RPe->SetLineColor(kBlue);
  h_theta_RPe->SetLineColor(kBlue);
  h_phi_RPe->SetLineColor(kBlue);
  h_E_RPe->SetLineColor(kBlue);
  h_eta_RPg->SetLineColor(kBlue);
  h_pt_RPg->SetLineColor(kBlue);
  h_p_RPg->SetLineColor(kBlue);
  h_theta_RPg->SetLineColor(kBlue);
  h_phi_RPg->SetLineColor(kBlue);
  h_E_RPg->SetLineColor(kBlue);
  h_t_RP->SetLineColor(kBlue);
  h_t_RPP->SetLineColor(kCyan);
  h_Q2_RP->SetLineColor(kBlue);
  h_xB_RP->SetLineColor(kBlue);
  h_y_RP->SetLineColor(kBlue);
  h_TPhi_RP->SetLineColor(kBlue);
  h_TPhi_RPP->SetLineColor(kCyan);
  h_dphi_RPeg->SetLineColor(kBlue);
  h_dtheta_RPeg->SetLineColor(kBlue);
  h_Emiss3_RP->SetLineColor(kBlue);
  h_Pmiss3_RP->SetLineColor(kBlue);
  h_Ptmiss3_RP->SetLineColor(kBlue);
  h_M2miss3_RP->SetLineColor(kBlue);
  h_M2miss2ep_RP->SetLineColor(kBlue);
  h_M2miss2eg_RP->SetLineColor(kBlue);
  h_QPQG_RP->SetLineColor(kBlue);
  h_Cone_RP->SetLineColor(kBlue);
  h_Cone_RPP->SetLineColor(kCyan);

  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Draw histograms
  gStyle->SetOptStat(00000000);

  // CANVAS 1: DVCS EVENT KINEMATICS - MC AND RECONSTRUCTED
  TCanvas* c1 = new TCanvas("c1","",1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  h_t_MC->SetMinimum(1);
  gPad->SetLogy(1);
  h_t_MC->GetXaxis()->SetTitle("|t| [(GeV/#it{c})^{2}]");
  h_t_MC->Draw();
  h_t_RP->Draw("same");
  h_t_RPP->Draw("same");
  TLegend* lC1P1 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC1P1->SetLineColorAlpha(kWhite,0);
  lC1P1->SetFillColorAlpha(kWhite,0);
  lC1P1->AddEntry(h_t_MC, "MC gen.", "l");
  lC1P1->AddEntry(h_t_RP, "Reco. (B0)", "l");
  lC1P1->AddEntry(h_t_RPP, "Reco. (RP)", "l");
  lC1P1->Draw();
  c1->cd(2);
  h_Q2_MC->GetXaxis()->SetTitle("Q^{2} [(GeV/#it{c})^{2}]");
  h_Q2_MC->Draw();
  h_Q2_RP->Draw("same");
  TLegend* lC1P2 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC1P2->SetLineColorAlpha(kWhite,0);
  lC1P2->SetFillColorAlpha(kWhite,0);
  lC1P2->AddEntry(h_Q2_MC, "MC gen.", "l");
  lC1P2->AddEntry(h_Q2_RP, "Reco", "l");
  lC1P2->Draw();
  c1->cd(3);
  //h_xB_MC->SetMaximum(1.1*h_xB_RP->GetMaximum());
  h_xB_MC->GetXaxis()->SetTitle("log(x_{B})");
  h_xB_MC->Draw();
  h_xB_RP->Draw("same");
  TLegend* lC1P3 = new TLegend(0.15, 0.7, 0.32, 0.87);
  lC1P3->SetLineColorAlpha(kWhite,0);
  lC1P3->SetFillColorAlpha(kWhite,0);
  lC1P3->AddEntry(h_xB_MC, "MC gen.", "l");
  lC1P3->AddEntry(h_xB_RP, "Reco", "l");
  lC1P3->Draw();
  c1->cd(4);
  h_y_MC->GetXaxis()->SetTitle("y");
  h_y_MC->Draw();
  h_y_RP->Draw("same");
  TLegend* lC1P4 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC1P4->SetLineColorAlpha(kWhite,0);
  lC1P4->SetFillColorAlpha(kWhite,0);
  lC1P4->AddEntry(h_y_MC, "MC gen.", "l");
  lC1P4->AddEntry(h_y_RP, "Reco", "l");
  lC1P4->Draw();
  c1->cd(5);
  //h_TPhi_MC->SetMaximum(1.1*h_TPhi_RP->GetMaximum());
  h_TPhi_MC->SetMinimum(1);
  gPad->SetLogy(1);
  h_TPhi_MC->GetXaxis()->SetTitle("#phi_{h} [rad]");
  h_TPhi_MC->Draw();
  h_TPhi_RP->Draw("same");
  h_TPhi_RPP->Draw("same");
  TLegend* lC1P5 = new TLegend(0.15, 0.7, 0.32, 0.87);
  lC1P5->SetLineColorAlpha(kWhite,0);
  lC1P5->SetFillColorAlpha(kWhite,0);
  lC1P5->AddEntry(h_TPhi_MC, "MC gen.", "l");
  lC1P5->AddEntry(h_TPhi_RP, "Reco. (B0)", "l");
  lC1P5->AddEntry(h_TPhi_RPP, "Reco. (RP)", "l");
  lC1P5->Draw();

  c1->Print("B0benchmark_temp010.pdf");
  c1->Close();

  // CANVAS 2: EXCLUSIVITY VARIABLES - MISSING KINEMATCS; FULL FINAL STATE
  TCanvas* c2 = new TCanvas("c2","",1200,800);
  c2->Divide(2,2);
  c2->cd(1);
  gPad->SetLogy();
  h_Emiss3_MC->GetXaxis()->SetTitle("E_{miss, e'p'#gamma} [GeV]");
  h_Emiss3_MC->SetMinimum(1);
  h_Emiss3_MC->Draw();
  h_Emiss3_MCA->Draw("same");
  h_Emiss3_RP->Draw("same");
  c2->cd(2);
  gPad->SetLogy();
  h_Pmiss3_MC->GetXaxis()->SetTitle("p_{miss, e'p'#gamma} [GeV/#it{c}]");
  h_Pmiss3_MC->SetMinimum(1);
  h_Pmiss3_MC->Draw();
  h_Pmiss3_MCA->Draw("same");
  h_Pmiss3_RP->Draw("same");
  c2->cd(3);
  gPad->SetLogy();
  h_Ptmiss3_MC->GetXaxis()->SetTitle("p_{T, miss, e'p'#gamma} [GeV/#it{c}]");
  h_Ptmiss3_MC->SetMinimum(1);
  h_Ptmiss3_MC->Draw();
  h_Ptmiss3_MCA->Draw("same");
  h_Ptmiss3_RP->Draw("same");
  c2->cd(4);
  gPad->SetLogy();
  h_M2miss3_MC->GetXaxis()->SetTitle("M^{2}_{miss, e'p'#gamma} [(GeV/#it{c}^{2})^{2}]");
  h_M2miss3_MC->SetMinimum(1);
  h_M2miss3_MC->Draw();
  h_M2miss3_MCA->Draw("same");
  h_M2miss3_RP->Draw("same");
  TLegend* lC2 = new TLegend(0.15,0.4,0.75,0.87);
  lC2->SetLineColorAlpha(kWhite,0);
  lC2->SetFillColorAlpha(kWhite,0);
  lC2->AddEntry(h_Emiss3_MC, "Generated MC", "l");
  lC2->AddEntry(h_Emiss3_MCA, "Associated MC", "l");
  lC2->AddEntry(h_Emiss3_RP, "Reconstructed", "l");
  //lC2->Draw();
    
  c2->Print("B0benchmark_temp020.pdf");
  c2->Close();

  // CANVAS 2A: EXCLUSIVITY VARIABLES - PARTIAL FINAL STATE
  TCanvas* c2a = new TCanvas("c2a","",1200,800);
  c2a->Divide(2,1);
  c2a->cd(1);
  gPad->SetLogy();
  h_M2miss2ep_MC->GetXaxis()->SetTitle("M^{2}_{miss, e'p'} [(GeV/#it{c}^{2})^{2}]");
  h_M2miss2ep_MC->SetMinimum(1);
  h_M2miss2ep_MC->Draw();
  h_M2miss2ep_MCA->Draw("same");
  h_M2miss2ep_RP->Draw("same");
  c2a->cd(2);
  gPad->SetLogy();
  h_M2miss2eg_MC->GetXaxis()->SetTitle("M^{2}_{miss, e'#gamma} [(GeV/#it{c}^{2})^{2}]");
  h_M2miss2eg_MC->SetMinimum(1);
  h_M2miss2eg_MC->Draw();
  h_M2miss2eg_MCA->Draw("same");
  h_M2miss2eg_RP->Draw("same");
  
  c2a->Print("B0benchmark_temp021.pdf");
  c2a->Close();

  // CANVAS 2B: EXCLUSIVITY VARIABLES - CONE ANGLE AND QP-QG PLANE ANGLE
  TCanvas* c2b = new TCanvas("c2b","",1200,800);
  c2b->Divide(2,1);
  c2b->cd(1);
  gPad->SetLogy();
  h_QPQG_MC->Draw();
  h_QPQG_MCA->Draw("same");
  h_QPQG_RP->Draw("same");
  TLegend* lC2BP1 = new TLegend(0.7,0.77,0.87,0.87);
  lC2BP1->SetLineColorAlpha(kWhite,0);
  lC2BP1->SetFillColorAlpha(kWhite,0);
  lC2BP1->AddEntry(h_QPQG_MC,"MC gen.","l");
  lC2BP1->AddEntry(h_QPQG_MCA,"MC assoc.","l");
  lC2BP1->AddEntry(h_QPQG_RP,"Reco.","l");
  lC2BP1->Draw();
  c2b->cd(2);
  gPad->SetLogy();
  h_Cone_MC->Draw();
  h_Cone_MC->SetMinimum(1);
  h_Cone_MCA->Draw("same");
  h_Cone_RP->Draw("same");
  h_Cone_RPP->Draw("same");
  TLegend* lC2BP2 = new TLegend(0.7,0.7,0.87,0.87);
  lC2BP2->SetLineColorAlpha(kWhite,0);
  lC2BP2->SetFillColorAlpha(kWhite,0);
  lC2BP2->AddEntry(h_Cone_MC,"MC gen.","l");
  lC2BP2->AddEntry(h_Cone_MCA,"MC assoc.","l");
  lC2BP2->AddEntry(h_Cone_RP,"Reco. (B0)","l");
  lC2BP2->AddEntry(h_Cone_RPP,"Reco. (RP)","l");
  lC2BP2->Draw();

  c2b->Print("B0benchmark_temp022.pdf");
  c2b->Close();

  // CANVAS 3: DETECTOR OCCUPANCIES - MC
  TCanvas* c3 = new TCanvas("c3","",1200,800);
  gPad->SetLogy();
  TH1F* h_deteta_MCe = (TH1F*)h_eta_MCe->Clone("h_deteta_MCe");
  TH1F* h_deteta_MCg = (TH1F*)h_eta_MCg->Clone("h_deteta_MCg");
  TH1F* h_deteta_MCp = (TH1F*)h_eta_MCp->Clone("h_deteta_MCp");
  h_deteta_MCe->SetLineColor(kRed);
  h_deteta_MCg->SetLineColor(kGreen);
  h_deteta_MCp->SetLineColor(kBlue);
  h_deteta_MCe->GetXaxis()->SetTitle("#eta_{MC}");
  h_deteta_MCe->Draw("hist");
  h_deteta_MCg->Draw("hist same");
  h_deteta_MCp->Draw("hist same");
  TLegend* lC3 = new TLegend(0.15,0.7,0.3,0.87);
  lC3->SetLineColorAlpha(kWhite,0);
  lC3->SetFillColorAlpha(kWhite,0);
  lC3->AddEntry(h_deteta_MCe, "e'", "l");
  lC3->AddEntry(h_deteta_MCg, "#gamma", "l");
  lC3->AddEntry(h_deteta_MCp, "p'", "l");
  lC3->Draw();
  
  c3->Print("B0benchmark_temp03.pdf");
  c3->Close();

  // CANVAS 3A: DETECTOR OCCUPANCIES - Reco
  TCanvas* c3a = new TCanvas("c3a","",1200,800);
  gPad->SetLogy();
  TH1F* h_deteta_RPe = (TH1F*)h_eta_RPe->Clone("h_deteta_RPe");
  TH1F* h_deteta_RPg = (TH1F*)h_eta_RPg->Clone("h_deteta_RPg");
  TH1F* h_deteta_RPp = (TH1F*)h_eta_RPp->Clone("h_deteta_RPp");
  TH1F* h_deteta_RPPp = (TH1F*)h_eta_RPPp->Clone("h_deteta_RPPp");
  h_deteta_RPe->SetLineColor(kRed);
  h_deteta_RPg->SetLineColor(kGreen);
  h_deteta_RPp->SetLineColor(kBlue);
  h_deteta_RPPp->SetLineColor(kCyan);
  h_deteta_RPe->GetXaxis()->SetTitle("#eta_{Reco.}");
  h_deteta_RPe->Draw("hist");
  h_deteta_RPg->Draw("hist same");
  h_deteta_RPp->Draw("hist same");
  h_deteta_RPPp->Draw("hist same");
  TLegend* lC3A = new TLegend(0.15,0.7,0.3,0.87);
  lC3A->SetLineColorAlpha(kWhite,0);
  lC3A->SetFillColorAlpha(kWhite,0);
  lC3A->AddEntry(h_deteta_RPe, "e'", "l");
  lC3A->AddEntry(h_deteta_RPg, "#gamma", "l");
  lC3A->AddEntry(h_deteta_RPp, "p' (B0)", "l");
  lC3A->AddEntry(h_deteta_RPPp, "p' (RP)", "l");
  lC3A->Draw();
  
  c3a->Print("B0benchmark_temp031.pdf");
  c3a->Close();

  // CANVAS 4: SCATTERED ELECTRON KINEMATICS
  TCanvas* c4 = new TCanvas("c4","",1200,800);
  c4->Divide(3,2);
  c4->cd(1);
  h_eta_MCe->SetMinimum(1);
  gPad->SetLogy(1);
  h_eta_MCe->GetXaxis()->SetTitle("#eta_{e'}");
  h_eta_MCe->Draw();
  h_eta_MCAe->Draw("same");
  h_eta_RPe->Draw("same");
  TLegend* lC4P1 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4P1->SetLineColorAlpha(kWhite,0);
  lC4P1->SetFillColorAlpha(kWhite,0);
  lC4P1->AddEntry(h_eta_MCe, "MC gen.", "l");
  lC4P1->AddEntry(h_eta_MCAe, "MC assoc.", "l");
  lC4P1->AddEntry(h_eta_RPe, "Reco", "l");
  lC4P1->Draw();
  c4->cd(2);
  h_theta_MCe->SetMinimum(1);
  gPad->SetLogy(1);
  h_theta_MCe->GetXaxis()->SetTitle("#theta_{e'} [rad]");
  h_theta_MCe->Draw();
  h_theta_MCAe->Draw("same");
  h_theta_RPe->Draw("same");
  TLegend* lC4P2 = new TLegend(0.12, 0.7, 0.29, 0.87);
  lC4P2->SetLineColorAlpha(kWhite,0);
  lC4P2->SetFillColorAlpha(kWhite,0);
  lC4P2->AddEntry(h_theta_MCe, "MC gen.", "l");
  lC4P2->AddEntry(h_theta_MCAe, "MC assoc.", "l");
  lC4P2->AddEntry(h_theta_RPe, "Reco", "l");
  lC4P2->Draw();
  c4->cd(3);
  h_phi_MCe->SetMinimum(1);
  gPad->SetLogy(1);
  h_phi_MCe->GetXaxis()->SetTitle("#phi_{e'} [rad]");
  h_phi_MCe->Draw();
  h_phi_MCAe->Draw("same");
  h_phi_RPe->Draw("same");
  TLegend* lC4P3 = new TLegend(0.7, 0.12, 0.87, 0.29);
  lC4P3->SetLineColorAlpha(kWhite,0);
  lC4P3->SetFillColorAlpha(kWhite,0);
  lC4P3->AddEntry(h_phi_MCe, "MC gen.", "l");
  lC4P3->AddEntry(h_phi_RPe, "Reco", "l");
  //lC4P3->Draw();
  c4->cd(4);
  h_p_MCe->SetMinimum(1);
  gPad->SetLogy(1);
  h_p_MCe->GetXaxis()->SetTitle("p_{e'} [GeV/#it{c}]");
  h_p_MCe->Draw();
  h_p_RPe->Draw("same");
  h_p_MCAe->Draw("same");
  TLegend* lC4P4 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4P4->SetLineColorAlpha(kWhite,0);
  lC4P4->SetFillColorAlpha(kWhite,0);
  lC4P4->AddEntry(h_p_MCe, "MC gen.", "l");
  lC4P4->AddEntry(h_p_MCAe, "MC assoc.", "l");
  lC4P4->AddEntry(h_p_RPe, "Reco", "l");
  lC4P4->Draw();
  c4->cd(5);
  h_pt_MCe->GetXaxis()->SetTitle("p_{T, e'} [GeV/#it{c}]");
  h_pt_MCe->Draw();
  h_pt_RPe->Draw("same");
  h_pt_MCAe->Draw("same");
  TLegend* lC4P5 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4P5->SetLineColorAlpha(kWhite,0);
  lC4P5->SetFillColorAlpha(kWhite,0);
  lC4P5->AddEntry(h_pt_MCe, "MC gen.", "l");
  lC4P5->AddEntry(h_pt_MCAe, "MC assoc.", "l");
  lC4P5->AddEntry(h_pt_RPe, "Reco", "l");
  lC4P5->Draw();
  c4->cd(6);
  gPad->SetLogy(1);
  h_E_MCe->GetXaxis()->SetTitle("E_{e'} [GeV]");
  h_E_MCe->Draw();
  h_E_MCAe->Draw("same");
  h_E_RPe->Draw("same");
  TLegend* lC4P6 = new TLegend(0.12, 0.7, 0.29, 0.87);
  lC4P6->SetLineColorAlpha(kWhite,0);
  lC4P6->SetFillColorAlpha(kWhite,0);
  lC4P6->AddEntry(h_E_MCe, "MC gen.", "l");
  lC4P6->AddEntry(h_E_MCAe, "MC assoc.", "l");
  lC4P6->AddEntry(h_E_RPe, "Reco", "l");
  lC4P6->Draw();
  
  c4->Print("B0benchmark_temp040.pdf");
  c4->Close();

  // CANVAS 4A: SCATTERED PROTON KINEMATICS
  TCanvas* c4a = new TCanvas("c4a","",1200,800);
  c4a->Divide(3,2);
  c4a->cd(1);
  h_eta_MCp->SetMinimum(1);
  gPad->SetLogy(1);
  h_eta_MCp->GetXaxis()->SetTitle("#eta_{p'}");
  h_eta_MCp->Draw();
  h_eta_MCAp->Draw("same");
  h_eta_RPp->Draw("same");
  h_eta_RPPp->Draw("same");
  TLegend* lC4AP1 = new TLegend(0.12, 0.7, 0.29, 0.87);
  lC4AP1->SetLineColorAlpha(kWhite,0);
  lC4AP1->SetFillColorAlpha(kWhite,0);
  lC4AP1->AddEntry(h_eta_MCp, "MC gen.", "l");
  lC4AP1->AddEntry(h_eta_MCAp, "MC assoc.", "l");
  lC4AP1->AddEntry(h_eta_RPp, "Reco. (B0)", "l");
  lC4AP1->AddEntry(h_eta_RPPp, "Reco. (RP)", "l");
  lC4AP1->Draw();
  c4a->cd(2);
  h_theta_MCp->SetMinimum(1);
  gPad->SetLogy(1);
  h_theta_MCp->GetXaxis()->SetTitle("#theta_{p'} [mrad]");
  h_theta_MCp->Draw();
  h_theta_MCAp->Draw("same");
  h_theta_RPp->Draw("same");
  h_theta_RPPp->Draw("same");
  TLegend* lC4AP2 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4AP2->SetLineColorAlpha(kWhite,0);
  lC4AP2->SetFillColorAlpha(kWhite,0);
  lC4AP2->AddEntry(h_theta_MCp, "MC gen.", "l");
  lC4AP2->AddEntry(h_theta_MCAp, "MC assoc.", "l");
  lC4AP2->AddEntry(h_theta_RPp, "Reco. (B0)", "l");
  lC4AP2->AddEntry(h_theta_RPPp, "Reco. (RP)", "l");
  lC4AP2->Draw();
  c4a->cd(3);
  h_phi_MCp->SetMinimum(1);
  gPad->SetLogy(1);
  h_phi_MCp->GetXaxis()->SetTitle("#phi_{p'} [rad]");
  h_phi_MCp->Draw();
  h_phi_MCAp->Draw("same");
  h_phi_RPp->Draw("same");
  h_phi_RPPp->Draw("same");
  TLegend* lC4AP3 = new TLegend(0.7, 0.12, 0.87, 0.29);
  lC4AP3->SetLineColorAlpha(kWhite,0);
  lC4AP3->SetFillColorAlpha(kWhite,0);
  lC4AP3->AddEntry(h_phi_MCp, "MC gen.", "l");
  lC4AP3->AddEntry(h_phi_RPp, "Reco. (B0)", "l");
  lC4AP3->AddEntry(h_phi_RPPp, "Reco. (RP)", "l");
  //lC4AP3->Draw();
  c4a->cd(4);
  h_p_MCp->SetMinimum(1);
  gPad->SetLogy(1);
  h_p_MCp->GetXaxis()->SetTitle("p_{p'} [GeV/#it{c}]");
  h_p_MCp->Draw();
  h_p_RPp->Draw("same");
  h_p_RPPp->Draw("same");
  h_p_MCAp->Draw("same");
  TLegend* lC4AP4 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4AP4->SetLineColorAlpha(kWhite,0);
  lC4AP4->SetFillColorAlpha(kWhite,0);
  lC4AP4->AddEntry(h_p_MCp, "MC gen.", "l");
  lC4AP4->AddEntry(h_p_MCAp, "MC assoc.", "l");
  lC4AP4->AddEntry(h_p_RPp, "Reco. (B0)", "l");
  lC4AP4->AddEntry(h_p_RPPp, "Reco. (RP)", "l");
  lC4AP4->Draw();
  c4a->cd(5);
  h_pt_MCp->SetMinimum(1);
  gPad->SetLogy(1);
  h_pt_MCp->GetXaxis()->SetTitle("p_{T, p'} [GeV/#it{c}]");
  h_pt_MCp->Draw();
  h_pt_RPp->Draw("same"); 
  h_pt_RPPp->Draw("same"); 
  h_pt_MCAp->Draw("same"); 
  TLegend* lC4AP5 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4AP5->SetLineColorAlpha(kWhite,0);
  lC4AP5->SetFillColorAlpha(kWhite,0);
  lC4AP5->AddEntry(h_pt_MCp, "MC gen.", "l");
  lC4AP5->AddEntry(h_pt_MCAp, "MC assoc.", "l");
  lC4AP5->AddEntry(h_pt_RPp, "Reco. (B0)", "l");
  lC4AP5->AddEntry(h_pt_RPPp, "Reco. (RP)", "l");
  lC4AP5->Draw();
  c4a->cd(6);
  gPad->SetLogy(1);
  h_E_MCp->GetXaxis()->SetTitle("E_{p'} [GeV]");
  h_E_MCp->Draw();
  h_E_MCAp->Draw("same");
  h_E_RPp->Draw("same");
  h_E_RPPp->Draw("same");
  TLegend* lC4AP6 = new TLegend(0.12, 0.7, 0.29, 0.87);
  lC4AP6->SetLineColorAlpha(kWhite,0);
  lC4AP6->SetFillColorAlpha(kWhite,0);
  lC4AP6->AddEntry(h_pt_MCp, "MC gen.", "l");
  lC4AP6->AddEntry(h_pt_MCAp, "MC assoc.", "l");
  lC4AP6->AddEntry(h_pt_RPp, "Reco. (B0)", "l");
  lC4AP6->AddEntry(h_pt_RPPp, "Reco. (RP)", "l");
  lC4AP6->Draw();

  c4a->Print("B0benchmark_temp041.pdf");
  c4a->Close();

  // CANVAS 4B: FINAL STATE PHOTON KINEMATICS
  TCanvas* c4b = new TCanvas("c4b","",1200,800);
  c4b->Divide(3,2);
  c4b->cd(1);
  h_eta_MCg->SetMinimum(1);
  //gPad->SetLogy(1);
  h_eta_MCg->GetXaxis()->SetTitle("#eta_{#gamma}");
  h_eta_MCg->Draw();
  h_eta_MCAg->Draw("same");
  h_eta_RPg->Draw("same");
  TLegend* lC4BP1 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4BP1->SetLineColorAlpha(kWhite,0);
  lC4BP1->SetFillColorAlpha(kWhite,0);
  lC4BP1->AddEntry(h_eta_MCg, "MC gen.", "l");
  lC4BP1->AddEntry(h_eta_MCAg, "MC assoc.", "l");
  lC4BP1->AddEntry(h_eta_RPg, "Reco", "l");
  lC4BP1->Draw();
  c4b->cd(2);
  h_theta_MCg->SetMinimum(1);
  gPad->SetLogy(1);
  h_theta_MCg->GetXaxis()->SetTitle("#theta_{#gamma} [rad]");
  h_theta_MCg->Draw();
  h_theta_MCAg->Draw("same");
  h_theta_RPg->Draw("same");
  TLegend* lC4BP2 = new TLegend(0.7, 0.12, 0.87, 0.29);
  lC4BP2->SetLineColorAlpha(kWhite,0);
  lC4BP2->SetFillColorAlpha(kWhite,0);
  lC4BP2->AddEntry(h_theta_MCg, "MC gen.", "l");
  lC4BP2->AddEntry(h_theta_MCAg, "MC assoc.", "l");
  lC4BP2->AddEntry(h_theta_RPg, "Reco", "l");
  lC4BP2->Draw();
  c4b->cd(3);
  h_phi_MCg->SetMinimum(1);
  gPad->SetLogy(1);
  h_phi_MCg->GetXaxis()->SetTitle("#phi_{#gamma} [rad]");
  h_phi_MCg->Draw();
  h_phi_MCAg->Draw("same");
  h_phi_RPg->Draw("same");
  TLegend* lC4BP3 = new TLegend(0.7, 0.12, 0.87, 0.29);
  lC4BP3->SetLineColorAlpha(kWhite,0);
  lC4BP3->SetFillColorAlpha(kWhite,0);
  lC4BP3->AddEntry(h_phi_MCg, "MC gen.", "l");
  lC4BP3->AddEntry(h_phi_RPg, "Reco", "l");
  //lC4BP3->Draw();
  c4b->cd(4);
  h_p_MCg->SetMinimum(1);
  gPad->SetLogy(1);
  h_p_MCg->GetXaxis()->SetTitle("p_{#gamma} [GeV/#it{c}]");
  h_p_MCg->Draw();
  h_p_MCAg->Draw("same");
  h_p_RPg->Draw("same");
  TLegend* lC4BP4 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4BP4->SetLineColorAlpha(kWhite,0);
  lC4BP4->SetFillColorAlpha(kWhite,0);
  lC4BP4->AddEntry(h_p_MCg, "MC gen.", "l");
  lC4BP4->AddEntry(h_p_MCAg, "MC assoc.", "l");
  lC4BP4->AddEntry(h_p_RPg, "Reco", "l");
  lC4BP4->Draw();
  c4b->cd(5);
  h_pt_MCg->SetMinimum(1);
  h_pt_MCg->GetXaxis()->SetTitle("p_{T, #gamma} [GeV/#it{c}]");
  h_pt_MCg->Draw();
  h_pt_MCAg->Draw("same");
  h_pt_RPg->Draw("same"); 
  TLegend* lC4BP5 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4BP5->SetLineColorAlpha(kWhite,0);
  lC4BP5->SetFillColorAlpha(kWhite,0);
  lC4BP5->AddEntry(h_pt_MCg, "MC gen.", "l");
  lC4BP5->AddEntry(h_pt_MCAg, "MC assoc.", "l");
  lC4BP5->AddEntry(h_pt_RPg, "Reco", "l");
  lC4BP5->Draw();
  c4b->cd(6);
  gPad->SetLogy(1);
  h_E_MCg->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
  h_E_MCg->Draw();
  h_E_MCAg->Draw("same");
  h_E_RPg->Draw("same");
  TLegend* lC4BP6 = new TLegend(0.7, 0.7, 0.87, 0.87);
  lC4BP6->SetLineColorAlpha(kWhite,0);
  lC4BP6->SetFillColorAlpha(kWhite,0);
  lC4BP6->AddEntry(h_E_MCg, "MC gen.", "l");
  lC4BP6->AddEntry(h_E_MCAg, "MC assoc.", "l");
  lC4BP6->AddEntry(h_E_RPg, "Reco", "l");
  lC4BP6->Draw();
  
  c4b->Print("B0benchmark_temp042.pdf");
  c4b->Close();

  // CANVAS 5: ANGULAR DISTANCES BETWEEN FINAL STATE PHOTON AND ELECTRON
  TCanvas* c5 = new TCanvas("c5","",1200,800);
  c5->Divide(2);
  c5->cd(1);
  h_eta_MCg->SetMinimum(1);
  //gPad->SetLogy(1);
  h_dphi_MCeg->GetXaxis()->SetTitle("#delta#phi_{e'#gamma} [rad]");
  h_dphi_MCeg->Draw();
  h_dphi_RPeg->Draw("same");  
  TLegend* lC5P1 = new TLegend(0.7, 0.75, 0.87, 0.87);
  lC5P1->SetLineColorAlpha(kWhite,0);
  lC5P1->SetFillColorAlpha(kWhite,0);
  lC5P1->AddEntry(h_dphi_MCeg, "MC gen.", "l");
  lC5P1->AddEntry(h_dphi_RPeg, "Reco", "l");
  lC5P1->Draw();
  c5->cd(2);
  h_eta_MCg->SetMinimum(1);
  //gPad->SetLogy(1);
  h_dtheta_MCeg->GetXaxis()->SetTitle("#delta#theta_{e'#gamma} [rad]");
  h_dtheta_MCeg->Draw();
  h_dtheta_RPeg->Draw("same");  
  TLegend* lC5P2 = new TLegend(0.7, 0.75, 0.87, 0.87);
  lC5P2->SetLineColorAlpha(kWhite,0);
  lC5P2->SetFillColorAlpha(kWhite,0);
  lC5P2->AddEntry(h_dtheta_MCeg, "MC gen.", "l");
  lC5P2->AddEntry(h_dtheta_RPeg, "Reco", "l");
  lC5P2->Draw();

  c5->Print("B0benchmark_temp05.pdf");
  c5->Close();

  // CANVAS 6: B0 TRACKER OCCUPANCIES AND ENERGY DEPOSITS
  TCanvas* c6 = new TCanvas("c6","",1200,800);
  c6->Divide(3,2);
  c6->cd(1);
  h_B0_occupancy_map_layer_0->SetTitle("");
  h_B0_occupancy_map_layer_0->GetYaxis()->SetTitleOffset(1.35);
  h_B0_occupancy_map_layer_0->Draw("col");
  c6->cd(2);
  h_B0_occupancy_map_layer_1->SetTitle("");
  h_B0_occupancy_map_layer_1->GetYaxis()->SetTitleOffset(1.35);
  h_B0_occupancy_map_layer_1->Draw("col");
  c6->cd(3);
  h_B0_hit_energy_deposit->Draw();
  c6->cd(4);
  h_B0_occupancy_map_layer_2->SetTitle("");
  h_B0_occupancy_map_layer_2->GetYaxis()->SetTitleOffset(1.35);
  h_B0_occupancy_map_layer_2->Draw("col");
  c6->cd(5);
  h_B0_occupancy_map_layer_3->SetTitle("");
  h_B0_occupancy_map_layer_3->GetYaxis()->SetTitleOffset(1.35);
  h_B0_occupancy_map_layer_3->Draw("col");

  c6->Print("B0benchmark_temp06.pdf");
  c6->Close();
 
  // CANVAS 7: B0 ECAL OCCUPANCY AND ENERGY DEPOSITS
  TCanvas* c7 = new TCanvas("c7","",1200,800);
  c7->Divide(2);
  c7->cd(1);
  h_B0_emcal_occupancy_map->GetYaxis()->SetTitleOffset(1.35);
  h_B0_emcal_occupancy_map->Draw();
  c7->cd(2)->SetLogy();
  h_B0_emcal_cluster_energy->Draw();

  c7->Print("B0benchmark_temp07.pdf");
  c7->Close();

  // CANVAS 8: ROMAN POT OCCUPANCIES
  TCanvas* c8 = new TCanvas("c8","",1200,800);
  h_rp_occupancy_map->Draw("colz");

  c8->Print("B0benchmark_temp08.pdf");
  c8->Close();

  // CANVAS 9 AND 9A: B0 P AND PT RESOLUTIONS
  TCanvas* c9 = new TCanvas("c9","",1200,900);
  h_b0_extracted_p_resolution->SetMinimum(0);
  h_b0_extracted_p_resolution->Draw();
  c9->Print("B0benchmark_temp09.pdf");
  c9->Close();

  TCanvas* c9a = new TCanvas("c9a","",1200,800);
  h_b0_extracted_pt_resolution->SetMinimum(0);
  h_b0_extracted_pt_resolution->Draw();
  c9a->Print("B0benchmark_temp091.pdf");
  c9a->Close();

  // CANVAS 10: RATIO BETWEEN ASSOCIATED MONTE CARLO AND RECONSTRUCTED MOMENTA
  TH1F* hAcc_p = (TH1F*)h_p_RPp->Clone("hAcc_p");
  hAcc_p->Divide(h_p_MCAp);
  hAcc_p->SetTitle(";p [GeV/#it{c}];#frac{Reco}{MC|Reco}");
  
  TCanvas* c10 = new TCanvas("c10","",1200,800);
  hAcc_p->Draw("pe");
  
  c10->Print("B0benchmark_temp100.pdf");
  c10->Close();

  // CANVAS 10A: 2D PLOT, ASSOCIATED MC VS RECONSTRUCTED PROTON MOMENTUM
  TCanvas* c10a = new TCanvas("c10a","",1200,800);
  h_p_2d->Draw();
  
  c10a->Print("B0benchmark_temp101.pdf");
  c10a->Close();
  
  // CANVAS 11: RATIO BETWEEN ASSOCIATED MONTE CARLO AND RECONSTRUCTED TRANSVERSE MOMENTA
  TH1F* hAcc_pt = (TH1F*)h_pt_RPp->Clone("hAcc_pt");
  hAcc_pt->Divide(h_pt_MCAp);  
  hAcc_pt->SetTitle(";p_{T} [GeV/#it{c}];#frac{Reco}{MC|Reco}");

  TCanvas* c11 = new TCanvas("c11","",1200,800);
  hAcc_pt->Draw("pe");

  c11->Print("B0benchmark_temp110.pdf");
  c11->Close();

  // CANVAS 11A: 2D PLOT, ASSOCIATED MC VS RECONSTRUCTED PROTON TRANSVERSE MOMENTUM
  TCanvas* c11a = new TCanvas("c11a","",1200,800);
  h_pt_2d->Draw();

  c11a->Print("B0benchmark_temp111.pdf");
  c11a->Close();

  // CANVAS 12: RATIO BETWEEN ASSOCIATED MONTE CARLO AND RECONSTRUCTED T
  TH1F* hAcc_t = (TH1F*)h_t_RP->Clone("hAcc_t");
  hAcc_t->Divide(h_t_MCA);
  hAcc_t->SetTitle(";t [GeV/#it{c}^{2}];#frac{Reco}{MC|Reco}");

  TCanvas* c12 = new TCanvas("c12","",1200,800);
  hAcc_t->Draw("pe");

  c12->Print("B0benchmark_temp120.pdf");
  c12->Close();
  
  // CANVAS 12A: 2D PLOT, ASSOCIATED MC VS RECONSTRUCTED T
  TCanvas* c12a = new TCanvas("c12a","",1200,800);
  h_t_2d->Draw();
  
  c12a->Print("B0benchmark_temp121.pdf");
  c12a->Close();

  // Combine PDFs into one and clean up
  std::cout<<"...Cleaning up files..."<<std::endl;
  TString filePlots = "$EIC_WORK_DIR/DVCS_Analysis/Plots/DVCSPlots_" + campaign + "_" + energy + "_" + setting + ".pdf";
  std::cout<<"Moving plots to "<<filePlots<<std::endl;
  TString pdfUniteCmd = "pdfunite B0benchmark*.pdf "+filePlots;
  gSystem->Exec(pdfUniteCmd);
  gSystem->Exec("rm B0benchmark_temp*.pdf");
  
  return;
}
