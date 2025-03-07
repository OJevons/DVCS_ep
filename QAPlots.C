using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

// Set histogram drawing options
void setHistOpts(TH1D* hist, TString type){
  // 1) Truth: kBlack, marker style 8 (filled circle), marker size 2
  // 2) Raw reco.: kBlue, marker style 71 (open circle, line thickness 2), marker size 2
  // 3) Corrected reco.: kGreen+3, marker style 33 (filled diamond), marker size 2

  if(type == "truth"){
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kBlack);
  }
  else if(type == "reco"){
    hist->SetLineColor(kBlue);
    hist->SetMarkerStyle(24);
    hist->SetMarkerColor(kBlue);
  }
  else if(type == "corr"){
    hist->SetLineColor(kGreen+3);
    hist->SetMarkerStyle(33);
    hist->SetMarkerColor(kGreen+3);
  }
  else{
    hist->SetLineColor(kBlack);
    hist->SetMarkerStyle(8);
    hist->SetMarkerColor(kBlack);
  }

  return;
}

//---------------------------------------------------------------------
// MAIN
//---------------------------------------------------------------------
void QAPlots(TString campaign, TString energy, TString setting){
  //---------------------------------------------------------------------
  // Get data file from DVCSAnalysis script
  //---------------------------------------------------------------------
  TString inFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_" + campaign + "_" + energy + "_QA.root";
  cout<<"Input file: "<<inFileName<<endl;
  TFile* inFile = new TFile(inFileName);

  //---------------------------------------------------------------------
  // Extract histograms
  //---------------------------------------------------------------------
  // Q2 
  TH1D* h_Q2_Truth = (TH1D*)inFile->Get("q2_truth");
  TH1D* h_Q2_Acc   = (TH1D*)inFile->Get("q2_acc");
  TH1D* h_Q2_Reco  = (TH1D*)inFile->Get("q2_reco");
  TH2D* h_Q2_Resp  = (TH2D*)inFile->Get("q2_resp");
  TH1D* h_Q2_Pur   = (TH1D*)inFile->Get("q2_pur");
  TH2D* h_dQ2vQ2   = (TH2D*)inFile->Get("dq2vq2");
  // x
  TH1D* h_xB_Truth = (TH1D*)inFile->Get("xb_truth");
  TH1D* h_xB_Acc   = (TH1D*)inFile->Get("xb_acc");
  TH1D* h_xB_Reco  = (TH1D*)inFile->Get("xb_reco");
  TH2D* h_xB_Resp  = (TH2D*)inFile->Get("xb_resp");
  TH1D* h_xB_Pur   = (TH1D*)inFile->Get("xb_pur");
  TH2D* h_dxBvxB   = (TH2D*)inFile->Get("dxbvxb");
  // y
  TH1D* h_y_Truth = (TH1D*)inFile->Get("y_truth");
  TH1D* h_y_Acc   = (TH1D*)inFile->Get("y_acc");
  TH1D* h_y_Reco  = (TH1D*)inFile->Get("y_reco");
  TH2D* h_y_Resp  = (TH2D*)inFile->Get("y_resp");
  TH1D* h_y_Pur   = (TH1D*)inFile->Get("y_pur");
  TH2D* h_dyvy    = (TH2D*)inFile->Get("dyvy");
  // Exclusive kinematic quantities
  TH1D* h_t_Truth  = (TH1D*)inFile->Get("t_truth");
  TH1D* h_t_B0Acc  = (TH1D*)inFile->Get("t_b0acc");
  TH1D* h_t_RPAcc  = (TH1D*)inFile->Get("t_rpacc");
  TH1D* h_t_B0Reco = (TH1D*)inFile->Get("t_b0reco");
  TH1D* h_t_RPReco = (TH1D*)inFile->Get("t_rpreco");
  TH2D* h_t_B0Resp = (TH2D*)inFile->Get("t_b0resp");
  TH2D* h_t_RPResp = (TH2D*)inFile->Get("t_rpresp");
  TH1D* h_t_B0Pur  = (TH1D*)inFile->Get("t_b0pur");
  TH1D* h_t_RPPur  = (TH1D*)inFile->Get("t_rppur");
  TH2D* h_B0dtvt   = (TH2D*)inFile->Get("b0dtvt");
  TH2D* h_RPdtvt   = (TH2D*)inFile->Get("rpdtvt");
  //Single particle kinematics - protons
  TH1D* h_theta_p_Truth  = (TH1D*)inFile->Get("theta_p_truth");
  TH1D* h_theta_p_B0Acc  = (TH1D*)inFile->Get("theta_p_b0acc");
  TH1D* h_theta_p_RPAcc  = (TH1D*)inFile->Get("theta_p_rpacc");
  TH1D* h_theta_p_B0Reco = (TH1D*)inFile->Get("theta_p_b0reco");
  TH1D* h_theta_p_RPReco = (TH1D*)inFile->Get("theta_p_rpreco");
  TH2D* h_theta_p_B0Resp = (TH2D*)inFile->Get("theta_p_b0resp");
  TH2D* h_theta_p_RPResp = (TH2D*)inFile->Get("theta_p_rpresp");
  TH1D* h_theta_p_B0Pur  = (TH1D*)inFile->Get("theta_p_b0pur");
  TH1D* h_theta_p_RPPur  = (TH1D*)inFile->Get("theta_p_rppur");
  TH1D* h_E_p_Truth  = (TH1D*)inFile->Get("E_p_truth");
  TH1D* h_E_p_B0Acc  = (TH1D*)inFile->Get("E_p_b0acc");
  TH1D* h_E_p_RPAcc  = (TH1D*)inFile->Get("E_p_rpacc");
  TH1D* h_E_p_B0Reco = (TH1D*)inFile->Get("E_p_b0reco");
  TH1D* h_E_p_RPReco = (TH1D*)inFile->Get("E_p_rpreco");
  TH2D* h_E_p_B0Resp = (TH2D*)inFile->Get("E_p_b0resp");
  TH2D* h_E_p_RPResp = (TH2D*)inFile->Get("E_p_rpresp");
  TH1D* h_E_p_B0Pur  = (TH1D*)inFile->Get("E_p_b0pur");
  TH1D* h_E_p_RPPur  = (TH1D*)inFile->Get("E_p_rppur");
  //Single particle kinematics - electrons
  TH1D* h_theta_e_Truth = (TH1D*)inFile->Get("theta_e_truth");
  TH1D* h_theta_e_Acc   = (TH1D*)inFile->Get("theta_e_acc");  
  TH1D* h_theta_e_Reco  = (TH1D*)inFile->Get("theta_e_reco"); 
  TH2D* h_theta_e_Resp  = (TH2D*)inFile->Get("theta_e_resp"); 
  TH1D* h_theta_e_Pur   = (TH1D*)inFile->Get("theta_e_pur");
  TH1D* h_E_e_Truth = (TH1D*)inFile->Get("E_e_truth");
  TH1D* h_E_e_Acc   = (TH1D*)inFile->Get("E_e_acc");
  TH1D* h_E_e_Reco  = (TH1D*)inFile->Get("E_e_reco");
  TH2D* h_E_e_Resp  = (TH2D*)inFile->Get("E_e_resp");
  TH1D* h_E_e_Pur   = (TH1D*)inFile->Get("E_e_pur");
  //Single particle kinematics - photons
  TH1D* h_theta_g_Truth = (TH1D*)inFile->Get("theta_g_truth");
  TH1D* h_theta_g_Acc   = (TH1D*)inFile->Get("theta_g_acc");
  TH1D* h_theta_g_Reco  = (TH1D*)inFile->Get("theta_g_reco");
  TH2D* h_theta_g_Resp  = (TH2D*)inFile->Get("theta_g_resp");
  TH1D* h_theta_g_Pur   = (TH1D*)inFile->Get("theta_g_pur");
  TH1D* h_E_g_Truth = (TH1D*)inFile->Get("E_g_truth");
  TH1D* h_E_g_Acc   = (TH1D*)inFile->Get("E_g_acc");
  TH1D* h_E_g_Reco  = (TH1D*)inFile->Get("E_g_reco");
  TH2D* h_E_g_Resp  = (TH2D*)inFile->Get("E_g_resp");
  TH1D* h_E_g_Pur   = (TH1D*)inFile->Get("E_g_pur");
  // HFS quantities
  TH1D* h_HFSSigma_Truth = (TH1D*)inFile->Get("hfssigma_truth");
  TH1D* h_HFSSigma_Acc   = (TH1D*)inFile->Get("hfssigma_acc");
  TH1D* h_HFSSigma_Reco  = (TH1D*)inFile->Get("hfssigma_reco");
  TH2D* h_HFSSigma_Resp  = (TH2D*)inFile->Get("hfssigma_resp");
  TH1D* h_HFSSigma_Pur   = (TH1D*)inFile->Get("hfssigma_pur");
  TH1D* h_HFSpT2_Truth = (TH1D*)inFile->Get("hfspt2_truth");
  TH1D* h_HFSpT2_Acc   = (TH1D*)inFile->Get("hfspt2_acc");
  TH1D* h_HFSpT2_Reco  = (TH1D*)inFile->Get("hfspt2_reco");
  TH2D* h_HFSpT2_Resp  = (TH2D*)inFile->Get("hfspt2_resp");
  TH1D* h_HFSpT2_Pur   = (TH1D*)inFile->Get("hfspt2_pur");
  TH1D* h_FullSigma_Truth = (TH1D*)inFile->Get("fullsigma_truth");
  TH1D* h_FullSigma_Acc   = (TH1D*)inFile->Get("fullsigma_acc");
  TH1D* h_FullSigma_Reco  = (TH1D*)inFile->Get("fullsigma_reco");
  TH2D* h_FullSigma_Resp  = (TH2D*)inFile->Get("fullsigma_resp");
  TH1D* h_FullSigma_Pur   = (TH1D*)inFile->Get("fullsigma_pur");
  // E/p for electron/photon
  TH1D* h_EOverp_e_Truth = (TH1D*)inFile->Get("eoverp_e_truth");
  TH1D* h_EOverp_e_Reco  = (TH1D*)inFile->Get("eoverp_e_reco");
  TH2D* h_EOverp_e_Resp  = (TH2D*)inFile->Get("eoverp_e_resp");
  TH1D* h_EOverp_g_Truth = (TH1D*)inFile->Get("eoverp_g_truth");
  TH1D* h_EOverp_g_Reco  = (TH1D*)inFile->Get("eoverp_g_reco");
  TH2D* h_EOverp_g_Resp  = (TH2D*)inFile->Get("eoverp_g_resp");

  // Duplicate reco. histograms for detector corrected
  TH1D* h_Q2_Corr        = (TH1D*)h_Q2_Reco->Clone("q2_corr");
  TH1D* h_xB_Corr        = (TH1D*)h_xB_Reco->Clone("xb_corr");
  TH1D* h_y_Corr         = (TH1D*)h_y_Reco->Clone("y_corr");
  TH1D* h_t_B0Corr       = (TH1D*)h_t_B0Reco->Clone("t_b0corr");
  TH1D* h_t_RPCorr       = (TH1D*)h_t_RPReco->Clone("t_rpcorr");
  TH1D* h_theta_p_B0Corr = (TH1D*)h_theta_p_B0Reco->Clone("theta_p_b0corr");
  TH1D* h_theta_p_RPCorr = (TH1D*)h_theta_p_RPReco->Clone("theta_p_rpcorr");
  TH1D* h_E_p_B0Corr     = (TH1D*)h_E_p_B0Reco->Clone("E_p_b0corr");
  TH1D* h_E_p_RPCorr     = (TH1D*)h_E_p_RPReco->Clone("E_p_rpcorr");
  TH1D* h_theta_e_Corr   = (TH1D*)h_theta_e_Reco->Clone("theta_e_corr"); 
  TH1D* h_E_e_Corr       = (TH1D*)h_E_e_Reco->Clone("E_e_corr");
  TH1D* h_theta_g_Corr   = (TH1D*)h_theta_g_Reco->Clone("theta_g_corr");
  TH1D* h_E_g_Corr       = (TH1D*)h_E_g_Reco->Clone("E_g_corr");
  TH1D* h_HFSSigma_Corr  = (TH1D*)h_HFSSigma_Reco->Clone("hfssigma_corr");
  TH1D* h_HFSpT2_Corr    = (TH1D*)h_HFSpT2_Reco->Clone("hfspt2_corr");
  TH1D* h_FullSigma_Corr = (TH1D*)h_FullSigma_Reco->Clone("fullsigma_corr");
  
  //--------------------------------------------------------------------------------------
  // Extract detector corrected histograms
  //--------------------------------------------------------------------------------------
  // 1. Divide MC associated by MC generated - gets efficiency
  h_Q2_Acc->Divide(h_Q2_Truth);
  h_xB_Acc->Divide(h_xB_Truth);
  h_y_Acc->Divide(h_y_Truth);
  h_t_B0Acc->Divide(h_t_Truth);
  h_t_RPAcc->Divide(h_t_Truth);
  h_theta_p_B0Acc->Divide(h_theta_p_Truth);
  h_theta_p_RPAcc->Divide(h_theta_p_Truth);
  h_E_p_B0Acc->Divide(h_E_p_Truth);
  h_E_p_RPAcc->Divide(h_E_p_Truth);
  h_theta_e_Acc->Divide(h_theta_e_Truth);
  h_E_e_Acc->Divide(h_E_e_Truth);
  h_theta_g_Acc->Divide(h_theta_g_Truth);
  h_E_g_Acc->Divide(h_E_g_Truth);
  h_HFSSigma_Acc->Divide(h_HFSSigma_Truth);
  h_HFSpT2_Acc->Divide(h_HFSpT2_Truth);
  h_FullSigma_Acc->Divide(h_FullSigma_Truth);
  
  // 2. Divide (cloned) reconstructed by efficiency
  h_Q2_Corr->Divide(h_Q2_Acc);
  h_xB_Corr->Divide(h_xB_Acc);
  h_y_Corr->Divide(h_y_Acc);
  h_t_B0Corr->Divide(h_t_B0Acc);
  h_t_RPCorr->Divide(h_t_RPAcc);
  h_theta_p_B0Corr->Divide(h_theta_p_B0Acc);
  h_theta_p_RPCorr->Divide(h_theta_p_RPAcc);
  h_E_p_B0Corr->Divide(h_E_p_B0Acc);
  h_E_p_RPCorr->Divide(h_E_p_RPAcc);
  h_theta_e_Corr->Divide(h_theta_e_Acc);
  h_E_e_Corr->Divide(h_E_e_Acc);
  h_theta_g_Corr->Divide(h_theta_g_Acc);
  h_E_g_Corr->Divide(h_E_g_Acc);
  h_HFSSigma_Corr->Divide(h_HFSSigma_Acc);
  h_HFSpT2_Corr->Divide(h_HFSpT2_Acc);
  h_FullSigma_Corr->Divide(h_FullSigma_Acc);

  //--------------------------------------------------------------------------------------
  // Set histogram draw options
  //--------------------------------------------------------------------------------------
  // Q2
  setHistOpts(h_Q2_Truth, "truth");
  setHistOpts(h_Q2_Reco, "reco");
  setHistOpts(h_Q2_Corr, "corr");
  setHistOpts(h_Q2_Pur, "reco");
  // x
  setHistOpts(h_xB_Truth, "truth");
  setHistOpts(h_xB_Reco, "reco");
  setHistOpts(h_xB_Corr, "corr");
  setHistOpts(h_xB_Pur, "reco");
  // y
  setHistOpts(h_y_Truth, "truth");
  setHistOpts(h_y_Reco, "reco");
  setHistOpts(h_y_Corr, "corr");
  setHistOpts(h_y_Pur, "reco");
  // Exclusive kinematic quantities
  setHistOpts(h_t_Truth, "truth");
  setHistOpts(h_t_B0Reco, "reco");
  setHistOpts(h_t_RPReco, "reco");
  setHistOpts(h_t_B0Corr, "corr");
  setHistOpts(h_t_RPCorr, "corr");
  setHistOpts(h_t_B0Pur, "reco");
  setHistOpts(h_t_RPPur, "reco");
  //Single particle kinematics - protons
  setHistOpts(h_theta_p_Truth, "truth");
  setHistOpts(h_theta_p_B0Reco, "reco");
  setHistOpts(h_theta_p_RPReco, "reco");
  setHistOpts(h_theta_p_B0Corr, "corr");
  setHistOpts(h_theta_p_RPCorr, "corr");
  setHistOpts(h_theta_p_B0Pur, "reco");
  setHistOpts(h_theta_p_RPPur, "reco");
  setHistOpts(h_E_p_Truth, "truth");
  setHistOpts(h_E_p_B0Reco, "reco");
  setHistOpts(h_E_p_RPReco, "reco");
  setHistOpts(h_E_p_B0Corr, "corr");
  setHistOpts(h_E_p_RPCorr, "corr");
  setHistOpts(h_E_p_B0Pur, "reco");
  setHistOpts(h_E_p_RPPur, "reco");
  //Single particle kinematics - electrons
  setHistOpts(h_theta_e_Truth, "truth");
  setHistOpts(h_theta_e_Reco, "reco");
  setHistOpts(h_theta_e_Corr, "corr");
  setHistOpts(h_theta_e_Pur, "reco");
  setHistOpts(h_E_e_Truth, "truth");
  setHistOpts(h_E_e_Reco, "reco");
  setHistOpts(h_E_e_Corr, "corr");
  setHistOpts(h_E_e_Pur, "reco");
  //Single particle kinematics - photons
  setHistOpts(h_theta_g_Truth, "truth");
  setHistOpts(h_theta_g_Reco, "reco");
  setHistOpts(h_theta_g_Corr, "corr");
  setHistOpts(h_theta_g_Pur, "reco");
  setHistOpts(h_E_g_Truth, "truth");
  setHistOpts(h_E_g_Reco, "reco");
  setHistOpts(h_E_g_Corr, "corr");
  setHistOpts(h_E_g_Pur, "reco");
  // HFS quantities
  setHistOpts(h_HFSSigma_Truth, "truth");
  setHistOpts(h_HFSSigma_Reco, "reco");
  setHistOpts(h_HFSSigma_Corr, "corr");
  setHistOpts(h_HFSSigma_Pur, "reco");
  setHistOpts(h_HFSpT2_Truth, "truth");
  setHistOpts(h_HFSpT2_Reco, "reco");
  setHistOpts(h_HFSpT2_Corr, "corr");
  setHistOpts(h_HFSpT2_Pur, "reco");
  setHistOpts(h_FullSigma_Truth, "truth");
  setHistOpts(h_FullSigma_Reco, "reco");
  setHistOpts(h_FullSigma_Corr, "corr");
  setHistOpts(h_FullSigma_Pur, "reco");
  // E/p for electron/photon
  setHistOpts(h_EOverp_e_Truth, "truth");
  setHistOpts(h_EOverp_e_Reco, "reco");
  setHistOpts(h_EOverp_g_Truth, "truth");
  setHistOpts(h_EOverp_g_Reco, "reco");
    
  //--------------------------------------------------------------------------------------
  // Draw histograms
  //--------------------------------------------------------------------------------------
  // Draw histograms
  gStyle->SetOptStat(00000000);
  // Set header strings
  TString sHead1 = "#it{#bf{ePIC} Simulation}";
  TString sHead2 = "#it{ep " + energy + " GeV}";
  // Create header text objects
  TLatex* tHead1 = new TLatex(0.10, 0.91, sHead1);
  tHead1->SetNDC();
  tHead1->SetTextSize(30);
  tHead1->SetTextFont(43);
  tHead1->SetTextColor(kBlack);
  TLatex* tHead2 = new TLatex(0.74, 0.91, sHead2);
  tHead2->SetNDC();
  tHead2->SetTextSize(30);
  tHead2->SetTextFont(43);
  tHead2->SetTextColor(kBlack);

  // CANVAS 1: Q2
  TCanvas* cQ2 = new TCanvas("cQ2","",1200,800);
  cQ2->Divide(2,2);
  cQ2->cd(1);
  h_Q2_Truth->GetXaxis()->SetTitle("Q^{2} [(GeV/#it{c}^{2})^{2}]");
  h_Q2_Truth->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Truth->Draw("ep");
  h_Q2_Reco->Draw("epsame");
  h_Q2_Corr->Draw("epsame");
  TLegend lQ2 = new TLegend(0.5,0.5,0.89,0.89);
  lQ2->SetHeader("Q^{2}; all electrons");
  lQ2->AddEntry(h_Q2_Truth,"MC truth","pl");
  lQ2->AddEntry(h_Q2_Reco,"Reco.","pl");
  lQ2->AddEntry(h_Q2_Corr,"Reco., #epsilon corr.","pl");
  cQ2->cd(2);
  h_Q2_Resp->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Resp->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Resp->GetXaxis()->SetMarkerSize(0.05);
  h_Q2_Resp->GetYaxis()->SetTitleSize(0.05);
  h_Q2_Resp->GetYaxis()->SetTitleOffset(0.85);
  h_Q2_Resp->GetYaxis()->SetMarkerSize(0.05);
  h_Q2_Resp->Draw("colz");
  cQ2->cd(3);
  h_Q2_Pur->GetXaxis()->SetTitleSize(0.05);
  h_Q2_Pur->GetXaxis()->SetTitleOffset(0.9);
  h_Q2_Pur->GetXaxis()->SetLabelSize(0.05);
  h_Q2_Pur->GetYaxis()->SetTitleSize(0.05);
  h_Q2_Pur->GetYaxis()->SetTitleOffset(0.88);
  h_Q2_Pur->Draw("ep");
  cQ2->cd(4);
  h_dQ2vQ2->GetXaxis()->SetTitleSize(0.05);
  h_dQ2vQ2->GetXaxis()->SetTitleOffset(0.9);
  h_dQ2vQ2->GetXaxis()->SetLabelSize(0.05);
  h_dQ2vQ2->GetYaxis()->SetTitleSize(0.05);
  h_dQ2vQ2->GetYaxis()->SetTitleOffset(0.8);
  h_dQ2vQ2->GetYaxis()->SetLabelSize(0.05);
  h_dQ2vQ2->Draw("colz");

  //cQ2->Print("DVCSQA_temp001.pdf");
  //cQ2->Close();

  TCanvas* cxB = new TCanvas("cxB","",1200,800);
  cxB->Divide(2,2);
  cxB->cd(1);
  h_xB_Truth->GetXaxis()->SetTitle("Q^{2} [(GeV/#it{c}^{2})^{2}]");
  h_xB_Truth->GetXaxis()->SetTitleSize(0.05);
  h_xB_Truth->GetXaxis()->SetTitleOffset(0.9);
  h_xB_Truth->Draw("ep");
  h_xB_Reco->Draw("epsame");
  h_xB_Corr->Draw("epsame");
  cxB->cd(2);
  h_xB_Resp->Draw("colz");
  cxB->cd(3);
  h_xB_Pur->Draw("p");
  cxB->cd(4);
  h_dxBvxB->Draw("colz");

  // Combine PDFs into one and clean up
  /*std::cout<<"...Cleaning up files..."<<std::endl;
  TString filePlots = "$EIC_WORK_DIR/DVCS_Analysis/Plots/DVCSPlots_" + campaign + "_" + energy + "_QA.pdf";
  std::cout<<"Moving plots to "<<filePlots<<std::endl;
  TString pdfUniteCmd = "pdfunite DVCSQA*.pdf "+filePlots;
  gSystem->Exec(pdfUniteCmd);
  gSystem->Exec("rm DVCSQA_temp*.pdf");*/
  
  return;
}
