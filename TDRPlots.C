using namespace std;

#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>


//---------------------------------------------------------------------
// MAIN
//---------------------------------------------------------------------
void TDRPlots(TString campaign, TString energy, TString setting){
  //---------------------------------------------------------------------
  // Get data file from DVCSAnalysis script
  //---------------------------------------------------------------------
  TString inFileName = "$EIC_WORK_DIR/DVCS_Analysis/RootFiles/ePIC_DVCS_" + campaign + "_" + energy + "_" + setting + ".root";
  //TString inFileName = "/scratch/oliver/ePIC_DVCS_" + campaign + "_" + energy + "_" + setting + ".root";
  cout<<"Input file: "<<inFileName<<endl;
  TFile* inFile = new TFile(inFileName);

  //---------------------------------------------------------------------
  // Extract histograms
  //---------------------------------------------------------------------
  // Particle kinematics - MC generated protons, electrons and photons
  TH1D* h_eta_MCp   = (TH1D*)inFile->Get("eta_MCp");
  TH1D* h_eta_MCe   = (TH1D*)inFile->Get("eta_MCe");
  TH1D* h_eta_MCg   = (TH1D*)inFile->Get("eta_MCg");
  // DVCS event kinematics - MC generated particles
  TH1D* h_t_MC    = (TH1D*)inFile->Get("t_MC");
  // Particle kinematics - Associated MC protons, electrons and photons
  TH1D* h_eta_MCAp   = (TH1D*)inFile->Get("eta_MCAp");
  TH1D* h_eta_MCAe   = (TH1D*)inFile->Get("eta_MCAe");
  TH1D* h_eta_MCAg   = (TH1D*)inFile->Get("eta_MCAg");
  // DVCS event kinematics - Associated MC particles
  TH1D* h_t_MCA    = (TH1D*)inFile->Get("t_MCA");
  // Particle kinematics - Reconstructed protons (B0 ONLY), electrons and photons
  TH1D* h_eta_RPp   = (TH1D*)inFile->Get("eta_RPp");
  TH1D* h_eta_RPe   = (TH1D*)inFile->Get("eta_RPe");
  TH1D* h_eta_RPg   = (TH1D*)inFile->Get("eta_RPg");
  // DVCS event kinematics - Reconstructed particles (BARREL/B0 ONLY)
  TH1D* h_t_RP    = (TH1D*)inFile->Get("t_RP");
  // Particle kinematics - Reconstructed protons (ROMAN POTS)
  TH1D* h_eta_RPPp   = (TH1D*)inFile->Get("eta_RPPp");
  // DVCS event kinematics - Reconstructed particles (ROMAN POTS ONLY)
  TH1D* h_t_RPP    = (TH1D*)inFile->Get("t_RPP");

  // 2D distibutions of kinematic variables (MC vs recon.)
  TH2D* h_p_2d = (TH2D*)inFile->Get("h_p_2d");
  TH2D* h_pt_2d = (TH2D*)inFile->Get("h_pt_2d");
  TH2D* h_t_2d = (TH2D*)inFile->Get("h_t_2d");
  TH2D* h_tRes_2d = (TH2D*)inFile->Get("h_tres_2d");
  TH1D* h_t_resolution = (TH1D*)inFile->Get("extracted_t_resolution");

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
  
  // Photon resolutions
  TH1D* h_PhotRes_E = (TH1D*)inFile->Get("photres_E");
  TH1D* h_PhotRes_theta = (TH1D*)inFile->Get("photres_theta");
  TH2D* h_PhotRes2D_theta = (TH2D*)inFile->Get("photres2d_theta");

  //--------------------------------------------------------------------------------------
  // Set histogram line colours: MC raw - black; MC associated - red; Reconstructed - blue
  //--------------------------------------------------------------------------------------
  // MC generated
  h_t_MC->SetLineColor(kBlack);
  h_Emiss3_MC->SetLineColor(kBlack);
  h_Pmiss3_MC->SetLineColor(kBlack);
  h_Ptmiss3_MC->SetLineColor(kBlack);
  h_M2miss3_MC->SetLineColor(kBlack);
  h_M2miss2ep_MC->SetLineColor(kBlack);
  h_M2miss2eg_MC->SetLineColor(kBlack);
  // Associated MC
  h_t_MCA->SetLineColor(kRed);
  h_Emiss3_MCA->SetLineColor(kRed);
  h_Pmiss3_MCA->SetLineColor(kRed);
  h_Ptmiss3_MCA->SetLineColor(kRed);
  h_M2miss3_MCA->SetLineColor(kRed);
  h_M2miss2ep_MCA->SetLineColor(kRed);
  h_M2miss2eg_MCA->SetLineColor(kRed);
  // Reconstructed
  h_t_RP->SetLineColor(kBlue);
  h_t_RP->SetMarkerColor(kBlue);
  h_t_RPP->SetLineColor(kCyan+1);
  h_t_RPP->SetMarkerColor(kCyan+1);
  h_Emiss3_RP->SetLineColor(kBlue);
  h_Pmiss3_RP->SetLineColor(kBlue);
  h_Ptmiss3_RP->SetLineColor(kBlue);
  h_M2miss3_RP->SetLineColor(kBlue);
  h_M2miss2ep_RP->SetLineColor(kBlue);
  h_M2miss2eg_RP->SetLineColor(kBlue);

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

  // CANVAS 1: T-DISTRIBUTION
  TCanvas* c1 = new TCanvas("c1","",1200,800);
  h_t_MC->SetMinimum(1);
  gPad->SetLogy(1);
  h_t_MC->GetXaxis()->SetTitle("|t| [GeV^{2}]");
  h_t_MC->GetYaxis()->SetTitle("Counts / 0.02 GeV^{2}");
  h_t_MC->SetLineWidth(2);
  h_t_MC->Draw();
  h_t_RP->SetMarkerStyle(20);
  h_t_RP->Draw("pesame");
  h_t_RPP->SetMarkerStyle(20);
  h_t_RPP->Draw("pesame");
  // Add header text
  tHead1->Draw("same");
  tHead2->Draw("same");
  // Add legend
  TLatex* tC1 = new TLatex(0.58, 0.83, "#splitline{#bf{EpIC} ep #rightarrow e'p'#gamma, Q^{2} #geq 1 GeV^{2}}{t_{RP} #leq 0.3 GeV^{2}, M_{miss}^{2} < 1 GeV^{2}}");
  tC1->SetNDC();
  tC1->SetTextSize(30);
  tC1->SetTextFont(43);
  tC1->SetTextColor(kBlack);
  tC1->Draw("same");
  TLegend* lC1 = new TLegend(0.57, 0.6, 0.8, 0.77);
  lC1->SetLineColorAlpha(kWhite,0);
  lC1->SetFillColorAlpha(kWhite,0);
  lC1->AddEntry(h_t_MC, "#bf{EpIC} MC gen.", "l");
  lC1->AddEntry(h_t_RP, "Reco. B0", "lp");
  lC1->AddEntry(h_t_RPP, "Reco. RP", "lp");
  lC1->Draw();
 
  //c1->Print("$EIC_WORK_DIR/DVCS_ANALYSIS/Plots/TDRPlots/DVCS_t.pdf");
  //c1->Close();

  // CANVAS 2: DETECTOR OCCUPANCIES (FULL REACTION)
  TCanvas* c2 = new TCanvas("c2","",1200,800);
  gPad->SetLogy();
  TH1F* h_deteta_MCe = (TH1F*)h_eta_MCe->Clone("h_deteta_MCe");
  TH1F* h_deteta_MCg = (TH1F*)h_eta_MCg->Clone("h_deteta_MCg");
  TH1F* h_deteta_MCp = (TH1F*)h_eta_MCp->Clone("h_deteta_MCp");
  TH1F* h_deteta_RPe = (TH1F*)h_eta_RPe->Clone("h_deteta_RPe");
  TH1F* h_deteta_RPg = (TH1F*)h_eta_RPg->Clone("h_deteta_RPg");
  TH1F* h_deteta_RPp = (TH1F*)h_eta_RPp->Clone("h_deteta_RPp");
  TH1F* h_deteta_RPPp = (TH1F*)h_eta_RPPp->Clone("h_deteta_RPPp");
  h_deteta_MCe->SetLineColor(kRed);
  h_deteta_MCe->SetLineWidth(2);
  h_deteta_MCg->SetLineColor(kGreen+1);
  h_deteta_MCg->SetLineWidth(2);
  h_deteta_MCp->SetLineColor(kBlue);
  h_deteta_MCp->SetLineWidth(2);
  h_deteta_RPe->SetLineColor(kRed);
  h_deteta_RPe->SetMarkerStyle(20);
  h_deteta_RPe->SetMarkerColor(kRed);
  h_deteta_RPg->SetLineColor(kGreen+1);
  h_deteta_RPg->SetMarkerColor(kGreen+1);
  h_deteta_RPg->SetMarkerStyle(20);
  h_deteta_RPp->SetLineColor(kBlue);
  h_deteta_RPp->SetMarkerStyle(20);
  h_deteta_RPp->SetMarkerColor(kBlue);
  h_deteta_RPPp->SetLineColor(kCyan+1);
  h_deteta_RPPp->SetMarkerStyle(20);
  h_deteta_RPPp->SetMarkerColor(kCyan+1);
  h_deteta_MCe->GetXaxis()->SetTitle("#eta");
  h_deteta_MCe->GetYaxis()->SetTitle("Counts");
  h_deteta_MCe->SetMaximum(10*h_deteta_MCe->GetMaximum());
  h_deteta_MCe->Draw("hist");
  h_deteta_MCg->Draw("hist same");
  h_deteta_MCp->Draw("hist same");
  h_deteta_RPe->Draw("pesame");
  h_deteta_RPg->Draw("pesame");
  h_deteta_RPp->Draw("pesame");
  h_deteta_RPPp->Draw("pesame");
  // Add header (and other) text
  tHead1->Draw("same");
  tHead2->Draw("same");
  TLatex* tC2t1 = new TLatex(0.52, 0.83, "#splitline{#theta_{p'} < 5.5 mrad for RP, 5.5 < #theta_{p'} < 20 mrad for B0}{Q^{2} #geq 1 GeV^{2}}");
  tC2t1->SetNDC();
  tC2t1->SetTextSize(23);
  tC2t1->SetTextFont(43);
  tC2t1->SetTextColor(kBlack);
  tC2t1->Draw("same");
  // Add legend
  TLegend* lC2 = new TLegend(0.12,0.6,0.55,0.88);
  lC2->SetLineColorAlpha(kWhite,0);
  lC2->SetLineWidth(0);
  lC2->SetFillStyle(0);
  lC2->SetFillColorAlpha(kWhite,0);
  lC2->SetHeader("#splitline{#bf{EpIC} ep #rightarrow e'p'#gamma, all identified tracks}{Lines: #bf{EpIC} generated, Points: Reconstructed}");
  //lC2->AddEntry((TObject*)0, "", "");
  lC2->AddEntry(h_deteta_RPe, "e'", "pl");
  lC2->AddEntry(h_deteta_RPg, "#gamma", "pl");
  lC2->AddEntry(h_deteta_RPp, "p' (in B0)", "pl");
  lC2->AddEntry(h_deteta_RPPp, "p' (in RP)", "pl");
  lC2->Draw();
  
  // CANVAS 2: DETECTOR OCCUPANCIES (SEPARATED BY SPECIES)
  TCanvas* c2a = new TCanvas("c2a","",1200,500);
  c2a->Divide(3,1,0,0);
  TLatex* tHeadLc2a = new TLatex(0.04, 0.91, sHead1);
  tHeadLc2a->SetNDC();
  tHeadLc2a->SetTextSize(25);
  tHeadLc2a->SetTextFont(43);
  tHeadLc2a->SetTextColor(kBlack);
  tHeadLc2a->Draw();
  TLatex* tHeadM1c2a = new TLatex(0.29, 0.91, "#bf{EpIC} ep #rightarrow e'p'#gamma");
  tHeadM1c2a->SetNDC();
  tHeadM1c2a->SetTextSize(25);
  tHeadM1c2a->SetTextFont(43);
  tHeadM1c2a->SetTextColor(kBlack);
  tHeadM1c2a->Draw("same");
  TLatex* tHeadM2c2a = new TLatex(0.59, 0.91, "Q^{2} #geq 1 GeV^{2}");
  tHeadM2c2a->SetNDC();
  tHeadM2c2a->SetTextSize(25);
  tHeadM2c2a->SetTextFont(43);
  tHeadM2c2a->SetTextColor(kBlack);
  tHeadM2c2a->Draw("same");
  TLatex* tHeadRc2a = new TLatex(0.83, 0.91, sHead2);
  tHeadRc2a->SetNDC();
  tHeadRc2a->SetTextSize(25);
  tHeadRc2a->SetTextFont(43);
  tHeadRc2a->SetTextColor(kBlack);
  tHeadRc2a->Draw("same");
  c2a->cd(1);
  gPad->SetLogy();
  h_eta_MCe->SetMaximum(10*h_eta_MCe->GetMaximum());
  h_eta_MCe->GetXaxis()->SetTitle("#eta");
  h_eta_MCe->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCe->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCe->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCe->GetXaxis()->SetLabelOffset(0.01);
  h_eta_MCe->GetXaxis()->CenterTitle();  
  h_eta_MCe->GetYaxis()->SetTitle("Counts");
  h_eta_MCe->GetYaxis()->SetTitleSize(0.05);
  h_eta_MCe->GetYaxis()->SetTitleOffset(0.95);
  h_eta_MCe->GetYaxis()->SetLabelSize(0.04);
  h_eta_MCe->GetYaxis()->SetLabelOffset(0.002);
  h_eta_MCe->SetLineColor(kBlack);
  h_eta_MCe->SetLineWidth(2);
  h_eta_RPe->SetLineColor(kBlack);
  h_eta_RPe->SetMarkerColor(kBlack);
  h_eta_RPe->SetMarkerStyle(20);
  h_eta_MCe->Draw();
  h_eta_RPe->Draw("pesame");
  // Add text
  TLatex* tHeadc2ap1 = new TLatex(0.48, 0.92, "e^{-}, p_{e'} #leq 11 GeV");
  tHeadc2ap1->SetNDC();
  tHeadc2ap1->SetTextSize(30);
  tHeadc2ap1->SetTextFont(43);
  tHeadc2ap1->SetTextColor(kBlack);
  tHeadc2ap1->Draw("same");
  // Add legend
  TLegend* lC2Ap1 = new TLegend(0.7,0.7,0.95,0.9);
  lC2Ap1->SetLineColorAlpha(kWhite,0);
  lC2Ap1->SetLineWidth(0);
  lC2Ap1->SetFillStyle(0);
  lC2Ap1->SetFillColorAlpha(kWhite,0);
  lC2Ap1->AddEntry(h_eta_MCe, "MC gen.", "pl");
  lC2Ap1->AddEntry(h_eta_RPe, "Reco.", "pl");
  lC2Ap1->Draw();
  c2a->cd(2);
  gPad->SetLogy();
  h_eta_MCg->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCg->GetXaxis()->SetTitle("#eta");
  h_eta_MCg->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCg->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCg->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCg->GetXaxis()->SetLabelOffset(0.01); 
  h_eta_MCg->SetLineColor(kRed);
  h_eta_MCg->SetLineWidth(2);
  h_eta_RPg->SetLineColor(kRed);
  h_eta_RPg->SetMarkerColor(kRed);
  h_eta_RPg->SetMarkerStyle(20);
  h_eta_MCg->Draw();
  h_eta_RPg->Draw("pesame");
  // Add text
  TLatex* tHeadc2ap2 = new TLatex(0.07, 0.92, "#gamma");
  tHeadc2ap2->SetNDC();
  tHeadc2ap2->SetTextSize(30);
  tHeadc2ap2->SetTextFont(43);
  tHeadc2ap2->SetTextColor(kBlack);
  tHeadc2ap2->Draw("same");
  c2a->cd(3);
  gPad->SetLogy();
  h_eta_MCp->SetMaximum(h_eta_MCe->GetMaximum()); // Match maxima across species
  h_eta_MCp->GetXaxis()->SetTitle("#eta");
  h_eta_MCp->GetXaxis()->SetTitleSize(0.06);
  h_eta_MCp->GetXaxis()->SetTitleOffset(0.8);
  h_eta_MCp->GetXaxis()->SetLabelSize(0.05);
  h_eta_MCp->GetXaxis()->SetLabelOffset(0.01);
  h_eta_MCp->SetLineColor(kBlue);
  h_eta_MCp->SetLineWidth(2);
  h_eta_RPp->SetLineColor(kViolet);
  h_eta_RPp->SetMarkerColor(kViolet);
  h_eta_RPp->SetMarkerStyle(20);
  h_eta_RPPp->SetLineColor(kCyan+1);
  h_eta_RPPp->SetMarkerColor(kCyan+1);
  h_eta_RPPp->SetMarkerStyle(20);
  h_eta_MCp->Draw();
  h_eta_RPp->Draw("pesame");
  h_eta_RPPp->Draw("pesame");
  // Add text
  TLatex* tHeadc2ap3 = new TLatex(0.07, 0.92, "p, p_{p'} #leq 110 GeV");
  tHeadc2ap3->SetNDC();
  tHeadc2ap3->SetTextSize(30);
  tHeadc2ap3->SetTextFont(43);
  tHeadc2ap3->SetTextColor(kBlack);
  tHeadc2ap3->Draw("same");
  // Add legend
  TLegend* lC2Ap3 = new TLegend(0.03,0.7,0.7,0.9);
  lC2Ap3->SetLineColorAlpha(kWhite,0);
  lC2Ap3->SetLineWidth(0);
  lC2Ap3->SetFillStyle(0);
  lC2Ap3->SetFillColorAlpha(kWhite,0);
  //lC2Ap3->AddEntry(h_eta_RPp, "B0 - 5.5 < #theta_{p'} < 20 mrad", "pl");
  lC2Ap3->AddEntry(h_eta_RPp, "Reco. B0", "pl");
  //lC2Ap3->AddEntry(h_eta_RPPp, "RP - #theta_{p'} < 5.5 mrad", "pl");
  lC2Ap3->AddEntry(h_eta_RPPp, "Reco. RP", "pl");
  lC2Ap3->Draw();
  
  
  //c2->Print("$EIC_WORK_DIR/DVCS_ANALYSIS/Plots/TDRPlots/DVCS_eta.pdf");
  //c2->Print("/scratch/DVCS_eta.pdf");
  //c2->Close();

  // CANVAS 3: PHOTON ANGLE RESOLUTION
  TCanvas* c3 = new TCanvas("c3","",1500,600);
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLogy();
  h_PhotRes_theta->SetLineColor(kBlack);
  h_PhotRes_theta->SetMarkerColor(kBlack);
  h_PhotRes_theta->SetMarkerStyle(20);
  h_PhotRes_theta->GetXaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetXaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->GetYaxis()->SetTitle("Counts / 5mrad");
  h_PhotRes_theta->GetYaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetYaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->Draw("pe");
  // Add header text
  TLatex* tHeadLc3 = new TLatex(0.1, 0.91, sHead1);
  tHeadLc3->SetNDC();
  tHeadLc3->SetTextSize(25);
  tHeadLc3->SetTextFont(43);
  tHeadLc3->SetTextColor(kBlack);
  tHeadLc3->Draw("same");
  c3->cd(2);
  h_PhotRes2D_theta->GetXaxis()->SetTitleSize(0.05);
  h_PhotRes2D_theta->GetXaxis()->SetTitleOffset(0.9);
  h_PhotRes2D_theta->GetYaxis()->SetTitleSize(0.05);
  h_PhotRes2D_theta->GetYaxis()->SetTitleOffset(0.9);
  h_PhotRes2D_theta->GetZaxis()->SetMaxDigits(2);
  h_PhotRes2D_theta->Draw("colz");
  // Add header text
  TLatex* tHeadRc3 = new TLatex(0.68, 0.91, sHead2);
  tHeadRc3->SetNDC();
  tHeadRc3->SetTextSize(25);
  tHeadRc3->SetTextFont(43);
  tHeadRc3->SetTextColor(kBlack);
  tHeadRc3->Draw("same");

  // CANVAS 3A: PHOTON ANGLE RESOLUTION (SUBFIGURE)
  TCanvas* c3a = new TCanvas("c3a","",1200,800);
  // Draw 1D plot as main figure
  h_PhotRes_theta->GetXaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetXaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->GetYaxis()->SetTitle("Counts / 5mrad");
  h_PhotRes_theta->GetYaxis()->SetTitleSize(0.05);
  h_PhotRes_theta->GetYaxis()->SetTitleOffset(0.9);
  h_PhotRes_theta->SetLineColor(kBlack);
  h_PhotRes_theta->SetMarkerColor(kBlack);
  h_PhotRes_theta->SetMarkerStyle(20);
  h_PhotRes_theta->Draw("pe");
  // Create subfigure TPad on top for 2D plot
  TPad* p3a = new TPad("p3a","",0.51,0.49,0.89,0.89);
  p3a->SetRightMargin(0.01);
  p3a->SetLeftMargin(0.15);
  p3a->SetTopMargin(0.02);
  p3a->SetBottomMargin(0.15);
  p3a->Draw("same");
  TLatex* tHeadLc3a = new TLatex(0.1, 0.94, sHead1);
  tHeadLc3a->SetNDC();
  tHeadLc3a->SetTextSize(30);
  tHeadLc3a->SetTextFont(43);
  tHeadLc3a->SetTextColor(kBlack);
  tHeadLc3a->Draw("same");
  tHead2->Draw("same");
  p3a->cd();
  h_PhotRes2D_theta->RebinX(2);
  h_PhotRes2D_theta->RebinY(8);
  h_PhotRes2D_theta->GetXaxis()->SetTitleSize(0.08);
  h_PhotRes2D_theta->GetXaxis()->SetTitleOffset(0.85);
  h_PhotRes2D_theta->GetXaxis()->SetLabelSize(0.06);
  h_PhotRes2D_theta->GetYaxis()->SetTitle("#theta_{#gamma}(Reco)-#theta_{#gamma}(MC) [rad]");
  h_PhotRes2D_theta->GetYaxis()->SetTitleSize(0.08);
  h_PhotRes2D_theta->GetYaxis()->SetTitleOffset(0.85);
  h_PhotRes2D_theta->GetYaxis()->SetLabelSize(0.06);
  h_PhotRes2D_theta->GetZaxis()->SetMaxDigits(2);
  h_PhotRes2D_theta->Draw("same");

  // CANVAS 4: T RESOLUTION (SUBFIGURE)
  TCanvas* c4 = new TCanvas("c4","",1200,800);
  // Draw 1D plot as main figure
  h_t_resolution->SetLineColor(kBlack);
  h_t_resolution->SetMarkerColor(kBlack);
  h_t_resolution->SetMarkerStyle(20);
  h_t_resolution->SetTitle("");
  h_t_resolution->GetXaxis()->SetTitleSize(0.05);
  h_t_resolution->GetXaxis()->SetTitleOffset(0.8);
  h_t_resolution->GetYaxis()->SetTitleSize(0.05);
  h_t_resolution->GetYaxis()->SetTitleOffset(0.9);
  h_t_resolution->Draw();
  // Add text
  tHeadLc3a->Draw("same");
  tHead2->Draw("same");
  // Create subfigure TPad on top for 2D plot
  TPad* p4 = new TPad("p4","",0.4,0.3,0.89,0.89);
  p4->SetRightMargin(0.01);
  p4->SetLeftMargin(0.15);
  p4->SetTopMargin(0.02);
  p4->SetBottomMargin(0.15);
  p4->Draw("same");
  p4->cd();
  h_tRes_2d->GetXaxis()->SetTitleSize(0.07);
  h_tRes_2d->GetXaxis()->SetTitleOffset(0.85);
  h_tRes_2d->GetXaxis()->SetLabelSize(0.05);
  h_tRes_2d->GetYaxis()->SetTitleSize(0.07);
  h_tRes_2d->GetYaxis()->SetTitleOffset(0.75);
  h_tRes_2d->GetYaxis()->SetLabelSize(0.05);
  h_tRes_2d->GetYaxis()->SetRangeUser(0,12);
  h_tRes_2d->Draw();

  return;
}
